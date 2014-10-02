
'''Definitions related to reaction equations and parsing of such equations'''

import re
import functools
from decimal import Decimal

from .expression.affine import Expression

@functools.total_ordering
class Compound(object):
    '''Represents a compound in a reaction equation

    A compound is a named entity in the reaction equations representing a
    chemical compound. A compound can represent a generalized chemical entity
    (e.g. polyphosphate) and the arguments can be used to instantiate a specific
    chemical entity (e.g. polyphosphate(3)) by passing a number as an argument
    or a partially specified entity by passing an expression (e.g. polyphosphate(n)).'''

    def __init__(self, name, compartment=None, arguments=()):
        self._name = str(name)
        self._compartment = None if compartment is None else str(compartment)
        self._arguments = tuple(arguments)

    @property
    def name(self):
        '''Name of compound'''
        return self._name

    @property
    def compartment(self):
        '''Compartment of compound'''
        return self._compartment

    @property
    def arguments(self):
        '''Expression argument for generalized compounds'''
        return self._arguments

    @property
    def id(self):
        '''Encode compound to unique identifier string

        >>> Compound('L-Glutamine', 'e').id
        'L_Glutamine_e'
        '''
        s = self._name + '_'.join(self._arguments)
        if self._compartment is not None:
            s += '_'+self._compartment
        if len(s) > 0 and s[0].isdigit():
            s = '_'+s
        return re.sub('[^\w]', '_', s)

    def translate(self, func):
        '''Translate compound name using given function

        >>> Compound('Pb').translate(lambda x: x.lower())
        Compound('pb')'''
        return self.__class__(func(self._name), self._compartment, self._arguments)

    def in_compartment(self, compartment):
        '''Return an instance of this compound in the specified compartment

        >>> Compound('H+').in_compartment('e')
        Compound('H+', 'e')
        '''
        return self.__class__(self._name, compartment, self._arguments)

    def __eq__(self, other):
        return (isinstance(other, Compound) and
            self._name == other._name and
            self._compartment == other._compartment and
            self._arguments == other._arguments)

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, Compound):
            return ((self._name, self._compartment, self._arguments) <
                    (other._name, other._compartment, other._arguments))
        return NotImplemented

    def __hash__(self):
        return (hash('Compound') ^ hash(self._name) ^
                hash(self._compartment) ^ hash(self._arguments))

    def __str__(self):
        '''String representation of compound

        >>> str(Compound('Phosphate'))
        'Phosphate'
        >>> str(Compound('Phosphate', 'e'))
        'Phosphate[e]'
        >>> str(Compound('Polyphosphate', None, [Expression('n')]))
        'Polyphosphate(n)'
        >>> str(Compound('Polyphosphate', 'p', [Expression('n')]))
        'Polyphosphate(n)[p]'
        '''
        s = self._name
        if len(self._arguments) > 0:
            s += '({})'.format(', '.join(str(a) for a in self._arguments))
        if self._compartment is not None:
            s += '[{}]'.format(self._compartment)
        return s

    def __repr__(self):
        def str_repr(*args):
            return 'Compound({})'.format(', '.join(repr(a) for a in args))

        if len(self._arguments) == 0:
            if self._compartment is None:
                return str_repr(self._name)
            return str_repr(self._name, self._compartment)
        return str_repr(self._name, self._compartment, self._arguments)

class Reaction(object):
    '''Reaction equation representation

    Each compound is associated with a stoichiometric value.'''

    Bidir = '<=>'
    Left = '<='
    Right = '=>'

    def __init__(self, direction, left, right):
        if direction not in (Reaction.Bidir, Reaction.Left, Reaction.Right):
            raise ValueError('Invalid direction: {}'.format(direction))
        self._direction = direction
        self._left = tuple(left)
        self._right = tuple(right)

    @property
    def direction(self):
        '''Direction of reaction equation'''
        return self._direction

    @property
    def left(self):
        '''Sequence of compounds on the left-hand side of the reaction equation'''
        return self._left

    @property
    def right(self):
        '''Sequence of compounds on the right-hand side of the reaction equation'''
        return self._right

    @property
    def compounds(self):
        '''Sequence of compounds on both sides of the reaction equation'''
        return self._left + self._right

    def normalized(self):
        '''Return normalized reaction

        The normalized reaction will have direction Bidir or Right.'''

        if self._direction == Reaction.Left:
            direction = Reaction.Right
            left = self._right
            right = self._left
        else:
            direction = self._direction
            left = self._left
            right = self._right

        return Reaction(direction, left, right)

    def translated_compounds(self, translate):
        '''Return reaction where compound names have been translated

        For each compound the translate function is called with the compound name
        and the returned value is used as the new compound name. A new reaction is
        returned with the substituted compound names.'''

        left = ((compound.translate(translate), count) for compound, count in self._left)
        right = ((compound.translate(translate), count) for compound, count in self._right)

        return Reaction(self._direction, left, right)

    def format(self, formatter=None):
        '''Format reaction as string using specified formatter'''

        if formatter is None:
            formatter = ModelSEED
        return formatter.format(self)

    def copy(self):
        '''Returns a distinct copy as a new Reaction object'''
        return self.__class__(self._direction, self._left, self._right)

    def __str__(self):
        return self.format()

    def __repr__(self):
        return 'Reaction({}, {}, {})'.format(repr(self._direction), repr(self._left), repr(self._right))

    def __eq__(self, other):
        '''Indicate equality of self and other'''
        return (self._direction == other._direction and
                self._left == other._left and
                self._right == other._right)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash('Reaction') ^ hash(self._direction) ^ hash(self._left) ^ hash(self._right)

class ParseError(Exception):
    '''Exception used to signal errors while parsing'''
    pass

class ModelSEED(object):
    '''Implementation of a parser for the textual representations of chemical
    reaction equations as used in ModelSEED.

    The representation does not seem to follow a published scheme so the
    following grammar has been reverse engineered from actual data sets. This
    parser is based on the derived grammer.::

        <reaction>     ::= <comp-list> ' ' <reaction-dir> ' ' <comp-list>
        <reaction-dir> ::= '<=' | '<=>' | '=>' | '?' | ''
        <comp-list>    ::= '' | <compound> | <compound> ' + ' <comp-list>
        <compound>     ::= <comp-count> ' ' <comp-spec> | <comp-spec>
        <comp-count>   ::= '(' <comp-number> ')' | <comp-number>
        <comp-number>  ::= <decimal>
        <comp-spec>    ::= '|' <comp-id> '|' | 'cdp' <cdp-id>
        <comp-id>      ::= <comp-name> '[' <comp-compart> ']' | <comp-name>
        <comp-compart> ::= <alpha>
        <comp-name>    ::= <any characters other than "|">
        <cdp-id>       ::= <five digits>   ; [sic]

    Note that the derived grammar is quite sloppy and could easily follow more
    strict rules that would make parsing easier. When converting parsed reactions
    back to string represetation using this module only a subset of the grammar
    is used.::

        <reaction>     ::= <comp-list> ' ' <reaction-dir> ' ' <comp-list>
        <reaction-dir> ::= '<=>' | '=>' | '?'
        <comp-list>    ::= '' | <compound> | <compound> ' + ' <comp-list>
        <compound>     ::= <comp-count> ' ' <comp-spec> | <comp-spec>
        <comp-count>   ::= '(' <decimal> ')'
        <comp-spec>    ::= '|' <comp-name> '|' | '|' <comp-name> '[' <comp-compart> ']' '|'
        <comp-compart> ::= <alpha>
        <comp-name>    ::= <any characters other than "|">
    '''

    @classmethod
    def parse(cls, s):
        '''Parse a reaction string'''

        tokens = list(cls._tokenize(s))
        direction = None
        for i, t in enumerate(tokens):
            if t in ('<=', '<=>', '=>', '?', ''):
                direction = t
                left = tokens[:i]
                right = tokens[i+1:]

        if direction is None:
            raise ParseError('Failed to parse reaction: {}'.format(tokens))

        if direction in ('?', ''):
            direction = Reaction.Bidir

        left = cls._parse_compound_list(left)
        right = cls._parse_compound_list(right)
        return Reaction(direction, left, right)

    @classmethod
    def format(cls, rx):
        '''Format reaction as string'''

        # Define helper functions
        def format_compound(compound, count):
            '''Format compound'''
            cpdspec = str(compound)
            if count != 1:
                return '({}) |{}|'.format(count, cpdspec)
            return '|{}|'.format(cpdspec)

        def format_compound_list(cmpds):
            '''Format compound list'''
            return ' + '.join(format_compound(compound, count) for compound, count in cmpds)

        return '{} {} {}'.format(format_compound_list(rx.left),
                                 '?' if rx.direction == '' else rx.direction,
                                 format_compound_list(rx.right))

    @classmethod
    def _tokenize(cls, s):
        '''Return tokens of reaction string'''

        s = s.lstrip()
        token = ''
        barquote = False
        for c in s:
            if c.isspace() and not barquote:
                yield token
                token = ''
                continue

            token += c

            if c == '|':
                barquote = not barquote

        if token != '':
            yield token

    @classmethod
    def _parse_compound_number(cls, number):
        '''Parse compound number

        Return plain int if possible, otherwise use Decimal.'''

        d = Decimal(number)
        if d % 1 == 0:
            return int(d)
        return d

    @classmethod
    def _parse_compound_count(cls, count):
        '''Parse compound count'''

        m = re.match(r'^\((.+)\)|(.+)$', count)
        if not m:
            raise ParseError('Unable to parse compound count: {}'.format(count))

        number = m.group(1) if m.group(1) is not None else m.group(2)
        return cls._parse_compound_number(number)

    @classmethod
    def _parse_compound_name(cls, name):
        '''Parse compound name'''

        m = re.match(r'^\|(.+)\||(cdp\d+)$', name)
        if not m:
            raise ParseError('Unable to parse compound name: {}'.format(name))

        return m.group(1) if m.group(1) is not None else m.group(2)

    @classmethod
    def _parse_compound(cls, cmpd):
        '''Parse compound'''

        count = 1
        if len(cmpd) == 2:
            count = cls._parse_compound_count(cmpd[0])
            name = cls._parse_compound_name(cmpd[1])
        elif len(cmpd) == 1:
            name = cls._parse_compound_name(cmpd[0])
        else:
            raise ParseError('Unexpected number of tokens in compound: {}'.format(cmpd))

        compartment = None
        m = re.match(r'^(.+?)\[(.)\]$', name)
        if m is not None:
            name = m.group(1)
            compartment = m.group(2)

        return Compound(name, compartment=compartment), count

    @classmethod
    def _parse_compound_list(cls, cmpds):
        '''Parse a list of compounds'''

        if len(cmpds) == 0:
            return

        cmpd = []
        for t in cmpds:
            if t == '+':
                yield cls._parse_compound(cmpd)
                cmpd = []
            else:
                cmpd.append(t)

        if len(cmpd) == 0:
            raise ParseError('Expected compound in compound list')

        yield cls._parse_compound(cmpd)

class SudenSimple(object):
    '''Parser for the simple reaction format in Suden model'''

    @classmethod
    def parse(cls, s):
        '''Parse a reaction string'''

        def parse_compound_list(s):
            for cpd in s.split('+'):
                if cpd == '':
                    continue
                count, spec = cpd.strip().split(' ')
                spec_split = spec.split('[')
                if len(spec_split) == 2:
                    # compartment
                    cpdid = spec_split[0]
                    comp = spec_split[1].rstrip(']')
                else:
                    cpdid = spec
                    comp = None
                d = Decimal(count)
                if d % 1 == 0:
                    d = int(d)
                yield Compound(cpdid, compartment=comp), d

        cpd_left, cpd_right = s.split('<=>')
        left = list(parse_compound_list(cpd_left))
        right = list(parse_compound_list(cpd_right))

        return Reaction('<=>', left, right)

class MetNet(object):
    '''Parser for the reaction format in MetNet model'''

    @classmethod
    def parse(cls, s):
        '''Parse a reaction string'''

        def parse_compound_list(s, global_comp):
            if s.strip() == '':
                return

            for cpd in s.strip().split('+'):
                cpdsplit = cpd.strip().split(') ')
                if len(cpdsplit) == 2:
                    count = Decimal(cpdsplit[0].lstrip('('))
                    if count % 1 == 0:
                        count = int(count)
                    cpdspec = cpdsplit[1]
                else:
                    count = 1
                    cpdspec = cpd.strip()

                if global_comp is None:
                    cpd_spec_split = cpdspec.split('[')
                    if len(cpd_spec_split) == 2:
                        comp = cpd_spec_split[1].rstrip(']')
                        cpdid = cpd_spec_split[0]
                    else:
                        cpdid = cpd_spec_split
                        comp = None
                else:
                    cpdid = cpdspec
                    comp = global_comp

                yield Compound(cpdid, compartment=comp), count

        # Split by colon for compartment information
        eq_split = s.split(':')

        if len(eq_split) == 2:
            global_comp, second = eq_split
            global_comp = global_comp.strip('[] ')
        else:
            second = eq_split[0]
            global_comp = None

        # Split by equation arrow
        direction = '=>'
        left_right = second.split('-->')
        if len(left_right) == 1:
            direction = '<=>'
            cpd_left, cpd_right = second.split('<==>')
        else:
            cpd_left, cpd_right = left_right

        # Remove incompatible characters from compound id
        def translate(cpdid):
            return 'cpd_' + cpdid.replace(',', '__')

        # Split by +
        left = ((compound.translate(translate), count) for compound, count in parse_compound_list(cpd_left, global_comp))
        right = ((compound.translate(translate), count) for compound, count in parse_compound_list(cpd_right, global_comp))

        return Reaction(direction, left, right)

class KEGG(object):
    '''Parser for the reaction format in KEGG'''

    @classmethod
    def parse(cls, s):
        '''Parse a KEGG reaction string'''

        def parse_count(s):
            m = re.match(r'^\((.+)\)$', s)
            if m is not None:
                s = m.group(1)

            m = re.match(r'^\d+$', s)
            if m is not None:
                return int(m.group(0))

            return Expression(s)

        def parse_compound(s):
            m = re.match(r'(.+)\((.+)\)', s)
            if m is not None:
                return Compound(m.group(1), arguments=[Expression(m.group(2))])
            return Compound(s)

        def parse_compound_list(s):
            for cpd in s.split(' + '):
                if cpd == '':
                    continue

                fields = cpd.strip().split(' ')
                if len(fields) > 2:
                    raise ParseError('Malformed compound specification: {}'.format(cpd))
                if len(fields) == 1:
                    count = 1
                    compound = parse_compound(fields[0])
                else:
                    count = parse_count(fields[0])
                    compound = parse_compound(fields[1])

                yield compound, count

        cpd_left, cpd_right = s.split('<=>')
        left = list(parse_compound_list(cpd_left.strip()))
        right = list(parse_compound_list(cpd_right.strip()))

        return Reaction('<=>', left, right)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
