
'''Definitions related to reaction equations and parsing of such equations'''

import re
from decimal import Decimal

class Reaction(object):
    '''Reaction equation representation'''

    def __init__(self, direction, left, right):
        self.direction = direction
        self.left = left
        self.right = right

    def normalized(self):
        '''Return normalized reaction

        >>> Reaction('<=', [('Au', 1, None)], [('Pb', 1, None)]).normalized()
        Reaction('=>', [('Pb', 1, None)], [('Au', 1, None)])
        '''

        if self.direction == '<=':
            direction = '=>'
            left = list(self.right)
            right = list(self.left)
        else:
            direction = self.direction
            left = list(self.left)
            right = list(self.right)

        return Reaction(direction, left, right)

    def translated_compounds(self, translate):
        '''Return reaction where compound names have been translated

        For each compound the translate function is called with the compound name
        and the returned value is used as the new compound name. A new reaction is
        returned with the substituted compound names.

        >>> rx = Reaction('=>', [('Pb', 1, None)], [('Au', 1, None)])
        >>> rx.translated_compounds(lambda name: name.lower())
        Reaction('=>', [('pb', 1, None)], [('au', 1, None)])
        '''

        def translate_compound_list(l):
            return [(translate(cpdname), count, comp) for cpdname, count, comp in l]
        left = translate_compound_list(self.left)
        right = translate_compound_list(self.right)

        return Reaction(self.direction, left, right)

    def format(self, formatter=None):
        '''Format reaction as string using specified formatter

        >>> Reaction('=>', [('Pb', 1, None)], [('Au', 1, None)]).format()
        '|Pb| => |Au|'
        '''

        if formatter is None:
            formatter = ModelSEED
        return formatter.format(self)

    def __str__(self):
        return self.format()

    def __repr__(self):
        return 'Reaction({}, {}, {})'.format(repr(self.direction), repr(self.left), repr(self.right))

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
        '''Parse a reaction string

        >>> ModelSEED.parse('|H2O| + |PPi| => (2) |Phosphate| + (2) |H+|')
        Reaction('=>', [('H2O', 1, None), ('PPi', 1, None)], [('Phosphate', 2, None), ('H+', 2, None)])

        >>> ModelSEED.parse('|H2| + (0.5) |O2| => |H2O|')
        Reaction('=>', [('H2', 1, None), ('O2', Decimal('0.5'), None)], [('H2O', 1, None)])
        '''

        tokens = list(cls._tokenize(s))
        direction = None
        for i, t in enumerate(tokens):
            if t in ('<=', '<=>', '=>', ''):
                direction = t
                left = tokens[:i]
                right = tokens[i+1:]

        if direction is None:
            raise ParseError('Failed to parse reaction: {}'.format(tokens))

        return Reaction(direction, list(cls._parse_compound_list(left)),
                        list(cls._parse_compound_list(right)))

    @classmethod
    def format(cls, rx):
        '''Format reaction as string

        >>> ModelSEED.format(Reaction('<=', [('H2O', 2, None)], [('H2', 2, None), ('O2', 1, None)]))
        '(2) |H2| + |O2| => (2) |H2O|'
        '''

        # Define helper functions
        def format_compound(cmpd):
            '''Format compound'''
            name, count, compartment = cmpd

            if compartment is not None:
                spec = '{}[{}]'.format(name, compartment)
            else:
                spec = name

            if count != 1:
                return '({}) |{}|'.format(count, spec)
            return '|{}|'.format(spec)

        def format_compound_list(cmpds):
            '''Format compound list'''
            return ' + '.join(format_compound(cmpd) for cmpd in cmpds)

        rx = rx.normalized()
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

        return (name, count, compartment)

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
        '''Parse a reaction string

        >>> SudenSimple.parse('1 H2O + 1 PPi <=> 2 Phosphate + 2 proton')
        Reaction('<=>', [('H2O', 1, None), ('PPi', 1, None)], [('Phosphate', 2, None), ('proton', 2, None)])

        >>> SudenSimple.parse('1 H2 + 0.5 O2 <=> 1 H2O')
        Reaction('<=>', [('H2', 1, None), ('O2', Decimal('0.5'), None)], [('H2O', 1, None)])
        '''
        def parse_compound_list(s):
            cpds = []
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
                cpds.append((cpdid, d, comp))
            return cpds

        cpd_left, cpd_right = s.split('<=>')
        left = parse_compound_list(cpd_left)
        right = parse_compound_list(cpd_right)

        return Reaction('<=>', left, right)

class MetNet(object):
    '''Parser for the reaction format in MetNet model'''

    @classmethod
    def parse(cls, s):
        '''Parse a reaction string

        >>> MetNet.parse('[c] : akg + ala-L <==> glu-L + pyr')
        Reaction('<=>', [('cpd_akg', 1, None), ('cpd_ala-L', 1, None)], [('cpd_glu-L', 1, None), ('cpd_pyr', 1, None)])

        >>> MetNet.parse('(2) ficytcc553[c] + so3[c] + h2o[c] --> (2) focytcc553[c] + so4[c] + (2) h[e]')
        Reaction('=>', [('cpd_ficytcc553', 2, None), ('cpd_so3', 1, None), ('cpd_h2o', 1, None)], [('cpd_focytcc553', 2, None), ('cpd_so4', 1, None), ('cpd_h', 2, 'e')])
        '''

        def parse_compound_list(s, global_comp):
            cpds = []

            if s.strip() == '':
                return []

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

                # By convention we set the cytosol compartment to None
                if comp == 'c':
                    comp = None

                cpds.append((cpdid, count, comp))

            return cpds

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
        left = [(translate(cpdid), count, comp) for cpdid, count, comp in parse_compound_list(cpd_left, global_comp)]
        right = [(translate(cpdid), count, comp) for cpdid, count, comp in parse_compound_list(cpd_right, global_comp)]

        return Reaction(direction, left, right)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
