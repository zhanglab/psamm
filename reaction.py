
'''Implementation of a parser for the textual representations of chemical
reaction equations as used in ModelSEED.

>>> parse('|H2O| + |PPi| => (2) |Phosphate| + (2) |H+|')
('=>', [('H2O', 1, None), ('PPi', 1, None)], [('Phosphate', 2, None), ('H+', 2, None)])

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

import re
from decimal import Decimal

# TODO This implementation is simply a hand-written parser. A proper parser
# could be generated using a parser generator and the grammar described above.

class ParseError(Exception):
    '''Exception used to signal errors while parsing'''
    pass

def tokenize(rx):
    '''Return tokens of reaction string'''

    rx = rx.lstrip()
    token = ''
    barquote = False
    for s in rx:
        if s.isspace() and not barquote:
            yield token
            token = ''
            continue

        token += s

        if s == '|':
            barquote = not barquote

    if token != '':
        yield token

def parse_compound_number(number):
    '''Parse compound number'''
    return Decimal(number)

def parse_compound_count(count):
    '''Parse compound count'''

    m = re.match(r'^\((.+)\)|(.+)$', count)
    if not m:
        raise ParseError('Unable to parse compound count: {}'.format(count))

    number = m.group(1) if m.group(1) is not None else m.group(2)
    return parse_compound_number(number)

def parse_compound_name(name):
    '''Parse compound name'''

    m = re.match(r'^\|(.+)\||(cdp\d+)$', name)
    if not m:
        raise ParseError('Unable to parse compound name: {}'.format(name))

    return m.group(1) if m.group(1) is not None else m.group(2)

def parse_compound(cmpd):
    '''Parse compound'''

    count = 1
    if len(cmpd) == 2:
        count = parse_compound_count(cmpd[0])
        name = parse_compound_name(cmpd[1])
    elif len(cmpd) == 1:
        name = parse_compound_name(cmpd[0])
    else:
        raise ParseError('Unexpected number of tokens in compound: {}'.format(cmpd))

    compartment = None
    m = re.match(r'^(.+?)\[(.)\]$', name)
    if m is not None:
        name = m.group(1)
        compartment = m.group(2)

    return (name, count, compartment)

def parse_compound_list(cmpds):
    '''Parse a list of compounds'''

    if len(cmpds) == 0:
        return

    cmpd = []
    for t in cmpds:
        if t == '+':
            yield parse_compound(cmpd)
            cmpd = []
        else:
            cmpd.append(t)

    if len(cmpd) == 0:
        raise ParseError('Expected compound in compound list')

    yield parse_compound(cmpd)

def parse(rx):
    '''Parse a reaction string'''

    tokens = list(tokenize(rx))
    direction = None
    for i, t in enumerate(tokens):
        if t in ('<=', '<=>', '=>', ''):
            direction = t
            left = tokens[:i]
            right = tokens[i+1:]

    if direction is None:
        raise ParseError('Failed to parse reaction: {}'.format(tokens))

    return (direction, list(parse_compound_list(left)),
            list(parse_compound_list(right)))

def normalize(rx):
    '''Normalize reaction by turning all left-directed reactions'''

    direction, left, right = rx

    if direction == '<=':
        direction = '=>'
        left = rx[2]
        right = rx[1]

    return (direction, left, right)

def translate_compounds(rx, translate):
    '''Translate compound names using translate function'''

    direction, left, right = rx

    return (direction, [(translate(name), count, compartment) for name, count, compartment in left],
            [(translate(name), count, compartment) for name, count, compartment in right])

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

def format(rx):
    '''Format parsed reaction using subset grammar'''
    direction, left, right = normalize(rx)
    return '{} {} {}'.format(format_compound_list(left),
                             '?' if direction == '' else direction,
                             format_compound_list(right))
