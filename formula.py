
'''Parser for chemical formulas'''

import pprint
import numbers
import re
from collections import defaultdict
from expression import Variable, Expression
import functools

class FormulaElement(object):
    def __add__(self, other):
        if isinstance(other, numbers.Number) and other == 0:
            return Formula({ self: 1 })
        elif isinstance(other, FormulaElement):
            if self == other:
                return Formula({ self: 2 })
            return Formula({ self: 1, other: 1 })
        elif isinstance(other, Formula):
            return other + self
        else:
            return NotImplemented

    def __radd__(self, other):
        if isinstance(other, numbers.Number) and other == 0:
            return self + other
        else:
            return NotImplemented

    def __mul__(self, other):
        if not isinstance(other, numbers.Number):
            return NotImplemented
        return Formula({ self: other })

    def __rmul__(self, other):
        return self * other

@functools.total_ordering
class Atom(FormulaElement):
    def __init__(self, symbol):
        self._symbol = symbol

    @property
    def symbol(self):
        return self._symbol

    def __hash__(self):
        return hash('Atom') ^ hash(self._symbol)

    def __eq__(self, other):
        return isinstance(other, Atom) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, Atom):
            return self._symbol < other._symbol
        return NotImplemented

    def __str__(self):
        return self._symbol

    def __repr__(self):
        return 'Atom({})'.format(repr(self._symbol))

class Radical(FormulaElement):
    def __init__(self, symbol):
        self._symbol = symbol

    @property
    def symbol(self):
        return self._symbol

    def __hash__(self):
        return hash('Radical') ^ hash(self._symbol)

    def __eq__(self, other):
        return isinstance(other, Radical) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        return self._symbol

    def __repr__(self):
        return 'Radical({})'.format(repr(self._symbol))

class Formula(object):
    def __init__(self, values={}):
        self._values = dict(values)

        self._variables = set()
        for element, value in self._values.iteritems():
            if callable(getattr(value, 'variables', None)):
                for var in value.variables():
                    self._variables.add(var)
            if callable(getattr(element, 'variables', None)):
                for var in element.variables():
                    self._variables.add(var)

    def substitute(self, **kwargs):
        result = self.__class__()
        for element, value in self._values.iteritems():
            if callable(getattr(value, 'substitute', None)):
                value = value.substitute(**kwargs)
                if isinstance(value, int) and value <= 0:
                    raise ValueError('Expression evaluated to non-positive number')
            if callable(getattr(element, 'substitute', None)):
                element = element.substitute(**kwargs)
            result += value * element
        return result

    def variables(self):
        return iter(self._variables)

    def is_variable(self):
        return len(self._variables) > 0

    def __call__(self, **kwargs):
        return self.substitute(**kwargs)

    def __str__(self):
        '''Return formula represented using Hill notation system'''
        def hill_sorted_elements(values):
            def element_sort_key(pair):
                element, value = pair
                if isinstance(element, Atom):
                    return 0, element.symbol
                elif isinstance(element, Radical):
                    return 1, element.symbol
                else:
                    return 2, None

            if Atom('C') in values:
                yield Atom('C'), values[Atom('C')]
                if Atom('H') in values:
                    yield Atom('H'), values[Atom('H')]
                for element, value in sorted(values.iteritems(), key=element_sort_key):
                    if element not in (Atom('C'), Atom('H')):
                        yield element, value
            else:
                for element, value in sorted(values.iteritems(), key=element_sort_key):
                    yield element, value

        s = ''
        for element, value in hill_sorted_elements(self._values):
            def grouped(element, value):
                return '({}){}'.format(element, value if value != 1 else '')
            def nongrouped(element, value):
                return '{}{}'.format(element, value if value != 1 else '')
            if isinstance(element, Radical):
                s += nongrouped(element, value) if len(element.symbol) == 1 else grouped(element, value)
            elif isinstance(element, Atom):
                s += nongrouped(element, value)
            else:
                s += grouped(element, value)
        return s

    def __repr__(self):
        return 'Formula({})'.format(pprint.pformat(self._values))

    def __add__(self, other):
        '''Sum of self and other

        >>> Formula({Atom('H'): 2, Atom('O'): 1}) + Formula({Atom('N'): 1, Atom('O'): 2})
        Formula({Atom('H'): 2, Atom('N'): 1, Atom('O'): 3})

        >>> Formula({Atom('H'): Expression.parse('n')}) + Formula({Atom('H'): Expression.parse('-n')})
        Formula({})'''

        if isinstance(other, FormulaElement):
            result = Formula(self._values)
            if other in self._values:
                result._values[other] += 1
                if isinstance(result._values[other], Expression):
                    result._values[other] = result._values[other].simplify()
            else:
                result._values[other] = 1
        elif isinstance(other, Formula):
            values = defaultdict(int)
            for f in (self._values, other._values):
                for key, value in f.items():
                    values[key] += value
                    if isinstance(values[key], Expression):
                        values[key] = values[key].simplify()
            result = Formula(values)
        else:
            return NotImplemented

        result._values = { key: value for key, value in result._values.items() if value != 0 }
        return result

    def __radd__(self, other):
        return self + other

    def __mul__(self, other):
        '''Multiply each count by other

        >>> Formula({Atom('H'): 2, Atom('O'): 1}) * 4
        Formula({Atom('H'): 8, Atom('O'): 4})'''

        if isinstance(other, numbers.Number):
            result = Formula({ key: value*other for key, value in self._values.items() })
        elif isinstance(other, Variable) or isinstance(other, Expression):
            result = Formula({ self: other.simplify() })
        else:
            return NotImplemented

        result._values = { key: value for key, value in result._values.items() if value != 0 }
        return result

    def __rmul__(self, other):
        '''Multiply each count by other

        >>> 2 * Formula({Atom('H'): 2, Atom('O'): 1})
        Formula({Atom('H'): 4, Atom('O'): 2})'''

        return self * other

    def __hash__(self):
        h = hash('Formula')
        for element, value in self._values.iteritems():
            h ^= hash(element) ^ hash(value)
        return h

    def __eq__(self, other):
        '''Equality of self and other

        >>> Formula({Atom('H'): 2, Atom('O'): 1}) == Formula({Atom('O'): 1, Atom('H'): 2})
        True
        >>> Formula({Atom('Au'): 1}) == Formula({Atom('Ag'): 1})
        False'''

        return isinstance(other, Formula) and self._values == other._values

    def __ne__(self, other):
        '''Not equality of self and other

        >>> Formula({Atom('Au'): 1}) != Formula({Atom('Ag'): 1})
        True'''

        return not self == other

    @classmethod
    def parse(cls, s):
        '''Parse a formula string (e.g. C6H10O2)

        >>> Formula.parse('H2O2')
        Formula({Atom('H'): 2, Atom('O'): 2})

        >>> Formula.parse('H2O')
        Formula({Atom('H'): 2, Atom('O'): 1})

        >>> Formula.parse('Zn')
        Formula({Atom('Zn'): 1})

        >>> Formula.parse('C2H5NO2')
        Formula({Atom('C'): 2, Atom('H'): 5, Atom('N'): 1, Atom('O'): 2})
        '''

        def tokenize(s):
            '''Tokenize formula string'''
            token = ''
            for c in s:
                if c == ' ':
                    continue
                elif c in '().':
                    if token != '':
                        yield token
                        token = ''
                    yield c
                elif c.isupper():
                    if token != '':
                        yield token
                    token = c
                else:
                    token += c
            if token != '':
                yield token

        def split_atom_string(a):
            '''Split atom string into name and count'''
            for i, c in enumerate(a):
                if c.isdigit():
                    return a[:i], int(a[i:])
            return a, 1

        def transform_subformula(form):
            '''Extract radical if subformula is a singleton with a radical'''
            if len(form._values) == 1:
                # A radical in a singleton subformula is interpreted as a
                # numbered radical.
                element, value = next(form._values.iteritems())
                if isinstance(element, Radical):
                    return Radical('{}{}'.format(element.symbol, value))
            return form

        stack = []
        formula = Formula()
        expect_group_spec = False
        for token in tokenize(s):
            if expect_group_spec:
                expect_group_spec = False
                if token.islower() or token.isdigit():
                    subformula = formula
                    formula = stack.pop()
                    if token.isdigit():
                        formula += int(token) * transform_subformula(subformula)
                    else:
                        expr = Expression.parse(token)
                        formula += expr * transform_subformula(subformula)
                    continue
                else:
                    subformula = formula
                    formula = stack.pop()
                    formula += transform_subformula(subformula)

            if token == '(':
                stack.append(formula)
                formula = Formula()
            elif token == ')':
                expect_group_spec = True
            else:
                symbol, value = split_atom_string(token)
                if symbol in 'RX':
                    formula += value * Radical(symbol)
                else:
                    formula += value * Atom(symbol)

        if expect_group_spec:
            subformula = formula
            formula = stack.pop()
            formula += transform_subformula(subformula)

        return formula

    @classmethod
    def balance(cls, lhs, rhs):
        '''Return formulas that need to be added to balance given formulas

        Given complete formulas for right side and left side of a reaction,
        calculate formulas for the missing compounds on both sides. Return
        as a left, right tuple.

        >>> Formula.balance(Formula.parse('H2O'), Formula.parse('OH'))
        (Formula({}), Formula({Atom('H'): 1}))

        >>> Formula.balance(Formula.parse('C3H6OH'), Formula.parse('CH6O2'))
        (Formula({Atom('O'): 1}), Formula({Atom('C'): 2, Atom('H'): 1}))

        >>> Formula.balance(Formula.parse('H2(CH2)n'), Formula.parse('CH3O(CH2)n'))
        (Formula({Atom('C'): 1, Atom('H'): 1, Atom('O'): 1}), Formula({}))'''

        def missing(formula, other):
            for element, value in formula._values.iteritems():
                if element not in other._values:
                    yield value*element
                # Check that
                elif value - other._values[element] > 0:
                    delta = value - other._values[element]
                    yield delta*element

        return sum(missing(rhs, lhs), Formula()), sum(missing(lhs, rhs), Formula())

if __name__ == '__main__':
    import doctest
    doctest.testmod()
