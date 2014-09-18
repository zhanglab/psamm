
'''Representations of affine expressions and variables

These classes can be used to represent affince expressions
and do manipulation and evalutaion with substitutions of
particular variables.'''

import re
import numbers
import functools
from collections import Counter

@functools.total_ordering
class Variable(object):
    '''Represents a variable in an expression

    Equality of variables is based on the symbol.'''

    def __init__(self, symbol):
        '''Create variable with given symbol

        >>> Variable('x')
        Variable('x')'''

        self._symbol = str(symbol)

    @property
    def symbol(self):
        '''Symbol of variable

        >>> Variable('x').symbol
        'x'
        '''
        return self._symbol

    def simplify(self):
        '''Return simplified expression

        The simplified form of a variable is always the
        variable itself.

        >>> Variable('x').simplify()
        Variable('x')'''
        return self

    def substitute(self, **kwargs):
        '''Return expression with variables substituted

        >>> Variable('x').substitute(x=567)
        567
        >>> Variable('x').substitute(y=42)
        Variable('x')
        >>> Variable('x').substitute(x=123, y=56)
        123'''

        if self._symbol in kwargs:
            return kwargs[self._symbol]
        return self

    def __add__(self, other):
        return Expression({ self: 1 }) + other

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return Expression({ self: 1 }) - other

    def __rsub__(self, other):
        return self - other

    def __mul__(self, other):
        return Expression({ self: 1 }) * other

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        return Expression({ self: 1 }) / other

    def __neg__(self):
        return Expression({ self: -1 })

    def __eq__(self, other):
        '''Check equality of variables

        >>> Variable('x') == Variable('x')
        True
        >>> Variable('x') == Variable('y')
        False
        >>> Variable('x') == 5
        False'''
        if isinstance(other, Expression):
            return other == self
        return isinstance(other, Variable) and self._symbol == other._symbol

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if isinstance(other, Variable):
            return self._symbol < other._symbol
        return NotImplemented

    def __hash__(self):
        return hash('Variable') ^ hash(self._symbol)

    def __repr__(self):
        return 'Variable({})'.format(repr(self._symbol))

class Expression(object):
    '''Represents an affine expression (e.g. 2x + 3y - z + 5)'''

    def __init__(self, variables={}, offset=0):
        '''Create expression of given variables and offset

        >>> Expression({ Variable('x'): 2 }, 3)
        <Expression '2x + 3'>
        >>> Expression({ Variable('x'): 1, Variable('y'): 1 })
        <Expression 'x + y'>'''
        self._variables = dict(variables)
        self._offset = offset

    def simplify(self):
        '''Return simplified expression

        If the expression is of the form 'x', the variable will be returned,
        and if the expression contains no variables, the offset will be returned
        as a number.

        >>> Expression({}, 5).simplify()
        5
        >>> Expression({ Variable('x'): 1 }).simplify()
        Variable('x')
        >>> Expression({ Variable('x'): 2 }).simplify()
        <Expression '2x'>'''
        result = self.__class__({ var: value for var, value in self._variables.iteritems() if value != 0 }, self._offset)
        if len(result._variables) == 0:
            return result._offset
        elif len(result._variables) == 1:
            var, value = next(result._variables.iteritems())
            if value == 1:
                return var
        return result

    def substitute(self, **kwargs):
        '''Return expression with variables substituted

        >>> Expression({ Variable('x'): 2 }, 1).substitute(x=2)
        5
        >>> Expression({ Variable('x'): 1 }).substitute(y=2)
        Variable('x')
        >>> Expression({ Variable('x'): 2 }, 1).substitute(y=2)
        <Expression '2x + 1'>
        >>> ex = Expression({ Variable('x'): 1, Variable('y'): 2 })
        >>> ex.substitute(y=Variable('x'))
        <Expression '3x'>'''

        expr = self.__class__()
        for var, value in self._variables.iteritems():
            expr += value * var.substitute(**kwargs)
        return (expr + self._offset).simplify()

    def variables(self):
        '''Return iterator of variables in expression

        >>> ex = Expression({ Variable('x'): 1, Variable('y'): 2 }, 3)
        >>> sorted(ex.variables())
        [Variable('x'), Variable('y')]'''
        return iter(self._variables)

    def __add__(self, other):
        '''Add expressions, variables or numbers

        >>> Expression({ Variable('x'): 1 }) + 1
        <Expression 'x + 1'>
        >>> Expression({ Variable('x'): 1 }) + Variable('y')
        <Expression 'x + y'>
        >>> Expression({ Variable('x'): 1 }) + Variable('x')
        <Expression '2x'>
        >>> Expression({ Variable('x'): 1 }) + Expression({ Variable('y'): 2 })
        <Expression 'x + 2y'>'''

        if isinstance(other, numbers.Number):
            return self.__class__(self._variables, self._offset + other)
        elif isinstance(other, Variable):
            return self + Expression({ other: 1 })
        elif isinstance(other, Expression):
            variables = Counter(self._variables)
            variables.update(other._variables)
            variables = { var: value for var, value in variables.iteritems() if value != 0 }
            return self.__class__(variables, self._offset + other._offset)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        '''Subtract expressions, variables or numbers

        >>> Expression({ Variable('x'): 1 }) - 1
        <Expression 'x - 1'>
        >>> Expression({ Variable('x'): 1 }) - Variable('y')
        <Expression 'x - y'>
        >>> Expression({ Variable('x'): 1 }) - Variable('x')
        <Expression '0'>
        >>> Expression({ Variable('x'): 1 }) - Expression({ Variable('y'): 2 })
        <Expression 'x - 2y'>'''

        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        '''Multiply by scalar

        >>> Expression({ Variable('x'): 1 }, 3) * 4
        <Expression '4x + 12'>'''
        if isinstance(other, numbers.Number):
            if other == 0:
                return self.__class__()
            return self.__class__({ var: value*other for var, value in self._variables.iteritems()}, self._offset*other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return self * -1

    def __eq__(self, other):
        '''Expression equality

        >>> Expression({ Variable('x'): 2, Variable('y'): 5 }) == Expression({ Variable('y'): 5, Variable('x'): 2 })
        True
        >>> Expression({ Variable('x'): 2, Variable('y'): 5 }) == Expression({ Variable('y'): 2, Variable('x'): -5 })
        False
        >>> Expression({ Variable('x'): 2, Variable('y'): 5 }) == Expression({ Variable('y'): 2, Variable('x'): -5, Variable('z'): 1 })
        False
        >>> Expression({ Variable('x'): 2, Variable('y'): 5 }, 4) == Expression({ Variable('y'): 5, Variable('x'): 2 })
        False
        >>> Expression({ Variable('x'): 1 }) == Variable('x')
        True
        >>> Expression({ Variable('x'): 2 }) == Variable('x')
        False
        >>> Expression({ Variable('x'): 1 }, 4) == Variable('x')
        False
        >>> Expression({ Variable('x'): 1, Variable('y'): 1 }) == Variable('x')
        False
        >>> Expression({}, 3) == 3
        True
        >>> Expression({}, -3) == 3
        False
        >>> Expression({ Variable('x'): 1 }, 3) == 3
        False
        '''
        if isinstance(other, Expression):
            return self._variables == other._variables and self._offset == other._offset
        elif isinstance(other, Variable):
            # Check that there is just one variable in the expression
            # with a coefficient of one.
            return (self._offset == 0 and len(self._variables) == 1 and
                    next(self._variables.iterkeys()) == other and
                    next(self._variables.itervalues()) == 1)
        elif isinstance(other, numbers.Number):
            return len(self._variables) == 0 and self._offset == other
        return False

    def __ne__(self, other):
        return not self == other

    def __str__(self):
        def all_terms():
            count_vars = 0
            for symbol, value in sorted((var.symbol, value) for var, value in self._variables.iteritems()):
                if value != 0:
                    count_vars += 1
                    yield symbol, value
            if self._offset != 0 or count_vars == 0:
                yield None, self._offset

        terms = []
        for i, spec in enumerate(all_terms()):
            symbol, value = spec
            if i == 0:
                # First term is special
                if symbol is None:
                    terms.append('{}'.format(value))
                elif abs(value) == 1:
                    terms.append(symbol if value > 0 else '-'+symbol)
                else:
                    terms.append('{}{}'.format(value, symbol))
            else:
                prefix = '+' if value >= 0 else '-'
                if symbol is None:
                    terms.append('{} {}'.format(prefix, abs(value)))
                elif abs(value) == 1:
                    terms.append('{} {}'.format(prefix, symbol))
                else:
                    terms.append('{} {}{}'.format(prefix, abs(value), symbol))
        return ' '.join(terms)

    def __repr__(self):
        return '<Expression \'{}\'>'.format(str(self))

    @classmethod
    def parse(cls, s):
        '''Parse expression string

        >>> Expression.parse('2x + 3')
        <Expression '2x + 3'>
        >>> Expression.parse('1')
        <Expression '1'>
        >>> Expression.parse('x + 4y + 2x + 0')
        <Expression '3x + 4y'>
        >>> Expression.parse('-2x1 + 5pi - 3x2')
        <Expression '5pi - 2x1 - 3x2'>'''

        expr = Expression()
        s = s.strip()
        if len(s) == 0:
            return expr

        if s[0] not in '+-':
            s = '+'+s
        for m in re.finditer(r'([+-])\s*(\d+)?([a-z]+\d*)?', s):
            term = 1 if m.group(1) == '+' else -1
            if m.group(2) is not None:
                term *= int(m.group(2))
            if m.group(3) is not None:
                term *= Variable(m.group(3))
            expr += term

        return expr


if __name__ == '__main__':
    import doctest
    doctest.testmod()
