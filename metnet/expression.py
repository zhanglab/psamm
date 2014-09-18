
'''Representations of affine expressions and variables

These classes can be used to represent affince expressions
and do manipulation and evalutaion with substitutions of
particular variables.'''

import re
import numbers
import functools
from collections import Counter, defaultdict

@functools.total_ordering
class Variable(object):
    '''Represents a variable in an expression

    Equality of variables is based on the symbol.'''

    def __init__(self, symbol):
        '''Create variable with given symbol

        Symbol must start with a letter or underscore but
        can contain numbers in other positions.

        >>> Variable('x')
        Variable('x')'''

        if not re.match(r'^[^\d\W]\w*\Z', symbol):
            raise ValueError('Invalid symbol {}'.format(symbol))
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
        123
        '''

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
        return -self + other

    def __mul__(self, other):
        return Expression({ self: 1 }) * other

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        return Expression({ self: 1 }) / other

    def __neg__(self):
        return Expression({ self: -1 })

    def __eq__(self, other):
        '''Check equality of variables'''
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

        self._variables = {}
        self._offset = offset

        for var, value in variables.iteritems():
            if not isinstance(var, Variable):
                raise ValueError('Not a variable: {}'.format(var))
            if value != 0:
                self._variables[var] = value

    @classmethod
    def parse(cls, s):
        return Expression(*cls._parse_string(s))

    @classmethod
    def _parse_string(cls, s):
        '''Parse expression string

        Variables must be valid variable symbols and
        coefficients must be integers.'''

        scanner = re.compile(r'''
            (\s+) |         # whitespace
            ([^\d\W]\w*) |  # variable
            (\d+) |         # number
            ([+-]) |        # sign
            (\Z) |          # end
            (.)             # error
        ''', re.DOTALL | re.VERBOSE)

        variables = defaultdict(int)
        offset = 0

        # Parse using four states:
        # 0: expect sign, variable, number or end (start state)
        # 1: expect sign or end
        # 2: expect variable or number
        # 3: expect sign, variable or end
        # All whitespace is ignored
        state = 0
        state_number = 1
        for match in re.finditer(scanner, s):
            whitespace, variable, number, sign, end, error = match.groups()
            if error is not None:
                raise ValueError('Invalid token in expression string: {}'.format(match.group(0)))
            elif whitespace is not None:
                continue
            elif variable is not None and state in (0, 2, 3):
                variables[Variable(variable)] += state_number
                state = 1
            elif sign is not None and state in (0, 1, 3):
                if state == 3:
                    offset += state_number
                state_number = 1 if sign == '+' else -1
                state = 2
            elif number is not None and state in (0, 2):
                state_number = state_number * int(number)
                state = 3
            elif end is not None and state in (0, 1, 3):
                if state == 3:
                    offset += state_number
            else:
                raise ValueError('Invalid token in expression string: {}'.format(match.group(0)))

        # Remove zero-coefficient elements
        variables = { var: value for var, value in variables.iteritems() if value != 0 }
        return variables, offset

    def simplify(self):
        '''Return simplified expression

        If the expression is of the form 'x', the variable will be returned,
        and if the expression contains no variables, the offset will be returned
        as a number.'''
        result = self.__class__(self._variables, self._offset)
        if len(result._variables) == 0:
            return result._offset
        elif len(result._variables) == 1 and result._offset == 0:
            var, value = next(result._variables.iteritems())
            if value == 1:
                return var
        return result

    def substitute(self, **kwargs):
        '''Return expression with variables substituted

        >>> Expression.parse('x + 2y').substitute(y=-3)
        <Expression 'x - 6'>
        >>> Expression.parse('x + 2y').substitute(y=Variable('z'))
        <Expression 'x + 2z'>
        '''
        expr = self.__class__()
        for var, value in self._variables.iteritems():
            expr += value * var.substitute(**kwargs)
        return (expr + self._offset).simplify()

    def variables(self):
        '''Return iterator of variables in expression'''
        return iter(self._variables)

    def __add__(self, other):
        '''Add expressions, variables or numbers'''

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
        '''Subtract expressions, variables or numbers'''
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        '''Multiply by scalar'''
        if isinstance(other, numbers.Number):
            if other == 0:
                return self.__class__()
            return self.__class__({ var: value*other for var, value in self._variables.iteritems()}, self._offset*other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __div__(self, other):
        '''Divide by scalar'''
        if isinstance(other, numbers.Number):
            return self.__class__({ var: value/other for var, value in self._variables.iteritems()}, self._offset/other)
        return NotImplemented

    def __neg__(self):
        return self * -1

    def __eq__(self, other):
        '''Expression equality'''
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


if __name__ == '__main__':
    import doctest
    doctest.testmod()
