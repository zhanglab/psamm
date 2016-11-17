# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Base objects for representation of LP problems.

A linear programming problem is built from a number of constraints and an
objective function. The objective function is a linear expression represented
by :class:`.Expression`. The constraints are represented by :class:`.Relation`,
created from a linear expression and a relation sense (equals, greater, less).

Expressions are built from variables defined in the :class:`.Problem` instance.
In addition, an expression can contain a :class:`.VariableSet` instead of a
single variable. This allows many similar expressions to be represented by one
:class:`.Expression` instance which means that the LP problem can be
constructed faster.
"""

from __future__ import unicode_literals

import math
import numbers
import operator
from collections import Counter, defaultdict
import abc
import enum

import six
from six import add_metaclass, iteritems, viewkeys, viewitems, text_type
from six.moves import range, reduce

_INF = float('inf')


class VariableSet(tuple):
    """A tuple used to represent sets of variables."""


class Product(tuple):
    """A tuple used to represent a variable product."""


class _RangedAccessor(object):
    """Accessor to value bounded by minimum/maximum."""
    def __init__(self, obj, prop):
        self._obj = obj
        self._prop = prop

    @property
    def value(self):
        """Value of property."""
        if self._prop.fget is None:
            raise AttributeError('Unable to read attribute')
        return self._prop.fget(self._obj)

    @value.setter
    def value(self, v):
        if self._prop.fset is None:
            raise AttributeError('Unable to write attribute')
        self._prop.fset(self._obj, v)

    @value.deleter
    def value(self):
        if self._prop.fdel is None:
            raise AttributeError('Unable to delete attribute')
        self._prop.fdel(self._obj)

    @property
    def min(self):
        """Minimum value."""
        if self._prop.fmin is None:
            return -_INF
        return self._prop.fmin(self._obj)

    @property
    def max(self):
        """Maximum value."""
        if self._prop.fmax is None:
            return _INF
        return self._prop.fmax(self._obj)


class RangedProperty(object):
    """Numeric property with minimum and maximum values.

    The value attribute is used to get/set the actual value of the propery.
    The min/max attributes are used to get the bounds. The range is not
    automatically enforced when the value is set.
    """
    def __init__(self, fget=None, fset=None, fdel=None,
                 fmin=None, fmax=None, doc=None):
        self.fget = fget
        self.fset = fset
        self.fdel = fdel
        self.fmin = fmin
        self.fmax = fmax

        if doc is None and fget is not None:
            doc = fget.__doc__
        self.__doc__ = doc

    def __get__(self, obj, objtype=None):
        if obj is None:
            return self
        return _RangedAccessor(obj, self)

    def __set__(self, obj, value):
        raise AttributeError('Unable to set attribute (try value property)')

    def __delete__(self, obj):
        raise AttributeError('Unable to delete attribute (try value property)')

    def getter(self, fget):
        return type(self)(
            fget, self.fset, self.fdel, self.fmin, self.fmax, self.__doc__)

    def setter(self, fset):
        return type(self)(
            self.fget, fset, self.fdel, self.fmin, self.fmax, self.__doc__)

    def deleter(self, fdel):
        return type(self)(
            self.fget, self.fset, fdel, self.fmin, self.fmax, self.__doc__)


def ranged_property(min=None, max=None):
    """Decorator for creating ranged property with fixed bounds."""
    min_value = -_INF if min is None else min
    max_value = _INF if max is None else max
    return lambda fget: RangedProperty(
        fget, fmin=lambda obj: min_value, fmax=lambda obj: max_value)


@six.python_2_unicode_compatible
class Expression(object):
    """Represents a linear expression

    The variables can be any hashable objects. If one or more variables are
    instead :class:`VariableSets <.VariableSet>`, this will be taken to
    represent a set of expressions separately using a different element of the
    :class:`.VariableSet`.

    >>> e = Expression({'x': 2, 'y': 3})
    >>> str(e)
    '2*x + 3*y'

    In order to provide a more natural syntax for creating
    :class:`Relations <.Relation>` the binary relation operators have been
    overloaded to return :class:`.Relation` instances.

    >>> rel = Expression({'x': 2}) >= Expression({'y': 3})
    >>> str(rel)
    '2*x - 3*y >= 0'
    """

    def __init__(self, variables={}, offset=0):
        self._variables = Counter(variables)
        self._offset = offset

    @property
    def offset(self):
        """Value of the offset"""
        return self._offset

    def variables(self):
        """Return immutable view of variables in expression."""
        return viewkeys(self._variables)

    def values(self):
        """Return immutable view of (variable, value)-pairs in expression."""
        return viewitems(self._variables)

    def value(self, variable):
        return self._variables.get(variable, 0)

    def value_sets(self):
        """Iterator of expression sets

        This will yield an iterator of (variable, value)-pairs for each
        expression in the expression set (each equivalent to values()). If none
        of the variables is a set variable then a single iterator will be
        yielded.
        """
        count = max(1 if not isinstance(var, VariableSet) else
                    len(var) for var in self._variables)

        def value_set(n):
            for variable, value in iteritems(self._variables):
                if isinstance(variable, VariableSet):
                    yield variable[n], value
                else:
                    yield variable, value
        for i in range(count):
            yield value_set(i)

    def __contains__(self, variable):
        return variable in self._variables

    def __add__(self, other):
        """Add expression with a number or another expression"""

        if isinstance(other, numbers.Number):
            if math.isinf(self._offset) or math.isinf(other):
                return self.__class__(offset=self._offset + other)
            return self.__class__(self._variables, self._offset + other)
        elif isinstance(other, self.__class__):
            if math.isinf(self._offset) or math.isinf(other._offset):
                return self.__class__(offset=self._offset + other._offset)
            variables = Counter(self._variables)
            variables.update(other._variables)
            return self.__class__(variables, self._offset + other._offset)
        return NotImplemented

    __radd__ = __add__

    def __iadd__(self, other):
        if isinstance(other, numbers.Number):
            if math.isinf(self._offset) or math.isinf(other):
                self._variables = {}
            self._offset += other
        elif isinstance(other, self.__class__):
            self._offset += other._offset
            if math.isinf(self._offset) or math.isinf(other._offset):
                self._variables = {}
            else:
                self._variables.update(other._variables)
        else:
            return NotImplemented

        return self

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
        self += -other
        return self

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            if math.isinf(other):
                return self.__class__(offset=float('nan'))

            return self.__class__(
                {var: value * other for var, value in
                 iteritems(self._variables)}, self._offset * other)
        elif isinstance(other, self.__class__):
            variables = defaultdict(int)
            for v1, value1 in iteritems(self._variables):
                for v2, value2 in iteritems(other._variables):
                    p1 = v1 if isinstance(v1, Product) else Product((v1,))
                    p2 = v2 if isinstance(v2, Product) else Product((v2,))
                    product = Product(sorted(p1 + p2))
                    variables[product] += value1 * value2

                if other._offset != 0:
                    variables[v1] += value1 * other._offset

            if self._offset != 0:
                for v2, value2 in iteritems(other._variables):
                    variables[v2] += value2 * self._offset

            offset = self._offset * other._offset
            return self.__class__(variables, offset=offset)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __imul__(self, other):
        if isinstance(other, numbers.Number):
            if math.isinf(other):
                self._variables = {}
                self._offset = float('nan')
            else:
                for var in self._variables:
                    self._variables[var] *= other
            self._offset *= other
        elif isinstance(other, self.__class__):
            variables = defaultdict(int)
            for v1, value1 in iteritems(self._variables):
                for v2, value2 in iteritems(other._variables):
                    p1 = v1 if isinstance(v1, Product) else Product((v1,))
                    p2 = v2 if isinstance(v2, Product) else Product((v2,))
                    product = Product(sorted(p1 + p2))
                    variables[product] += value1 * value2

                if other._offset != 0:
                    variables[v1] += value1 * other._offset

            if self._offset != 0:
                for v2, value2 in iteritems(other._variables):
                    variables[v2] += value2 * self._offset

            self._offset *= other._offset
            self._variables = dict(variables)
        else:
            return NotImplemented

        return self

    def __neg__(self):
        return self.__class__(
            {var: -value for var, value in iteritems(self._variables)},
            -self._offset)

    def __pow__(self, other):
        if isinstance(other, int):
            if other < 0:
                raise ValueError('Exponent must be positive')
            elif other == 0:
                return self.__class__(offset=1)

            r = self.__class__(self._variables, self._offset)
            for i in range(other - 1):
                r.__imul__(self)

            return r
        else:
            return NotImplemented

    def __rpow__(self, other):
        return pow(other, self)

    def __ipow__(self, other):
        if isinstance(other, int):
            if other < 0:
                raise ValueError('Exponent must be positive')
            elif other == 0:
                self._variables = {}
                self._offset = 1
            else:
                tmp = self.__class__(self._variables, self._offset)
                for i in range(other - 1):
                    self.__imul__(tmp)
        else:
            return NotImplemented

        return self

    def __eq__(self, other):
        """Return equality relation (equation): self == other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(RelationSense.Equals, self - other)

    def __ge__(self, other):
        """Return greater-than relation (inequality): self >= other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(RelationSense.Greater, self - other)

    def __le__(self, other):
        """Return less-than relation (inequality): self <= other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(RelationSense.Less, self - other)

    def __gt__(self, other):
        """Return strictly greater-than relation (inequality): self > other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(RelationSense.StrictlyGreater, self - other)

    def __lt__(self, other):
        """Return strictly less-than relation (inequality): self < other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(RelationSense.StrictlyLess, self - other)

    def __str__(self):
        """Return string representation of expression"""

        def all_terms():
            count_vars = 0
            for name, value in sorted(
                    iteritems(self._variables),
                    key=lambda p: (isinstance(p[0], Product), p[0])):
                if value != 0:
                    count_vars += 1
                    if isinstance(name, VariableSet):
                        yield '<set>', value
                    elif isinstance(name, Product):
                        yield '*'.join(text_type(n) for n in name), value
                    else:
                        yield name, value
            if self._offset != 0 or count_vars == 0:
                yield None, self._offset

        terms = []
        for i, (name, value) in enumerate(all_terms()):
            if i == 0:
                # First term is special
                if name is None:
                    terms.append('{}'.format(value))
                elif abs(value) == 1:
                    terms.append(
                        text_type(name) if value > 0 else '-{}'.format(name))
                else:
                    terms.append('{}*{}'.format(value, name))
            else:
                prefix = '+' if value >= 0 else '-'
                if name is None:
                    terms.append('{} {}'.format(prefix, abs(value)))
                elif abs(value) == 1:
                    terms.append('{} {}'.format(prefix, name))
                else:
                    terms.append('{} {}*{}'.format(prefix, abs(value), name))
        return ' '.join(terms)

    def __repr__(self):
        return str('<{} {}>').format(
            self.__class__.__name__, repr(str(self)))


@enum.unique
class RelationSense(enum.Enum):
    Equals = '=='
    Greater = '>='
    Less = '<='
    StrictlyGreater = '>'
    StrictlyLess = '<'


@six.python_2_unicode_compatible
class Relation(object):
    """Represents a binary relation (equation or inequality)

    Relations can be equalities or inequalities. All relations
    of this type can be represented as a left-hand side expression
    and the type of relation. In this representation, the right-hand
    side is always zero.
    """

    def __init__(self, sense, expression):
        self._sense = sense
        self._expression = expression

    @property
    def sense(self):
        """Type of relation (equality or inequality)

        Can be one of Equal, Greater or Less, or one of the
        strict relations, StrictlyGreater or StrictlyLess.
        """
        return self._sense

    @property
    def expression(self):
        """Left-hand side expression"""
        return self._expression

    def __str__(self):
        """Convert relation to string representation"""
        var_expr = Expression(dict(self._expression.values()))
        return '{} {} {}'.format(
            str(var_expr), self._sense.value, -self._expression.offset)

    def __repr__(self):
        return str('<{} {}>').format(self.__class__.__name__, repr(str(self)))


@enum.unique
class ObjectiveSense(enum.Enum):
    """Enumeration of objective sense values"""

    Minimize = -1
    """Minimize objective function"""

    Maximize = 1
    """Maximize objective function"""


@enum.unique
class VariableType(enum.Enum):
    """Enumeration of variable types"""

    Continuous = 'C'
    """Continuous variable type"""

    Integer = 'I'
    """Integer variable type"""

    Binary = 'B'
    """Binary variable type (0 or 1)"""


@add_metaclass(abc.ABCMeta)
class Solver(object):
    """Factory for LP Problem instances"""

    @abc.abstractmethod
    def create_problem(self):
        """Create a new :class:`.Problem` instance"""


@add_metaclass(abc.ABCMeta)
class Constraint(object):
    """Represents a constraint within an LP Problem"""

    @abc.abstractmethod
    def delete(self):
        """Remove constraint from Problem instance"""


class VariableNamespace(object):
    def __init__(self, problem, **kwargs):
        self._problem = problem
        self._define_kwargs = kwargs

    def define(self, names, **kwargs):
        define_kwargs = dict(self._define_kwargs)
        define_kwargs.update(kwargs)
        self._problem.define(
            *((self, name) for name in names), **define_kwargs)

    def __call__(self, name):
        return self._problem.var((self, name))

    def set(self, names):
        return self._problem.set((self, name) for name in names)

    def sum(self, names):
        return Expression({(self, name): 1 for name in names})

    def expr(self, items):
        return Expression({(self, name): value for name, value in items})

    def value(self, name):
        return self._problem.result.get_value((self, name))

    def __repr__(self):
        return '<{} of {} ({})>'.format(
            self.__class__.__name__, repr(self._problem), id(self))


@add_metaclass(abc.ABCMeta)
class Problem(object):
    """Representation of LP Problem instance

    Variable names in the problem can be any hashable object. It is the
    responsibility of the solver interface to translate the object into a
    unique string if required by the underlying LP solver.
    """

    @abc.abstractmethod
    def define(self, *names, **kwargs):
        """Define a variable in the problem"""

    @abc.abstractmethod
    def has_variable(self, name):
        """Check whether a variable is defined in the problem."""

    def namespace(self, names=None, **kwargs):
        ns = VariableNamespace(self, **kwargs)
        if names is not None:
            ns.define(names)
        return ns

    def var(self, name):
        """Return variable as an :class:`.Expression`."""
        if not self.has_variable(name):
            raise ValueError('Undefined variable: {}'.format(name))
        return Expression({name: 1})

    def expr(self, values, offset=0):
        """Return the given dictionary of values as an :class:`.Expression`."""
        if isinstance(values, dict):
            for name in values:
                if not self.has_variable(name):
                    raise ValueError('Undefined variable: {}'.format(name))
            return Expression(values, offset=offset)

        if not self.has_variable(values):
            raise ValueError('Undefined variable: {}'.format(values))
        return Expression({values: 1}, offset=offset)

    def set(self, names):
        """Return variables as a set expression.

        This returns an :class:`.Expression` containing a
        :class:`.VariableSet`.
        """
        names = tuple(names)
        if any(not self.has_variable(name) for name in names):
            raise ValueError('Undefined variables: {}'.format(
                set(names) - set(self._variables)))
        return Expression({VariableSet(names): 1})

    def _check_relation(self, relation):
        """Check whether the given relation is valid.

        Return false normally, or true if the relation is determined to be
        tautologically true. Raises ``ValueError`` if the relation is
        tautologically false. This includes relations given as a ``False`` bool
        value or certain types of relations that have an infinity offset.
        """
        if isinstance(relation, bool):
            if not relation:
                raise ValueError('Unsatisfiable relation added')
            return True

        if relation.sense in (
                RelationSense.StrictlyGreater, RelationSense.StrictlyLess):
            raise ValueError(
                'Strict relations are invalid in LP-problems:'
                ' {}'.format(relation))

        if relation.expression.offset == -_INF:
            if relation.sense == RelationSense.Less:
                return True
            else:
                raise ValueError('Unsatisfiable relation added: {}'.format(
                    relation))
        elif relation.expression.offset == _INF:
            if relation.sense == RelationSense.Greater:
                return True
            else:
                raise ValueError('Unsatisfiable relation added: {}'.format(
                    relation))

        return False

    @abc.abstractmethod
    def add_linear_constraints(self, *relations):
        """Add constraints to the problem.

        Each constraint is given as a :class:`.Relation`, and the expression
        in that relation can be a set expression. Returns a sequence of
        :class:`Constraints <.Constraint>`.
        """

    @abc.abstractmethod
    def set_objective(self, expression):
        """Set objective of the problem to the given :class:`.Expression`."""

    set_linear_objective = set_objective
    """Set objective of the problem.

    .. deprecated:: 0.19
       Use :meth:`set_objective` instead.
    """

    @abc.abstractmethod
    def set_objective_sense(self, sense):
        """Set type of problem (minimize or maximize)"""

    @abc.abstractmethod
    def solve(self, sense=None):
        """Solve problem and return result"""

    @abc.abstractproperty
    def result(self):
        """Result of solved problem"""

    @ranged_property(min=None, max=None)
    def feasibility_tolerance(self):
        raise NotImplementedError('Feasiblity tolerance not available')

    @ranged_property(min=None, max=None)
    def optimality_tolerance(self):
        raise NotImplementedError('Optimality tolerance not available')

    @ranged_property(min=None, max=None)
    def integrality_tolerance(self):
        raise NotImplementedError('Integrality tolerance not available')


class InvalidResultError(Exception):
    """Raised when a result that has been invalidated is accessed"""

    def __init__(self, msg=None):
        if msg is None:
            msg = ('Previous result is no longer valid as problem'
                   ' has been solved again')
        super(InvalidResultError, self).__init__(msg)


class Result(object):
    """Result of solving an LP problem

    The result is tied to the solver instance and is valid at least until the
    problem is solved again. If the problem has been solved again an
    :class:`InvalidResultError` may be raised.
    """

    @abc.abstractproperty
    def success(self):
        """Whether solution was optimal"""

    def __nonzero__(self):
        """Whether solution was optimal"""
        return self.success

    __bool__ = __nonzero__

    @abc.abstractproperty
    def status(self):
        """String indicating the status of the problem result"""

    @abc.abstractproperty
    def unbounded(self):
        """Whether solution is unbounded"""

    @abc.abstractmethod
    def _get_value(self, var):
        """Return the solution value of a single variable."""

    @abc.abstractmethod
    def _has_variable(self, var):
        """Return whether variable exists in the solution."""

    def _evaluate_expression(self, expr):
        """Evaluate an :class:`.Expression` using :meth:`_get_value`."""
        total = expr.offset
        for var, value in expr.values():
            if not isinstance(var, Product):
                total += self._get_value(var) * value
            else:
                total += reduce(
                    operator.mul, (self._get_value(v) for v in var), value)
        return total

    def get_value(self, expression):
        """Get value of variable or expression in result

        Expression can be an object defined as a name in the problem, in which
        case the corresponding value is simply returned. If expression is an
        actual :class:`.Expression` object, it will be evaluated using the
        values from the result.
        """
        if isinstance(expression, Expression):
            return self._evaluate_expression(expression)
        elif not self._has_variable(expression):
            raise ValueError('Unknown expression: {}'.format(expression))

        return self._get_value(expression)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
