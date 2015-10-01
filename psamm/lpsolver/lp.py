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

import numbers
from collections import Counter
import abc

from six import add_metaclass, iteritems
from six.moves import range


class VariableSet(tuple):
    """A tuple used to represent sets of variables"""


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
        """Iterator of variables in expression"""
        return iter(self._variables)

    def values(self):
        """Iterator of (variable, value)-pairs in expression"""
        return iteritems(self._variables)

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
            return self.__class__(self._variables, self._offset + other)
        elif isinstance(other, self.__class__):
            variables = Counter(self._variables)
            variables.update(other._variables)
            return self.__class__(variables, self._offset + other._offset)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        return self.__class__(
            {var: value*other for var, value in iteritems(self._variables)},
            self._offset*other)

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return self * -1

    def __eq__(self, other):
        """Return equality relation (equation): self == other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(Relation.Equals, self - other)

    def __ge__(self, other):
        """Return greater-than relation (inequality): self >= other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(Relation.Greater, self - other)

    def __le__(self, other):
        """Return less-than relation (inequality): self <= other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(Relation.Less, self - other)

    def __gt__(self, other):
        """Return strictly greater-than relation (inequality): self > other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(Relation.StrictlyGreater, self - other)

    def __lt__(self, other):
        """Return strictly less-than relation (inequality): self < other

        This method is overloaded so that relations can be
        formed using a natural syntax.
        """

        return Relation(Relation.StrictlyLess, self - other)

    def __str__(self):
        """Return string representation of expression"""

        def all_terms():
            count_vars = 0
            for name, value in sorted(iteritems(self._variables)):
                if value != 0:
                    count_vars += 1
                    if isinstance(name, VariableSet):
                        yield '<set>', value
                    else:
                        yield name, value
            if self._offset != 0 or count_vars == 0:
                yield None, self._offset

        terms = []
        for i, spec in enumerate(all_terms()):
            name, value = spec
            if i == 0:
                # First term is special
                if name is None:
                    terms.append('{}'.format(value))
                elif abs(value) == 1:
                    terms.append(name if value > 0 else '-{}'.format(name))
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
        return '<{} {}>'.format(self.__class__.__name__, repr(str(self)))


class Relation(object):
    """Represents a binary relation (equation or inequality)

    Relations can be equalities or inequalities. All relations
    of this type can be represented as a left-hand side expression
    and the type of relation. In this representation, the right-hand
    side is always zero.
    """

    Equals = object()
    Greater = object()
    Less = object()
    StrictlyGreater = object()
    StrictlyLess = object()

    SYMBOL = {
        Equals: '==',
        Greater: '>=',
        Less: '<=',
        StrictlyGreater: '>',
        StrictlyLess: '<'
    }

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
        return '{} {} 0'.format(
            str(self._expression), Relation.SYMBOL[self._sense])

    def __repr__(self):
        return '<{} {}>'.format(self.__class__.__name__, repr(str(self)))


class ObjectiveSense(object):
    """Enumeration of objective sense values"""

    Minimize = object()
    """Minimize objective function"""

    Maximize = object()
    """Maximize objective function"""


class VariableType(object):
    """Enumeration of variable types"""

    Continuous = object()
    """Continuous variable type"""

    Integer = object()
    """Integer variable type"""

    Binary = object()
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

        Raises `ValueError` if the relation is invalid. This method also
        accepts relation given as `bool` and will raise an error if the value
        is `False`.
        """
        if isinstance(relation, bool):
            if not relation:
                raise ValueError('Unsatisfiable relation added')
        else:
            if relation.sense in (
                    Relation.StrictlyGreater, Relation.StrictlyLess):
                raise ValueError(
                    'Strict relations are invalid in LP-problems:'
                    ' {}'.format(relation))

    @abc.abstractmethod
    def add_linear_constraints(self, *relations):
        """Add constraints to the problem.

        Each constraint is given as a :class:`.Relation`, and the expression
        in that relation can be a set expression. Returns a sequence of
        :class:`Constraints <.Constraint>`.
        """

    @abc.abstractmethod
    def set_linear_objective(self, expression):
        """Set linear objective of the problem to the given
        :class:`.Expression`.
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

    @abc.abstractmethod
    def get_value(self, expression):
        """Get value of variable or expression in result

        Expression can be an object defined as a name in the problem, in which
        case the corresponding value is simply returned. If expression is an
        actual :class:`.Expression` object, it will be evaluated using the
        values from the result.
        """


if __name__ == '__main__':
    import doctest
    doctest.testmod()
