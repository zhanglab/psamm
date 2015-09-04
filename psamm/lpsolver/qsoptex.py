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

"""Linear programming solver using QSopt_ex."""

from __future__ import absolute_import

from itertools import repeat, count
import numbers

from six import iteritems
from six.moves import zip
import qsoptex

from .lp import Solver as BaseSolver
from .lp import Constraint as BaseConstraint
from .lp import Problem as BaseProblem
from .lp import Result as BaseResult
from .lp import (Expression, Relation, ObjectiveSense, VariableType,
                 InvalidResultError)


class Solver(BaseSolver):
    """Represents an LP solver using QSopt_ex"""

    def create_problem(self, **kwargs):
        """Create a new LP-problem using the solver"""
        return Problem(**kwargs)


class Problem(BaseProblem):
    """Represents an LP-problem of a qsoptex.Solver"""

    CONSTR_SENSE_MAP = {
        Relation.Equals: qsoptex.ConstraintSense.EQUAL,
        Relation.Greater: qsoptex.ConstraintSense.GREATER,
        Relation.Less: qsoptex.ConstraintSense.LESS
    }

    def __init__(self, **kwargs):
        self._p = qsoptex.ExactProblem()
        self._p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)

        self._variables = {}
        self._var_names = ('x'+str(i) for i in count(1))
        self._constr_names = ('c'+str(i) for i in count(1))

        self._result = None

    @property
    def qsoptex(self):
        """The underlying qsoptex.ExactProblem object"""
        return self._p

    def define(self, *names, **kwargs):
        """Define variable in the problem

        Variables must be defined before they can be accessed by var() or
        set(). This function takes keyword arguments lower and upper to define
        the bounds of the variable (default: -inf to inf). The keyword argument
        types can be used to select the type of the variable (Only Continuous
        is suported).
        """
        names = tuple(names)
        lower = kwargs.get('lower', None)
        upper = kwargs.get('upper', None)
        vartype = kwargs.get('types', None)

        # Repeat values if a scalar is given
        if lower is None or isinstance(lower, numbers.Number):
            lower = repeat(lower, len(names))
        if upper is None or isinstance(upper, numbers.Number):
            upper = repeat(upper, len(names))
        if vartype is None or vartype in (
                VariableType.Continuous, VariableType.Binary,
                VariableType.Integer):
            vartype = repeat(vartype, len(names))

        lp_names = tuple(next(self._var_names) for name in names)

        # Assign default values
        vartype = (VariableType.Continuous if value is None else value
                   for value in vartype)

        self._variables.update(zip(names, lp_names))
        for name, lower, upper, t in zip(lp_names, lower, upper, vartype):
            if t != VariableType.Continuous:
                raise ValueError(
                    'Solver does not support non-continuous types')
            self._p.add_variable(0, lower, upper, name)

    def has_variable(self, name):
        """Check whether variable is defined in the model."""
        return name in self._variables

    def _add_constraints(self, relation):
        """Add the given relation as one or more constraints

        Return a list of the names of the constraints added.
        """
        expression = relation.expression
        names = []
        for value_set in expression.value_sets():
            values = ((self._variables[variable], value)
                      for variable, value in value_set)
            constr_name = next(self._constr_names)
            sense = self.CONSTR_SENSE_MAP[relation.sense]
            self._p.add_linear_constraint(
                sense=sense, values=values, rhs=-expression.offset,
                name=constr_name)
            names.append(constr_name)

        return names

    def add_linear_constraints(self, *relations):
        """Add constraints to the problem

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.
        """
        constraints = []

        for relation in relations:
            self._check_relation(relation)
            if isinstance(relation, bool):
                # A bool in place of a relation is accepted to mean
                # a relation that does not involve any variables and
                # has therefore been evaluated to a truth-value (e.g
                # '0 == 0' or '2 >= 3').
                constraints.append(Constraint(self, None))
            else:
                for name in self._add_constraints(relation):
                    constraints.append(Constraint(self, name))

        return constraints

    def set_linear_objective(self, expression):
        """Set linear objective of problem"""

        if isinstance(expression, numbers.Number):
            # Allow expressions with no variables as objective,
            # represented as a number
            expression = Expression()

        self._p.set_linear_objective(
            (lp_name, expression.value(var))
            for var, lp_name in iteritems(self._variables))

    def set_objective_sense(self, sense):
        """Set type of problem (maximize or minimize)"""
        if sense == ObjectiveSense.Minimize:
            self._p.set_objective_sense(qsoptex.ObjectiveSense.MINIMIZE)
        elif sense == ObjectiveSense.Maximize:
            self._p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
        else:
            raise ValueError('Invalid objective sense')

    def solve(self, sense=None):
        """Solve problem"""
        if sense is not None:
            self.set_objective_sense(sense)
        self._p.solve()

        self._result = Result(self)
        return self._result

    @property
    def result(self):
        return self._result


class Constraint(BaseConstraint):
    """Represents a constraint in a qsoptex.Problem"""

    def __init__(self, prob, name):
        self._prob = prob
        self._name = name

    def delete(self):
        if self._name is not None:
            self._prob._p.delete_linear_constraint(self._name)


class Result(BaseResult):
    """Represents the solution to a qsoptex.Problem

    This object will be returned from the Problem.solve() method or by
    accessing the Problem.result property after solving a problem. This
    class should not be instantiated manually.

    Result will evaluate to a boolean according to the success of the
    solution, so checking the truth value of the result will immediately
    indicate whether solving was successful.
    """

    def __init__(self, prob):
        self._problem = prob

    def _check_valid(self):
        if self._problem.result != self:
            raise InvalidResultError()

    @property
    def success(self):
        """Return boolean indicating whether a solution was found"""
        self._check_valid()
        return self._problem._p.get_status() == 1

    @property
    def status(self):
        """Return string indicating the error encountered on failure"""
        self._check_valid()
        return 'Status: {}'.format(self._problem._p.get_status())

    def get_value(self, expression):
        """Return value of expression"""

        self._check_valid()
        if isinstance(expression, Expression):
            return sum(self._problem._p.get_value(
                self._problem._variables[var])*value
                for var, value in expression.values())
        elif expression not in self._problem._variables:
            raise ValueError('Unknown expression: {}'.format(expression))
        return self._problem._p.get_value(self._problem._variables[expression])
