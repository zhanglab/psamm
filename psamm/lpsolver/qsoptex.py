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

"""Linear programming solver using QSopt_ex"""

from __future__ import absolute_import

from itertools import repeat, count, izip
import numbers

import qsoptex

from .lp import Solver as BaseSolver
from .lp import Problem as BaseProblem
from .lp import Result as BaseResult
from .lp import (VariableSet, Expression, Relation,
                    ObjectiveSense, VariableType, InvalidResultError)


class Solver(BaseSolver):
    """Represents an LP solver using QSopt_ex"""

    def create_problem(self, **kwargs):
        """Create a new LP-problem using the solver"""
        return Problem(**kwargs)


class Problem(BaseProblem):
    """Represents an LP-problem of a qsoptex.Solver"""

    def __init__(self, **kwargs):
        self._p = qsoptex.ExactProblem()
        self._p.set_param(qsoptex.Parameter.SIMPLEX_DISPLAY, 1)

        self._variables = {}
        self._var_names = ('x'+str(i) for i in count(1))

        self._result = None

    @property
    def qsoptex(self):
        """The underlying qsoptex.ExactProblem object"""
        return self._p

    def define(self, *names, **kwargs):
        """Define variable in the problem

        Variables must be defined before they can be accessed by var() or set().
        This function takes keyword arguments lower and upper to define the
        bounds of the variable (default: -inf to inf). The keyword argument types
        can be used to select the type of the variable (Only Continuous is suported).
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
        if vartype is None or vartype in (VariableType.Continuous, VariableType.Binary,
                                            VariableType.Integer):
            vartype = repeat(vartype, len(names))

        lp_names = tuple(next(self._var_names) for name in names)

        # Assign default values
        vartype = (VariableType.Continuous if value is None else value for value in vartype)

        self._variables.update(izip(names, lp_names))
        for name, lower, upper, t in izip(lp_names, lower, upper, vartype):
            if t != VariableType.Continuous:
                raise ValueError('Solver does not support non-continuous types')
            self._p.add_variable(0, lower, upper, name)

    def var(self, name):
        """Return the variable as an expression"""
        if name not in self._variables:
            raise ValueError('Undefined variable: {}'.format(name))
        return Expression({ name: 1 })

    def set(self, names):
        """Return the set of variables as an expression"""
        names = tuple(names)
        if any(name not in self._variables for name in names):
            raise ValueError('Undefined variables: {}'.format(set(names) - set(self._variables)))
        return Expression({ VariableSet(names): 1 })

    def add_linear_constraints(self, *relations):
        """Add constraints to the problem

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.
        """
        for relation in relations:
            if isinstance(relation, bool):
                # A bool in place of a relation is accepted to mean
                # a relation that does not involve any variables and
                # has therefore been evaluated to a truth-value (e.g
                # '0 == 0' or '2 >= 3').
                if not relation:
                    raise ValueError('Unsatisfiable relation added')
            else:
                if relation.sense in (Relation.StrictlyGreater, Relation.StrictlyLess):
                    raise ValueError('Strict relations are invalid in LP-problems: {}'.format(relation))

                expression = relation.expression
                pairs = []
                for value_set in expression.value_sets():
                    values = ((self._variables[variable], value) for variable, value in value_set)
                    self._p.add_linear_constraint(relation.sense, values, -expression.offset)

    def set_linear_objective(self, expression):
        """Set linear objective of problem"""

        if isinstance(expression, numbers.Number):
            # Allow expressions with no variables as objective,
            # represented as a number
            expression = Expression()

        self._p.set_linear_objective((lp_name, expression.value(var)) for var, lp_name in self._variables.iteritems())

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
            return sum(self._problem._p.get_value(self._problem._variables[var])*value for var, value in expression.values())
        elif expression not in self._problem._variables:
            raise ValueError('Unknown expression: {}'.format(expression))
        return self._problem._p.get_value(self._problem._variables[expression])
