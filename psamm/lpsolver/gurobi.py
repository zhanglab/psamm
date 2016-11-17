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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Linear programming solver using Gurobi."""

from __future__ import absolute_import

import logging
from itertools import repeat, count
import numbers

from six import raise_from
from six.moves import zip

import gurobipy

from .lp import Solver as BaseSolver
from .lp import Constraint as BaseConstraint
from .lp import Problem as BaseProblem
from .lp import Result as BaseResult
from .lp import (Expression, Product, RelationSense, ObjectiveSense,
                 VariableType, InvalidResultError, ranged_property)

# Module-level logging
logger = logging.getLogger(__name__)


class Solver(BaseSolver):
    """Represents an LP-solver using Gurobi."""

    def create_problem(self, **kwargs):
        """Create a new LP-problem using the solver."""
        return Problem(**kwargs)


class Problem(BaseProblem):
    """Represents an LP-problem of a gurobi.Solver."""

    VARTYPE_MAP = {
        VariableType.Continuous: gurobipy.GRB.CONTINUOUS,
        VariableType.Binary: gurobipy.GRB.BINARY,
        VariableType.Integer: gurobipy.GRB.INTEGER
    }

    CONSTR_SENSE_MAP = {
        RelationSense.Equals: gurobipy.GRB.EQUAL,
        RelationSense.Greater: gurobipy.GRB.GREATER_EQUAL,
        RelationSense.Less: gurobipy.GRB.LESS_EQUAL
    }

    OBJ_SENSE_MAP = {
        ObjectiveSense.Maximize: gurobipy.GRB.MAXIMIZE,
        ObjectiveSense.Minimize: gurobipy.GRB.MINIMIZE
    }

    def __init__(self, **kwargs):
        self._p = gurobipy.Model()

        # TODO Logging should use the Python logging system like every other
        # part of PSAMM. Unfortunately, Gurobi does not make this easy as the
        # logging output is directed to stdout/(stderr?) and/or a file. For
        # now, let's just disable logging to console so we don't mess up any
        # of our output. This can be done with the LogToConsole parameter which
        # disables output to console EXCEPT that when this parameter is set
        # a log message is generated explaining that the parameter was
        # changed... We are fortunate enough that this does not happen when
        # using the OutputFlag parameter but instead we lose the log file.
        self._p.params.OutputFlag = 0

        # Set tolerances. By default, we decrease to 1e-9.
        self.feasibility_tolerance.value = kwargs.get(
            'feasibility_tolerance', 1e-9)
        self.optimality_tolerance.value = kwargs.get(
            'optimality_tolerance', 1e-9)

        # Set number of threads
        if 'threads' in kwargs:
            logger.info('Setting threads to {!r}'.format(kwargs['threads']))
            self._p.params.Threads = kwargs['threads']

        self._variables = {}
        self._var_names = ('x'+str(i) for i in count(1))
        self._constr_names = ('c'+str(i) for i in count(1))

        self._result = None

    @property
    def gurobi(self):
        """The underlying Gurobi Model object."""
        return self._p

    def define(self, *names, **kwargs):
        """Define variable in the problem.

        Variables must be defined before they can be accessed by var() or
        set(). This function takes keyword arguments lower and upper to define
        the bounds of the variable (default: -inf to inf). The keyword argument
        types can be used to select the type of the variable (Continuous
        (default), Binary or Integer). Setting any variables different than
        Continuous will turn the problem into an MILP problem.
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
        lower = (-gurobipy.GRB.INFINITY if value is None else float(value)
                 for value in lower)
        upper = (gurobipy.GRB.INFINITY if value is None else float(value)
                 for value in upper)
        vartype = tuple(VariableType.Continuous if value is None else value
                        for value in vartype)

        self._variables.update(zip(names, lp_names))
        for name, lb, ub, vt in zip(lp_names, lower, upper, vartype):
            self._p.addVar(name=name, lb=lb, ub=ub, vtype=self.VARTYPE_MAP[vt])

        self._p.update()

    def has_variable(self, name):
        """Check whether variable is defined in the model."""
        return name in self._variables

    def _grb_expr_from_value_set(self, value_set):
        linear = []
        quad = []

        for var, val in value_set:
            if not isinstance(var, Product):
                linear.append((
                    float(val), self._p.getVarByName(self._variables[var])))
            else:
                if len(var) > 2:
                    raise ValueError('Invalid objective: {}'.format(var))
                quad.append((
                    float(val), self._p.getVarByName(self._variables[var[0]]),
                    self._p.getVarByName(self._variables[var[1]])))

        if len(quad) > 0:
            expr = gurobipy.QuadExpr()
            if len(linear) > 0:
                expr.addTerms(*zip(*linear))
            expr.addTerms(*zip(*quad))
        else:
            expr = gurobipy.LinExpr(linear)

        return expr

    def _add_constraints(self, relation):
        """Add the given relation as one or more constraints.

        Return a list of the names of the constraints added.
        """
        names = []
        expression = relation.expression
        for value_set in expression.value_sets():
            name = next(self._constr_names)
            self._p.addConstr(
                self._grb_expr_from_value_set(value_set),
                self.CONSTR_SENSE_MAP[relation.sense],
                -float(expression.offset), name)
            names.append(name)

        self._p.update()

        return names

    def add_linear_constraints(self, *relations):
        """Add constraints to the problem.

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.
        """
        constraints = []

        for relation in relations:
            if self._check_relation(relation):
                constraints.append(Constraint(self, None))
            else:
                for name in self._add_constraints(relation):
                    constraints.append(Constraint(self, name))

        return constraints

    def set_objective(self, expression):
        """Set linear objective of problem."""

        if isinstance(expression, numbers.Number):
            # Allow expressions with no variables as objective,
            # represented as a number
            expression = Expression()

        self._p.setObjective(
            self._grb_expr_from_value_set(expression.values()))

    set_linear_objective = set_objective
    """Set objective of the problem.

    .. deprecated:: 0.19
       Use :meth:`set_objective` instead.
    """

    def set_objective_sense(self, sense):
        """Set type of problem (maximize or minimize)."""

        if sense not in (ObjectiveSense.Minimize, ObjectiveSense.Maximize):
            raise ValueError('Invalid objective sense')

        self._p.ModelSense = self.OBJ_SENSE_MAP[sense]

    def solve(self, sense=None):
        """Solve problem."""
        if sense is not None:
            self.set_objective_sense(sense)
        self._p.optimize()

        self._result = Result(self)
        return self._result

    @property
    def result(self):
        return self._result

    @ranged_property(min=1e-9, max=1e-2)
    def feasibility_tolerance(self):
        """Feasibility tolerance."""
        return self._p.params.FeasibilityTol

    @feasibility_tolerance.setter
    def feasibility_tolerance(self, value):
        logger.info('Setting feasibility tolerance to {!r}'.format(value))
        try:
            self._p.params.FeasibilityTol = value
        except gurobipy.GurobiError as e:
            raise_from(ValueError(e.message), e)

    @ranged_property(min=1e-9, max=1e-2)
    def optimality_tolerance(self):
        """Optimality tolerance."""
        return self._p.params.OptimalityTol

    @optimality_tolerance.setter
    def optimality_tolerance(self, value):
        logger.info('Setting optimality tolerance to {!r}'.format(value))
        try:
            self._p.params.OptimalityTol = value
        except gurobipy.GurobiError as e:
            raise_from(ValueError(e.message), e)

    @ranged_property(min=1e-9, max=1e-1)
    def integrality_tolerance(self):
        """Integrality tolerance."""
        return self._p.params.IntFeasTol

    @integrality_tolerance.setter
    def integrality_tolerance(self, value):
        logger.info('Setting integrality tolerance to {!r}'.format(value))
        try:
            self._p.params.IntFeasTol = value
        except gurobipy.GurobiError as e:
            raise_from(ValueError(e.message), e)


class Constraint(BaseConstraint):
    """Represents a constraint in a gurobi.Problem."""

    def __init__(self, prob, name):
        self._prob = prob
        self._name = name

    def delete(self):
        if self._name is not None:
            self._prob._p.remove(self._prob._p.getConstrByName(self._name))


class Result(BaseResult):
    """Represents the solution to a gurobi.Problem.

    This object will be returned from the gurobi.Problem.solve() method or by
    accessing the gurobi.Problem.result property after solving a problem. This
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
        """Return boolean indicating whether a solution was found."""
        self._check_valid()
        return self._problem._p.Status == gurobipy.GRB.OPTIMAL

    @property
    def unbounded(self):
        """Whether solution is unbounded"""
        self._check_valid()

        status = self._problem._p.Status
        if (status == gurobipy.GRB.INF_OR_UNBD and
                self._problem._p.params.DualReductions):
            # Disable dual reductions to obtain a definitve answer
            self._problem._p.params.DualReductions = 0
            try:
                self._problem._p.optimize()
            finally:
                self._problem._p.params.DualReductions = 1

            status = self._problem._p.Status
        return status == gurobipy.GRB.UNBOUNDED

    @property
    def status(self):
        """Return string indicating the error encountered on failure."""
        self._check_valid()
        return str(self._problem._p.Status)

    def _get_value(self, var):
        """Return value of variable in solution."""
        return self._problem._p.getVarByName(self._problem._variables[var]).x

    def _has_variable(self, var):
        """Whether variable exists in the solution."""
        return self._problem.has_variable(var)

    def get_value(self, expression):
        """Return value of expression."""

        self._check_valid()
        return super(Result, self).get_value(expression)
