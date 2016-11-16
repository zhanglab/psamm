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

"""Linear programming solver using GLPK."""

from __future__ import absolute_import

import logging
from itertools import repeat, count
import numbers

from six.moves import zip

import swiglpk

from .lp import Solver as BaseSolver
from .lp import Constraint as BaseConstraint
from .lp import Problem as BaseProblem
from .lp import Result as BaseResult
from .lp import (Expression, RelationSense, ObjectiveSense, VariableType,
                 InvalidResultError, ranged_property)

# Module-level logging
logger = logging.getLogger(__name__)

# Logger specific to log messages from GLPK library
_glpk_logger = logging.getLogger('glpk')


# Redirect GLPK terminal output to logger
def _term_hook(s):
    _glpk_logger.debug(s.rstrip())


swiglpk.glp_term_hook(_term_hook)


_INF = float('inf')


class GLPKError(Exception):
    """Error from calling GLPK library."""


class Solver(BaseSolver):
    """Represents an LP-solver using Gurobi."""

    def __init__(self):
        super(Solver, self).__init__()
        logger.warn('Support for GLPK solver is experimental!')

    def create_problem(self, **kwargs):
        """Create a new LP-problem using the solver."""
        return Problem(**kwargs)


class Problem(BaseProblem):
    """Represents an LP-problem of a gurobi.Solver."""

    VARTYPE_MAP = {
        VariableType.Continuous: swiglpk.GLP_CV,
        VariableType.Binary: swiglpk.GLP_BV,
        VariableType.Integer: swiglpk.GLP_IV
    }

    OBJ_SENSE_MAP = {
        ObjectiveSense.Maximize: swiglpk.GLP_MAX,
        ObjectiveSense.Minimize: swiglpk.GLP_MIN
    }

    def __init__(self, **kwargs):
        self._p = swiglpk.glp_create_prob()

        self._variables = {}
        self._constraints = {}

        self._do_presolve = True

        # Initialize simplex tolerances to default values
        parm = swiglpk.glp_smcp()
        swiglpk.glp_init_smcp(parm)
        self._feasibility_tolerance = parm.tol_bnd
        self._optimality_tolerance = parm.tol_dj

        # Initialize mip tolerance to default value
        parm = swiglpk.glp_iocp()
        swiglpk.glp_init_iocp(parm)
        self._integrality_tolerance = parm.tol_int

        self._result = None

    def __del__(self):
        swiglpk.glp_delete_prob(self._p)

    @property
    def glpk(self):
        """The underlying GLPK (SWIG) object."""
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

        # Assign default values
        vartype = tuple(VariableType.Continuous if value is None else value
                        for value in vartype)

        if len(names) == 0:
            return

        var_indices = count(swiglpk.glp_add_cols(self._p, len(names)))

        for i, name, lb, ub, vt in zip(
                var_indices, names, lower, upper, vartype):
            self._variables[name] = i

            lb = None if lb == -_INF else lb
            ub = None if ub == _INF else ub

            if lb is None and ub is None:
                swiglpk.glp_set_col_bnds(self._p, i, swiglpk.GLP_FR, 0, 0)
            elif lb is None:
                swiglpk.glp_set_col_bnds(
                    self._p, i, swiglpk.GLP_UP, 0, float(ub))
            elif ub is None:
                swiglpk.glp_set_col_bnds(
                    self._p, i, swiglpk.GLP_LO, float(lb), 0)
            elif lb == ub:
                swiglpk.glp_set_col_bnds(
                    self._p, i, swiglpk.GLP_FX, float(lb), 0)
            else:
                swiglpk.glp_set_col_bnds(
                    self._p, i, swiglpk.GLP_DB, float(lb), float(ub))

            if vt != VariableType.Continuous:
                swiglpk.glp_set_col_kind(self._p, i, self.VARTYPE_MAP[vt])

        self._do_presolve = True

    def has_variable(self, name):
        """Check whether variable is defined in the model."""
        return name in self._variables

    def _add_constraints(self, relation):
        """Add the given relation as one or more constraints.

        Return a list of the names of the constraints added.
        """

        expression = relation.expression
        constr_count = sum(True for _ in expression.value_sets())
        if constr_count == 0:
            return []

        row_indices = count(swiglpk.glp_add_rows(self._p, constr_count))

        names = []
        for i, value_set in zip(row_indices, expression.value_sets()):
            value_set = list(value_set)
            var_indices = swiglpk.intArray(1 + len(value_set))
            var_values = swiglpk.doubleArray(1 + len(value_set))
            for j, (variable, coeff) in enumerate(value_set):
                var_indices[1 + j] = self._variables[variable]
                var_values[1 + j] = float(coeff)

            swiglpk.glp_set_mat_row(
                self._p, i, len(value_set), var_indices, var_values)

            if relation.sense == RelationSense.Greater:
                swiglpk.glp_set_row_bnds(
                    self._p, i, swiglpk.GLP_LO, -float(expression.offset), 0)
            elif relation.sense == RelationSense.Less:
                swiglpk.glp_set_row_bnds(
                    self._p, i, swiglpk.GLP_UP, 0, -float(expression.offset))
            else:
                swiglpk.glp_set_row_bnds(
                    self._p, i, swiglpk.GLP_FX, -float(expression.offset), 0)

            names.append(i)

        self._do_presolve = True

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
        """Set objective of problem."""

        if isinstance(expression, numbers.Number):
            # Allow expressions with no variables as objective,
            # represented as a number
            expression = Expression(offset=expression)

        # Clear previous objective
        for i in range(swiglpk.glp_get_num_cols(self._p)):
            swiglpk.glp_set_obj_coef(self._p, 1 + i, 0)

        for variable, value in expression.values():
            var_index = self._variables[variable]
            swiglpk.glp_set_obj_coef(self._p, var_index, value)

        swiglpk.glp_set_obj_coef(self._p, 0, expression.offset)

    set_linear_objective = set_objective
    """Set objective of the problem.

    .. deprecated:: 0.19
       Use :meth:`set_objective` instead.
    """

    def set_objective_sense(self, sense):
        """Set type of problem (maximize or minimize)."""

        if sense not in (ObjectiveSense.Minimize, ObjectiveSense.Maximize):
            raise ValueError('Invalid objective sense')

        swiglpk.glp_set_obj_dir(self._p, self.OBJ_SENSE_MAP[sense])

    def solve(self, sense=None):
        """Solve problem."""
        if sense is not None:
            self.set_objective_sense(sense)

        parm = swiglpk.glp_smcp()
        swiglpk.glp_init_smcp(parm)
        if self._do_presolve:
            parm.presolve = swiglpk.GLP_ON
            self._do_presolve = False
        else:
            parm.presolve = swiglpk.GLP_OFF

        parm.tol_bnd = self._feasibility_tolerance
        parm.tol_dj = self._optimality_tolerance

        logger.debug('Solving using glp_simplex()')
        r = swiglpk.glp_simplex(self._p, parm)
        if r in (swiglpk.GLP_ENOPFS, swiglpk.GLP_ENODFS):
            self._result = Result(self, r)
            return self._result
        elif r != 0:
            raise GLPKError('glp_simplex: {}'.format(r))

        if swiglpk.glp_get_num_int(self._p) == 0:
            self._result = Result(self)
        else:
            # The intopt MILP solver needs an optimal solution to the LP
            # relaxation. Therefore, glp_simplex has to run before glp_intopt
            # for MILP problems.
            logger.debug('Solving using glp_intopt()')
            parm = swiglpk.glp_iocp()
            swiglpk.glp_init_iocp(parm)
            parm.tol_int = self._integrality_tolerance
            r = swiglpk.glp_intopt(self._p, parm)
            if r != 0:
                raise GLPKError('glp_intopt: {}'.format(r))

            self._result = MIPResult(self)

        return self._result

    @property
    def result(self):
        return self._result

    @ranged_property(min=None, max=None)
    def feasibility_tolerance(self):
        """Feasibility tolerance."""
        return self._feasibility_tolerance

    @feasibility_tolerance.setter
    def feasibility_tolerance(self, value):
        logger.info('Setting feasibility tolerance to {!r}'.format(value))
        self._feasibility_tolerance = value

    @ranged_property(min=None, max=None)
    def optimality_tolerance(self):
        """Optimality tolerance."""
        return self._optimality_tolerance

    @optimality_tolerance.setter
    def optimality_tolerance(self, value):
        logger.info('Setting optimality tolerance to {!r}'.format(value))
        self._optimality_tolerance = value

    @ranged_property(min=None, max=None)
    def integrality_tolerance(self):
        """Integrality tolerance."""
        return self._integrality_tolerance

    @integrality_tolerance.setter
    def integrality_tolerance(self, value):
        logger.info('Setting integrality tolerance to {!r}'.format(value))
        self._integrality_tolerance = value


class Constraint(BaseConstraint):
    """Represents a constraint in a gurobi.Problem."""

    def __init__(self, prob, name):
        self._prob = prob
        self._name = name

    def delete(self):
        if self._name is not None:
            swiglpk.glp_set_row_bnds(
                self._prob._p, self._name, swiglpk.GLP_FR, 0, 0)


class Result(BaseResult):
    """Represents the solution to a gurobi.Problem.

    This object will be returned from the gurobi.Problem.solve() method or by
    accessing the gurobi.Problem.result property after solving a problem. This
    class should not be instantiated manually.

    Result will evaluate to a boolean according to the success of the
    solution, so checking the truth value of the result will immediately
    indicate whether solving was successful.
    """

    def __init__(self, prob, ret_val=0):
        self._problem = prob
        self._ret_val = ret_val

    def _check_valid(self):
        if self._problem.result != self:
            raise InvalidResultError()

    @property
    def success(self):
        """Return boolean indicating whether a solution was found."""
        self._check_valid()
        if self._ret_val != 0:
            return False
        return swiglpk.glp_get_status(self._problem._p) == swiglpk.GLP_OPT

    @property
    def unbounded(self):
        """Whether solution is unbounded"""
        self._check_valid()
        if self._ret_val != 0:
            return self._ret_val == swiglpk.GLP_ENODFS
        return swiglpk.glp_get_status(self._problem._p) == swiglpk.GLP_UNBND

    @property
    def status(self):
        """Return string indicating the error encountered on failure."""
        self._check_valid()
        if self._ret_val == swiglpk.GLP_ENOPFS:
            return 'No primal feasible solution'
        elif self._ret_val == swiglpk.GLP_ENODFS:
            return 'No dual feasible solution'
        return str(swiglpk.glp_get_status(self._problem._p))

    def _get_var_value(self, variable):
        self._check_valid()
        if variable not in self._problem._variables:
            raise ValueError('Unknown variable: {}'.format(variable))
        return swiglpk.glp_get_col_prim(
            self._problem._p, self._problem._variables[variable])

    def get_value(self, expression):
        """Return value of expression."""

        self._check_valid()
        if isinstance(expression, Expression):
            return sum(self._get_var_value(var) * value
                       for var, value in expression.values())
        return self._get_var_value(expression)


class MIPResult(Result):
    """Specialization of Result for MIP problems."""

    @property
    def success(self):
        self._check_valid()
        return swiglpk.glp_mip_status(self._problem._p) == swiglpk.GLP_OPT

    @property
    def status(self):
        self._check_valid()
        return str(swiglpk.glp_mip_status(self._problem._p))

    def _get_var_value(self, variable):
        self._check_valid()
        if variable not in self._problem._variables:
            raise ValueError('Unknown variable: {}'.format(variable))
        return swiglpk.glp_mip_col_val(
            self._problem._p, self._problem._variables[variable])
