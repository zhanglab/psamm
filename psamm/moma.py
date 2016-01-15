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
# Copyright 2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Implementation of Minimization of Metabolic Adjustments (MOMA)."""

import logging
from itertools import product

from six import iteritems

from psamm.lpsolver import lp

logger = logging.getLogger(__name__)


class MOMAError(Exception):
    """Error indicating an error solving MOMA."""


class ConstraintGroup(object):
    def __init__(self, moma, *args):
        self._moma = moma
        self._constrs = []
        if len(args) > 0:
            self.add(*args)

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        self.delete()

    def add(self, *args):
        self._constrs.extend(self._moma._prob.add_linear_constraints(*args))

    def delete(self):
        self._moma._remove_constr.extend(self._constrs)


class MOMAProblem(object):
    def __init__(self, model, solver):
        self._prob = solver.create_problem()
        self._model = model

        self._remove_constr = []

        self._l1_constraints = None

        # Define flux variables
        for reaction_id in self._model.reactions:
            lower, upper = self._model.limits[reaction_id]
            self._prob.define(('v_wt', reaction_id), lower=lower, upper=upper)
            self._prob.define(('v', reaction_id), lower=lower, upper=upper)

        # Define constraints
        mass_balance = self.constraints()
        massbalance_lhs = {
            spec: 0 for spec in product(model.compounds, ('wt', 'mod'))}
        for (compound, reaction_id), value in iteritems(self._model.matrix):
            v_wt = self._prob.var(('v_wt', reaction_id))
            v = self._prob.var(('v', reaction_id))
            massbalance_lhs[compound, 'wt'] += v_wt * value
            massbalance_lhs[compound, 'mod'] += v * value
        for compound, lhs in iteritems(massbalance_lhs):
            mass_balance.add(lhs == 0)

    @property
    def prob(self):
        return self._prob

    def constraints(self, *args):
        return ConstraintGroup(self, *args)

    def _adjustment_reactions(self):
        for reaction_id in self._model.reactions:
            if not self._model.is_exchange(reaction_id):
                yield reaction_id

    def _solve(self, sense=None):
        # Remove temporary constraints
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()

        try:
            return self._prob.solve(sense=sense)
        finally:
            self._remove_constr = []

    def _solve_fba(self, objective):
        v_wt = self._prob.var(('v_wt', objective))
        self._prob.set_objective(v_wt)

        result = self._solve(lp.ObjectiveSense.Maximize)
        if not result:
            raise MOMAError('Unable to solve initial FBA: {}'.format(
                result.status))

        return result.get_value(v_wt)

    def minimize_l1(self, objective):
        wt_obj = self._solve_fba(objective)

        if self._l1_constraints is None:
            reactions = set(self._adjustment_reactions())
            self._prob.define(*(('z', r) for r in reactions), lower=0)
            z = self._prob.set(('z', r) for r in reactions)
            v_wt = self._prob.set(('v_wt', r) for r in reactions)
            v = self._prob.set(('v', r) for r in reactions)

            self._l1_constraints = self.constraints(
                z >= v_wt - v, v_wt - v >= -z)

        obj_expr = self._prob.expr({
            ('z', r): 1 for r in self._adjustment_reactions()})
        self._prob.set_objective(obj_expr)

        v_wt_obj = self._prob.var(('v_wt', objective))
        with self.constraints(v_wt_obj == wt_obj):
            result = self._solve(lp.ObjectiveSense.Minimize)

        if not result:
            raise MOMAError('Unable to solve L1 MOMA: {}'.format(
                result.status))

    def get_flux(self, reaction):
        return self._prob.result.get_value(('v', reaction))

    def get_flux_var(self, reaction):
        return self._prob.var(('v', reaction))
