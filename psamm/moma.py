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

        self._v_wt = v_wt = self._prob.namespace()
        self._v = v = self._prob.namespace()
        self._z = None

        # Define flux variables
        for reaction_id in self._model.reactions:
            lower, upper = self._model.limits[reaction_id]
            v_wt.define([reaction_id], lower=lower, upper=upper)
            v.define([reaction_id], lower=lower, upper=upper)

        # Define constraints
        mass_balance = self.constraints()
        massbalance_lhs = {
            spec: 0 for spec in product(model.compounds, ('wt', 'mod'))}
        for (compound, reaction_id), value in iteritems(self._model.matrix):
            massbalance_lhs[compound, 'wt'] += v_wt(reaction_id) * value
            massbalance_lhs[compound, 'mod'] += v(reaction_id) * value
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
        self._prob.set_objective(self._v_wt(objective))

        result = self._solve(lp.ObjectiveSense.Maximize)
        if not result:
            raise MOMAError('Unable to solve initial FBA: {}'.format(
                result.status))

        return result.get_value(self._v_wt(objective))

    def minimize_l1(self, objective):
        wt_obj = self._solve_fba(objective)

        if self._z is None:
            self._z = z = self._prob.namespace()
            reactions = set(self._adjustment_reactions())
            z.define(reactions, lower=0)

            zs = z.set(reactions)
            vs_wt = self._v_wt.set(reactions)
            vs = self._v.set(reactions)

            self.constraints(zs >= vs_wt - vs, vs_wt - vs >= -zs)
        else:
            z = self._z

        self._prob.set_objective(z.sum(self._adjustment_reactions()))

        v_wt_obj = self._v_wt(objective)
        with self.constraints(v_wt_obj == wt_obj):
            result = self._solve(lp.ObjectiveSense.Minimize)

        if not result:
            raise MOMAError('Unable to solve L1 MOMA: {}'.format(
                result.status))

    def minimize_l2(self, objective):
        wt_obj = self._solve_fba(objective)

        obj_expr = 0
        for reaction in self._adjustment_reactions():
            v_wt = self._v_wt(reaction)
            v = self._v(reaction)
            obj_expr += (v_wt - v)**2

        self._prob.set_objective(obj_expr)

        v_wt_obj = self._v_wt(objective)
        with self.constraints(v_wt_obj == wt_obj):
            result = self._solve(lp.ObjectiveSense.Minimize)

        if not result:
            raise MOMAError('Unable to solve L2 MOMA: {}'.format(
                result.status))

    def get_flux(self, reaction):
        return self._prob.result.get_value(self._v(reaction))

    def get_flux_var(self, reaction):
        return self._v(reaction)
