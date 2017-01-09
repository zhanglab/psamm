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
        self._prob._cp.set_log_stream(None)
        self._prob._cp.set_error_stream(None)
        self._prob._cp.set_warning_stream(None)
        self._prob._cp.set_results_stream(None)
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

    # _solve is this function
    # solve exists inside the generic class
    def _solve(self, sense=None):
        # Remove temporary constraints
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()
            #print("popping constraints")

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
        return result

    def get_fba_flux(self, objective):
        flux_result = self._solve_fba(objective)
        # Creat a dictionary of all the reaction ids and their fluxes.
        # {'reaction1': 1000, 'reaction2': 2000, ...}
        fba_fluxes = {}

        # key = 'reaction1'
        for key in self._model.reactions:
            fba_fluxes[key] = flux_result.get_value(self._v_wt(key))
        return fba_fluxes


    def get_fba_biomass(self, objective):
        # Run the FBA
        flux_result = self._solve_fba(objective)
        # Pull out the value
        return flux_result.get_value(self._v_wt(objective))

    # Implemetnation of the LP2 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same as the wildtype.
    # This helps avoid the assumption that an organism will perform optimally directly after
    # removing a gene.
    def minimize_lp2(self, objective, wt_fluxes):
        # These are all of our non eachange reactions
        reactions = set(self._adjustment_reactions())
        # First times in the loop
        if self._z is None:
            self._z = z = self._prob.namespace()
            z.define(reactions, lower=0)
        # After the first time, we can just use the same z object as before
        else:
            z = self._z

        # Get the v model and all the reactions 'variables' associated with it
        v = self._v

        # Create a constraints object for this model
        constr = self.constraints()

        # Loop through all the flux reactions
        for f_reaction, f_value in iteritems(wt_fluxes):
            # Makes sure the reaction we are trying to put a constraint on is a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such that the differenct between the wildtype flux
                # is similar to the knockout flux.
                constr.add(z(f_reaction) >= f_value - v(f_reaction), f_value - v(f_reaction) >= -z(f_reaction))

        # If we minimize the sum of the z vector then we will minimize the |vs_wt - vs| from above
        self._prob.set_objective(z.sum(self._adjustment_reactions()))

        # The result of minimizing z
        result = self._solve(lp.ObjectiveSense.Minimize)


        # Need to make sure we catch any problems that occur in the solver
        if not result:
            print("ERRROR")
            raise MOMAError('Unable to solve LP2 MOMA: {}'.format(
                result.status))

        #print("Gene knockout flux: " + str(self.get_flux(objective)))

        # Print out the flux
        #print("Gene knockout flux: " + str(self.get_flux(objective)))

        # We need to go back and manually delete all the constaints
        constr.delete()

    def minimize_l1(self, objective, wt_obj, wt_fluxes):
        #wt_result = self._solve_fba(objective)
        #wt_obj = wt_result.get_value(self._v(objective))
        # These are all of our non eachange reactions
        reactions = set(self._adjustment_reactions())
        # First times in the loop
        if self._z is None:
            self._z = z = self._prob.namespace()
            z.define(reactions, lower=0)
        # After the first time, we can just use the same z object as before
        else:
            z = self._z

        # Get the v model and all the reactions 'variables' associated with it
        v = self._v

        v_wt = self._v_wt

        # Create a constraints object for this model
        constr = self.constraints()

        # Loop through all the flux reactions
        for f_reaction, f_value in iteritems(wt_fluxes):
            # Makes sure the reaction we are trying to put a constraint on is a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such that the differenct between the wildtype flux
                # is similar to the knockout flux.
                constr.add(z(f_reaction) >= v_wt(f_reaction) - v(f_reaction), v_wt(f_reaction) - v(f_reaction) >= -z(f_reaction))

        self._prob.set_objective(z.sum(self._adjustment_reactions()))

        v_wt_obj = self._v_wt(objective)
        with self.constraints(v_wt_obj == wt_obj):
            result = self._solve(lp.ObjectiveSense.Minimize)

        #print("Gene knockout flux: " + str(self.get_flux(objective)))

        if not result:
            raise MOMAError('Unable to solve L1 MOMA: {}'.format(
                result.status))

    #def minimize_l1(self, objective, wt_obj):
    #    #wt_result = self._solve_fba(objective)
    #    #wt_obj = wt_result.get_value(self._v(objective))
#
#        if self._z is None:
#            self._z = z = self._prob.namespace()
#            reactions = set(self._adjustment_reactions())
#            z.define(reactions, lower=0)
#
#            zs = z.set(reactions)
#            vs_wt = self._v_wt.set(reactions)
#            vs = self._v.set(reactions)
#
#            self.constraints(zs >= vs_wt - vs, vs_wt - vs >= -zs)
#        else:
#            z = self._z
#
#        self._prob.set_objective(z.sum(self._adjustment_reactions()))
#
#        v_wt_obj = self._v_wt(objective)
#        with self.constraints(v_wt_obj == wt_obj):
#            result = self._solve(lp.ObjectiveSense.Minimize)
#
#        #print("Gene knockout flux: " + str(self.get_flux(objective)))
#
#        if not result:
#            raise MOMAError('Unable to solve L1 MOMA: {}'.format(
#                result.status))

    # Implemetnation of the QLP2 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same as the wildtype.
    # This helps avoid the assumption that an organism will perform optimally directly after
    # removing a gene.
    def minimize_qlp2(self, objective, wt_fluxes):


        # These are all of our non eachange reactions
        reactions = set(self._adjustment_reactions())

        # Get the v model and all the reactions 'variables' associated with it
        v = self._v

        # Create a constraints object for this model
        constr = self.constraints()

        obj_expr = 0
        # Loop through all the flux reactions
        for f_reaction, f_value in iteritems(wt_fluxes):
            # Makes sure the reaction we are trying to put a constraint on is a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such that the differenct between the wildtype flux
                # is similar to the knockout flux.
                obj_expr += (f_value - v(f_reaction))**2

        # If we minimize the sum of the z vector then we will minimize the |vs_wt - vs| from above
        #self._prob.set_objective(z.sum(self._adjustment_reactions()))
        self._prob.set_objective(obj_expr)

        # The result of minimizing z
        result = self._solve(lp.ObjectiveSense.Minimize)

        # We need to go back and manually delete all the constaints
        constr.delete()

        # Need to make sure we catch any problems that occur in the solver
        if not result:
            print("ERROR")
            raise MOMAError('Unable to solve QLP2 MOMA: {}'.format(result.status))
            print("Gene knockout flux: " + str(self.get_flux(objective)))

        # Print out the flux
        #print("Gene knockout flux: " + str(self.get_flux(objective)))

    def minimize_l2(self, objective):
            wt_result = self._solve_fba(objective)
            wt_obj = wt_result.get_value(self._v_wt(objective))

            obj_expr = 0
            for reaction in self._adjustment_reactions():
                v_wt = self._v_wt(reaction)
                v = self._v(reaction)
                obj_expr += (v_wt - v)**2

            self._prob.set_objective(obj_expr)

            v_wt_obj = self._v_wt(objective)
            with self.constraints(v_wt_obj == wt_obj):
                result = self._solve(lp.ObjectiveSense.Minimize)



            #print("Gene knockout flux: " + str(self.get_flux(objective)))


    def get_flux(self, reaction):
        return self._prob.result.get_value(self._v(reaction))

    def get_flux_var(self, reaction):
        return self._v(reaction)
