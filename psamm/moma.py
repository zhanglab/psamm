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
    # Returns the problem that we are working on
    def prob(self):
        return self._prob

    # Creates a constraints object
    def constraints(self, *args):
        return ConstraintGroup(self, *args)

    # Returns a generator of all the non exchange reactions in the model
    def _adjustment_reactions(self):
        for reaction_id in self._model.reactions:
            if not self._model.is_exchange(reaction_id):
                yield reaction_id

    # Solve the linear programming problem
    # Sense tells the solver to minimize or maximize the result
    def _solve(self, sense=None):
        # Remove temporary constraints
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()

        # Try to solve the problem
        try:
            return self._prob.solve(sense=sense)
        # Reset the remove constraints list before the next solve
        finally:
            self._remove_constr = []

    # Solves for the standard FBA of the problem
    # Objective is the biomass
    def _solve_fba(self, objective):
        self._prob.set_objective(self._v_wt(objective))

        # Solve and store the result
        result = self._solve(lp.ObjectiveSense.Maximize)

        # If no solution was found an error is raised and told to the user
        if not result:
            raise MOMAError('Unable to solve initial FBA: {}'.format(
                result.status))

        # Return the result object
        return result

    # Returns all the flux values for each of the reactions. This is used in
    # MOMA LP2 and QLP2 to minimize the change in the flux flow following a
    # gene deletion.
    def get_fba_flux(self, objective):
        # Get the result of the wild type FBA problem
        flux_result = self._solve_fba(objective)

        # Creat a dictionary of all the reaction ids and their fluxes.
        # {'reaction1': 1000, 'reaction2': 2000, ...}
        fba_fluxes = {}

        # Place all the flux values in a dictionary
        for key in self._model.reactions:
            fba_fluxes[key] = flux_result.get_value(self._v_wt(key))

        # Return the dictionary with all of the fluxes associated with the
        # reaction id
        # {'reaction1': 1000, 'reaction2': 2000, ...}
        return fba_fluxes

    # Solves for the FBA biomass. Uses the result returned by _solve_fba and
    # finds the value of the objective (biomass) function.
    def get_fba_biomass(self, objective):
        # Run FBA and store the result
        flux_result = self._solve_fba(objective)
        # Pull out the biomass value
        return flux_result.get_value(self._v_wt(objective))

    # Implemetnation of the LP2 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same
    # as the wildtype. This helps avoid the assumption that an organism will
    # perform optimally directly after removing a gene.
    # This might not be the optimal minimization solution for the wild type
    # and the knockout vectors.
    # Minimizes sum of |wild type - knockout|
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
            # Makes sure the reaction we are trying to put a constraint on is
            # a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such
                # that the difference between the wildtype flux is similar
                # to the knockout flux.
                constr.add(z(f_reaction) >= f_value - v(f_reaction),
                    f_value - v(f_reaction) >= -z(f_reaction))

        # If we minimize the sum of the z vector then we will minimize
        # the |vs_wt - vs| from above
        self._prob.set_objective(z.sum(self._adjustment_reactions()))

        # Store the result of minimizing the sum of z
        result = self._solve(lp.ObjectiveSense.Minimize)

        # Need to make sure we catch any problems that occur in the solver
        if not result:
            raise MOMAError('Unable to solve LP2 MOMA: {}'.format(
                result.status))

        # Go back and manually delete all the constaints
        constr.delete()

    # Implemetnation of the LP3 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same
    # as the wildtype. This helps avoid the assumption that an organism will
    # perform optimally directly after removing a gene.
    # Creates the constraint that the we select the optimal flux vector that
    # is closest to the wildtype.
    # Minimizes sum of |wild type - knockout|
    def minimize_l1(self, objective, wt_obj, wt_fluxes):
        # These are all of our non eachange reactions
        reactions = set(self._adjustment_reactions())
        # First times in the loop
        if self._z is None:
            self._z = z = self._prob.namespace()
            z.define(reactions, lower=0)
        # After the first time, we can just use the same z object as before
        else:
            z = self._z

        # Get the v (knockout) and v_wt (wild type) model and all the
        # reactions 'variables' associated with it
        v = self._v
        v_wt = self._v_wt

        # Create a constraints object for this model
        constr = self.constraints()

        # Loop through all the flux reactions
        for f_reaction, f_value in iteritems(wt_fluxes):
            # Makes sure the reaction we are trying to put a constraint on
            # is a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such
                # that the difference between the wildtype flux
                # is similar to the knockout flux.
                constr.add(z(f_reaction) >= v_wt(f_reaction) - v(f_reaction),
                    v_wt(f_reaction) - v(f_reaction) >= -z(f_reaction))

        # Set the objective for the sum of z
        self._prob.set_objective(z.sum(self._adjustment_reactions()))

        # Set a constraint that says the biomass flux of the wild type is
        # the same as the FBA analysis
        v_wt_obj = self._v_wt(objective)
        with self.constraints(v_wt_obj == wt_obj):
            result = self._solve(lp.ObjectiveSense.Minimize)

        # Raise any errors
        if not result:
            raise MOMAError('Unable to solve L1 MOMA: {}'.format(
                result.status))
    # Implemetnation of the QLP2 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same
    # as the wildtype. This helps avoid the assumption that an organism will
    # perform optimally directly after removing a gene.
    # Minimizes sum of (wild type - knockout)^2
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
            # Makes sure the reaction we are trying to put a constraint on is
            # a non exchange reaction
            if f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such
                # that the differenct between the wildtype flux
                # is similar to the knockout flux.
                obj_expr += (f_value - v(f_reaction))**2

        # Set the objective
        self._prob.set_objective(obj_expr)

        # Minimize the squared distance between the wild type flux values and
        # the knockout values
        result = self._solve(lp.ObjectiveSense.Minimize)

        # We need to go back and manually delete all the constaints
        constr.delete()

        # Need to make sure we catch any problems that occur in the solver
        if not result:
            raise MOMAError('Unable to solve QLP2 MOMA: {}'.format(
                result.status))

    # Implemetnation of the QLP3 MOMA algorithm
    # We try to maximize biomass, by keeping the fluxes relitively the same
    # as the wildtype. This helps avoid the assumption that an organism will
    # perform optimally directly after removing a gene.
    # Minimizes sum of (wild type - knockout)^2
    def minimize_l2(self, objective):
            # Solve using FBA
            wt_result = self._solve_fba(objective)

            # Extract the result of the FBA biomass
            wt_obj = wt_result.get_value(self._v_wt(objective))

            obj_expr = 0

            # Loop through all the non exchange reactions and add a constraint
            # to reduce the squared difference between the wild type and
            # the knockout
            for reaction in self._adjustment_reactions():
                v_wt = self._v_wt(reaction)
                v = self._v(reaction)
                obj_expr += (v_wt - v)**2

            self._prob.set_objective(obj_expr)

            # Set a constraint that says the biomass flux of the wild type is
            # the same as the FBA analysis
            v_wt_obj = self._v_wt(objective)
            with self.constraints(v_wt_obj == wt_obj):
                result = self._solve(lp.ObjectiveSense.Minimize)

            # Need to make sure we catch any problems that occur in the solver
            if not result:
                raise MOMAError('Unable to solve QLP3 MOMA: {}'.format(
                    result.status))

    # Get the flux of the knockout model
    def get_flux(self, reaction):
        return self._prob.result.get_value(self._v(reaction))

    # Get the variable object for a specific reaction
    def get_flux_var(self, reaction):
        return self._v(reaction)
