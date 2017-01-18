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
# Copyright 2016-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016-2017  Brian Bishop <brian_bishop@my.uri.edu>

"""Implementation of Minimization of Metabolic Adjustments (MOMA)."""

import logging
from itertools import product

from six import iteritems, text_type, raise_from

from psamm.lpsolver import lp

logger = logging.getLogger(__name__)


class MOMAError(Exception):
    """Error indicating an error solving MOMA."""


class ConstraintGroup(object):
    """Constraints that will be imposed on the model when solving.

    Args:
        moma: MOMAProblem object for the proposed constraints.
        *args: The constraints that are imposed on the model.
    """
    def __init__(self, moma, *args):
        self._moma = moma
        self._constrs = []
        if len(args) > 0:
            self.add(*args)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.delete()

    def add(self, *args):
        """Add a constraints to the model."""
        self._constrs.extend(self._moma._prob.add_linear_constraints(*args))

    def delete(self):
        """Set up the constraints to get deleted on the next solve."""
        self._moma._remove_constr.extend(self._constrs)


class MOMAProblem(object):
    """Model as a flux optimization problem with minimal flux redistribution.

    Create a representation of the model as an LP optimization problem with
    steady state assumption and a minimal redistribution of metabolic fluxes
    with respect to the wild type configuration.

    The problem can be solved using any of the four MOMA variants described in
    [Segre02]_ and [Mo09]_. MOMA is formulated to avoid the FBA assumption that
    that growth efficiency has evolved to an optimal point directly following
    model perturbation. MOMA finds the optimal solution for a model with
    minimal flux redistribution with respect to the wild type flux
    configuration.

    MOMA is implemented with two variations of a linear optimization problem
    (:meth:`.lin_moma` and :meth:`.lin_moma2`) and two variations of a
    quadratic optimization problem (:meth:`.moma` and :meth:`.moma2`).
    Further information on these methods can be found within their respective
    documentation.

    The problem can be modified and solved as many times as needed. The flux
    of a reaction can be obtained after solving using :meth:`.get_flux`.

    Args:
        model: :class:`MetabolicModel` to solve.
        solver: LP solver instance to use.
    """
    def __init__(self, model, solver):
        self._prob = solver.create_problem()
        self._model = model

        self._remove_constr = []

        self._v_wt = v_wt = self._prob.namespace(name='v_wt')
        self._v = v = self._prob.namespace(name='v')
        self._z = z = self._prob.namespace(name='z')
        self._z_diff = z_diff = self._prob.namespace(name='z_diff')

        # Define flux variables
        for reaction_id in self._model.reactions:
            lower, upper = self._model.limits[reaction_id]
            v_wt.define([reaction_id], lower=lower, upper=upper)
            v.define([reaction_id], lower=lower, upper=upper)
            z.define([reaction_id], lower=0, upper=max(abs(lower), abs(upper)))

            flux_range = float(upper) - float(lower)
            z_diff.define([reaction_id], lower=0, upper=flux_range)

        # Define constraints
        mass_balance = self.constraints()
        massbalance_lhs = {
            spec: 0 for spec in product(model.compounds, ('wt', 'mod'))}
        for (compound, reaction_id), value in iteritems(self._model.matrix):
            massbalance_lhs[compound, 'wt'] += v_wt[reaction_id] * value
            massbalance_lhs[compound, 'mod'] += v[reaction_id] * value
        for compound, lhs in iteritems(massbalance_lhs):
            mass_balance.add(lhs == 0)

    @property
    def prob(self):
        """Return the underlying LP problem."""
        return self._prob

    def constraints(self, *args):
        """Return a constraint object."""
        return ConstraintGroup(self, *args)

    def _adjustment_reactions(self):
        """Yield all the non exchange reactions in the model."""
        for reaction_id in self._model.reactions:
            if not self._model.is_exchange(reaction_id):
                yield reaction_id

    def _solve(self, sense=None):
        """Remove old constraints and then solve the current problem.

        Args:
            sense: Minimize or maximize the objective.
                (:class:`.lp.ObjectiveSense)

        Returns:
            The Result object for the solved LP problem
        """
        # Remove the constraints from the last run
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()

        try:
            return self._prob.solve(sense=sense)
        except lp.SolverError as e:
            raise_from(MOMAError(text_type(e)), e)
        finally:
            self._remove_constr = []

    def solve_fba(self, objective):
        """Solve the wild type problem using FBA.

        Args:
            objective: The objective reaction to be maximized.

        Returns:
            The LP Result object for the solved FBA problem.
        """
        self._prob.set_objective(self._v_wt[objective])
        return self._solve(lp.ObjectiveSense.Maximize)

    def get_fba_flux(self, objective):
        """Return a dictionary of all the fluxes solved by FBA.

        Dictionary of fluxes is used in :meth:`.lin_moma` and :meth:`.moma`
        to minimize changes in the flux distributions following model
        perturbation.

        Args:
            objective: The objective reaction that is maximized.

        Returns:
            Dictionary of fluxes for each reaction in the model.
        """
        flux_result = self.solve_fba(objective)
        fba_fluxes = {}

        # Place all the flux values in a dictionary
        for key in self._model.reactions:
            fba_fluxes[key] = flux_result.get_value(self._v_wt[key])
        return fba_fluxes

    def get_minimal_fba_flux(self, objective):
        """Find the FBA solution that minimizes all the flux values.

        Maximize the objective flux then minimize all other fluxes
        while keeping the objective flux at the maximum.

        Args:
            objective: The objective reaction that is maximized.

        Returns:
            A dictionary of all the reactions and their minimized fluxes.
        """
        if self._z is None:
            self._define_z_vars()

        # Define constraints
        vs_wt = self._v_wt.set(self._model.reactions)
        zs = self._z.set(self._model.reactions)

        wt_obj_flux = self.get_fba_obj_flux(objective)

        with self.constraints() as constr:
            constr.add(
                zs >= vs_wt, vs_wt >= -zs,
                self._v_wt[objective] >= wt_obj_flux)
            self._prob.set_objective(self._z.sum(self._model.reactions))
            result = self._solve(lp.ObjectiveSense.Minimize)

        fba_fluxes = {}
        for key in self._model.reactions:
            fba_fluxes[key] = result.get_value(self._v_wt[key])
        return fba_fluxes

    def minimize_qlp2(self, objective, wt_fluxes):
        """Implemetnation of the MOMA algorithm.

        We try to maximize biomass, by keeping the fluxes relitively the same
        as the wildtype. This helps avoid the assumption that an organism will
        perform optimally directly after removing a gene.

        Minimizes sum of (wild type - knockout)^2

        Args:
            objective: Biomass reaction for the model.
            wt_fluxes: All the wt_fluxes found by FBA in the form of a
                dictionary. Use :meth: _get_fba_flux(objective) to return
                a dictionary.
        """
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
        """Return the knockout flux for a specific reaction."""
        return self._prob.result.get_value(self._v[reaction])

    def get_flux_var(self, reaction):
        """Return the LP variable for a specific reaction."""
        return self._v[reaction]
