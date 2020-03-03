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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

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
        """Add constraints to the model."""
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

    def get_fba_obj_flux(self, objective):
        """Return the maximum objective flux solved by FBA."""
        flux_result = self.solve_fba(objective)
        return flux_result.get_value(self._v_wt[objective])

    def lin_moma(self, wt_fluxes):
        """Minimize the redistribution of fluxes using a linear objective.

        The change in flux distribution is mimimized by minimizing the sum
        of the absolute values of the differences of wild type FBA solution
        and the knockout strain flux solution.

        This formulation bases the solution on the wild type fluxes that
        are specified by the user. If these wild type fluxes were calculated
        using FBA, then an arbitrary flux vector that optimizes the objective
        function is used. See [Segre`_02] for more information.

        Args:
            wt_fluxes: Dictionary of all the wild type fluxes. Use
                :meth:`.get_fba_flux(objective)` to return a dictionary of
                fluxes found by FBA.
        """
        reactions = set(self._adjustment_reactions())

        z_diff = self._z_diff
        v = self._v

        with self.constraints() as constr:
            for f_reaction, f_value in iteritems(wt_fluxes):
                if f_reaction in reactions:
                    # Add the constraint that finds the optimal solution, such
                    # that the difference between the wildtype flux is similar
                    # to the knockout flux.
                    constr.add(
                        z_diff[f_reaction] >= f_value - v[f_reaction],
                        f_value - v[f_reaction] >= -z_diff[f_reaction])

            # If we minimize the sum of the z vector then we will minimize
            # the |vs_wt - vs| from above
            self._prob.set_objective(z_diff.sum(reactions))

            self._solve(lp.ObjectiveSense.Minimize)

    def lin_moma2(self, objective, wt_obj):
        """Find the smallest redistribution vector using a linear objective.

        The change in flux distribution is mimimized by minimizing the sum
        of the absolute values of the differences of wild type FBA solution
        and the knockout strain flux solution.

        Creates the constraint that the we select the optimal flux vector that
        is closest to the wildtype. This might still return an arbitrary flux
        vector the maximizes the objective function.

        Args:
            objective: Objective reaction for the model.
            wt_obj: The flux value for your wild type objective reactions.
                Can either use an expiremental value or on determined by FBA
                by using :meth:`.get_fba_obj_flux(objective)`.
        """
        reactions = set(self._adjustment_reactions())

        z_diff = self._z_diff
        v = self._v
        v_wt = self._v_wt

        with self.constraints() as constr:
            for f_reaction in reactions:
                # Add the constraint that finds the optimal solution, such
                # that the difference between the wildtype flux
                # is similar to the knockout flux.
                constr.add(
                    z_diff[f_reaction] >= v_wt[f_reaction] - v[f_reaction],
                    v_wt[f_reaction] - v[f_reaction] >= -z_diff[f_reaction])

            # If we minimize the sum of the z vector then we will minimize
            # the |v_wt - v| from above
            self._prob.set_objective(z_diff.sum(reactions))

            constr.add(self._v_wt[objective] >= wt_obj)

            self._solve(lp.ObjectiveSense.Minimize)

    def moma(self, wt_fluxes):
        """Minimize the redistribution of fluxes using Euclidean distance.

        Minimizing the redistribution of fluxes using a quadratic objective
        function. The distance is minimized by minimizing the sum of
        (wild type - knockout)^2.

        Args:
            wt_fluxes: Dictionary of all the wild type fluxes that will be
                used to find a close MOMA solution. Fluxes can be expiremental
                or calculated using :meth: get_fba_flux(objective).
        """
        reactions = set(self._adjustment_reactions())
        v = self._v

        obj_expr = 0
        for f_reaction, f_value in iteritems(wt_fluxes):
            if f_reaction in reactions:
                # Minimize the Euclidean distance between the two vectors
                obj_expr += (f_value - v[f_reaction])**2

        self._prob.set_objective(obj_expr)
        self._solve(lp.ObjectiveSense.Minimize)

    def moma2(self, objective, wt_obj):
        """Find the smallest redistribution vector using Euclidean distance.

        Minimizing the redistribution of fluxes using a quadratic objective
        function. The distance is minimized by minimizing the sum of
        (wild type - knockout)^2.

        Creates the constraint that the we select the optimal flux vector that
        is closest to the wildtype. This might still return an arbitrary flux
        vector the maximizes the objective function.

        Args:
            objective: Objective reaction for the model.
            wt_obj: The flux value for your wild type objective reactions.
                Can either use an expiremental value or on determined by FBA
                by using :meth:`.get_fba_obj_flux(objective)`.
        """
        obj_expr = 0
        for reaction in self._adjustment_reactions():
            v_wt = self._v_wt[reaction]
            v = self._v[reaction]
            obj_expr += (v_wt - v)**2

        self._prob.set_objective(obj_expr)

        with self.constraints(self._v_wt[objective] >= wt_obj):
            self._solve(lp.ObjectiveSense.Minimize)

    def get_flux(self, reaction):
        """Return the knockout flux for a specific reaction."""
        return self._prob.result.get_value(self._v[reaction])

    def get_flux_var(self, reaction):
        """Return the LP variable for a specific reaction."""
        return self._v[reaction]
