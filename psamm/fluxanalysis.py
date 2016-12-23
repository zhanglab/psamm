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

"""Implementation of Flux Balance Analysis."""

from __future__ import unicode_literals

import logging
import random

from .lpsolver import lp

from six import iteritems

# Module-level logging
logger = logging.getLogger(__name__)

_INF = float('inf')


def _get_fba_problem(model, tfba, solver):
    """Convenience function for returning the right FBA problem instance"""
    p = FluxBalanceProblem(model, solver)
    if tfba:
        p.add_thermodynamic()
    return p


class FluxBalanceError(Exception):
    """Error indicating that a flux balance cannot be solved."""

    def __init__(self, *args, **kwargs):
        self._result = kwargs.pop('result')
        super(FluxBalanceError, self).__init__(*args, **kwargs)

    @property
    def result(self):
        return self._result


class FluxBalanceProblem(object):
    """Model as a flux optimization problem with steady state assumption.

    Create a representation of the model as an LP optimization problem with
    steady state assumption, i.e. the concentrations of compounds are always
    zero.

    The problem can be modified and solved as many times as needed. The flux
    of a reaction can be obtained after solving using :meth:`.get_flux`.

    Args:
        model: :class:`MetabolicModel` to solve.
        solver: LP solver instance to use.
    """

    def __init__(self, model, solver):
        self._prob = solver.create_problem()
        self._model = model

        self._z = None
        self._v = v = self._prob.namespace()

        self._has_minimization_vars = False

        # We keep track of temporary constraints from the last optimization so
        # that we can remove them during the next call to solve(). This is
        # necessary because removing the constraint immediately after solving
        # a model can invalidate the solution in some solvers.
        self._temp_constr = []
        self._remove_constr = []

        # Define flux variables
        for reaction_id in self._model.reactions:
            lower, upper = self._model.limits[reaction_id]
            v.define([reaction_id], lower=lower, upper=upper)

        # Define constraints
        massbalance_lhs = {compound: 0 for compound in model.compounds}
        for spec, value in iteritems(self._model.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value
        for compound, lhs in iteritems(massbalance_lhs):
            self._prob.add_linear_constraints(lhs == 0)

    @property
    def prob(self):
        """Return the underlying LP problem.

        This can be used to add additional constraints on the problem. Calling
        solve on the underlying problem is not guaranteed to work correctly,
        instead use the methods on this object that solves the problem or
        make a subclass with a method that calls :meth:`_solve`.
        """
        return self._prob

    def add_thermodynamic(self, em=1000):
        """Apply thermodynamic constraints to the model.

        Adding these constraints restricts the solution space to only
        contain solutions that have no internal loops [Schilling00]_. This is
        solved as a MILP problem as described in [Muller13]_. The time to solve
        a problem with thermodynamic constraints is usually much longer than a
        normal FBA problem.

        The ``em`` parameter is the upper bound on the delta mu reaction
        variables. This parameter has to be balanced based on the model size
        since setting the value too low can result in the correct solutions
        being infeasible and setting the value too high can result in
        numerical instability which again makes the correct solutions
        infeasible. The default value should work in all cases as long as the
        model is not unusually large.
        """

        internal = set(r for r in self._model.reactions
                       if not self._model.is_exchange(r))

        # Reaction fluxes
        v = self._v

        # Indicator variable
        alpha = self._prob.namespace(internal, types=lp.VariableType.Binary)

        # Delta mu is the stoichiometrically weighted sum of the compound mus.
        dmu = self._prob.namespace(internal)

        for reaction_id in self._model.reactions:
            if not self._model.is_exchange(reaction_id):
                flux = v(reaction_id)
                alpha_r = alpha(reaction_id)
                dmu_r = dmu(reaction_id)

                lower, upper = self._model.limits[reaction_id]

                # Constrain the reaction to a direction determined by alpha
                # and contrain the delta mu to a value in [-em; -1] if
                # alpha is one, otherwise in [1; em].
                self._prob.add_linear_constraints(
                    flux >= lower * (1 - alpha_r),
                    flux <= upper * alpha_r,
                    dmu_r >= -em * alpha_r + (1 - alpha_r),
                    dmu_r <= em * (1 - alpha_r) - alpha_r)

        # Define mu variables
        mu = self._prob.namespace(self._model.compounds)

        tdbalance_lhs = {reaction_id: 0
                         for reaction_id in self._model.reactions}
        for spec, value in iteritems(self._model.matrix):
            compound, reaction_id = spec
            if not self._model.is_exchange(reaction_id):
                tdbalance_lhs[reaction_id] += mu(compound) * value
        for reaction_id, lhs in iteritems(tdbalance_lhs):
            if not self._model.is_exchange(reaction_id):
                self._prob.add_linear_constraints(lhs == dmu(reaction_id))

    def maximize(self, reaction):
        """Solve the model by maximizing the given reaction.

        If reaction is a dictionary object, each entry is interpreted as a
        weight on the objective for that reaction (non-existent reaction will
        have zero weight).
        """

        self._prob.set_objective(self.flux_expr(reaction))
        self._solve()

    def flux_bound(self, reaction, direction):
        """Return the flux bound of the reaction.

        Direction must be a positive number to obtain the upper bound or a
        negative number to obtain the lower bound. A value of inf or -inf is
        returned if the problem is unbounded.
        """
        try:
            self.maximize({reaction: direction})
        except FluxBalanceError as e:
            if not e.result.unbounded:
                raise
            return direction * _INF
        else:
            return self.get_flux(reaction)

    def _add_minimization_vars(self):
        """Add variables and constraints for L1 norm minimization."""

        self._z = self._prob.namespace(self._model.reactions, lower=0)

        # Define constraints
        v = self._v.set(self._model.reactions)
        z = self._z.set(self._model.reactions)

        self._prob.add_linear_constraints(z >= v, v >= -z)

    def minimize_l1(self, weights={}):
        """Solve the model by minimizing the L1 norm of the fluxes.

        If the weights dictionary is given, the weighted L1 norm if minimized
        instead. The dictionary contains the weights of each reaction
        (default 1).
        """

        if self._z is None:
            self._add_minimization_vars()

        objective = self._z.expr(
            (reaction_id, -weights.get(reaction_id, 1))
            for reaction_id in self._model.reactions)
        self._prob.set_objective(objective)

        self._solve()

    def max_min_l1(self, reaction, weights={}):
        """Maximize flux of reaction then minimize the L1 norm.

        During minimization the given reaction will be fixed at the maximum
        obtained from the first solution. If reaction is a dictionary object,
        each entry is interpreted as a weight on the objective for that
        reaction (non-existent reaction will have zero weight).
        """

        self.maximize(reaction)

        if isinstance(reaction, dict):
            reactions = list(reaction)
        else:
            reactions = [reaction]

        # Save flux values before modifying the LP problem
        fluxes = {r: self.get_flux(r) for r in reactions}

        # Add constraints on the maximized reactions
        for r in reactions:
            flux_var = self.get_flux_var(r)
            c, = self._prob.add_linear_constraints(flux_var == fluxes[r])
            self._temp_constr.append(c)

        self.minimize_l1(weights)

    def check_constraints(self):
        """Optimize without objective to check that solution is possible.

        Raises :class:`FluxBalanceError` if no flux solution is possible.
        """
        self._prob.set_objective(0)
        self._solve()

    def _solve(self):
        """Solve the problem with the current objective."""

        # Remove temporary constraints
        while len(self._remove_constr) > 0:
            self._remove_constr.pop().delete()

        try:
            result = self._prob.solve(lp.ObjectiveSense.Maximize)
            if not result:
                raise FluxBalanceError('Non-optimal solution: {}'.format(
                    result.status), result=result)
        finally:
            # Set temporary constraints to be removed on next solve call
            self._remove_constr = self._temp_constr
            self._temp_constr = []

    def get_flux_var(self, reaction):
        """Get LP variable representing the reaction flux."""
        return self._v(reaction)

    def flux_expr(self, reaction):
        """Get LP expression representing the reaction flux."""
        if isinstance(reaction, dict):
            return self._v.expr(iteritems(reaction))
        return self._v(reaction)

    def get_flux(self, reaction):
        """Get resulting flux value for reaction."""
        return self._prob.result.get_value(self._v(reaction))


def flux_balance(model, reaction, tfba, solver):
    """Run flux balance analysis on the given model.

    Yields the reaction id and flux value for each reaction in the model.

    This is a convenience function for sertting up and running the
    FluxBalanceProblem. If the FBA is solved for more than one parameter
    it is recommended to setup and reuse the FluxBalanceProblem manually
    for a speed up.

    This is an implementation of flux balance analysis (FBA) as described in
    [Orth10]_ and [Fell86]_.

    Args:
        model: MetabolicModel to solve.
        reaction: Reaction to maximize. If a dict is given, this instead
            represents the objective function weights on each reaction.
        tfba: If True enable thermodynamic constraints.
        solver: LP solver instance to use.

    Returns:
        Iterator over reaction ID and reaction flux pairs.
    """

    fba = _get_fba_problem(model, tfba, solver)
    fba.maximize(reaction)
    for reaction in model.reactions:
        yield reaction, fba.get_flux(reaction)


def flux_variability(model, reactions, fixed, tfba, solver):
    """Find the variability of each reaction while fixing certain fluxes.

    Yields the reaction id, and a tuple of minimum and maximum value for each
    of the given reactions. The fixed reactions are given in a dictionary as
    a reaction id to value mapping.

    This is an implementation of flux variability analysis (FVA) as described
    in [Mahadevan03]_.

    Args:
        model: MetabolicModel to solve.
        reactions: Reactions on which to report variablity.
        fixed: dict of additional lower bounds on reaction fluxes.
        tfba: If True enable thermodynamic constraints.
        solver: LP solver instance to use.

    Returns:
        Iterator over pairs of reaction ID and bounds. Bounds are returned as
        pairs of lower and upper values.
    """

    fba = _get_fba_problem(model, tfba, solver)

    for reaction_id, value in iteritems(fixed):
        flux = fba.get_flux_var(reaction_id)
        fba.prob.add_linear_constraints(flux >= value)

    def min_max_solve(reaction_id):
        for direction in (-1, 1):
            yield fba.flux_bound(reaction_id, direction)

    # Solve for each reaction
    for reaction_id in reactions:
        yield reaction_id, tuple(min_max_solve(reaction_id))


def flux_minimization(model, fixed, solver, weights={}):
    """Minimize flux of all reactions while keeping certain fluxes fixed.

    The fixed reactions are given in a dictionary as reaction id
    to value mapping. The weighted L1-norm of the fluxes is minimized.

    Args:
        model: MetabolicModel to solve.
        fixed: dict of additional lower bounds on reaction fluxes.
        solver: LP solver instance to use.
        weights: dict of weights on the L1-norm terms.

    Returns:
        An iterator of reaction ID and reaction flux pairs.
    """

    fba = FluxBalanceProblem(model, solver)

    for reaction_id, value in iteritems(fixed):
        flux = fba.get_flux_var(reaction_id)
        fba.prob.add_linear_constraints(flux >= value)

    fba.minimize_l1()

    return ((reaction_id, fba.get_flux(reaction_id))
            for reaction_id in model.reactions)


def flux_randomization(model, threshold, tfba, solver):
    """Find a random flux solution on the boundary of the solution space.

    The reactions in the threshold dictionary are constrained with the
    associated lower bound.

    Args:
        model: MetabolicModel to solve.
        threshold: dict of additional lower bounds on reaction fluxes.
        tfba: If True enable thermodynamic constraints.
        solver: LP solver instance to use.

    Returns:
        An iterator of reaction ID and reaction flux pairs.
    """

    optimize = {}
    for reaction_id in model.reactions:
        if model.is_reversible(reaction_id):
            optimize[reaction_id] = 2*random.random() - 1.0
        else:
            optimize[reaction_id] = random.random()

    fba = _get_fba_problem(model, tfba, solver)
    for reaction_id, value in iteritems(threshold):
        fba.prob.add_linear_constraints(fba.get_flux_var(reaction_id) >= value)

    fba.maximize(optimize)
    for reaction_id in model.reactions:
        yield reaction_id, fba.get_flux(reaction_id)


def consistency_check(model, subset, epsilon, tfba, solver):
    """Check that reaction subset of model is consistent using FBA.

    Yields all reactions that are *not* flux consistent. A reaction is
    consistent if there is at least one flux solution to the model that both
    respects the model constraints and also allows the reaction in question to
    have non-zero flux.

    This can be determined by running FBA on each reaction in turn
    and checking whether the flux in the solution is non-zero. Since FBA
    only tries to maximize the flux (and the flux can be negative for
    reversible reactions), we have to try to both maximize and minimize
    the flux. An optimization to this method is implemented such that if
    checking one reaction results in flux in another unchecked reaction,
    that reaction will immediately be marked flux consistent.

    Args:
        model: MetabolicModel to check for consistency.
        subset: Subset of model reactions to check.
        epsilon: The threshold at which the flux is considered non-zero.
        tfba: If True enable thermodynamic constraints.
        solver: LP solver instance to use.

    Returns:
        An iterator of flux inconsistent reactions in the subset.
    """

    fba = _get_fba_problem(model, tfba, solver)

    subset = set(subset)
    while len(subset) > 0:
        reaction = next(iter(subset))

        logger.info('{} left, checking {}...'.format(len(subset), reaction))

        fba.maximize(reaction)
        subset = set(reaction_id for reaction_id in subset
                     if abs(fba.get_flux(reaction_id)) <= epsilon)
        if reaction not in subset:
            continue
        elif model.is_reversible(reaction):
            fba.maximize({reaction: -1})
            subset = set(reaction_id for reaction_id in subset
                         if abs(fba.get_flux(reaction_id)) <= epsilon)
            if reaction not in subset:
                continue

        logger.info('{} not consistent!'.format(reaction))

        yield reaction
        subset.remove(reaction)
