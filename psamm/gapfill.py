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

"""Identify blocked metabolites and possible reconstructions.

This implements a variant of the algorithms described in [Kumar07]_.
"""

import logging

from six import iteritems

from .lpsolver import lp

logger = logging.getLogger(__name__)


class GapFillError(Exception):
    """Indicates an error while running GapFind/GapFill"""


def _find_integer_tolerance(epsilon, v_max, min_tol):
    """Find appropriate integer tolerance for gap-filling problems."""
    int_tol = min(epsilon / (10 * v_max), 0.1)
    min_tol = max(1e-10, min_tol)
    if int_tol < min_tol:
        eps_lower = min_tol * 10 * v_max
        logger.warning(
            'When the maximum flux is {}, it is recommended that'
            ' epsilon > {} to avoid numerical issues with this'
            ' solver. Results may be incorrect with'
            ' the current settings!'.format(v_max, eps_lower))
        return min_tol

    return int_tol


def gapfind(model, solver, epsilon=0.001, v_max=1000):
    """Identify compounds in the model that cannot be produced.

    Yields all compounds that cannot be produced. This method
    assumes implicit sinks for all compounds in the model so
    the only factor that influences whether a compound can be
    produced is the presence of the compounds needed to produce it.

    Epsilon indicates the threshold amount of reaction flux for the products
    to be considered non-blocked. V_max indicates the maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.
    """
    prob = solver.create_problem()

    # Set integrality tolerance such that w constraints are correct
    min_tol = prob.integrality_tolerance.min
    int_tol = _find_integer_tolerance(epsilon, v_max, min_tol)
    if int_tol < prob.integrality_tolerance.value:
        prob.integrality_tolerance.value = int_tol

    # Define flux variables
    v = prob.namespace()
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        v.define([reaction_id], lower=lower, upper=upper)

    # Define constraints on production of metabolites in reaction
    w = prob.namespace(types=lp.VariableType.Binary)
    binary_cons_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        if value != 0:
            w.define([spec])
            w_var = w(spec)

            lower, upper = (float(x) for x in model.limits[reaction_id])
            if value > 0:
                dv = v(reaction_id)
            else:
                dv = -v(reaction_id)
                lower, upper = -upper, -lower

            prob.add_linear_constraints(
                dv <= upper * w_var,
                dv >= epsilon + (lower - epsilon) * (1 - w_var))

            binary_cons_lhs[compound] += w_var

    xp = prob.namespace(model.compounds, types=lp.VariableType.Binary)
    objective = xp.sum(model.compounds)
    prob.set_objective(objective)

    for compound, lhs in iteritems(binary_cons_lhs):
        prob.add_linear_constraints(lhs >= xp(compound))

    # Define mass balance constraints
    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        massbalance_lhs[compound] += v(reaction_id) * value
    for compound, lhs in iteritems(massbalance_lhs):
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise GapFillError('Non-optimal solution: {}'.format(result.status))

    for compound in model.compounds:
        if result.get_value(xp(compound)) < 0.5:
            yield compound


def gapfill(
        model, core, blocked, exclude, solver, epsilon=0.001, v_max=1000,
        weights={}):
    """Find a set of reactions to add such that no compounds are blocked.

    Returns two iterators: first an iterator of reactions not in
    core, that were added to resolve the model. Second, an
    iterator of reactions in core that had flux bounds expanded (i.e.
    irreversible reactions become reversible). Similarly to
    GapFind, this method assumes implicit sinks for all compounds in
    the model so the only factor that influences whether a compound
    can be produced is the presence of the compounds needed to produce
    it. This means that the resulting model will not necessarily be
    flux consistent.

    Core indicates the core set of reactions in the model. GapFill will
    minimize the number of added reactions that are not in core. Blocked
    indicates the set of compounds to be resolved. Exclude is a set of
    reactions that cannot be changed used for gap-filling. Epsilon indicates
    the threshold amount of a compound produced for it to not be considered
    blocked. V_max indicates the maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.
    """
    prob = solver.create_problem()

    # Set integrality tolerance such that w constraints are correct
    min_tol = prob.integrality_tolerance.min
    int_tol = _find_integer_tolerance(epsilon, v_max, min_tol)
    if int_tol < prob.integrality_tolerance.value:
        prob.integrality_tolerance.value = int_tol

    # Define flux variables
    v = prob.namespace(model.reactions, lower=-v_max, upper=v_max)

    # Add binary indicator variables
    database_reactions = set(model.reactions).difference(core, exclude)
    ym = prob.namespace(model.reactions, types=lp.VariableType.Binary)
    yd = prob.namespace(database_reactions, types=lp.VariableType.Binary)

    objective = ym.expr(
        (rxnid, weights.get(rxnid, 1)) for rxnid in model.reactions)
    objective += yd.expr(
        (rxnid, weights.get(rxnid, 1)) for rxnid in database_reactions)
    prob.set_objective(objective)

    # Add constraints on all reactions
    for reaction_id in model.reactions:
        lower, upper = (float(x) for x in model.limits[reaction_id])

        if reaction_id in exclude:
            prob.add_linear_constraints(
                upper >= v(reaction_id), v(reaction_id) >= lower)
        else:
            # Allow flux bounds to expand up to v_max with penalty
            delta_lower = min(0, -v_max - lower)
            delta_upper = max(0, v_max - upper)
            prob.add_linear_constraints(
                v(reaction_id) >= lower + ym(reaction_id) * delta_lower,
                v(reaction_id) <= upper + ym(reaction_id) * delta_upper)

    # Add constraints on database reactions
    for reaction_id in database_reactions:
        lower, upper = model.limits[reaction_id]
        prob.add_linear_constraints(
            v(reaction_id) >= yd(reaction_id) * -v_max,
            v(reaction_id) <= yd(reaction_id) * v_max)

    # Define constraints on production of blocked metabolites in reaction
    w = prob.namespace(types=lp.VariableType.Binary)
    binary_cons_lhs = {compound: 0 for compound in blocked}
    for (compound, reaction_id), value in iteritems(model.matrix):
        if reaction_id not in exclude and compound in blocked and value != 0:
            w.define([(compound, reaction_id)])
            w_var = w((compound, reaction_id))

            dv = v(reaction_id) if value > 0 else -v(reaction_id)
            prob.add_linear_constraints(
                dv <= v_max * w_var,
                dv >= epsilon + (-v_max - epsilon) * (1 - w_var))

            binary_cons_lhs[compound] += w_var

    for compound, lhs in iteritems(binary_cons_lhs):
        prob.add_linear_constraints(lhs >= 1)

    # Define mass balance constraints
    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for (compound, reaction_id), value in iteritems(model.matrix):
        if reaction_id not in exclude:
            massbalance_lhs[compound] += v(reaction_id) * value
    for compound, lhs in iteritems(massbalance_lhs):
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Minimize)
    if not result:
        raise GapFillError('Non-optimal solution: {}'.format(result.status))

    def added_iter():
        for reaction_id in database_reactions:
            if yd.value(reaction_id) > 0.5:
                yield reaction_id

    def no_bounds_iter():
        for reaction_id in model.reactions:
            if ym.value(reaction_id) > 0.5:
                yield reaction_id

    return added_iter(), no_bounds_iter()
