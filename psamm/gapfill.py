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

from six import iteritems

from .lpsolver import lp


class GapFillError(Exception):
    """Indicates an error while running GapFind/GapFill"""


def gapfind(model, solver, epsilon=1e-5, v_max=1000):
    """Identify compounds in the model that cannot be produced.

    Yields all compounds that cannot be produced. This method
    assumes implicit sinks for all compounds in the model so
    the only factor that influences whether a compound can be
    produced is the presence of the compounds needed to produce it.

    Epsilon indicates the threshold amount of a compound produced
    for it to not be considered blocked. V_max indicates the
    maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.
    """
    prob = solver.create_problem()

    # Define flux variables
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        prob.define(('v', reaction_id), lower=lower, upper=upper)

    # Define constraints on production of metabolites in reaction
    binary_cons_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        if value != 0 and (reaction_id in model.reversible or value > 0):
            prob.define(('w', reaction_id, compound),
                        types=lp.VariableType.Binary)

            w = prob.var(('w', reaction_id, compound))
            sv = float(value) * prob.var(('v', reaction_id))

            prob.add_linear_constraints(sv <= v_max*w)
            if model.is_reversible(reaction_id):
                prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            else:
                prob.add_linear_constraints(sv >= epsilon*w)

            binary_cons_lhs[compound] += w

    prob.define(*(('xp', compound) for compound in model.compounds),
                types=lp.VariableType.Binary)

    objective = prob.expr(
        {('xp', compound): 1 for compound in model.compounds})
    prob.set_linear_objective(objective)

    for compound, lhs in iteritems(binary_cons_lhs):
        prob.add_linear_constraints(lhs >= prob.var(('xp', compound)))

    # Define mass balance constraints
    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        massbalance_lhs[compound] += prob.var(('v', reaction_id)) * value
    for compound, lhs in iteritems(massbalance_lhs):
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise GapFillError('Non-optimal solution: {}'.format(result.status))

    for compound in model.compounds:
        if result.get_value(('xp', compound)) == 0:
            yield compound


def gapfill(model, core, blocked, solver, epsilon=1e-5, v_max=1000):
    """Find a set of reactions to add such that no compounds are blocked.

    Returns two iterators: first an iterator of reactions not in
    core, that were added to resolve the model. Second, an
    iterator of reaction in core, that were reversed. Similarly to
    GapFind, this method assumes implicit sinks for all compounds in
    the model so the only factor that influences whether a compound
    can be produced is the presence of the compounds needed to produce
    it. This means that the resulting model will not necessarily be
    flux consistent.

    Core indicates the core set of reactions in the model. GapFill will
    minimize the number of added reactions that are not in core. Blocked
    indicates the set of compounds to be resolved. Epsilon indicates the
    threshold amount of a compound produced for it to not be considered
    blocked. V_max indicates the maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.
    """
    prob = solver.create_problem()

    # Define flux variables
    prob.define(*(('v', reaction_id) for reaction_id in model.reactions),
                lower=-v_max, upper=v_max)

    # Add binary indicator variables
    database_reactions = set(model.reactions).difference(core)
    prob.define(*(('ym', reaction_id) for reaction_id in core),
                types=lp.VariableType.Binary)
    prob.define(*(('yd', reaction_id) for reaction_id in database_reactions),
                types=lp.VariableType.Binary)

    objective = prob.expr(
        {('ym', reaction_id): 1 for reaction_id in core})
    objective += prob.expr(
        {('yd', reaction_id): 1 for reaction_id in database_reactions})
    prob.set_linear_objective(objective)

    # Add constraints on core reactions
    for reaction_id in core:
        v = prob.var(('v', reaction_id))
        if model.is_reversible(reaction_id):
            prob.add_linear_constraints(v >= model.limits[reaction_id].lower)
        else:
            ym = prob.var(('ym', reaction_id))
            prob.add_linear_constraints(v >= -v_max*ym)
        prob.add_linear_constraints(v <= model.limits[reaction_id].upper)

    # Add constraints on database reactions
    for reaction_id in database_reactions:
        v = prob.var(('v', reaction_id))
        yd = prob.var(('yd', reaction_id))
        prob.add_linear_constraints(v >= yd * model.limits[reaction_id].lower)
        prob.add_linear_constraints(v <= yd * model.limits[reaction_id].upper)

    # Define constraints on production of blocked metabolites in reaction
    binary_cons_lhs = {compound: 0 for compound in blocked}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        if compound in blocked and value != 0:
            prob.define(('w', reaction_id, compound),
                        types=lp.VariableType.Binary)

            w = prob.var(('w', reaction_id, compound))
            sv = float(value) * prob.var(('v', reaction_id))

            prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            prob.add_linear_constraints(sv <= v_max*w)

            if reaction_id in model.reversible or value > 0:
                binary_cons_lhs[compound] += w

    for compound, lhs in iteritems(binary_cons_lhs):
        if compound in blocked:
            prob.add_linear_constraints(lhs >= 1)

    # Define mass balance constraints
    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        massbalance_lhs[compound] += prob.var(('v', reaction_id)) * value
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
            if result.get_value(('yd', reaction_id)) > 0:
                yield reaction_id

    def reversed_iter():
        for reaction_id in core:
            if result.get_value(('ym', reaction_id)) > 0:
                yield reaction_id

    return added_iter(), reversed_iter()
