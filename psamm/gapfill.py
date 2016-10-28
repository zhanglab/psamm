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


def gapfind(model, solver, epsilon=0.001, v_max=1000):
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
    v = prob.namespace()
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        v.define([reaction_id], lower=lower, upper=upper)

    # Define constraints on production of metabolites in reaction
    w = prob.namespace(types=lp.VariableType.Binary)
    binary_cons_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, reaction_id = spec
        if value != 0 and (reaction_id in model.reversible or value > 0):
            w.define([spec])
            w_var = w(spec)
            sv = v(reaction_id) * float(value)

            prob.add_linear_constraints(sv <= v_max * w_var)
            if model.is_reversible(reaction_id):
                prob.add_linear_constraints(
                    sv >= epsilon - v_max * (1 - w_var))
            else:
                prob.add_linear_constraints(sv >= epsilon * w_var)

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


def gapfill(model, core, blocked, solver, epsilon=0.001, v_max=1000):
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
    v = prob.namespace(model.reactions, lower=-v_max, upper=v_max)

    # Add binary indicator variables
    database_reactions = set(model.reactions).difference(core)
    ym = prob.namespace(core, types=lp.VariableType.Binary)
    yd = prob.namespace(database_reactions, types=lp.VariableType.Binary)

    objective = ym.sum(core) + yd.sum(database_reactions)
    prob.set_objective(objective)

    # Add constraints on core reactions
    for r in core:
        if model.is_reversible(r):
            prob.add_linear_constraints(v(r) >= model.limits[r].lower)
        else:
            prob.add_linear_constraints(v(r) >= -v_max * ym(r))
        prob.add_linear_constraints(v(r) <= model.limits[r].upper)

    # Add constraints on database reactions
    for r in database_reactions:
        prob.add_linear_constraints(v(r) >= yd(r) * model.limits[r].lower)
        prob.add_linear_constraints(v(r) <= yd(r) * model.limits[r].upper)

    # Define constraints on production of blocked metabolites in reaction
    w = prob.namespace(types=lp.VariableType.Binary)
    yn = prob.namespace(types=lp.VariableType.Binary)
    binary_cons_lhs = {compound: 0 for compound in blocked}
    for (compound, reaction_id), value in iteritems(model.matrix):
        if compound in blocked and value != 0:
            w.define([(compound, reaction_id)])
            w_var = w((compound, reaction_id))

            sv = float(value) * v(reaction_id)
            prob.add_linear_constraints(
                sv >= epsilon - v_max * (1 - w_var))
            prob.add_linear_constraints(sv <= v_max * w_var)

            if reaction_id in model.reversible or value > 0:
                binary_cons_lhs[compound] += w_var
            elif reaction_id in core:
                # In this case, we need to perform a logical AND on the w and
                # ym variables. This is done by introducing another helper
                # variable, yn.
                yn.define([(reaction_id, compound)])
                yn_var = yn((reaction_id, compound))
                prob.add_linear_constraints(
                    2 * yn_var <= w_var + ym(reaction_id))
                binary_cons_lhs[compound] += yn_var

    for compound, lhs in iteritems(binary_cons_lhs):
        if compound in blocked:
            prob.add_linear_constraints(lhs >= 1)

    # Define mass balance constraints
    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for (compound, reaction_id), value in iteritems(model.matrix):
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

    def reversed_iter():
        for reaction_id in core:
            if ym.value(reaction_id) > 0.5:
                yield reaction_id

    return added_iter(), reversed_iter()
