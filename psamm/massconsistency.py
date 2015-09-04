# -*- coding: utf-8 -*-
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

"""Mass consistency analysis of metabolic databases

A stoichiometric matrix, S, is said to be mass-consistent if
S^Tm = 0 has a positive solution (m_i > 0). This corresponds to
assigning a positive mass to each compound in the stoichiometric
matrix and having each reaction preserve mass. Exchange reactions
will have to be excluded from this check, as they are not able to
preserve mass (by definition). In addition some databases may contain
pseudo-compounds (e.g. "photon") that also has to be excluded.
"""

from .lpsolver import lp

from six import iteritems


class MassConsistencyError(Exception):
    """Indicates an error while checking for mass consistency"""


def _non_localized_compounds(database):
    """Return set of non-localized compounds in database"""
    return set(c.in_compartment(None) for c in database.compounds)


def is_consistent(database, solver, exchange=set(), zeromass=set()):
    """Try to assign a positive mass to each compound

    Return True if successful. The masses are simply constrained by m_i > 1 and
    finding a solution under these conditions proves that the database is mass
    consistent.
    """

    prob = solver.create_problem()
    compound_set = _non_localized_compounds(database)

    # Define mass variables
    for compound in compound_set:
        if compound.name not in zeromass:
            prob.define(('m', compound), lower=1)
        else:
            prob.define(('m', compound), lower=0, upper=0)

    prob.set_linear_objective(prob.expr(
        {('m', compound): 1 for compound in compound_set}))

    # Define constraints
    massbalance_lhs = {reaction: 0 for reaction in database.reactions}
    for spec, value in iteritems(database.matrix):
        compound, reaction = spec
        mass = prob.var(('m', compound.in_compartment(None)))
        massbalance_lhs[reaction] += mass * value
    for reaction, lhs in iteritems(massbalance_lhs):
        if reaction not in exchange:
            prob.add_linear_constraints(lhs == 0)

    result = prob.solve(lp.ObjectiveSense.Minimize)
    return result.success


def check_reaction_consistency(database, solver, exchange=set(),
                               checked=set(), zeromass=set(), weights={}):
    """Check inconsistent reactions by minimizing mass residuals

    Return a reaction iterable, and compound iterable. The reaction iterable
    yields reaction ids and mass residuals. The compound iterable yields
    compound ids and mass assignments.

    Each compound is assigned a mass of at least one, and the masses are
    balanced using the stoichiometric matrix. In addition, each reaction has a
    residual mass that is included in the mass balance equations. The L1-norm
    of the residuals is minimized. Reactions in the checked set are assumed to
    have been manually checked and therefore have the residual fixed at zero.
    """

    # Create Flux balance problem
    prob = solver.create_problem()
    compound_set = _non_localized_compounds(database)

    # Define mass variables
    for compound in compound_set:
        if compound.name not in zeromass:
            prob.define(('m', compound), lower=1)
        else:
            prob.define(('m', compound), lower=0, upper=0)

    # Define residual mass variables and objective constriants
    prob.define(*(('z', reaction_id) for reaction_id in database.reactions),
                lower=0)
    prob.define(*(('r', reaction_id) for reaction_id in database.reactions))

    objective = prob.expr(
        {('z', reaction_id): weights.get(reaction_id, 1)
         for reaction_id in database.reactions})
    prob.set_linear_objective(objective)

    r = prob.set(('r', reaction_id) for reaction_id in database.reactions)
    z = prob.set(('z', reaction_id) for reaction_id in database.reactions)
    prob.add_linear_constraints(z >= r, r >= -z)

    massbalance_lhs = {reaction_id: 0 for reaction_id in database.reactions}
    for spec, value in iteritems(database.matrix):
        compound, reaction_id = spec
        mass_var = prob.var(('m', compound.in_compartment(None)))
        massbalance_lhs[reaction_id] += value * mass_var
    for reaction_id, lhs in iteritems(massbalance_lhs):
        if reaction_id not in exchange:
            if reaction_id not in checked:
                residual = prob.var(('r', reaction_id))
                prob.add_linear_constraints(lhs + residual == 0)
            else:
                prob.add_linear_constraints(lhs == 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Minimize)
    if not result:
        raise MassConsistencyError('Non-optimal solution: {}'.format(
            result.status))

    def iterate_reactions():
        for reaction_id in database.reactions:
            residual = result.get_value(('r', reaction_id))
            yield reaction_id, residual

    def iterate_compounds():
        for compound in compound_set:
            yield compound, result.get_value(('m', compound))

    return iterate_reactions(), iterate_compounds()


def check_compound_consistency(database, solver, exchange=set(),
                               zeromass=set()):
    """Yield each compound in the database with assigned mass

    Each compound will be assigned a mass and the number of compounds having a
    positive mass will be approximately maximized.

    This is an implementation of the solution originally proposed by
    [Gevorgyan08]_  but using the new method proposed by [Thiele14]_ to avoid
    MILP constraints. This is similar to the way Fastcore avoids MILP
    contraints.
    """

    # Create mass balance problem
    prob = solver.create_problem()
    compound_set = _non_localized_compounds(database)

    # Define mass variables
    for compound in compound_set:
        if compound.name not in zeromass:
            prob.define(('m', compound), lower=0)
        else:
            prob.define(('m', compound), lower=0, upper=0)

    # Define z variables
    prob.define(*(('z', compound) for compound in compound_set),
                lower=0, upper=1)
    prob.set_linear_objective(prob.expr(
        {('z', compound): 1 for compound in compound_set}))

    z = prob.set(('z', compound) for compound in compound_set)
    m = prob.set(('m', compound) for compound in compound_set)
    prob.add_linear_constraints(m >= z)

    massbalance_lhs = {reaction_id: 0 for reaction_id in database.reactions}
    for spec, value in iteritems(database.matrix):
        compound, reaction_id = spec
        mass_var = prob.var(('m', compound.in_compartment(None)))
        massbalance_lhs[reaction_id] += value * mass_var
    for reaction_id, lhs in iteritems(massbalance_lhs):
        if reaction_id not in exchange:
            prob.add_linear_constraints(lhs == 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise MassConsistencyError('Non-optimal solution: {}'.format(
            result.status))

    for compound in compound_set:
        yield compound, result.get_value(('m', compound))
