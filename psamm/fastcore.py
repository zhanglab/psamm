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

"""Fastcore module implementing the fastcore algorithm

This is an implementation of the algorithms described in [Vlassis14]_.
"""

import logging

from .lpsolver import lp
from .fluxanalysis import flux_balance
from .metabolicmodel import FlipableModelView

from six import iteritems

# Module-level logging
logger = logging.getLogger(__name__)


def support(fluxiter, threshold=None):
    """Yield reactions in the support set of the fluxes

    The fluxiter argument must be an iterable returning reaction, flux-pairs.
    The support set is the set of reactions that have an absolute flux value at
    or above a certain threshold. If threshold is None the mathematical
    definition (v != 0) is used.
    """

    if threshold is None:
        return (rxnid for rxnid, v in fluxiter if v != 0)
    return (rxnid for rxnid, v in fluxiter if abs(v) >= threshold)


def support_positive(fluxiter, threshold=None):
    """Yield reactions in the support set of the fluxes (only positive fluxes)

    The fluxiter argument must be an iterable returning reaction, flux-pairs.
    The support set is the set of reactions that have a flux value at or above
    a certain threshold. If threshold is None the mathematical definition
    (v > 0) is used.
    """

    if threshold is None:
        return (rxnid for rxnid, v in fluxiter if v > 0)
    return (rxnid for rxnid, v in fluxiter if v >= threshold)


class FastcoreError(Exception):
    """Indicates an error while running Fastcore"""


def lp7(model, reaction_subset, epsilon, solver):
    """Approximately maximize the number of reaction with flux above epsilon

    This is similar to FBA but approximately maximizing the number of reactions
    in subset with flux > epsilon, instead of just maximizing the flux of one
    particular reaction. LP7 prefers "flux splitting" over
    "flux concentrating".
    """

    # Create LP-7 problem of Fastcore
    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reactions:
        lower, upper = model.limits[rxnid]
        prob.define(('v', rxnid), lower=lower, upper=upper)

    # Define z variables
    prob.define(*(('z', rxnid) for rxnid in reaction_subset),
                lower=0, upper=epsilon)
    prob.set_linear_objective(prob.expr(
        {('z', rxnid): 1 for rxnid in reaction_subset}))
    v = prob.set(('v', rxnid) for rxnid in reaction_subset)
    z = prob.set(('z', rxnid) for rxnid in reaction_subset)
    prob.add_linear_constraints(v >= z)

    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, rxnid = spec
        massbalance_lhs[compound] += prob.var(('v', rxnid)) * value
    prob.add_linear_constraints(
        *(lhs == 0 for compound, lhs in iteritems(massbalance_lhs)))

    # Solve
    result = prob.solve(lp.ObjectiveSense.Maximize)
    if not result:
        raise FastcoreError('Non-optimal solution: {}'.format(result.status))

    for rxnid in sorted(model.reactions):
        yield rxnid, result.get_value(('v', rxnid))


def lp10(model, subset_k, subset_p, epsilon, scaling, solver, weights={}):
    """Force reactions in K above epsilon while minimizing support of P

    This program forces reactions in subset K to attain flux > epsilon
    while minimizing the sum of absolute flux values for reactions
    in subset P (L1-regularization).
    """

    if len(subset_k) == 0:
        return

    # Create LP-10 problem of Fastcore
    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reactions:
        lower, upper = model.limits[rxnid]
        if rxnid in subset_k:
            lower = max(lower, epsilon)
        prob.define(('v', rxnid), lower=lower*scaling, upper=upper*scaling)

    # Define z variables
    prob.define(*(('z', rxnid) for rxnid in subset_p), lower=0)
    prob.set_linear_objective(prob.expr(
        {('z', rxnid): weights.get(rxnid, 1) for rxnid in subset_p}))

    z = prob.set(('z', rxnid) for rxnid in subset_p)
    v = prob.set(('v', rxnid) for rxnid in subset_p)
    prob.add_linear_constraints(z >= v, v >= -z)

    massbalance_lhs = {compound: 0 for compound in model.compounds}
    for spec, value in iteritems(model.matrix):
        compound, rxnid = spec
        massbalance_lhs[compound] += prob.var(('v', rxnid)) * value
    prob.add_linear_constraints(
        *(lhs == 0 for compound, lhs in iteritems(massbalance_lhs)))

    # Solve
    result = prob.solve(lp.ObjectiveSense.Minimize)
    if not result:
        raise FastcoreError('Non-optimal solution: {}'.format(result.status))

    for reaction_id in model.reactions:
        yield reaction_id, result.get_value(('v', reaction_id))


def fastcc(model, epsilon, solver):
    """Check consistency of model reactions

    Yields all reactions in the model that are not part
    of the consistent subset.
    """

    reaction_set = set(model.reactions)
    subset = reaction_set.difference(model.reversible)

    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    consistent_subset = set(
        support(lp7(model, subset, epsilon, solver), epsilon))

    logger.debug('|A| = {}, A = {}'.format(
        len(consistent_subset), consistent_subset))

    for reaction in subset - consistent_subset:
        # Inconsistent reaction
        yield reaction

    # Check remaining reactions
    subset = (reaction_set - subset) - consistent_subset

    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    # Wrap model in flipable proxy so reactions can be flipped
    model = FlipableModelView(model)

    flipped = False
    singleton = False
    while len(subset) > 0:
        if singleton:
            reaction = next(iter(subset))
            subset_i = {reaction}

            logger.debug('LP3 on {}'.format(subset_i))
            supp = support(flux_balance(
                model, reaction, tfba=False, solver=solver), epsilon)
        else:
            subset_i = subset

            logger.debug('LP7 on {}'.format(subset_i))
            supp = support(lp7(model, subset_i, epsilon, solver), epsilon)
        consistent_subset.update(supp)

        logger.debug('|A| = {}, A = {}'.format(
            len(consistent_subset), consistent_subset))

        if not subset.isdisjoint(consistent_subset):
            subset -= consistent_subset
            logger.debug('|J| = {}, J = {}'.format(len(subset), subset))
            flipped = False
        else:
            # TODO: irreversible reactions are taken care of before the
            # loop so at this point all reactions in subset_i are reversble(?).
            subset_rev_i = subset_i & model.reversible
            if flipped or len(subset_rev_i) == 0:
                flipped = False
                if singleton:
                    subset -= subset_rev_i
                    for reaction in subset_rev_i:
                        yield reaction
                else:
                    singleton = True
            else:
                model.flip(subset_rev_i)
                flipped = True
                logger.debug('Flip')


def fastcc_is_consistent(model, epsilon, solver):
    """Quickly check whether model is consistent

    Returns true if the model is consistent. If it is only necessary to know
    whether a model is consistent, this function is faster as it will return
    the result as soon as it finds a single inconsistent reaction.
    """

    for reaction in fastcc(model, epsilon, solver):
        return False
    return True


def fastcc_consistent_subset(model, epsilon, solver):
    """Return consistent subset of model

    The largest consistent subset is returned as
    a set of reaction names.
    """

    reaction_set = set(model.reactions)
    return reaction_set.difference(fastcc(model, epsilon, solver))


def find_sparse_mode(model, core, additional, epsilon, scaling, solver,
                     weights={}):
    """Find a sparse mode containing reactions of the core subset

    Return an iterator of the support of a sparse mode that contains as many
    reactions from core as possible, and as few reactions from additional as
    possible (approximately). A dictionary of weights can be supplied which
    gives further penalties for including specific additional reactions.
    """

    if len(core) == 0:
        return iter(())

    supp = support_positive(lp7(model, core, epsilon, solver), epsilon)
    k = core.intersection(supp)
    if len(k) == 0:
        return iter(())

    return support(lp10(model, k, additional, epsilon, solver=solver,
                        scaling=scaling, weights=weights), epsilon)


def fastcore(model, core, epsilon, solver, scaling=1e8, weights={}):
    """Find a flux consistent subnetwork containing the core subset

    The result will contain the core subset and as few of the additional
    reactions as possible.
    """

    consistent_subset = set()
    reaction_set = set(model.reactions)

    subset = core - model.reversible
    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    penalty_set = reaction_set - core
    logger.debug('|P| = {}, P = {}'.format(len(penalty_set), penalty_set))

    mode = set(find_sparse_mode(
        model, subset, penalty_set, epsilon, solver=solver, scaling=scaling,
        weights=weights))
    if not subset.issubset(mode):
        raise FastcoreError('Inconsistent irreversible core reactions:'
                            ' {}'.format(subset - mode))

    consistent_subset |= mode
    logger.debug('|A| = {}, A = {}'.format(
        len(consistent_subset), consistent_subset))

    subset = core - mode
    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    # Wrap model in flipable proxy so reactions can be flipped
    model = FlipableModelView(model)

    flipped = False
    singleton = False
    while len(subset) > 0:
        penalty_set -= consistent_subset
        if singleton:
            subset_i = set((next(iter(subset)),))
        else:
            subset_i = subset

        mode = find_sparse_mode(
            model, subset_i, penalty_set, epsilon, scaling=scaling,
            solver=solver, weights=weights)
        consistent_subset.update(mode)
        logger.debug('|A| = {}, A = {}'.format(
            len(consistent_subset), consistent_subset))

        if not subset.isdisjoint(consistent_subset):
            logger.debug('Subset improved {} -> {}'.format(
                len(subset), len(subset - consistent_subset)))
            subset -= consistent_subset
            logger.debug('|J| = {}, J = {}'.format(len(subset), subset))
            flipped = False
        else:
            logger.debug('Nothing found, changing state...')
            subset_rev_i = subset_i & model.reversible
            if flipped or len(subset_rev_i) == 0:
                if singleton:
                    raise FastcoreError('Global network inconsistent:'
                                        ' {}'.format(subset_rev_i))

                logger.debug('Going to non-flipped, singleton state...')
                singleton = True
                flipped = False
            else:
                model.flip(subset_rev_i)
                flipped = True
                logger.debug('Going to flipped state... {}'.format(
                    subset_rev_i))

    return consistent_subset
