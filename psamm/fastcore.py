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

from .fluxanalysis import FluxBalanceProblem

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


class FastcoreProblem(FluxBalanceProblem):
    def __init__(self, *args, **kwargs):
        self._epsilon = kwargs.pop('epsilon', 1e-5)
        super(FastcoreProblem, self).__init__(*args, **kwargs)
        self._has_maximization_vars = False
        self._flipped = set()

    def _add_maximization_vars(self):
        if self._has_maximization_vars:
            return

        # Define maximization variables
        self._prob.define(*(('zl', rxnid) for rxnid in self._model.reactions),
                          lower=0, upper=self._epsilon)

        self._has_maximization_vars = True

    def lp7(self, reaction_subset):
        """Approximately maximize the number of reaction with flux.

        This is similar to FBA but approximately maximizing the number of
        reactions in subset with flux > epsilon, instead of just maximizing the
        flux of one particular reaction. LP7 prefers "flux splitting" over
        "flux concentrating".
        """

        self._add_maximization_vars()

        positive = set(reaction_subset) - self._flipped
        negative = set(reaction_subset) & self._flipped

        v = self._prob.set(('v', rxnid) for rxnid in positive)
        zl = self._prob.set(('zl', rxnid) for rxnid in positive)
        cs = self._prob.add_linear_constraints(v >= zl)
        self._temp_constr.extend(cs)

        v = self._prob.set(('v', rxnid) for rxnid in negative)
        zl = self._prob.set(('zl', rxnid) for rxnid in negative)
        cs = self._prob.add_linear_constraints(v <= -zl)
        self._temp_constr.extend(cs)

        self._prob.set_linear_objective(self._prob.expr(
            {('zl', rxnid): 1 for rxnid in reaction_subset}))

        self._solve()

    def lp10(self, subset_k, subset_p, weights={}):
        """Force reactions in K above epsilon while minimizing support of P.

        This program forces reactions in subset K to attain flux > epsilon
        while minimizing the sum of absolute flux values for reactions
        in subset P (L1-regularization).
        """

        self._add_minimization_vars()

        positive = set(subset_k) - self._flipped
        negative = set(subset_k) & self._flipped

        v = self._prob.set(('v', rxnid) for rxnid in positive)
        cs = self._prob.add_linear_constraints(v >= self._epsilon)
        self._temp_constr.extend(cs)

        v = self._prob.set(('v', rxnid) for rxnid in negative)
        cs = self._prob.add_linear_constraints(v <= -self._epsilon)
        self._temp_constr.extend(cs)

        self._prob.set_linear_objective(self._prob.expr(
            {('z', rxnid): -weights.get(rxnid, 1) for rxnid in subset_p}))

        self._solve()

    def find_sparse_mode(self, core, additional, scaling, weights={}):
        """Find a sparse mode containing reactions of the core subset.

        Return an iterator of the support of a sparse mode that contains as
        many reactions from core as possible, and as few reactions from
        additional as possible (approximately). A dictionary of weights can be
        supplied which gives further penalties for including specific
        additional reactions.
        """

        if len(core) == 0:
            return

        self.lp7(core)
        k = set()
        for reaction_id in core:
            flux = self.get_flux(reaction_id)
            if self.is_flipped(reaction_id):
                flux *= -1
            if flux >= self._epsilon:
                k.add(reaction_id)

        if len(k) == 0:
            return

        self.lp10(k, additional, weights)
        for reaction_id in self._model.reactions:
            flux = self.get_flux(reaction_id)
            if abs(flux) >= self._epsilon / scaling:
                yield reaction_id

    def flip(self, reactions):
        for reaction in reactions:
            if reaction in self._flipped:
                self._flipped.remove(reaction)
            else:
                self._flipped.add(reaction)

    def is_flipped(self, reaction):
        return reaction in self._flipped


def fastcc(model, epsilon, solver):
    """Check consistency of model reactions

    Yields all reactions in the model that are not part
    of the consistent subset.
    """

    reaction_set = set(model.reactions)
    subset = set(reaction_id for reaction_id in reaction_set
                 if model.limits[reaction_id].lower >= 0)

    logger.info('Checking {} irreversible reactions...'.format(len(subset)))
    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    p = FastcoreProblem(model, solver, epsilon=epsilon)
    p.lp7(subset)

    consistent_subset = set(
        reaction_id for reaction_id in model.reactions
        if abs(p.get_flux(reaction_id)) >= 0.999 * epsilon)

    logger.debug('|A| = {}, A = {}'.format(
        len(consistent_subset), consistent_subset))

    for reaction in subset - consistent_subset:
        # Inconsistent reaction
        yield reaction

    # Check remaining reactions
    subset = (reaction_set - subset) - consistent_subset

    logger.info('Checking reversible reactions...')
    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    flipped = False
    singleton = False
    while len(subset) > 0:
        logger.info('{} reversible reactions left to check...'.format(
            len(subset)))
        if singleton:
            reaction = next(iter(subset))
            subset_i = {reaction}

            logger.debug('LP3 on {}'.format(subset_i))
            p.maximize({reaction: -1 if p.is_flipped(reaction) else 1})
        else:
            subset_i = subset

            logger.debug('LP7 on {}'.format(subset_i))
            p.lp7(subset_i)

        consistent_subset.update(
            reaction_id for reaction_id in subset
            if abs(p.get_flux(reaction_id) >= 0.999 * epsilon))

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
                        logger.info('Inconsistent: {}'.format(reaction))
                        yield reaction
                else:
                    singleton = True
            else:
                p.flip(subset_rev_i)
                flipped = True
                logger.info('Flipped {} reactions'.format(len(subset_rev_i)))


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


def fastcore(model, core, epsilon, solver, scaling=1e5, weights={}):
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

    p = FastcoreProblem(model, solver, epsilon=epsilon)
    mode = set(p.find_sparse_mode(subset, penalty_set, scaling, weights))
    if not subset.issubset(mode):
        raise FastcoreError('Inconsistent irreversible core reactions:'
                            ' {}'.format(subset - mode))

    consistent_subset |= mode
    logger.debug('|A| = {}, A = {}'.format(
        len(consistent_subset), consistent_subset))

    subset = core - mode
    logger.debug('|J| = {}, J = {}'.format(len(subset), subset))

    flipped = False
    singleton = False
    while len(subset) > 0:
        penalty_set -= consistent_subset
        if singleton:
            subset_i = set((next(iter(subset)),))
        else:
            subset_i = subset

        mode = set(p.find_sparse_mode(subset_i, penalty_set, scaling, weights))
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
                p.flip(subset_rev_i)
                flipped = True
                logger.debug('Flipped {} reactions'.format(
                    len(subset_rev_i)))

    return consistent_subset
