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
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2016  Chao Liu <lcddzyx@gmail.com>

from __future__ import unicode_literals

import random
import logging

from . import fluxanalysis
from .expression import boolean

from six import iteritems, string_types

logger = logging.getLogger(__name__)


class ReactionDeletionStrategy(object):
    """Deleting reactions strategy class.

    When initializing instances of this class,
    :func:`.get_exchange_reactions` can be useful if exchange reactions
    are used as the test set.
    """

    def __init__(self, model, reaction_set=None):
        self._model = model
        if reaction_set is None:
            self._test_set = set(self._model.reactions)
        else:
            self._test_set = set(reaction_set)
        self._reactions = set(self._test_set)

    @property
    def entities(self):
        return self._reactions

    def iter_tests(self):
        while len(self._test_set) > 0:
            reaction = random.choice(tuple(self._test_set))
            self._test_set.remove(reaction)
            yield reaction, {reaction}

    def delete(self, entity, deleted_reactions):
        pass


class GeneDeletionStrategy(object):
    """Deleting genes strategy class.

    When initializing instances of this class, :func:`get_gene_associations`
    can be called to obtain the gene association dict from the model.
    """

    def __init__(self, model, gene_assoc):
        self._model = model
        self._genes = set()
        self._gene_assoc = dict(gene_assoc)

        for reaction, assoc in iteritems(self._gene_assoc):
            self._genes.update(v.symbol for v in assoc.variables)

        self._reactions = set(self._model.reactions)
        self._test_set = set(self._genes)

    @property
    def entities(self):
        return set(self._genes)

    def iter_tests(self):
        while len(self._test_set) > 0:
            gene = random.choice(tuple(self._test_set))
            self._test_set.remove(gene)

            deleted_reactions = set()
            self._new_gene_assoc = {}
            for reaction in self._reactions:
                if reaction not in self._gene_assoc:
                    continue
                assoc = self._gene_assoc[reaction]
                if boolean.Variable(gene) in assoc.variables:
                    new_assoc = assoc.substitute(
                        lambda v: v if v.symbol != gene else False)
                    if new_assoc.has_value() and not new_assoc.value:
                        deleted_reactions.add(reaction)
                    else:
                        self._new_gene_assoc[reaction] = new_assoc
                else:
                    self._new_gene_assoc[reaction] = assoc

            yield gene, deleted_reactions

    def delete(self, entity, deleted_reactions):
        self._reactions -= deleted_reactions
        self._gene_assoc = self._new_gene_assoc


def get_gene_associations(model):
    """Create gene association for class :class:`.GeneDeletionStrategy`.

    Return a dict mapping reaction IDs to
    :class:`psamm.expression.boolean.Expression` objects,
    representing relationships between reactions and related genes. This helper
    function should be called when creating :class:`.GeneDeletionStrategy`
    objects.

    Args:
        model: :class:`psamm.datasource.native.NativeModel`.
    """

    for reaction in model.parse_reactions():
        assoc = None
        if reaction.genes is None:
            continue
        elif isinstance(reaction.genes, string_types):
            assoc = boolean.Expression(reaction.genes)
        else:
            variables = [boolean.Variable(g) for g in reaction.genes]
            assoc = boolean.Expression(boolean.And(*variables))
        yield reaction.id, assoc


def get_exchange_reactions(model):
    """Yield IDs of all exchange reactions from model.

    This helper function would be useful when creating
    :class:`.ReactionDeletionStrategy` objects.

    Args:
        model: :class:`psamm.metabolicmodel.MetabolicModel`.
    """
    for reaction_id in model.reactions:
        if model.is_exchange(reaction_id):
            yield reaction_id


def random_sparse(strategy, prob, obj_reaction, flux_threshold):
    """Find a random minimal network of model reactions.

    Given a reaction to optimize and a threshold, delete entities randomly
    until the flux of the reaction to optimize falls under the threshold.
    Keep deleting until no more entities can be deleted. It works
    with two strategies: deleting reactions or deleting genes (reactions
    related to certain genes).

    Args:
        strategy: :class:`.ReactionDeletionStrategy` or
            :class:`.GeneDeletionStrategy`.
        prob: :class:`psamm.fluxanalysis.FluxBalanceProblem`.
        obj_reaction: objective reactions to optimize.
        flux_threshold: threshold of max reaction flux.
    """

    essential = set()
    deleted = set()
    for entity, deleted_reactions in strategy.iter_tests():
        if obj_reaction in deleted_reactions:
            logger.info(
                'Marking entity {} as essential because the objective'
                ' reaction depends on this entity...'.format(entity))
            essential.add(entity)
            continue

        if len(deleted_reactions) == 0:
            logger.info(
                'No reactions were removed when entity {}'
                ' was deleted'.format(entity))
            deleted.add(entity)
            strategy.delete(entity, deleted_reactions)
            continue

        logger.info('Deleted reactions: {}'.format(
            ', '.join(deleted_reactions)))

        constr = []
        for r in deleted_reactions:
            flux_var = prob.get_flux_var(r)
            c, = prob.prob.add_linear_constraints(flux_var == 0)
            constr.append(c)

        logger.info('Trying FBA without reactions {}...'.format(
            ', '.join(deleted_reactions)))

        try:
            prob.maximize(obj_reaction)
        except fluxanalysis.FluxBalanceError:
            logger.info(
                'FBA is infeasible, marking {} as essential'.format(
                    entity))
            for c in constr:
                c.delete()
            essential.add(entity)
            continue

        logger.debug('Reaction {} has flux {}'.format(
            obj_reaction, prob.get_flux(obj_reaction)))

        if prob.get_flux(obj_reaction) < flux_threshold:
            for c in constr:
                c.delete()
            essential.add(entity)
            logger.info('Entity {} was essential'.format(
                entity))
        else:
            deleted.add(entity)
            strategy.delete(entity, deleted_reactions)
            logger.info('Entity {} was deleted'.format(entity))

    return essential, deleted
