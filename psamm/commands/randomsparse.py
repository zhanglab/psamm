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

import logging

from ..command import (Command, MetabolicMixin, LoopRemovalMixin,
                       ObjectiveMixin, SolverCommandMixin)
from .. import fluxanalysis, util
from .. import randomsparse

logger = logging.getLogger(__name__)


class RandomSparseNetworkCommand(MetabolicMixin, LoopRemovalMixin,
                                 ObjectiveMixin, SolverCommandMixin, Command):
    """Find a random minimal network of model reactions.

    Given a reaction to optimize and a threshold, delete reactions randomly
    until the flux of the reaction to optimize falls under the threshold.
    Keep deleting reactions until no more reactions can be deleted. By default
    this uses standard FBA (not tFBA). Since the internal fluxes are irrelevant
    the FBA and tFBA are equivalent for this purpose.

    The threshold can be specified as an absolute flux (e.g. '1.23') or a
    relative flux of the full model flux (e.g. '40.5%').
    """

    _supported_loop_removal = ['none', 'tfba']

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            'threshold', help='Threshold of max reaction flux',
            type=util.MaybeRelative)
        parser.add_argument(
            '--type', help='Type of deletion to perform',
            choices=['reactions', 'exchange', 'genes'], type=str,
            required=True)
        super(RandomSparseNetworkCommand, cls).init_parser(parser)

    def run(self):
        reaction = self._get_objective()
        if not self._mm.has_reaction(reaction):
            self.fail(
                'Specified reaction is not in model: {}'.format(reaction))

        loop_removal = self._get_loop_removal_option()
        enable_tfba = loop_removal == 'tfba'
        if not enable_tfba:
            solver = self._get_solver()
        else:
            solver = self._get_solver(integer=True)

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        if enable_tfba:
            p.add_thermodynamic()

        try:
            p.maximize(reaction)
        except fluxanalysis.FluxBalanceError as e:
            self.report_flux_balance_error(e)

        threshold = self._args.threshold
        if threshold.relative:
            threshold.reference = p.get_flux(reaction)

        flux_threshold = float(threshold)

        logger.info('Flux threshold for {} is {}'.format(
            reaction, flux_threshold))

        if self._args.type == 'exchange':
            strategy = randomsparse.ReactionDeletionStrategy(
                self._mm, randomsparse.get_exchange_reactions(self._mm))
            entity = 'reactions'
        elif self._args.type == 'reactions':
            strategy = randomsparse.ReactionDeletionStrategy(
                self._mm)
            entity = 'reactions'
        else:
            strategy = randomsparse.GeneDeletionStrategy(
                self._mm, randomsparse.get_gene_associations(self._model))
            entity = 'genes'

        essential, deleted = randomsparse.random_sparse(
            strategy, p, reaction, flux_threshold)

        logger.info('Essential {}: {}/{}'.format(
            entity, len(essential), len(strategy.entities)))
        logger.info('Deleted {}: {}/{}'.format(
            entity, len(deleted), len(strategy.entities)))

        for entity_id in sorted(strategy.entities):
            value = 0 if entity_id in deleted else 1
            print('{}\t{}'.format(entity_id, value))
