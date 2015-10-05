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

import time
import random
import logging

from ..command import Command, SolverCommandMixin, CommandError
from .. import fluxanalysis, util

logger = logging.getLogger(__name__)


class RandomSparseNetworkCommand(SolverCommandMixin, Command):
    """Find a random minimal network of model reactions.

    Given a reaction to optimize and a threshold, delete reactions randomly
    until the flux of the reaction to optimize falls under the threshold.
    Keep deleting reactions until no more reactions can be deleted. By default
    this uses standard FBA (not tFBA). Since the internal fluxes are irrelevant
    the FBA and tFBA are equivalent for this purpose.

    The threshold can be specified as an absolute flux (e.g. '1.23') or a
    relative flux of the full model flux (e.g. '40.5%').
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--tfba', help='Enable thermodynamic constraints on FBA',
            action='store_true', default=False)
        parser.add_argument('--objective', help='Reaction flux to maximize')
        parser.add_argument(
            'threshold', help='Threshold of max reaction flux',
            type=util.MaybeRelative)
        parser.add_argument(
            '--exchange', help='Only analyze the exchange reaction in model',
            action='store_true')
        super(RandomSparseNetworkCommand, cls).init_parser(parser)

    def run(self):
        if self._args.objective is not None:
            reaction = self._args.objective
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        if not self._args.tfba:
            solver = self._get_solver()
        else:
            solver = self._get_solver(integer=True)

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        if self._args.tfba:
            p.add_thermodynamic()

        threshold = self._args.threshold
        if threshold.relative:
            p.maximize(reaction)
            threshold.reference = p.get_flux(reaction)

        flux_threshold = float(threshold)

        logger.info('Flux threshold for {} is {}'.format(
            reaction, flux_threshold))

        essential = {reaction}
        deleted = set()
        if self._args.exchange:
            reactions = set()
            for reaction_id in self._mm.reactions:
                if self._mm.is_exchange(reaction_id):
                    reactions.add(reaction_id)
        else:
            reactions = set(self._mm.reactions)

        test_set = reactions - essential

        start_time = time.time()

        while len(test_set) > 0:
            testing_reaction = random.sample(test_set, 1)[0]
            test_set.remove(testing_reaction)

            flux_var = p.get_flux_var(testing_reaction)
            c, = p.prob.add_linear_constraints(flux_var == 0)

            logger.info('Trying FBA without reaction {}...'.format(
                testing_reaction))

            try:
                p.maximize(reaction)
            except fluxanalysis.FluxBalanceError:
                logger.info(
                    'FBA is infeasible, marking {} as essential'.format(
                        testing_reaction))
                c.delete()
                essential.add(testing_reaction)
                continue

            logger.debug('Reaction {} has flux {}'.format(
                reaction, p.get_flux(reaction)))

            if p.get_flux(reaction) < flux_threshold:
                c.delete()
                essential.add(testing_reaction)
                logger.info('Reaction {} was essential'.format(
                    testing_reaction))
            else:
                deleted.add(testing_reaction)
                logger.info('Reaction {} was deleted'.format(testing_reaction))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        for reaction_id in sorted(reactions):
            value = 0 if reaction_id in deleted else 1
            print('{}\t{}'.format(reaction_id, value))
