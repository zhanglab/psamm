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

from __future__ import unicode_literals

import time
import logging
from itertools import product

from six import text_type

from ..command import (Command, MetabolicMixin, LoopRemovalMixin,
                       SolverCommandMixin, ParallelTaskMixin)
from .. import fluxanalysis, fastcore

logger = logging.getLogger(__name__)


class FluxConsistencyCommand(MetabolicMixin, LoopRemovalMixin,
                             SolverCommandMixin, ParallelTaskMixin, Command):
    """Check that reactions are flux consistent in a model.

    A reaction is flux consistent if there exists any steady-state flux
    solution where the flux of the given reaction is non-zero. The
    bounds on the exchange reactions can be removed when performing the
    consistency check by providing the ``--unrestricted`` option.
    """

    _supported_loop_removal = ['none', 'tfba']

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--fastcore', help='Enable use of Fastcore algorithm',
            action='store_true')
        parser.add_argument(
            '--reduce-lp',
            help='Try to reduce the number of LP problems to solve',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float, help='Flux threshold',
            default=1e-5)
        parser.add_argument(
            '--unrestricted', action='store_true',
            help='Remove limits on exchange reactions before checking')
        super(FluxConsistencyCommand, cls).init_parser(parser)

    def run(self):
        """Run flux consistency check command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        epsilon = self._args.epsilon

        if self._args.unrestricted:
            # Allow all exchange reactions with no flux limits
            for reaction in self._mm.reactions:
                if self._mm.is_exchange(reaction):
                    del self._mm.limits[reaction].bounds

        loop_removal = self._get_loop_removal_option()
        enable_tfba = loop_removal == 'tfba'
        enable_fastcore = self._args.fastcore

        if enable_tfba and enable_fastcore:
            self.argument_error(
                'Using Fastcore with thermodynamic constraints'
                ' is not supported!')
        start_time = time.time()

        if enable_fastcore:
            solver = self._get_solver()
            try:
                inconsistent = set(fastcore.fastcc(
                    self._mm, epsilon, solver=solver))
            except fluxanalysis.FluxBalanceError as e:
                self.report_flux_balance_error(e)
        else:
            if enable_tfba:
                solver = self._get_solver(integer=True)
            else:
                solver = self._get_solver()

            if self._args.reduce_lp:
                logger.info('Running with reduced number of LP problems.')
                try:
                    inconsistent = set(
                        fluxanalysis.consistency_check(
                            self._mm, self._mm.reactions, epsilon,
                            tfba=enable_tfba, solver=solver))
                except fluxanalysis.FluxBalanceError as e:
                    self.report_flux_balance_error(e)
            else:
                logger.info('Using flux bounds to determine consistency.')
                try:
                    inconsistent = set(self._run_fva_fluxcheck(
                        self._mm, solver, enable_tfba, epsilon))
                except FluxCheckFVATaskError:
                    self.report_flux_balance_error()

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        # Count the number of reactions that are fixed at zero. While these
        # reactions are still inconsistent, they are inconsistent because they
        # have been explicitly disabled.
        disabled_exchange = 0
        disabled_internal = 0

        count_exchange = 0
        total_exchange = 0

        count_internal = 0
        total_internal = 0

        # Print result
        for reaction in sorted(self._mm.reactions):
            disabled = self._mm.limits[reaction].bounds == (0, 0)

            if self._mm.is_exchange(reaction):
                total_exchange += 1
                count_exchange += int(reaction in inconsistent)
                disabled_exchange += int(disabled)
            else:
                total_internal += 1
                count_internal += int(reaction in inconsistent)
                disabled_internal += int(disabled)

            if reaction in inconsistent:
                rx = self._mm.get_reaction(reaction)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}'.format(reaction, rxt))

        logger.info('Model has {}/{} inconsistent internal reactions'
                    ' ({} disabled by user)'.format(
                        count_internal, total_internal, disabled_internal))
        logger.info('Model has {}/{} inconsistent exchange reactions'
                    ' ({} disabled by user)'.format(
                        count_exchange, total_exchange, disabled_exchange))

    def _run_fva_fluxcheck(self, model, solver, enable_tfba, epsilon):
        handler_args = model, solver, enable_tfba
        executor = self._create_executor(
            FluxCheckFVATaskHandler, handler_args, cpus_per_worker=2)

        results = {}
        with executor:
            for (reaction_id, direction), value in executor.imap_unordered(
                    product(model.reactions, (1, -1)), 16):

                if reaction_id not in results:
                    results[reaction_id] = value
                    continue

                other_value = results[reaction_id]
                if direction == -1:
                    bounds = value, other_value
                else:
                    bounds = other_value, value

                lower, upper = bounds
                if abs(lower) < epsilon and abs(upper) < epsilon:
                    yield reaction_id

        executor.join()


class FluxCheckFVATaskError(Exception):
    """Error raised from parallel flux check task on failure."""


class FluxCheckFVATaskHandler(object):
    def __init__(self, model, solver, enable_tfba):
        self._problem = fluxanalysis.FluxBalanceProblem(model, solver)
        if enable_tfba:
            self._problem.add_thermodynamic()

    def handle_task(self, reaction_id, direction):
        try:
            return self._problem.flux_bound(reaction_id, direction)
        except fluxanalysis.FluxBalanceError as e:
            # FluxBalanceError is not picklable. Reraise as picklable
            # exception.
            raise FluxCheckFVATaskError(text_type(e))
