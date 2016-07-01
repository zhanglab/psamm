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

from ..command import SolverCommandMixin, MetabolicMixin, Command
from .. import fluxanalysis

logger = logging.getLogger(__name__)


class FluxBalanceCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Run flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument(
            '--all-reactions', help='Show all reaction fluxes',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float,
            help='Threshold for non-zero reaction fluxes', default=1e-5)
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        super(FluxBalanceCommand, cls).init_parser(parser)

    def run(self):
        """Run flux analysis command."""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        # Reaction genes information
        reaction_genes = {}
        for reaction in self._model.parse_reactions():
            if 'genes' in reaction.properties:
                reaction_genes[reaction.id] = reaction.properties['genes']

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                self.argument_error('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            self.fail(
                'Specified reaction is not in model: {}'.format(reaction))

        logger.info('Using {} as objective'.format(reaction))

        loop_removal = self._args.loop_removal
        if loop_removal == 'none':
            logger.info('Loop removal disabled; spurious loops are allowed')
            result = self.run_fba(reaction)
        elif loop_removal == 'l1min':
            logger.info('Loop removal using L1 minimization')
            result = self.run_fba_minimized(reaction)
        elif loop_removal == 'tfba':
            logger.info('Loop removal using thermodynamic constraints')
            result = self.run_tfba(reaction)
        else:
            self.argument_error(
                'Invalid loop constraint mode: {}'.format(loop_removal))

        optimum = None
        total_reactions = 0
        nonzero_reactions = 0
        for reaction_id, flux in sorted(result):
            total_reactions += 1
            if abs(flux) > self._args.epsilon:
                nonzero_reactions += 1

            if abs(flux) > self._args.epsilon or self._args.all_reactions:
                rx = self._mm.get_reaction(reaction_id)
                rx_trans = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                genes = reaction_genes.get(reaction_id, '')
                print('{}\t{}\t{}\t{}'.format(
                    reaction_id, flux, rx_trans, genes))

            # Remember flux of requested reaction
            if reaction_id == reaction:
                optimum = flux

        logger.info('Objective flux: {}'.format(optimum))
        logger.info('Reactions at zero flux: {}/{}'.format(
            total_reactions - nonzero_reactions, total_reactions))

    def run_fba(self, reaction):
        """Run standard FBA on model."""
        solver = self._get_solver()
        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        start_time = time.time()

        p.maximize(reaction)

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        for reaction_id in self._mm.reactions:
            yield reaction_id, p.get_flux(reaction_id)

    def run_fba_minimized(self, reaction):
        """Run normal FBA and flux minimization on model."""

        epsilon = self._args.epsilon
        solver = self._get_solver()

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        start_time = time.time()

        # Maximize reaction flux
        p.maximize(reaction)
        fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

        # Run flux minimization
        flux_var = p.get_flux_var(reaction)
        p.prob.add_linear_constraints(flux_var == p.get_flux(reaction))
        p.minimize_l1()

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        count = 0
        for reaction_id in self._mm.reactions:
            flux = p.get_flux(reaction_id)
            if abs(flux - fluxes[reaction_id]) > epsilon:
                count += 1
            yield reaction_id, flux
        logger.info('Minimized reactions: {}'.format(count))

    def run_tfba(self, reaction):
        """Run FBA and tFBA on model."""

        solver = self._get_solver(integer=True)
        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        start_time = time.time()

        p.add_thermodynamic()
        p.maximize(reaction)

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        for reaction_id in self._mm.reactions:
            yield reaction_id, p.get_flux(reaction_id)
