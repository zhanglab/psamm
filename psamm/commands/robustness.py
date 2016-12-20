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

from ..command import (Command, MetabolicMixin, LoopRemovalMixin,
                       ObjectiveMixin, SolverCommandMixin, ParallelTaskMixin)
from .. import fluxanalysis

from six.moves import range

logger = logging.getLogger(__name__)


class RobustnessCommand(MetabolicMixin, LoopRemovalMixin, ObjectiveMixin,
                        SolverCommandMixin, ParallelTaskMixin, Command):
    """Run robustness analysis on the model.

    Given a reaction to maximize and a reaction to vary,
    the robustness analysis will run FBA while fixing the
    reaction to vary at each iteration. The reaction will
    be fixed at the specified number of steps between the
    minimum and maximum flux value specified in the model.
    """

    _supported_loop_removal = ['none', 'tfba', 'l1min']

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--steps', metavar='N', type=int, default=10,
            help='Number of flux value steps for varying reaction')
        parser.add_argument(
            '--minimum', metavar='V', type=float,
            help='Minumum flux value of varying reacton')
        parser.add_argument(
            '--maximum', metavar='V', type=float,
            help='Maximum flux value of varying reacton')
        parser.add_argument(
            '--all-reaction-fluxes',
            help='Print reaction flux for all model reactions',
            action='store_true')
        parser.add_argument('varying', help='Reaction to vary')
        super(RobustnessCommand, cls).init_parser(parser)

    def run(self):
        """Run robustness command."""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        reaction = self._get_objective()
        if not self._mm.has_reaction(reaction):
            self.fail('Specified biomass reaction is not in model: {}'.format(
                reaction))

        varying_reaction = self._args.varying
        if not self._mm.has_reaction(varying_reaction):
            self.fail('Specified varying reaction is not in model: {}'.format(
                varying_reaction))

        steps = self._args.steps
        if steps <= 0:
            self.argument_error('Invalid number of steps: {}\n'.format(steps))

        loop_removal = self._get_loop_removal_option()
        if loop_removal == 'tfba':
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)
        if loop_removal == 'tfba':
            p.add_thermodynamic()

        try:
            p.check_constraints()
        except fluxanalysis.FluxBalanceError as e:
            self.report_flux_balance_error(e)

        # Determine minimum and maximum flux for varying reaction
        if self._args.maximum is None:
            p.maximize(varying_reaction)
            flux_max = p.get_flux(varying_reaction)
        else:
            flux_max = self._args.maximum

        if self._args.minimum is None:
            p.maximize({varying_reaction: -1})
            flux_min = p.get_flux(varying_reaction)
        else:
            flux_min = self._args.minimum

        if flux_min > flux_max:
            self.argument_error('Invalid flux range: {}, {}\n'.format(
                flux_min, flux_max))

        logger.info('Varying {} in {} steps between {} and {}'.format(
            varying_reaction, steps, flux_min, flux_max))

        start_time = time.time()

        handler_args = (
            self._mm, solver, loop_removal, self._args.all_reaction_fluxes)
        executor = self._create_executor(
            RobustnessTaskHandler, handler_args, cpus_per_worker=2)

        def iter_tasks():
            for i in range(steps):
                fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
                constraint = varying_reaction, fixed_flux
                yield constraint, reaction

        # Run FBA on model at different fixed flux values
        with executor:
            for task, result in executor.imap_unordered(iter_tasks(), 16):
                (varying_reaction, fixed_flux), _ = task
                if result is None:
                    logger.warning('No solution found for {} at {}'.format(
                        varying_reaction, fixed_flux))
                elif self._args.all_reaction_fluxes:
                    for other_reaction in self._mm.reactions:
                        print('{}\t{}\t{}'.format(
                            other_reaction, fixed_flux,
                            result[other_reaction]))
                else:
                    print('{}\t{}'.format(fixed_flux, result))

        executor.join()

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))


class RobustnessTaskHandler(object):
    def __init__(self, model, solver, loop_removal, all_reactions):
        self._problem = fluxanalysis.FluxBalanceProblem(model, solver)

        if loop_removal == 'none':
            self._run_fba = self._problem.maximize
        elif loop_removal == 'l1min':
            self._run_fba = self._problem.max_min_l1
        elif loop_removal == 'tfba':
            self._problem.add_thermodynamic()
            self._run_fba = self._problem.maximize

        self._reactions = list(model.reactions) if all_reactions else None

    def handle_task(self, constraint, reaction):
        varying_reaction, fixed_flux = constraint
        flux_var = self._problem.get_flux_var(varying_reaction)
        c, = self._problem.prob.add_linear_constraints(flux_var == fixed_flux)

        try:
            self._run_fba(reaction)
            if self._reactions is not None:
                return {reaction: self._problem.get_flux(reaction)
                        for reaction in self._reactions}
            else:
                return self._problem.get_flux(reaction)
        except fluxanalysis.FluxBalanceError:
            return None
        finally:
            c.delete()
