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

import logging

from ..command import Command, SolverCommandMixin, CommandError
from .. import fluxanalysis

from six.moves import range

logger = logging.getLogger(__name__)


class RobustnessCommand(SolverCommandMixin, Command):
    """Run robustness analysis on the model.

    Given a reaction to maximize and a reaction to vary, the robustness
    analysis will run FBA while fixing the reaction to vary at each iteration.
    The reaction will be fixed at the specified number of steps between the
    minimum and maximum flux value.

    The loop removal option can be used to apply thermodynamic constraints or
    flux minimization. The fluxes option can be used to select which reaction
    fluxes are calculated and displayed.
    """

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
            '--loop-removal', help='Select type of loop removal constraints',
            choices=['none', 'tfba', 'l1min'], default='none')
        parser.add_argument(
            '--reaction', help='Reaction to maximize', nargs='?')
        parser.add_argument(
            '--fluxes', help='Which reaction fluxes to calculate',
            choices=['objective', 'all', 'range'], default='objective')
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

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        varying_reaction = self._args.varying
        if not self._mm.has_reaction(varying_reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                varying_reaction))

        steps = self._args.steps
        if steps <= 0:
            raise CommandError('Invalid number of steps: {}\n'.format(steps))

        loop_removal = self._args.loop_removal

        # Automatically change loop removal in cases where the choice of
        # flux output is incompatible.
        if self._args.fluxes == 'objective' and loop_removal != 'none':
            logger.warning(
                'Loop removal has no effect when only the objective is'
                ' calculated; setting loop removal to "none".')
            loop_removal = 'none'

        if self._args.fluxes == 'range' and loop_removal == 'l1min':
            logger.warning(
                'Flux minimization has no effect when calculating flux ranges;'
                ' setting loop removal to "none".')
            loop_removal = 'none'

        if loop_removal == 'tfba':
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        if loop_removal == 'none':
            logger.info('Loop removal disabled; spurious loops are allowed')
            run_fba = p.maximize
        elif loop_removal == 'l1min':
            logger.info('Loop removal using L1 minimization')
            run_fba = p.max_min_l1
        elif loop_removal == 'tfba':
            logger.info('Loop removal using thermodynamic constraints')
            p.add_thermodynamic()
            run_fba = p.maximize
        else:
            raise CommandError('Invalid loop constraint mode: {}'.format(
                loop_removal))

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
            raise CommandError('Invalid flux range: {}, {}\n'.format(
                flux_min, flux_max))

        logger.info('Varying {} in {} steps between {} and {}'.format(
            varying_reaction, steps, flux_min, flux_max))

        # Run FBA on model at different fixed flux values
        for i in range(steps):
            fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
            flux_var = p.get_flux_var(varying_reaction)
            c, = p.prob.add_linear_constraints(flux_var == fixed_flux)

            try:
                if self._args.fluxes != 'range':
                    self._run_flux_iteration(p, run_fba, reaction, fixed_flux)
                else:
                    self._run_range_iteration(p, reaction, fixed_flux)
            except fluxanalysis.FluxBalanceError:
                logger.warning('No solution found for {} at {}'.format(
                    varying_reaction, fixed_flux))
            finally:
                c.delete()

    def _run_flux_iteration(self, problem, function, reaction, fixed_flux):
        function(reaction)

        if self._args.fluxes == 'objective':
            print('{}\t{}'.format(fixed_flux, problem.get_flux(reaction)))
        else:
            for other_reaction in self._mm.reactions:
                print('{}\t{}\t{}'.format(
                    other_reaction, fixed_flux,
                    problem.get_flux(other_reaction)))

    def _run_range_iteration(self, problem, reaction, fixed_flux):
        problem.maximize(reaction)
        optimum = problem.get_flux(reaction)
        flux_var = problem.get_flux_var(reaction)
        c, = problem.prob.add_linear_constraints(flux_var == optimum)

        try:
            for range_reaction in self._mm.reactions:
                max_flux = None
                try:
                    problem.maximize(range_reaction)
                    max_flux = problem.get_flux(range_reaction)
                except fluxanalysis.FluxBalanceError:
                    logger.warning(
                        'No solution found when determining maximum flux'
                        ' for {}'.format(range_reaction))

                min_flux = None
                try:
                    problem.maximize({range_reaction: -1})
                    min_flux = problem.get_flux(range_reaction)
                except fluxanalysis.FluxBalanceError:
                    logger.warning(
                        'No solution found when determining minimum flux'
                        ' for {}'.format(range_reaction))

                print('{}\t{}\t{}\t{}'.format(
                    range_reaction, fixed_flux, min_flux, max_flux))
        finally:
            c.delete()
