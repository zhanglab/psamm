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

from six import iteritems

from ..command import SolverCommandMixin, Command, CommandError
from .. import fluxanalysis

logger = logging.getLogger(__name__)


class FluxBalanceCommand(SolverCommandMixin, Command):
    """Run flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--no-tfba', help='Disable thermodynamic constraints on FBA',
            action='store_true')
        parser.add_argument(
            '--epsilon', type=float, help='Threshold for flux minimization',
            default=1e-5)
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

        if self._args.reaction is not None:
            reaction = self._args.reaction
        else:
            reaction = self._model.get_biomass_reaction()
            if reaction is None:
                raise CommandError('The biomass reaction was not specified')

        if not self._mm.has_reaction(reaction):
            raise CommandError('Specified reaction is not in model: {}'.format(
                reaction))

        if self._args.no_tfba:
            result = self.run_fba_minimized(reaction)
        else:
            result = self.run_tfba(reaction)

        optimum = None
        for reaction_id, fba_flux, flux in sorted(result):
            rx = self._mm.get_reaction(reaction_id)
            print('{}\t{}\t{}\t{}'.format(
                reaction_id, fba_flux, flux,
                rx.translated_compounds(lambda x: compound_name.get(x, x))))
            # Remember flux of requested reaction
            if reaction_id == reaction:
                optimum = flux

        logger.info('Maximum flux: {}'.format(optimum))

    def run_fba_minimized(self, reaction):
        """Run normal FBA and flux minimization on model."""

        epsilon = self._args.epsilon

        solver = self._get_solver()
        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        # Maximize reaction flux
        p.maximize(reaction)
        fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

        # Run flux minimization
        flux_var = p.get_flux_var(reaction)
        p.prob.add_linear_constraints(flux_var == p.get_flux(reaction))
        p.minimize_l1()
        count = 0
        for reaction_id in self._mm.reactions:
            flux = p.get_flux(reaction_id)
            if abs(flux - fluxes[reaction_id]) > epsilon:
                count += 1
            yield reaction_id, fluxes[reaction_id], flux
        logger.info('Minimized reactions: {}'.format(count))

    def run_tfba(self, reaction):
        """Run FBA and tFBA on model."""

        solver = self._get_solver(integer=True)
        p = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        p.maximize(reaction)
        fba_fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

        p.add_thermodynamic()
        p.maximize(reaction)
        fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

        for reaction_id, flux in iteritems(fluxes):
            yield reaction_id, fba_fluxes[reaction_id], flux
