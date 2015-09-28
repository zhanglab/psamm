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

import time
import logging

from ..command import Command, SolverCommandMixin, CommandError
from .. import fluxanalysis

logger = logging.getLogger(__name__)


class FluxVariabilityCommand(SolverCommandMixin, Command):
    """Run flux variablity analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--tfba', help='Enable thermodynamic constraints on FVA',
            action='store_true')
        parser.add_argument('reaction', help='Reaction to maximize', nargs='?')
        super(FluxVariabilityCommand, cls).init_parser(parser)

    def run(self):
        """Run flux variability command"""

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

        enable_tfba = self._args.tfba
        if enable_tfba:
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        start_time = time.time()

        fba_fluxes = dict(fluxanalysis.flux_balance(
            self._mm, reaction, tfba=False, solver=solver))
        optimum = fba_fluxes[reaction]

        flux_bounds = fluxanalysis.flux_variability(
            self._mm, sorted(self._mm.reactions), {reaction: optimum},
            tfba=enable_tfba, solver=solver)
        for reaction_id, bounds in flux_bounds:
            rx = self._mm.get_reaction(reaction_id)
            rxt = rx.translated_compounds(lambda x: compound_name.get(x, x))
            print('{}\t{}\t{}\t{}'.format(
                reaction_id, bounds[0], bounds[1], rxt))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))
