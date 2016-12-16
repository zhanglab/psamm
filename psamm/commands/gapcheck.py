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
# Copyright 2014-2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

from __future__ import unicode_literals

import logging

from ..command import Command, MetabolicMixin, SolverCommandMixin
from ..gapfill import gapfind, GapFillError

logger = logging.getLogger(__name__)


class GapCheckCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Check for compound production gaps in model."""
    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--epsilon', type=float, default=1e-5,
            help='Threshold for compound production')
        super(GapCheckCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        solver = self._get_solver(integer=True)
        extracellular_comp = self._model.extracellular_compartment
        epsilon = self._args.epsilon
        v_max = float(self._model.default_flux_limit)

        # Run GapFind on model
        logger.info('Searching for blocked compounds...')
        result = gapfind(self._mm, solver=solver, epsilon=epsilon, v_max=v_max)

        try:
            blocked = set(compound for compound in result
                          if compound.compartment != extracellular_comp)
        except GapFillError as e:
            self._log_epsilon_and_fail(epsilon, e)

        if len(blocked) > 0:
            logger.info('Blocked compounds: {}'.format(len(blocked)))
            for compound in sorted(blocked):
                name = compound_name.get(compound.name, compound.name)
                print('{}\t{}'.format(compound, name))
        else:
            logger.info('No blocked compounds found')

    def _log_epsilon_and_fail(self, epsilon, exc):
        msg = ('Finding blocked compounds failed with epsilon set to {}. Try'
               ' lowering the epsilon value to reduce artifical constraints'
               ' on the model.'.format(epsilon))
        self.fail(msg, exc)
