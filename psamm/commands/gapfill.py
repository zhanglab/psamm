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

from ..command import Command, SolverCommandMixin
from ..gapfill import gapfind, gapfill

logger = logging.getLogger(__name__)


class GapFillCommand(SolverCommandMixin, Command):
    """Run the GapFind and GapFill algorithms on the model."""

    def run(self):
        """Run GapFill command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

        model_compartments = set(self._mm.compartments)

        solver = self._get_solver(integer=True)

        # Run GapFind on model
        logger.info('Searching for blocked compounds')
        blocked = set(compound for compound in gapfind(self._mm, solver=solver)
                      if compound.compartment is not 'e')
        if len(blocked) > 0:
            logger.info('Blocked compounds')
            for compound in blocked:
                print(compound.translate(lambda x: compound_name.get(x, x)))

        if len(blocked) > 0:
            # Add exchange and transport reactions to database
            model_complete = self._mm.copy()
            logger.info('Adding database, exchange and transport reactions')
            model_complete.add_all_database_reactions(model_compartments)
            model_complete.add_all_exchange_reactions()
            model_complete.add_all_transport_reactions()

            logger.info('Searching for reactions to fill gaps')
            added_reactions, reversed_reactions = gapfill(
                model_complete, self._mm.reactions, blocked, solver=solver)

            for rxnid in added_reactions:
                rx = model_complete.get_reaction(rxnid)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(rxnid, 'Add', rxt))

            for rxnid in reversed_reactions:
                rx = model_complete.get_reaction(rxnid)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(rxnid, 'Reverse', rxt))
        else:
            logger.info('No blocked compounds found')
