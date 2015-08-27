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

import math
import logging

from ..command import Command

logger = logging.getLogger(__name__)


class ChargeBalanceCommand(Command):
    """Check whether compound charge is balanced.

    Balanced reactions are those reactions where the total charge
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    def run(self):
        """Run charge balance command"""

        # Load compound information
        compound_name = {}
        compound_charge = {}
        for compound in self._model.parse_compounds():
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)
            if hasattr(compound, 'charge') and compound.charge is not None:
                compound_charge[compound.id] = compound.charge

        # Create a set of known charge-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def reaction_charges(reaction_id):
            for compound, value in self._mm.get_reaction_values(reaction_id):
                charge = compound_charge.get(compound.name, float('nan'))
                yield charge * float(value)

        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction in sorted(self._mm.reactions):
            if reaction in exchange:
                continue

            count += 1
            charge = sum(reaction_charges(reaction))
            if math.isnan(charge):
                logger.debug('Not checking reaction {};'
                             ' missing charge'.format(reaction))
                unchecked += 1
            elif charge != 0:
                unbalanced += 1
                rx = self._mm.get_reaction(reaction)
                rxt = rx.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction, charge, rxt))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing charge: {}/{}'.format(
            unchecked, count))
