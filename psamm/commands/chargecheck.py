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

import math
import logging

from ..command import Command, FilePrefixAppendAction
from ..balancecheck import charge_balance

logger = logging.getLogger(__name__)


class ChargeBalanceCommand(Command):
    """Check whether compound charge is balanced.

    Balanced reactions are those reactions where the total charge
    is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', action=FilePrefixAppendAction,
            type=str, default=[], help='Exclude reaction from balance check')
        parser.add_argument(
            '--epsilon', metavar='epsilon', type=float, default=1e-6,
            help='Threshold for charge imbalance to be considered zero'
        )
        super(ChargeBalanceCommand, cls).init_parser(parser)

    def run(self):
        """Run charge balance command"""

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)

        # Create a set of excluded reactions
        exclude = set(self._args.exclude)
        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction, charge in charge_balance(self._model):
            count += 1

            if reaction.id in exclude or reaction.equation is None:
                continue

            if math.isnan(charge):
                logger.debug('Not checking reaction {};'
                             ' missing charge'.format(reaction.id))
                unchecked += 1
            elif abs(charge) > self._args.epsilon:
                unbalanced += 1
                rxt = reaction.equation.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction.id, charge, rxt))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing charge: {}/{}'.format(
            unchecked, count))
        logger.info('Reactions excluded from check: {}/{}'.format(
            len(exclude), count))
