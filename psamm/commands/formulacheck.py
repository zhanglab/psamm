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

import operator
import logging

from ..command import Command, FilePrefixAppendAction
from ..formula import Formula

logger = logging.getLogger(__name__)


class FormulaBalanceCommand(Command):
    """Check whether reactions in the model are elementally balanced.

    Balanced reactions are those reactions where the number of elements
    (atoms) is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', action=FilePrefixAppendAction,
            type=str, default=[], help='Exclude reaction from balance check')
        super(FormulaBalanceCommand, cls).init_parser(parser)

    def run(self):
        """Run formula balance command"""

        # Mapping from compound id to formula
        compound_formula = {}
        for compound in self._model.parse_compounds():
            if compound.formula is not None:
                try:
                    f = Formula.parse(compound.formula).flattened()
                    compound_formula[compound.id] = f
                except ValueError:
                    logger.warning(
                        'Error parsing formula for compound {}: {}'.format(
                            compound.id, compound.formula), exc_info=True)

        # Create a set of excluded reactions
        exclude = set(self._args.exclude)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula[compound.name]

        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction in self._model.parse_reactions():
            count += 1

            if reaction in exclude or reaction.equation is None:
                continue

            # Skip reaction if any compounds have undefined formula
            eq = reaction.equation
            for compound, _ in eq.compounds:
                if compound.name not in compound_formula:
                    unchecked += 1
                    break
            else:
                left_form = reduce(
                    operator.or_, multiply_formula(eq.left), Formula())
                right_form = reduce(
                    operator.or_, multiply_formula(eq.right), Formula())

                if right_form != left_form:
                    unbalanced += 1
                    right_missing, left_missing = Formula.balance(
                        right_form, left_form)

                    print('{}\t{}\t{}\t{}\t{}'.format(
                        reaction.id, left_form, right_form,
                        left_missing, right_missing))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing formula: {}/{}'.format(
            unchecked, count))
        logger.info('Reactions excluded from check: {}/{}'.format(
            len(exclude), count))
