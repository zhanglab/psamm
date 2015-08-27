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

import operator
import logging

from ..command import Command
from ..formula import Formula

logger = logging.getLogger(__name__)


class FormulaBalanceCommand(Command):
    """Check whether reactions in the model are elementally balanced.

    Balanced reactions are those reactions where the number of elements
    (atoms) is consistent on the left and right side of the reaction equation.
    Reactions that are not balanced will be printed out.
    """

    def run(self):
        """Run formula balance command"""

        # Mapping from compound id to formula
        compound_formula = {}
        for compound in self._model.parse_compounds():
            # Only cpd11632 (Photon) is allowed to have an empty formula.
            if compound.formula is None:
                if compound.id == 'cpd11632':
                    compound_formula[compound.id] = Formula()
            else:
                try:
                    f = Formula.parse(compound.formula).flattened()
                    compound_formula[compound.id] = f
                except ValueError as e:
                    logger.warning(
                        'Error parsing {}: {}'.format(compound.formula, e))

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula.get(compound.name, Formula())

        count = 0
        unbalanced = 0
        unchecked = 0
        for reaction in self._mm.reactions:
            if reaction in exchange:
                continue

            count += 1
            rx = self._mm.get_reaction(reaction)

            # Skip reaction if any compounds have undefined formula
            for compound, _ in rx.compounds:
                if compound.name not in compound_formula:
                    unchecked += 1
                    break
            else:
                left_form = reduce(
                    operator.or_, multiply_formula(rx.left), Formula())
                right_form = reduce(
                    operator.or_, multiply_formula(rx.right), Formula())

                if right_form != left_form:
                    unbalanced += 1
                    right_missing, left_missing = Formula.balance(
                        right_form, left_form)

                    print('{}\t{}\t{}\t{}\t{}'.format(
                        reaction, left_form, right_form,
                        left_missing, right_missing))

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing formula: {}/{}'.format(
            unchecked, count))
