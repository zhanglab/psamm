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
# Copyright 2016  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016  Chao Liu <lcddzyx@gmail.com>

from __future__ import unicode_literals

import operator
import logging

from six.moves import reduce

from .formula import Formula

logger = logging.getLogger(__name__)


def reaction_charge(reaction, compound_charge):
    """Calculate the overall charge for the specified reaction.

    Args:
        reaction: :class:`psamm.reaction.Reaction`.
        compound_charge: a map from each compound to charge values.
    """

    charge_sum = 0.0
    for compound, value in reaction.compounds:
        charge = compound_charge.get(compound.name, float('nan'))
        charge_sum += charge * float(value)
    return charge_sum


def charge_balance(model):
    """Calculate the overall charge for all reactions in the model.

    Yield (reaction, charge) pairs.

    Args:
        model: :class:`psamm.datasource.native.NativeModel`.
    """

    compound_charge = {}
    for compound in model.parse_compounds():
        if hasattr(compound, 'charge') and compound.charge is not None:
            compound_charge[compound.id] = compound.charge

    for reaction in model.parse_reactions():
        charge = reaction_charge(reaction.equation, compound_charge)
        yield reaction, charge


def reaction_formula(reaction, compound_formula):
    """Calculate formula compositions for both sides of the specified reaction.

    If the compounds in the reaction all have formula, then calculate and
    return the chemical compositions for both sides, otherwise return `None`.

    Args:
        reaction: :class:`psamm.reaction.Reaction`.
        compound_formula: a map from compound id to formula.
    """

    def multiply_formula(compound_list):
        for compound, count in compound_list:
            yield count * compound_formula[compound.name]

    for compound, _ in reaction.compounds:
        if compound.name not in compound_formula:
            return None
    else:
        left_form = reduce(
            operator.or_, multiply_formula(reaction.left), Formula())
        right_form = reduce(
            operator.or_, multiply_formula(reaction.right), Formula())
    return left_form, right_form


def formula_balance(model):
    """Calculate formula compositions for each reaction.

    Call :func:`reaction_formula` for each reaction.
    Yield (reaction, result) pairs, where result has two formula compositions
    or `None`.

    Args:
        model: :class:`psamm.datasource.native.NativeModel`.
    """

    # Mapping from compound id to formula
    compound_formula = {}
    for compound in model.parse_compounds():
        if compound.formula is not None:
            try:
                f = Formula.parse(compound.formula).flattened()
                compound_formula[compound.id] = f
            except ValueError:
                logger.warning(
                    'Error parsing formula for compound {}: {}'.format(
                        compound.id, compound.formula), exc_info=True)

    for reaction in model.parse_reactions():
        yield reaction, reaction_formula(reaction.equation, compound_formula)
