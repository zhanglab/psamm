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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import logging

from ..command import Command, FilePrefixAppendAction, convert_to_unicode
from ..formula import Formula
from ..balancecheck import formula_balance

logger = logging.getLogger(__name__)

def formula_check(self):
    # Create a set of excluded reactions
    if hasattr(self, '_args'):
        exclude = set(self._args.exclude)
    else:
        exclude = set()
    count = 0
    unbalanced = 0
    unchecked = 0
    form_list = []
    form_dict = {}
    for reaction, result in formula_balance(self._model):
        count += 1

        if reaction.id in exclude or reaction.equation is None:
            continue

        if result is None:
            unchecked += 1
            continue

        left_form, right_form = result
        if right_form != left_form:
            unbalanced += 1
            right_missing, left_missing = Formula.balance(
                right_form, left_form)
            form_list.append(reaction.id)
            form_dict[reaction.id]=(left_form, left_missing, right_form, right_missing)

            print('{}\t{}\t{}\t{}\t{}'.format(
                reaction.id, left_form, right_form,
                left_missing, right_missing))
    return(unbalanced, count, unchecked, exclude, form_list, form_dict)

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
            type=convert_to_unicode, default=[],
            help='Exclude reaction from balance check')
        super(FormulaBalanceCommand, cls).init_parser(parser)

    def run(self):
        """Run formula balance command"""
        unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(self)

        logger.info('Unbalanced reactions: {}/{}'.format(unbalanced, count))
        logger.info('Unchecked reactions due to missing formula: {}/{}'.format(
            unchecked, count))
        logger.info('Reactions excluded from check: {}/{}'.format(
            len(exclude), count))
