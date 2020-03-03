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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import logging

from six import itervalues

from ..command import Command

logger = logging.getLogger(__name__)


def reaction_signature(eq, direction=False, stoichiometry=False):
    """Return unique signature object for :class:`Reaction`.

    Signature objects are hashable, and compare equal only if the reactions
    are considered the same according to the specified rules.

    Args:
        direction: Include reaction directionality when considering equality.
        stoichiometry: Include stoichiometry when considering equality.
    """
    def compounds_sig(compounds):
        if stoichiometry:
            return tuple(sorted(compounds))
        else:
            return tuple(sorted(compound for compound, _ in compounds))

    left = compounds_sig(eq.left)
    right = compounds_sig(eq.right)

    if left < right:
        reaction_sig = left, right
        direction_sig = eq.direction
    else:
        reaction_sig = right, left
        direction_sig = eq.direction.flipped()

    if direction:
        return reaction_sig, direction_sig
    return reaction_sig


class DuplicatesCheck(Command):
    """Check for duplicated reactions in the model.

    This command reports any reactions in the model that appear more than
    once. Stoichiometry and reaction directionality is by default disregarded
    when checking for duplicates but can be enabled using the options.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--compare-direction', action='store_true',
            help='Take reaction directionality into consideration.')
        parser.add_argument(
            '--compare-stoichiometry', action='store_true',
            help='Take stoichiometry into consideration.')
        super(DuplicatesCheck, cls).init_parser(parser)

    def run(self):
        """Run check for duplicates"""

        # Create dictonary of signatures
        database_signatures = {}
        for entry in self._model.reactions:
            signature = reaction_signature(
                entry.equation, direction=self._args.compare_direction,
                stoichiometry=self._args.compare_stoichiometry)
            database_signatures.setdefault(signature, set()).add(
                (entry.id, entry.equation, entry.filemark))

        for reaction_set in itervalues(database_signatures):
            if len(reaction_set) > 1:
                print('Found {} duplicate reactions:'.format(
                    len(reaction_set)))
                for reaction, equation, filemark in reaction_set:
                    result = ' - {}: {}'.format(reaction, equation)
                    if filemark is not None:
                        result += ' (found in {})'.format(filemark)
                    print(result)
