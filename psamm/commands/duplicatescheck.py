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

# File originally based on import from fastgapfill.py @ commit db72df89a7c0007d2392584e94385987e374164e

from __future__ import unicode_literals

import argparse
import logging

from ..command import (Command)

logger = logging.getLogger(__name__)




# This function is from Keith's script, "reaction_signature.py". Eventually it should be moved to a tools module or a stand-alone module.
def reaction_signature(eq):
    return (tuple(sorted(str(compound) for compound, _ in eq.left)),
            tuple(sorted(str(compound) for compound, _ in eq.right)))




class DuplicatesCheck(Command):
    """Run duplicate checking algorithm on model. This finds any reactions that are duplicates. It does not remove them however."""

# Assumptions: Reactions are balenced correctly.

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--epsilon', type=float, default=1e-5,
            help='Threshold for compound production')
        super(DuplicatesCheck, cls).init_parser(parser)

    def run(self):
        """Run DuplicatesCheck command"""

        database_signatures = {}
        # Reading reaction entries from model. Model is passed into self, but in not yet parsed.
        for entry in self._model.parse_reactions():
            signature = reaction_signature(entry.equation)
            database_signatures.setdefault(signature, set()).add((entry.id,entry.equation))

        ## R_48: A C B -> X Y Z
        ## R_56: A C B -> X Y Z

        ## A C B -> X Y Z
        ## Z Y X -> A B C

        # ABC->XYZ
        # XYZ->ABC

        # ABC,l -> XYZ,r
        # XYZ,l -> ABC,r

        # ABC,l -> XYZ,r
        # ABC,r -> XYZ,l

        # sort both sides
        # label each side list, such that it is a tuple of list+{left-side reaction|right-side reaction}
        # sort the two lists
        #loss = arrangment on each side

        for reaction_set in database_signatures.itervalues():
            if len(reaction_set) > 1:
                print('Signature compounds identical for {} reactions:'.format(
                    len(reaction_set)))
                for reaction, equation in reaction_set:
                    print(' - {}: {}'.format(
                        reaction, equation))


    ###
