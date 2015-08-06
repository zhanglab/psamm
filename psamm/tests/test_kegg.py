#!/usr/bin/env python
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

import unittest

from psamm.datasource import kegg
from psamm.reaction import Reaction, Compound
from psamm.expression.affine import Expression


class TestKEGG(unittest.TestCase):
    def test_kegg_parse(self):
        r = kegg.parse_reaction('C00013 + C00001 <=> 2 C00009')
        self.assertEqual(r, Reaction(Reaction.Bidir, [(Compound('C00013'), 1), (Compound('C00001'), 1)],
                                      [(Compound('C00009'), 2)]))

    def test_kegg_parse_with_count_expression(self):
        r = kegg.parse_reaction('2n C00404 + n C00001 <=> (n+1) C02174')
        self.assertEqual(r, Reaction(Reaction.Bidir, [(Compound('C00404'), Expression('2n')),
                                                       (Compound('C00001'), Expression('n'))],
                                      [(Compound('C02174'), Expression('n+1'))]))

    def test_kegg_parse_with_compound_argument(self):
        r = kegg.parse_reaction('C00039(n) <=> C00013 + C00039(n+1)')
        self.assertEqual(r, Reaction(Reaction.Bidir, [(Compound('C00039', arguments=[Expression('n')]), 1)],
                                      [(Compound('C00013'), 1), (Compound('C00039', arguments=[Expression('n+1')]), 1)]))


if __name__ == '__main__':
    unittest.main()
