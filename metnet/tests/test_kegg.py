#!/usr/bin/env python

import unittest

from metnet.datasource import kegg
from metnet.reaction import Reaction, Compound
from metnet.expression.affine import Expression

class TestKEGG(unittest.TestCase):
    def test_kegg_parse(self):
        r = kegg.parse_reaction('C00013 + C00001 <=> 2 C00009')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00013'), 1), (Compound('C00001'), 1)],
                                      [(Compound('C00009'), 2)]))

    def test_kegg_parse_with_count_expression(self):
        r = kegg.parse_reaction('2n C00404 + n C00001 <=> (n+1) C02174')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00404'), Expression('2n')),
                                                       (Compound('C00001'), Expression('n'))],
                                      [(Compound('C02174'), Expression('n+1'))]))

    def test_kegg_parse_with_compound_argument(self):
        r = kegg.parse_reaction('C00039(n) <=> C00013 + C00039(n+1)')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00039', arguments=[Expression('n')]), 1)],
                                      [(Compound('C00013'), 1), (Compound('C00039', arguments=[Expression('n+1')]), 1)]))


if __name__ == '__main__':
    unittest.main()
