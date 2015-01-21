#!/usr/bin/env python

import unittest
from decimal import Decimal

from metnet.reaction import Reaction, Compound
from metnet.datasource import modelseed

class TestModelSEED(unittest.TestCase):
    def test_modelseed_parse(self):
        r = modelseed.parse_reaction('|H2O| + |PPi| => (2) |Phosphate| + (2) |H+|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2O'), 1), (Compound('PPi'), 1)],
                                      [(Compound('Phosphate'), 2), (Compound('H+'), 2)]))

    def test_modelseed_parse_with_decimal(self):
        r = modelseed.parse_reaction('|H2| + (0.5) |O2| => |H2O|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2'), 1), (Compound('O2'), Decimal('0.5'))],
                                      [(Compound('H2O'), 1)]))

    def test_modelseed_parse_with_compartment(self):
        r = modelseed.parse_reaction('(2) |H2| + |O2| => (2) |H2O[e]|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2'), 2), (Compound('O2'), 1)],
                                      [(Compound('H2O', compartment='e'), 2)]))

    def test_modelseed_parse_with_multichar_compartment(self):
        r = modelseed.parse_reaction('(2) |H2[C_c]| + |O2[C_c]| => (2) |H2O[C_e]|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2', compartment='C_c'), 2),
                                                        (Compound('O2', compartment='C_c'), 1)],
                                      [(Compound('H2O', compartment='C_e'), 2)]))

    def test_modelseed_str(self):
        r = Reaction(Reaction.Left, [(Compound('H2O'), 2)], [(Compound('H2'), 2), (Compound('O2'), 1)])
        self.assertEquals(str(r), '(2) |H2O| <= (2) |H2| + |O2|')


if __name__ == '__main__':
    unittest.main()
