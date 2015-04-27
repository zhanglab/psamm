
"""Test miscellaneous reaction parsers"""

import unittest
from decimal import Decimal

from psamm.reaction import Reaction, Compound
from psamm.datasource import misc

class TestSudenSimple(unittest.TestCase):
    def test_sudensimple_parse(self):
        r = misc.parse_sudensimple_reaction('1 H2O + 1 PPi <=> 2 Phosphate + 2 proton')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('H2O'), 1), (Compound('PPi'), 1)],
                                      [(Compound('Phosphate'), 2), (Compound('proton'), 2)]))

    def test_sudensimple_parse_with_implicit_count(self):
        r = misc.parse_sudensimple_reaction('H2O + PPi <=> 2 Phosphate + 2 proton')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('H2O'), 1), (Compound('PPi'), 1)],
                                      [(Compound('Phosphate'), 2), (Compound('proton'), 2)]))

    def test_sudensimple_parse_with_decimal(self):
        r = misc.parse_sudensimple_reaction('1 H2 + 0.5 O2 <=> 1 H2O')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('H2'), 1), (Compound('O2'), Decimal('0.5'))],
                                      [(Compound('H2O'), 1)]))

class TestMetNet(unittest.TestCase):
    def test_metnet_parse_with_global_compartment(self):
        r = misc.parse_metnet_reaction('[c] : akg + ala-L <==> glu-L + pyr')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('akg', 'c'), 1), (Compound('ala-L', 'c'), 1)],
                                      [(Compound('glu-L', 'c'), 1), (Compound('pyr', 'c'), 1)]))

    def test_metnet_parse_with_local_compartment(self):
        r =  misc.parse_metnet_reaction('(2) ficytcc553[c] + so3[c] + h2o[c] --> (2) focytcc553[c] + so4[c] + (2) h[e]')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('ficytcc553', 'c'), 2), (Compound('so3', 'c'), 1),
                                                       (Compound('h2o', 'c'), 1)],
                                      [(Compound('focytcc553', 'c'), 2), (Compound('so4', 'c'), 1),
                                       (Compound('h', 'e'), 2)]))

    def test_metnet_parse_global_with_colon_in_name(self):
        r = misc.parse_metnet_reaction('[c] : fdxr-4:2 + h + nadp <==> fdxo-4:2 + nadph')
        self.assertEquals(r, Reaction(Reaction.Bidir,
                                        [(Compound('fdxr-4:2', 'c'), 1), (Compound('h', 'c'), 1),
                                         (Compound('nadp', 'c'), 1)],
                                        [(Compound('fdxo-4:2', 'c'), 1), (Compound('nadph', 'c'), 1)]))

    def test_metnet_parse_local_with_colon_in_name(self):
        r = misc.parse_metnet_reaction('fdxr-4:2[c] + h[c] + nadp[c] <==> fdxo-4:2[c] + nadph[c]')
        self.assertEquals(r, Reaction(Reaction.Bidir,
                                        [(Compound('fdxr-4:2', 'c'), 1), (Compound('h', 'c'), 1),
                                         (Compound('nadp', 'c'), 1)],
                                        [(Compound('fdxo-4:2', 'c'), 1), (Compound('nadph', 'c'), 1)]))

    def test_metnet_parse_with_numeric_id(self):
        r = misc.parse_metnet_reaction('[c] : 3pg + atp <==> 13dpg + adp')
        self.assertEquals(r, Reaction(Reaction.Bidir,
                                      [(Compound('3pg', 'c'), 1), (Compound('atp', 'c'), 1)],
                                      [(Compound('13dpg', 'c'), 1), (Compound('adp', 'c'), 1)]))

if __name__ == '__main__':
    unittest.main()
