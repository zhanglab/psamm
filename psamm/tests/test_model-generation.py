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
# Copyright 2021 Christopher Powers <c-11060@uri.edu>

import unittest
import sys
from psamm import mapmaker
from psamm import generate_model

class TestGenerateDatabase(unittest.TestCase):
    def test_Compound_Download(self):
        cpd = ['C00001']
        cpd_out = generate_model._download_kegg_entries(".", cpd, generate_model.CompoundEntry, context=None)
        self.assertTrue(len(list(cpd_out)) == 1)

    def test_Compound_Contents(self):
        cpd = ['C00001']
        cpd_out = generate_model._download_kegg_entries(".", cpd, generate_model.CompoundEntry, context=None)
        cpd_out = list(cpd_out)
        self.assertEqual(list(cpd_out)[0].id, "C00001")
        self.assertEqual(list(cpd_out)[0].name, "H2O")
        self.assertEqual(list(cpd_out)[0].formula, "H2O")
        self.assertEqual(list(cpd_out)[0].mol_weight, 18.0153)
        self.assertEqual(list(cpd_out)[0].chebi, "15377")
        self.assertEqual(list(cpd_out)[0].charge, 0)

    def test_Reaction_Download(self):
        rxn = ['R00351']
        rxn_out = generate_model._download_kegg_entries(".", rxn, generate_model.ReactionEntry, context=None)
        self.assertTrue(len(list(rxn_out)) == 1)

    def test_Reaction_Contents(self):
        rxn = ['R00351']
        rxn_out = generate_model._download_kegg_entries(".", rxn, generate_model.ReactionEntry, context=None)
        rxn_out = list(rxn_out)
        self.assertEqual(rxn_out[0].id, "R00351")
        self.assertEqual(rxn_out[0].name, "acetyl-CoA:oxaloacetate C-acetyltransferase (thioester-hydrolysing)")
        self.assertEqual(rxn_out[0].definition, "Citrate + CoA <=> Acetyl-CoA + H2O + Oxaloacetate")
        self.assertEqual(rxn_out[0].equation, "C00158 + C00010 <=> C00024 + C00001 + C00036")
        self.assertEqual(list(rxn_out[0].enzymes), ["2.3.3.1","2.3.3.3","2.3.3.16"])
        self.assertIsNone(rxn_out[0].tcdb_family)
        self.assertIsNone(rxn_out[0].substrates)
        self.assertIsNone(rxn_out[0].genes)
        self.assertTrue(len(list(rxn_out[0].pathways))==8)
        self.assertEqual(list(rxn_out[0].pathways)[0], ('rn00020',  'Citrate cycle (TCA cycle)'))
