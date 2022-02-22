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
import os
from psamm import mapmaker
from psamm import generate_model

class TestGenerateDatabase(unittest.TestCase):
    def test_Compound_Download(self):
        # Test that the download of compounds works
        cpd = ['C00001']
        cpd_out = generate_model._download_kegg_entries(".", cpd, generate_model.CompoundEntry, context=None)
        self.assertTrue(len(list(cpd_out)) == 1)

    def test_Compound_Contents(self):
        # Test that the downloaded compound contains the relevant information
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
        # Test that the download of reactions works
        rxn = ['R00351']
        rxn_out = generate_model._download_kegg_entries(".", rxn, generate_model.ReactionEntry, context=None)
        self.assertTrue(len(list(rxn_out)) == 1)

    def test_Reaction_Contents(self):
        # Test that the downloaded reaction contains the relevant information
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

    def test_Model_Reactions(self):
        # Tests the conversion of the reactions object to a dictionary format
        rxn = ['R00351']
        rxn_out = []
        for entry in generate_model._download_kegg_entries(".", rxn, generate_model.ReactionEntry, context=None):
            rxn_out.append(entry)
        rxn_model = generate_model.model_reactions(rxn_out)
        rxn_model = list(rxn_model)
        self.assertEqual(rxn_out[0].id, rxn_model[0]['id'])
        self.assertEqual(rxn_out[0].name, rxn_model[0]['name'])
        self.assertEqual(rxn_out[0].definition, rxn_model[0]['KEGG_definition'])
        self.assertEqual(rxn_out[0].equation, rxn_model[0]['equation'])
        self.assertEqual(list(rxn_out[0].enzymes), rxn_model[0]['enzymes'])
        path = []
        for i in list(rxn_out[0].pathways):
            path.append(i[1])
        self.assertEqual(len(list(rxn_out[0].pathways)), len(rxn_model[0]['pathways']))
        self.assertEqual(path, rxn_model[0]['pathways'])

    def test_logfile(self):
        # Test that the logfile gets written
        rxn = ['R00351']
        rxn_out = generate_model._download_kegg_entries(".", rxn, generate_model.ReactionEntry, context=None)
        self.assertTrue(os.path.exists("./log.tsv"))

    def test_parseOrtho(self):
        # Test the ability to parse the defaul eggnog output
        test_egg = "gene1\ta\tb\tc\td\te\tf\tg\th\ti\t2.3.3.1\tK01647\tk\tl\tR00351\tm\tn\to\tp\tq\n"
        outfile = open("testin.tsv", "w")
        outfile.write(test_egg)
        outfile.close()
        asso = generate_model.parse_orthology("testin.tsv", "R", None)
        print(asso)
        self.assertTrue(len(asso)==1)
        self.assertTrue(asso["R00351"]==["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "KO", None)
        self.assertTrue(len(asso)==1)
        self.assertTrue(asso["K01647"]==["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "EC", None)
        self.assertTrue(len(asso)==1)
        self.assertTrue(asso["2.3.3.1"]==["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "R", 2)
        self.assertTrue(len(asso)==1)
        self.assertTrue(asso["a"]==["gene1"])
        os.remove("testin.tsv")
