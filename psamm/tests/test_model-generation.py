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
# Copyright 2022 Jason Vailionis <jason_vailionis@uri.edu>

import unittest
import sys
import os
from psamm import mapmaker
from psamm import generate_model
from psamm import generate_biomass
from pkg_resources import resource_filename
from psamm.datasource.native import NativeModel, ModelReader
from psamm.datasource.reaction import Reaction, Compound, Direction
from collections import defaultdict
import re
import yaml
import pandas as pd
import numpy as np
from io import StringIO


if sys.version_info.minor < 6:
    raise unittest.SkipTest("not compatible with py3.5")

class TestGenerateTransporters(unittest.TestCase):

    def test_parse_orthology(self):
        test_egg = "gene1\ta\tb\tc\td\te\tf\tg\th\ti\t2.3.3.1\tK01647\tk"
        test_egg += "\tl\tR00351\tm\tn\t2.A.49.5.2\tp\tq\n"
        outfile = open("testin.tsv", "w")
        outfile.write(test_egg)
        outfile.close()
        asso = generate_model.parse_orthology("testin.tsv", "tcdb", None)
        self.assertTrue(len(asso) == 1)
        self.assertTrue(asso["2.A.49.5.2"] == ["gene1"])

    def test_transporters(self):
        asso = {"R00351": ["gene1"], "R00375": ["gene2"]}
        generate_model.create_model_api(".", asso, False, False, "c")
        substrate = resource_filename("psamm",
                                      "external-data/tcdb_substrates.tsv")
        family = resource_filename("psamm",
                                   "external-data/tcdb_families.tsv")

        # read in the model and build a dictionary of Chebi IDs
        mr = ModelReader.reader_from_path(".""")
        nm = mr.create_model()
        mm = nm.create_metabolic_model()
        chebi_dict = defaultdict(lambda: [])
        for cpd in nm.compounds:
            if "ChEBI" in cpd.__dict__["_properties"]:
                chebi = re.split(" ", cpd.__dict__["_properties"]["ChEBI"])
                for i in chebi:
                    chebi_dict["CHEBI:{}".format(i)].append(cpd.id)
        # Read in the reaction substrates
        tp_sub_dict = defaultdict(lambda: [])
        with open(substrate, "r") as infile:
            for line in infile:
                line = line.rstrip()
                listall = re.split("\t", line)
                substrates = re.split(r"\|", listall[1])
                sub_out = []
                for i in substrates:
                    sub_out.append(re.split(";", i)[0])
                tp_sub_dict[listall[0]] = sub_out

        # read in the reaction families
        tp_fam_dict = defaultdict(lambda: "")
        with open(family, "r", encoding="utf8") as infile:
            for line in infile:
                line = line.rstrip()
                listall = re.split("\t", line)
                tp_fam_dict[listall[0]] = listall[1]
        tp_asso = {"2.A.11.1.6": ["gene2"], "1.A.8.9.20": ["gene3"]}
        generate_model.gen_transporters(".", tp_asso, tp_sub_dict,
                                        tp_fam_dict, chebi_dict, "c", "e")
        rxn = open("transporters.yaml", "r").read()
        self.assertIn("TP_2.A.11.1.6", rxn)
        self.assertIn("[c] <=> [e]", rxn)
        self.assertIn("CHEBI:34982", rxn)
        self.assertIn("CHEBI:3308", rxn)
        self.assertIn("CHEBI:50744", rxn)
        self.assertIn("C00001_diff", rxn)
        self.assertIn("[c] <=> [e]", rxn)
        self.assertIn("The Major Intrinsic Protein (MIP) Family", rxn)
        self.assertIn("CHEBI:15377", rxn)
        self.assertIn("C00001[c] <=> C00001[e]", rxn)
        # check different compartments
        generate_model.gen_transporters(".", tp_asso, tp_sub_dict,
                                        tp_fam_dict, chebi_dict, "d", "f")
        rxn = open("transporters.yaml", "r").read()
        self.assertIn("C00001[d] <=> C00001[f]", rxn)
        self.assertIn("[d] <=> [f]", rxn)
        os.remove("log.tsv")
        os.remove("reactions.yaml")
        os.remove("reactions_generic.yaml")
        os.remove("model.yaml")
        os.remove("compounds.yaml")
        os.remove("compounds_generic.yaml")
        os.remove("gene-association.tsv")
        os.remove("gene-association_generic.tsv")
        os.remove("model_def.tsv")
        os.remove("transporters.yaml")
        os.remove("transporter_log.tsv")


class TestGenerateDatabase(unittest.TestCase):
    """
    Test cases largely designed to test download of information
    from KEGG through biopython and generate a reaction database
    properly. Note that if these start to fail, one of the avenues
    to check for errors is the composition of the ec, ko, reaction, or
    compound representation in KEGG.
    """


    def test_overall(self):
        # Check ability to create expected output files
        asso = {"R00351": ["gene1"], "R00375": ["gene2"]}
        generate_model.create_model_api(".", asso, False, False, "c")
        self.assertTrue(os.path.exists("log.tsv"))
        self.assertTrue(os.path.exists("reactions.yaml"))
        self.assertTrue(os.path.exists("compounds.yaml"))
        self.assertTrue(os.path.exists("model.yaml"))
        self.assertTrue(os.path.exists("gene-association.tsv"))
        self.assertTrue(os.path.exists("gene-association_generic.tsv"))
        self.assertTrue(os.path.exists("reactions_generic.yaml"))
        self.assertTrue(os.path.exists("compounds_generic.yaml"))
        # Check contents of the created files
        mr = ModelReader.reader_from_path(".")
        nm = mr.create_model()
        mm = nm.create_metabolic_model()
        self.assertIn("R00351", set(mm.reactions))
        self.assertNotIn("R00375", set(mm.reactions))
        self.assertIn(Compound("C00010", "c"), set(mm.compounds))
        self.assertNotIn(Compound("C00039", "c"), set(mm.compounds))
        self.assertIn("2.3.3.1", nm.reactions["R00351"].
                      __dict__["_properties"]["enzymes"])
        self.assertIn("Glyoxylate and dicarboxylate metabolism",
                      nm.reactions["R00351"].
                      __dict__["_properties"]["pathways"])
        self.assertIn("K01647", nm.reactions["R00351"].
                      __dict__["_properties"]["orthology"])
        self.assertEqual("C6H5O7", nm.compounds["C00158"].
                         __dict__["_properties"]["formula"])
        self.assertEqual("16947", nm.compounds["C00158"].
                         __dict__["_properties"]["ChEBI"])
        self.assertEqual(-3, nm.compounds["C00158"].
                         __dict__["_properties"]["charge"])
        # Check for relevant entries in the generic and other files
        cpd = open("compounds_generic.yaml", "r").read()
        self.assertIn("C00039", cpd)
        rxn = open("reactions_generic.yaml", "r").read()
        self.assertIn("R00375", rxn)
        ga = open("gene-association.tsv", "r").read()
        self.assertIn("R00351", ga)
        ga = open("gene-association_generic.tsv", "r").read()
        self.assertIn("R00375", ga)
        os.remove("log.tsv")
        os.remove("reactions.yaml")
        os.remove("reactions_generic.yaml")
        os.remove("model.yaml")
        os.remove("compounds.yaml")
        os.remove("compounds_generic.yaml")
        os.remove("gene-association.tsv")
        os.remove("gene-association_generic.tsv")
        os.remove("model_def.tsv")

    def test_Rhea(self):
        # Test that the rhea flag properly captures charge of atp
        rhea_db = generate_model.RheaDb(resource_filename("psamm",
                                                          "external-data/"
                                                          "chebi_pH7_3_"
                                                          "mapping.tsv"))
        cpd = ["C00002"]
        cpd_out = generate_model._download_kegg_entries(".", cpd, None,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = list(cpd_out)
        nonRhea = list(generate_model.model_compounds(cpd_out))
        self.assertTrue(nonRhea[0]["charge"] == 0)
        cpd_out = generate_model._download_kegg_entries(".", cpd, rhea_db,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = list(cpd_out)
        Rhea = list(generate_model.model_compounds(cpd_out))
        self.assertTrue(Rhea[0]["charge"] == -4)
        os.remove("log.tsv")

    def test_EC_download(self):
        # Test when EC has one reaction
        rxn_mapping = {"2.3.3.1": ["Gene1"]}
        ec = generate_model.parse_rxns_from_EC(rxn_mapping, ".", False)
        self.assertTrue(len(ec) == 1)
        self.assertTrue(len(ec["R00351"]) == 1)
        self.assertTrue("R00351" in ec)
        self.assertTrue(ec["R00351"] == ["Gene1"])
        # Test when EC has multiple reactions
        rxn_mapping = {"4.2.1.3": ["Gene1"]}
        ec = generate_model.parse_rxns_from_EC(rxn_mapping, ".", False)
        self.assertTrue(len(ec) == 3)
        self.assertTrue("R01324" in ec)
        self.assertTrue("R01325" in ec)
        self.assertTrue("R01900" in ec)
        self.assertTrue(ec["R01324"] == ["Gene1"])
        self.assertTrue(ec["R01325"] == ["Gene1"])
        self.assertTrue(ec["R01900"] == ["Gene1"])
        # Test for multiple genes
        rxn_mapping = {"2.3.3.1": ["Gene1", "Gene2"]}
        ec = generate_model.parse_rxns_from_EC(rxn_mapping, ".", False)
        self.assertTrue(len(ec) == 1)
        self.assertTrue(len(ec["R00351"]) == 2)
        self.assertTrue("R00351" in ec)
        self.assertTrue(ec["R00351"] == ["Gene1", "Gene2"])
        os.remove("log.tsv")

    def test_KO_download(self):
        # Test when EC has one reaction
        rxn_mapping = {"K01647": ["Gene1"]}
        ko = generate_model.parse_rxns_from_KO(rxn_mapping, ".", False)
        self.assertTrue(len(ko) == 1)
        self.assertTrue(len(ko["R00351"]) == 1)
        self.assertTrue("R00351" in ko)
        self.assertTrue(ko["R00351"] == ["Gene1"])
        # Test when EC has multiple reactions
        rxn_mapping = {"K01681": ["Gene1"]}
        ko = generate_model.parse_rxns_from_KO(rxn_mapping, ".", False)
        self.assertTrue(len(ko) == 3)
        self.assertTrue("R01324" in ko)
        self.assertTrue("R01325" in ko)
        self.assertTrue("R01900" in ko)
        self.assertTrue(ko["R01324"] == ["Gene1"])
        self.assertTrue(ko["R01325"] == ["Gene1"])
        self.assertTrue(ko["R01900"] == ["Gene1"])
        # Test for multiple genes
        rxn_mapping = {"K01647": ["Gene1", "Gene2"]}
        ko = generate_model.parse_rxns_from_KO(rxn_mapping, ".", False)
        self.assertTrue(len(ko) == 1)
        self.assertTrue(len(ko["R00351"]) == 2)
        self.assertTrue("R00351" in ko)
        self.assertTrue(ko["R00351"] == ["Gene1", "Gene2"])
        os.remove("log.tsv")

    def test_model_compounds(self):
        # Tests that compounds are properly sorted into generic compounds
        # with the proper attributes.
        cpd = ["C02987", "C00001"]
        cpd_out = generate_model._download_kegg_entries(".", cpd, None,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = list(cpd_out)
        generic = generate_model.check_generic_compounds(cpd_out)
        generic_out = generate_model._download_kegg_entries(".", generic, None,
                                                            generate_model.
                                                            CompoundEntry,
                                                            context=None)
        generic_out = generate_model.model_generic_compounds(list(generic_out))
        generic_out = list(generic_out)
        self.assertTrue(len(cpd_out) == 2)
        self.assertTrue(len(generic_out) == 1)
        self.assertTrue("C02987" == cpd_out[0].id)
        self.assertTrue("C00001" == cpd_out[1].id)
        self.assertTrue("C02987" == generic_out[0]["id"])
        self.assertTrue("L-Glutamyl-tRNA(Glu)" == generic_out[0]["name"])
        self.assertTrue("C20H28N6O13PR(C5H8O6PR)n" ==
                        generic_out[0]["formula"])
        self.assertTrue("29157" == generic_out[0]["ChEBI"])
        cpd_out = list(generate_model.model_compounds(cpd_out))
        self.assertTrue(len(cpd_out) == 1)
        self.assertTrue("C00001" == cpd_out[0]["id"])
        self.assertTrue("H2O" == cpd_out[0]["name"])
        self.assertTrue("H2O" == cpd_out[0]["formula"])
        self.assertTrue("15377" == cpd_out[0]["ChEBI"])

    def test_generic_compoundID(self):
        # Test that the download of compounds works
        cpd = ["C02987", "C00001"]
        cpd_out = generate_model._download_kegg_entries(".", cpd, None,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = generate_model.check_generic_compounds(list(cpd_out))
        cpd_out = list(cpd_out)
        self.assertTrue(len(cpd_out) == 1)
        self.assertTrue("C02987" in cpd_out)
        self.assertTrue("C00001" not in cpd_out)
        os.remove("log.tsv")

    def test_Compound_Download(self):
        # Test that the download of compounds works
        cpd = ["C00001"]
        cpd_out = generate_model._download_kegg_entries(".", cpd, None,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = list(cpd_out)
        self.assertTrue(len(cpd_out) == 1)
        os.remove("log.tsv")

    def test_rxn_clean(self):
        # Test teh function that reformats stoichiometry
        rxn = ["R04347"]
        rxn_out = generate_model._download_kegg_entries(".", rxn, None,
                                                        generate_model.
                                                        ReactionEntry,
                                                        context=None)
        rxn_out, gen= generate_model.clean_reaction_equations(list(rxn_out))
        self.assertTrue("R04347" in gen)

    def test_Compound_Contents(self):
        # Test that the downloaded compound contains the relevant information
        cpd = ["C00001"]
        cpd_out = generate_model._download_kegg_entries(".", cpd, None,
                                                        generate_model.
                                                        CompoundEntry,
                                                        context=None)
        cpd_out = list(cpd_out)
        self.assertTrue(len(cpd_out) == 1)
        self.assertEqual(list(cpd_out)[0].id, "C00001")
        self.assertEqual(list(cpd_out)[0].name, "H2O")
        self.assertEqual(list(cpd_out)[0].formula, "H2O")
        self.assertEqual(list(cpd_out)[0].mol_weight, 18.0153)
        self.assertEqual(list(cpd_out)[0].chebi, "15377")
        self.assertEqual(list(cpd_out)[0].charge, 0)
        os.remove("log.tsv")

    def test_Reaction_Download(self):
        # Test that the download of reactions works
        rxn = ["R00351"]
        rxn_out = generate_model._download_kegg_entries(".", rxn, None,
                                                        generate_model.
                                                        ReactionEntry,
                                                        context=None)
        self.assertTrue(len(list(rxn_out)) == 1)
        os.remove("log.tsv")

    def test_Reaction_Contents(self):
        # Test that the downloaded reaction contains the relevant information
        rxn = ["R00351"]
        rxn_out = generate_model._download_kegg_entries(".", rxn, None,
                                                        generate_model.
                                                        ReactionEntry,
                                                        context=None)
        rxn_out = list(rxn_out)
        self.assertEqual(rxn_out[0].id, "R00351")
        self.assertEqual(rxn_out[0].name, "acetyl-CoA:oxaloacetate "
                         "C-acetyltransferase (thioester-hydrolysing)")
        self.assertEqual(rxn_out[0].definition, "Citrate + CoA <=> Acetyl-CoA"
                         " + H2O + Oxaloacetate")
        self.assertEqual(rxn_out[0].equation, "C00158 + C00010 <=> C00024 "
                         "+ C00001 + C00036")
        self.assertEqual(list(rxn_out[0].enzymes),
                         ["2.3.3.1", "2.3.3.3", "2.3.3.16"])
        self.assertIsNone(rxn_out[0].tcdb_family)
        self.assertIsNone(rxn_out[0].substrates)
        self.assertIsNone(rxn_out[0].genes)
        self.assertTrue(len(list(rxn_out[0].pathways)) == 8)
        self.assertEqual(list(rxn_out[0].pathways)[0],
                         ("rn00020", "Citrate cycle (TCA cycle)"))
        os.remove("log.tsv")

    def test_Model_Reactions(self):
        # Tests the conversion of the reactions object to a dictionary format
        rxn = ["R00351"]
        rxn_out = []
        for entry in generate_model._download_kegg_entries(".", rxn, None,
                                                           generate_model.
                                                           ReactionEntry,
                                                           context=None):
            rxn_out.append(entry)
        rxn_model = generate_model.model_reactions(rxn_out)
        rxn_model = list(rxn_model)
        self.assertEqual(rxn_out[0].id, rxn_model[0]["id"])
        self.assertEqual(rxn_out[0].name, rxn_model[0]["name"])
        self.assertEqual(rxn_out[0].definition, rxn_model[0]
                         ["KEGG_definition"])
        self.assertEqual(rxn_out[0].equation, rxn_model[0]["equation"])
        self.assertEqual(list(rxn_out[0].enzymes), rxn_model[0]["enzymes"])
        path = []
        for i in list(rxn_out[0].pathways):
            path.append(i[1])
        self.assertEqual(len(list(rxn_out[0].pathways)), len(rxn_model[0]
                         ["pathways"]))
        self.assertEqual(path, rxn_model[0]["pathways"])
        os.remove("log.tsv")

    def test_parseOrtho(self):
        # Test the ability to parse the defaul eggnog output
        test_egg = "gene1\ta\tb\tc\td\te\tf\tg\th\ti\t2.3.3.1\t"
        test_egg += "K01647\tk\tl\tR00351\tm\tn\to\tp\tq\n"
        outfile = open("testin.tsv", "w")
        outfile.write(test_egg)
        outfile.close()
        asso = generate_model.parse_orthology("testin.tsv", "R", None)
        self.assertTrue(len(asso) == 1)
        self.assertTrue(asso["R00351"] == ["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "KO", None)
        self.assertTrue(len(asso) == 1)
        self.assertTrue(asso["K01647"] == ["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "EC", None)
        self.assertTrue(len(asso) == 1)
        self.assertTrue(asso["2.3.3.1"] == ["gene1"])
        asso = generate_model.parse_orthology("testin.tsv", "R", 2)
        self.assertTrue(len(asso) == 1)
        self.assertTrue(asso["a"] == ["gene1"])
        os.remove("testin.tsv")


class TestGenerateBiomass(unittest.TestCase):

    def test_cpd_database_load(self):
        df = generate_biomass.load_compound_data()
        self.assertTrue(df.shape == (79, 7))
        config = generate_biomass.make_template(df)
        self.assertTrue(config.shape == (79, 3))
        config.custom_id = list(range(79))  # config w/ custom ids from 0 to 78
        config.to_csv("config_file.csv", index=False)
        modified_df = generate_biomass.apply_config(df, "config_file.csv")
        self.assertTrue(modified_df.shape == (79, 8))
        self.assertTrue(modified_df.loc["biomass", "id"] == 0)
        self.assertTrue(modified_df.loc["C02554", "id"] == 78)
        os.remove("config_file.csv")

    def test_model_parsing(self):
        # THIS will include the way to build a fake model
        # Then, test the load_model, check_missing_cpds, fix_cell_compartment,
        # and update_model_yaml methods
        yaml_args = {"default_flow_style": False,
                     "sort_keys": False,
                     "encoding": "utf-8",
                     "allow_unicode": True,
                     "width": float("inf")}
        test_model_dict = {"default_compartment": None, "biomass": None,
                           "compounds": [{"include": "compounds.yaml"}],
                           "reactions": [{"include": "reactions.yaml"}],
                           "model": [{"include": "model_def.tsv"}]}
        tmp_model_path = "tmp_model.yaml"
        df = generate_biomass.load_compound_data()
        with open(tmp_model_path, "w") as f:
            yaml.dump(test_model_dict, f, **yaml_args)
        model_dict = generate_biomass.load_model(tmp_model_path)
        self.assertTrue(test_model_dict == model_dict)
        # FAKE compounds.yaml which contains one of the 79 biomass compounds
        # Testing if generate_biomass.yaml has the other 78 compounds
        with open("compounds.yaml", "w") as f:
            f.write("- id: C02554")
        generate_biomass.check_missing_cpds(
            tmp_model_path, df, None, yaml_args)
        # Read in the biomass_compounds.yaml and make sure theres 78 cpds
        # matching the first 78 in the df
        with open("biomass_compounds.yaml", "r") as f:
            matches = re.findall("- id:.*\n", f.read())
            cpds = [x.rstrip("\n").replace("- id: ", "") for x in matches]
        self.assertTrue(len(cpds) == 78)
        self.assertTrue(list(df.id)[:-1] == cpds)
        # Check fix_cell_compartment()
        model_dict["default_compartment"] = None
        self.assertTrue(
            generate_biomass.fix_cell_compartment(model_dict) == "c")
        model_dict["default_compartment"] = "x"
        self.assertTrue(
            generate_biomass.fix_cell_compartment(model_dict) == "x")
        del model_dict["default_compartment"]
        self.assertTrue(
            generate_biomass.fix_cell_compartment(model_dict) == "c")
        # Check update_model_yaml()
        test_model_dict["compounds"].append(
            {"include": "./biomass_compounds.yaml"})
        test_model_dict["reactions"].append(
            {"include": "./biomass_reactions.yaml"})
        test_model_dict["biomass"] = "biomass_name"
        test_model_dict["default_compartment"] = "c"
        generate_biomass.update_model_yaml(
            tmp_model_path, model_dict, "c", "biomass_name", yaml_args)
        updated_model_dict = generate_biomass.load_model(tmp_model_path)
        self.assertTrue(test_model_dict == updated_model_dict)
        os.remove("biomass_compounds.yaml")
        os.remove(tmp_model_path)
        os.remove("compounds.yaml")

    def test_calc_stoichiometry(self):
        df = pd.DataFrame({"mol_weight": [1, 2, 3, 3]})
        out = generate_biomass.calc_stoichiometry_df(df, [2, 2, 3, 3])
        self.assertTrue(out.shape == (4, 5))
        self.assertTrue(list(out.columns) == ["mol_weight", "counts",
                                              "composition", "mass",
                                              "stoichiometry"])
        np.testing.assert_almost_equal(out.stoichiometry,
                                       [0.0833333, 0.0833333, 0.125, 0.125])

    def test_dna_entry(self):
        from Bio import SeqIO
        with StringIO(">seq1\nAAGGGT\n>seq2\nAGCT") as seqfile:
            genome = {seq.name: seq.upper()
                      for seq in SeqIO.parse(seqfile, "fasta")}
        df = generate_biomass.load_compound_data()
        dna_output = generate_biomass.count_DNA(genome, df, "c")
        print(dna_output)
        proper_output = {"id": "dna_met", "name": "DNA",
                         "equation": "(0.613065) C00131[c] + (0.40871) "
                         "C00459[c] + (0.81742) C00286[c] + (0.204355) "
                         "C00458[c] + (4.0871) C00002[c] + (4.0871) C00001[c] "
                         "=> dna[c] + (4.0871) C00008[c] + (4.0871) C00009[c] "
                         "+ (2.04355) C00013[c]", "pathways": ["Biomass"]}
        self.assertTrue(dna_output == proper_output)

    def test_rna_entry(self):
        from Bio import SeqIO
        df = generate_biomass.load_compound_data()
        with StringIO(">seq1\nAAGGGT\n>seq2\nAGCT") as fna:
            genome = {seq.name: seq.upper()
                      for seq in SeqIO.parse(fna, "fasta")}
        with open("temp.gff", "w") as gff:
            gff.write("seq1\tProtein Homology\tCDS\t2\t6\t.\t+\t0\ttest\n" +
                      "seq2\tProtein Homology\tCDS\t2\t4\t.\t-\t0\ttest")
        rna_output = generate_biomass.count_RNA(genome, "temp.gff", df, "c")
        proper_output = {"id": "rna_met", "name": "RNA",
                         "equation": "(4.248569) C00002[c] + (0.326813) "
                         "C00075[c] + (1.307252) C00044[c] + (3.921756) "
                         "C00001[c] => rna[c] + (3.921756) C00008[c] + "
                         "(3.921756) C00009[c] + (1.960878) C00013[c]",
                         "pathways": ["Biomass"]}
        self.assertTrue(rna_output == proper_output)
        os.remove("temp.gff")

    def test_prot_entry(self):
        from Bio import SeqIO
        df = generate_biomass.load_compound_data()
        with StringIO(">seq1\nLEW\n>seq2\nMIG") as faa:
            proteome = {seq.name: seq.upper()
                        for seq in SeqIO.parse(faa, "fasta")}
        prot_output = generate_biomass.count_Prot(proteome, df, "c")
        proper_output = {"id": "pro_met", "name": "Protein",
                         "equation": "(1.194786) C02987[c] + (1.194786) "
                         "C02412[c] + (1.194786) C03127[c] + (1.194786) "
                         "C02047[c] + (1.194786) C02430[c] + (1.194786) "
                         "C03512[c] + (14.337431) C00044[c] + (14.337431) "
                         "C00001[c] => (1.194786) C01641[c] + (1.194786) "
                         "C01642[c] + (1.194786) C01644[c] + (1.194786) "
                         "C01645[c] + (1.194786) C01647[c] + (1.194786) "
                         "C01652[c] + protein[c] + (14.337431) C00035[c] + "
                         "(14.337431) C00009[c]", "pathways": ["Biomass"]}
        self.assertTrue(prot_output == proper_output)

    def test_bio_entry(self):
        df = generate_biomass.load_compound_data()
        bio_output = generate_biomass.return_biomass_rxn(df, "biomass", "c")
        proper_output = {"id": "biomass", "name": "Biomass",
                         "equation": "dna[c] + rna[c] + protein[c] => "
                         "biomass[c]", "pathways": ["Biomass"]}
        self.assertTrue(bio_output == proper_output)

    def test_sink_entry(self):
        df = generate_biomass.load_compound_data()
        sink_output = generate_biomass.return_bio_sink_rxn(df, "c")
        proper_output = {"id": "sink_biomass", "name": "Biomass accumulation",
                         "equation": "biomass[c] =>", "pathways": ["Biomass"]}
        self.assertTrue(sink_output == proper_output)

    def test_tRNA_entries(self):
        df = generate_biomass.load_compound_data()
        trna_output = generate_biomass.return_trna_rxns(df, "c")
        self.assertTrue(len(trna_output) == 20)
        proper_output_lastaa = {"id": "R03665",
                                "name": "L-Valine:tRNAVal ligase "
                                "(AMP-forming)",
                                "equation": "C00002[c] + C00183[c] + C01653[c]"
                                " <=> C00020[c] + C00013[c] + C02554[c]",
                                "enzyme": "['6.1.1.9']",
                                "pathways": ["Biomass"]}
        self.assertTrue(trna_output[19] == proper_output_lastaa)
