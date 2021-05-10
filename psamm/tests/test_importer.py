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
# Copyright 2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import os
import tempfile
import shutil
import json
import unittest
from decimal import Decimal

import numpy as np
from scipy.io import savemat, loadmat
from scipy.sparse.csc import csc_matrix

from six import BytesIO

from psamm import importer
from psamm.importer import ModelLoadError
from psamm.datasource.native import NativeModel
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry,
                                    DictCompartmentEntry as CompartmentEntry)
from psamm.datasource.reaction import parse_reaction
from psamm.reaction import Reaction, Compound
from psamm.formula import Formula
from psamm.expression import boolean

from psamm.importers.sbml import StrictImporter as SBMLStrictImporter
from psamm.importers.sbml import NonstrictImporter as SBMLNonstrictImporter
from psamm.importers.cobrajson import Importer as CobraJSONImporter
from psamm.importers.matlab import Importer as MatlabImporter


class TestImporterBaseClass(unittest.TestCase):
    def setUp(self):
        self.importer = importer.Importer()

    def test_try_parse_formula(self):
        result = self.importer._try_parse_formula('cpd_1', ' ')
        self.assertIsNone(result)

        # Should not return parsed result!
        result = self.importer._try_parse_formula('cpd_1', 'CO2')
        self.assertEqual(result, 'CO2')

        result = self.importer._try_parse_formula('cpd_1', 'XYZ...\u03c0')
        self.assertEqual(result, 'XYZ...\u03c0')

    def test_try_parse_reaction(self):
        result = self.importer._try_parse_reaction('rxn_1', 'A => (3) B')
        self.assertIsInstance(result, Reaction)

        with self.assertRaises(importer.ParseError):
            self.importer._try_parse_reaction('rxn_1', '+ + ==> \u03c0 5')

    def test_try_parse_gene_association(self):
        result = self.importer._try_parse_gene_association(
            'rxn_1', 'g1 and (g2 or g3)')
        self.assertIsInstance(result, boolean.Expression)

        result = self.importer._try_parse_gene_association('rxn_1', ' ')
        self.assertIsNone(result)

        result = self.importer._try_parse_gene_association('rxn_1', '123 and')
        self.assertEqual(result, '123 and')


class TestImporterCheckAndWriteModel(unittest.TestCase):
    def setUp(self):
        self.dest = tempfile.mkdtemp()

        self.model = NativeModel({
            'id': 'test_mode',
            'name': 'Test model',
        })

        # Compounds
        self.model.compounds.add_entry(CompoundEntry({
            'id': 'cpd_1',
            'formula': Formula.parse('CO2')
        }))
        self.model.compounds.add_entry(CompoundEntry({
            'id': 'cpd_2',
        }))
        self.model.compounds.add_entry(CompoundEntry({
            'id': 'cpd_3',
        }))

        # Compartments
        self.model.compartments.add_entry(CompartmentEntry({
            'id': 'c'
        }))
        self.model.compartments.add_entry(CompartmentEntry({
            'id': 'e'
        }))

        # Reactions
        self.model.reactions.add_entry(ReactionEntry({
            'id': 'rxn_1',
            'equation': parse_reaction('(2) cpd_1[c] <=> cpd_2[e] + cpd_3[c]'),
            'genes': boolean.Expression('g_1 and (g_2 or g_3)'),
            'subsystem': 'Some subsystem'
        }))
        self.model.reactions.add_entry(ReactionEntry({
            'id': 'rxn_2',
            'equation': parse_reaction(
                '(0.234) cpd_1[c] + (1.34) cpd_2[c] => cpd_3[c]'),
            'subsystem': 'Biomass'
        }))

        self.model.biomass_reaction = 'rxn_2'
        self.model.limits['rxn_1'] = 'rxn_1', Decimal(-42), Decimal(1000)
        self.model.exchange[Compound('cpd_2', 'e')] = (
            Compound('cpd_2', 'e'), 'EX_cpd_2', -10, 0)

    def tearDown(self):
        shutil.rmtree(self.dest)

    def test_get_default_compartment(self):
        importer.get_default_compartment(self.model)

    def test_detect_best_flux_limit(self):
        flux_limit = importer.detect_best_flux_limit(self.model)
        self.assertEqual(flux_limit, 1000)

    def test_infer_compartment_entries(self):
        self.model.compartments.clear()
        importer.infer_compartment_entries(self.model)

        self.assertEqual(len(self.model.compartments), 2)
        self.assertEqual({c.id for c in self.model.compartments}, {'c', 'e'})

    def test_get_output_limit(self):
        self.assertIsNone(importer._get_output_limit(1000, 1000))
        self.assertEqual(importer._get_output_limit(-5, 0), -5)

    def test_generate_limit_items(self):
        item = dict(importer._generate_limit_items(None, None))
        self.assertEqual(item, {})

        item = dict(importer._generate_limit_items(None, 400))
        self.assertEqual(item, {'upper': 400})

        item = dict(importer._generate_limit_items(-56, None))
        self.assertEqual(item, {'lower': -56})

        item = dict(importer._generate_limit_items(-224, 0))
        self.assertEqual(item, {'lower': -224, 'upper': 0})

        item = dict(importer._generate_limit_items(-5, -5))
        self.assertEqual(item, {'fixed': -5})

    def test_count_genes(self):
        self.assertEqual(importer.count_genes(self.model), 3)

    def test_write_yaml_model(self):
        importer.write_yaml_model(self.model, self.dest)


class TestSBMLImporter(unittest.TestCase):
    def setUp(self):
        self.dest = tempfile.mkdtemp()
        with open(os.path.join(self.dest, 'model.sbml'), 'wb') as f:
            f.write('''<?xml version="1.0" encoding="UTF-8"?>
    <sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
          xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
          xmlns:html="http://www.w3.org/1999/xhtml"
          level="3" version="1"
          fbc:required="false">
     <model id="test_model" name="Test model">
      <notes>
       <body xmlns="http://www.w3.org/1999/xhtml">
        <p>This is model information intended to be seen by humans.</p>
       </body>
      </notes>
      <listOfCompartments>
       <compartment id="C_c" name="cell" constant="true"/>
       <compartment id="C_b" name="boundary" constant="true"/>
      </listOfCompartments>
      <listOfSpecies>
       <species id="M_Glucose_LPAREN_c_RPAREN_"
                metaid="meta_M_Glucose_LPAREN_c_RPAREN_" name="Glucose"
                compartment="C_c" constant="false" boundaryCondition="false"
                hasOnlySubstanceUnits="false" fbc:charge="0"
                fbc:chemicalFormula="C6H12O6">
        <notes>
         <body xmlns="http://www.w3.org/1999/xhtml">
          <p>This is compound information intended to be seen by humans.</p>
         </body>
        </notes>
       </species>
       <species id="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_" name="Glucose-6-P"
                compartment="C_c" constant="false" boundaryCondition="false"
                hasOnlySubstanceUnits="false" fbc:charge="-2"
                fbc:chemicalFormula="C6H11O9P"/>
       <species id="M_H2O_LPAREN_c_RPAREN_" name="H2O" compartment="C_c"
                constant="false" boundaryCondition="false"
                hasOnlySubstanceUnits="false" fbc:charge="0"
                fbc:chemicalFormula="H2O"/>
       <species id="M_Phosphate_LPAREN_c_RPAREN_" name="Phosphate"
                compartment="C_c" constant="false" boundaryCondition="false"
                hasOnlySubstanceUnits="false" fbc:charge="-2"
                fbc:chemicalFormula="HO4P"/>
       <species id="M_Biomass" name="Biomass" compartment="C_b" constant="false"
                boundaryCondition="true" hasOnlySubstanceUnits="false"/>
      </listOfSpecies>
      <listOfParameters>
       <parameter constant="true" id="P_lower_G6Pase" value="-10"/>
       <parameter constant="true" id="P_lower_zero" value="0"/>
       <parameter constant="true" id="P_upper_default" value="1000"/>
      </listOfParameters>
      <listOfReactions>
       <reaction id="R_G6Pase" metaid="meta_R_G6Pase" reversible="true"
                 fast="false" fbc:lowerFluxBound="P_lower_G6Pase"
                 fbc:upperFluxBound="P_upper_default">
        <listOfReactants>
         <speciesReference species="M_Glucose_LPAREN_c_RPAREN_" stoichiometry="2"
                           constant="true"/>
         <speciesReference species="M_Phosphate_LPAREN_c_RPAREN_" stoichiometry="2"
                           constant="true"/>
        </listOfReactants>
        <listOfProducts>
         <speciesReference species="M_H2O_LPAREN_c_RPAREN_" stoichiometry="2"
                           constant="true"/>
         <speciesReference species="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_"
                           stoichiometry="2" constant="true"/>
        </listOfProducts>
        <notes>
         <body xmlns="http://www.w3.org/1999/xhtml">
          <p>This is reaction information intended to be seen by humans.</p>
         </body>
        </notes>
        <annotation>
         <RDF xmlns="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
              xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
          <Description about="#meta_R_G6Pase">
           <bqbiol:isVersionOf>
            <Bag>
             <li resource="http://identifiers.org/ec-code/3.1.3.9"/>
            </Bag>
           </bqbiol:isVersionOf>
          </Description>
         </RDF>
        </annotation>
       </reaction>
       <reaction id="R_Biomass" reversible="false" fast="false"
                 fbc:lowerFluxBound="P_lower_zero"
                 fbc:upperFluxBound="P_upper_default">
        <listOfReactants>
         <speciesReference species="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_"
                           stoichiometry="0.56" constant="true"/>
         <speciesReference species="M_Glucose_LPAREN_c_RPAREN_"
                           stoichiometry="0.88" constant="true"/>
        </listOfReactants>
        <listOfProducts>
         <speciesReference species="M_Biomass" stoichiometry="1" constant="true"/>
        </listOfProducts>
       </reaction>
      </listOfReactions>
      <fbc:listOfObjectives fbc:activeObjective="obj1">
       <fbc:objective fbc:id="obj1" fbc:name="Objective 1" fbc:type="maximize">
        <fbc:listOfFluxObjectives>
         <fbc:fluxObjective fbc:reaction="R_Biomass" fbc:coefficient="1"/>
        </fbc:listOfFluxObjectives>
       </fbc:objective>
      </fbc:listOfObjectives>
     </model>
    </sbml>'''.encode('utf-8'))

    def tearDown(self):
        shutil.rmtree(self.dest)

    def test_create_strict_importer(self):
        importer = SBMLStrictImporter()
        importer.help()
        importer.import_model(os.path.join(self.dest, 'model.sbml'))

    def test_strict_importer_main(self):
        output_dir = os.path.join(self.dest, 'output')
        os.mkdir(output_dir)
        importer.main(importer_class=SBMLStrictImporter, args=[
            '--source', os.path.join(self.dest, 'model.sbml'),
            '--dest', output_dir
        ])

    def test_create_nonstrict_importer(self):
        importer = SBMLNonstrictImporter()
        importer.help()
        importer.import_model(os.path.join(self.dest, 'model.sbml'))

    def test_nonstrict_importer_main(self):
        output_dir = os.path.join(self.dest, 'output')
        os.mkdir(output_dir)
        importer.main(importer_class=SBMLNonstrictImporter, args=[
            '--source', os.path.join(self.dest, 'model.sbml'),
            '--dest', output_dir
        ])


class TestCobraJSONImporter(unittest.TestCase):
    def setUp(self):
        self.dest = tempfile.mkdtemp()
        with open(os.path.join(self.dest, 'model.json'), 'w') as f:
            json.dump({
                'id': 'test_model',
                'compartments': {
                    'c': 'Cytosol',
                    'e': 'Extracellular',
                },
                'metabolites': [
                    {
                        'id': 'cpd_1',
                        'name': 'Compound 1',
                        'formula': 'CO2',
                        'compartment': 'c'
                    },
                    {
                        'id': 'cpd_2',
                        'name': 'Compound 2',
                        'charge': -2,
                        'compartment': 'e'
                    },
                    {
                        'id': 'cpd_3',
                        'compartment': 'c'
                    }
                ],
                'reactions': [
                    {
                        'id': 'rxn_1',
                        'metabolites': {
                            'cpd_1': -1,
                            'cpd_2': 2,
                            'cpd_3': 1,
                        },
                        'lower_bound': -20,
                        'upper_bound': 100,
                        'gene_reaction_rule': 'gene_1 and (gene_2 or gene_3)'
                    }
                ]
            }, f)

    def tearDown(self):
        shutil.rmtree(self.dest)

    def test_create_importer(self):
        importer = CobraJSONImporter()
        importer.help()
        importer.import_model(os.path.join(self.dest, 'model.json'))

    def test_cobra_json_main(self):
        output_dir = os.path.join(self.dest, 'output')
        os.mkdir(output_dir)
        importer.main(importer_class=CobraJSONImporter, args=[
            '--source', os.path.join(self.dest, 'model.json'),
            '--dest', output_dir
        ])

    def test_cobra_json_main_bigg(self):
        def mock_urlopen(url):
            return open(os.path.join(self.dest, 'model.json'), 'rb')

        output_dir = os.path.join(self.dest, 'output')
        os.mkdir(output_dir)
        importer.main_bigg(
            args=['mock_model', '--dest', output_dir], urlopen=mock_urlopen)


class TestBiGGImportMain(unittest.TestCase):
    def test_list_models(self):
        def mock_urlopen(url):
            return BytesIO(b'''
{"results_count": 4,
 "results": [
   {"organism": "Escherichia coli str. K-12 substr. MG1655",
    "metabolite_count": 72,
    "bigg_id": "e_coli_core",
    "gene_count": 137,
    "reaction_count": 95},
   {"organism": "Homo sapiens",
    "metabolite_count": 342,
    "bigg_id": "iAB_RBC_283",
    "gene_count": 346,
    "reaction_count": 469},
   {"organism": "Escherichia coli str. K-12 substr. MG1655",
    "metabolite_count": 1668,
    "bigg_id": "iAF1260",
    "gene_count": 1261,
    "reaction_count": 2382}]}''')

        importer.main_bigg(['list'], urlopen=mock_urlopen)


class TestMatlabImporter(unittest.TestCase):
    def setUp(self):
        self.dest = tempfile.mkdtemp()
        dt = np.dtype([
            ('comps', 'O'),
            ('compNames', 'O'),
            ('rxns', 'O'),
            ('rxnNames', 'O'),
            ('rxnGeneMat', 'O'),
            ('ub', 'O'),
            ('lb', 'O'),
            ('S', 'O'),
            ('grRules', 'O'),
            ('subSystems', 'O'),
            ('mets', 'O'),
            ('metFormulas', 'O'),
            ('metNames', 'O'),
            ('metCharges', 'O'),
            ('metNotes', 'O'),
            ('c', 'O'),
        ])
        model1 = np.ndarray((1, 1), dtype=dt)
        # compartments
        model1['comps'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['comps'][0, 0][0][0] = np.array(['c'])
        model1['comps'][0, 0][1][0] = np.array(['e'])
        model1['comps'][0, 0][2][0] = np.array(['b'])
        model1['compNames'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['compNames'][0, 0][0][0] = np.array(['Cell'])
        model1['compNames'][0, 0][1][0] = np.array(['Extracellular'])
        model1['compNames'][0, 0][2][0] = np.array(['Boundary'])
        # reactions
        model1['rxns'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['rxns'][0, 0][0][0] = np.array(['rxn1'])
        model1['rxns'][0, 0][1][0] = np.array(['rxn2'])
        model1['rxns'][0, 0][2][0] = np.array(['rxn3'])
        model1['rxnNames'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['rxnNames'][0, 0][0][0] = np.array(['rxn1'])
        model1['rxnNames'][0, 0][1][0] = np.array([])
        model1['rxnNames'][0, 0][2][0] = np.array(['rxn3'])
        model1['rxnGeneMat'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['rxnGeneMat'][0, 0][0][0] = np.array([])
        model1['rxnGeneMat'][0, 0][1][0] = np.array([0])
        model1['rxnGeneMat'][0, 0][2][0] = np.array([0])
        model1['ub'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['ub'][0, 0][0][0] = np.array([1000])
        model1['ub'][0, 0][1][0] = np.array([0])
        model1['ub'][0, 0][2][0] = np.array([1000])
        model1['lb'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['lb'][0, 0][0][0] = np.array([0])
        model1['lb'][0, 0][1][0] = np.array([-1000])
        model1['lb'][0, 0][2][0] = np.array([-1000])
        # equation matrix
        data = np.array([-1, -2, 1, -1, 2, -1, 1])
        row = np.array([0, 1, 2, 2, 3, 3, 4])
        col = np.array([0, 0, 0, 1, 1, 2, 2])
        model1['S'][0, 0] = csc_matrix((data, (row, col)), shape=(5, 3))
        # genes
        model1['grRules'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['grRules'][0, 0][0][0] = np.array(['g1 and g2'])
        model1['grRules'][0, 0][1][0] = np.array(['g3 or g4'])
        model1['grRules'][0, 0][2][0] = np.array([])
        # subsystems
        model1['subSystems'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['subSystems'][0, 0][0][0] = np.array(['system1'])
        model1['subSystems'][0, 0][1][0] = np.ndarray((1, 1), 'O')
        model1['subSystems'][0, 0][1][0][0][0] = np.array(['system2'])
        model1['subSystems'][0, 0][2][0] = np.array(['system3'])
        # compounds
        model1['mets'][0, 0] = np.ndarray((5, 1), dtype='O')
        model1['mets'][0, 0][0][0] = np.array(['A[e]'])
        model1['mets'][0, 0][1][0] = np.array(['B_c'])
        model1['mets'][0, 0][2][0] = np.array(['C(c)'])
        model1['mets'][0, 0][3][0] = np.array(['D[c]'])
        model1['mets'][0, 0][4][0] = np.array(['E001e'])
        model1['metFormulas'][0, 0] = np.ndarray((5, 1), dtype='O')
        model1['metFormulas'][0, 0][0][0] = np.array(['H2O'])
        model1['metFormulas'][0, 0][1][0] = np.array(['C5H10O2'])
        model1['metFormulas'][0, 0][2][0] = np.array(['H2SPO4'])
        model1['metFormulas'][0, 0][3][0] = np.array(['CH4N2'])
        model1['metFormulas'][0, 0][4][0] = np.array(['CO2'])
        model1['metNames'][0, 0] = np.ndarray((5, 1), dtype='O')
        model1['metNames'][0, 0][0][0] = np.array(['A'])
        model1['metNames'][0, 0][1][0] = np.array(['B'])
        model1['metNames'][0, 0][2][0] = np.array(['C'])
        model1['metNames'][0, 0][3][0] = np.array(['D'])
        model1['metNames'][0, 0][4][0] = np.array(['E'])
        model1['metCharges'][0, 0] = np.ndarray((5, 1), dtype='O')
        model1['metCharges'][0, 0][0][0] = 1
        model1['metCharges'][0, 0][1][0] = 2
        model1['metCharges'][0, 0][2][0] = -1.0
        model1['metCharges'][0, 0][3][0] = 0
        model1['metCharges'][0, 0][4][0] = 0
        model1['metNotes'][0, 0] = np.ndarray((5, 1), dtype='O')
        model1['metNotes'][0, 0][0][0] = np.array(['A'])
        model1['metNotes'][0, 0][1][0] = np.array(['B'])
        model1['metNotes'][0, 0][2][0] = np.array([])
        model1['metNotes'][0, 0][3][0] = np.array(['D'])
        model1['metNotes'][0, 0][4][0] = np.array(['E'])
        # biomass reaction
        model1['c'][0, 0] = np.ndarray((3, 1), dtype='O')
        model1['c'][0, 0][0][0] = 0
        model1['c'][0, 0][1][0] = 0
        model1['c'][0, 0][2][0] = 1

        os.mkdir(os.path.join(self.dest, 'models'))
        savemat(os.path.join(self.dest, 'model1.mat'), {'model1': model1})
        savemat(os.path.join(self.dest, 'models', 'model1.mat'),
                {'model1': model1})
        model2 = model1.copy()
        model2['c'][0, 0][1][0] = 1
        savemat(os.path.join(self.dest, 'models', 'model2.mat'),
                {'model1': model1, 'model2': model2})

    def tearDown(self):
        shutil.rmtree(self.dest)

    def test_create_importer(self):
        importer = MatlabImporter()
        importer.help()
        importer.import_model(os.path.join(self.dest, 'model1.mat'))
        importer.import_model(self.dest)
        importer.import_model(os.path.join(self.dest, 'models', 'model2.mat'),
                              'model2')
        with self.assertRaises(ModelLoadError):
            importer.import_model(os.path.join(self.dest, 'models'))
        with self.assertRaises(ModelLoadError):
            importer.import_model(os.path.join(
                self.dest, 'models', 'model2.mat'))

    def test_matlab_importer_main(self):
        output_dir = os.path.join(self.dest, 'output')
        os.mkdir(output_dir)
        importer.main(importer_class=MatlabImporter,
                      args=[
                          '--source', os.path.join(self.dest, 'model1.mat'),
                          '--dest', output_dir,
                      ])
