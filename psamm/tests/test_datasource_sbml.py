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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import unittest

from psamm.datasource import sbml, native, entry
from psamm.datasource.reaction import parse_reaction
from psamm.reaction import Reaction, Compound, Direction

from decimal import Decimal
from fractions import Fraction
from six import BytesIO, itervalues


class TestSBMLDatabaseL1V2(unittest.TestCase):
    """Test parsing of a simple level 1 version 2 SBML file"""

    def setUp(self):
        self.doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level1"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="1" version="2">
 <model name="Test model">
  <listOfCompartments>
   <compartment name="boundary"/>
   <compartment name="cell"/>
  </listOfCompartments>
  <listOfSpecies>
   <species name="Glucose" compartment="cell" initialAmount="1"
            charge="0"/>
   <species name="Glucose_6_P" compartment="cell" initialAmount="1"
            charge="-2"/>
   <species name="H2O" compartment="cell" initialAmount="1" charge="0"/>
   <species name="Phosphate" compartment="cell" initialAmount="1"
            boundaryCondition="false"/>
   <species name="Biomass" compartment="boundary" initialAmount="1"
            boundaryCondition="true"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction name="G6Pase" reversible="true">
    <listOfReactants>
     <speciesReference species="Glucose" stoichiometry="2"/>
     <speciesReference species="Phosphate" stoichiometry="2"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="H2O" stoichiometry="2"/>
     <speciesReference species="Glucose_6_P" stoichiometry="2"/>
    </listOfProducts>
    <notes>
     <html:p>Glucose 6-phosphatase</html:p>
    </notes>
   </reaction>
   <reaction name="Biomass" reversible="false">
    <listOfReactants>
     <speciesReference species="Glucose_6_P" stoichiometry="56"
                       denominator="100"/>
     <speciesReference species="Glucose" stoichiometry="88"
                       denominator="100"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="Biomass"/>
    </listOfProducts>
   </reaction>
  </listOfReactions>
 </model>
</sbml>'''.encode('utf-8'))

    def test_model_name(self):
        reader = sbml.SBMLReader(self.doc)
        self.assertEqual(reader.name, 'Test model')

    def test_compartment_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 2)
        self.assertEqual(compartments['cell'].id, 'cell')
        self.assertEqual(compartments['cell'].name, 'cell')
        self.assertEqual(compartments['boundary'].id, 'boundary')
        self.assertEqual(compartments['boundary'].name, 'boundary')

    def test_compartment_exists_with_ignore_boundary(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=True)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 1)
        self.assertEqual(compartments['cell'].id, 'cell')
        self.assertEqual(compartments['cell'].name, 'cell')

    def test_compounds_exist(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        species = {entry.id: entry for entry in reader.species}
        self.assertEqual(len(species), 5)

        self.assertEqual(species['Glucose'].id, 'Glucose')
        self.assertEqual(species['Glucose'].name, 'Glucose')
        self.assertEqual(species['Glucose'].compartment, 'cell')
        self.assertFalse(species['Glucose'].boundary)
        self.assertEqual(species['Glucose'].charge, 0)

        self.assertEqual(species['Glucose_6_P'].id, 'Glucose_6_P')
        self.assertEqual(species['Glucose_6_P'].name, 'Glucose_6_P')
        self.assertEqual(species['Glucose_6_P'].compartment, 'cell')
        self.assertFalse(species['Glucose_6_P'].boundary)
        self.assertEqual(species['Glucose_6_P'].charge, -2)

        self.assertEqual(species['H2O'].id, 'H2O')
        self.assertEqual(species['H2O'].name, 'H2O')
        self.assertEqual(species['H2O'].compartment, 'cell')
        self.assertFalse(species['H2O'].boundary)
        self.assertEqual(species['H2O'].charge, 0)

        self.assertFalse(species['Phosphate'].boundary)
        self.assertTrue(species['Biomass'].boundary)

    def test_g6pase_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Both,
                                   [(Compound('Glucose', 'cell'), 2),
                                    (Compound('Phosphate', 'cell'), 2)],
                                   [(Compound('H2O', 'cell'), 2),
                                    (Compound('Glucose_6_P', 'cell'), 2)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_biomass_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Forward, [
            (Compound('Glucose_6_P', 'cell'), Fraction(56, 100)),
            (Compound('Glucose', 'cell'), Fraction(88, 100))
        ], [
            (Compound('Biomass', 'boundary'), 1)
        ])
        self.assertEqual(reaction.equation, actual_equation)

    def test_reaction_xml_notes(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('G6Pase')
        notes = reaction.xml_notes

        notes_tags = list(notes)
        self.assertEqual(len(notes_tags), 1)
        self.assertEqual(notes_tags[0].tag, '{http://www.w3.org/1999/xhtml}p')
        self.assertEqual(notes_tags[0].text, 'Glucose 6-phosphatase')

    def test_objective_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = list(reader.objectives)
        self.assertEqual(len(objectives), 0)
        self.assertIsNone(reader.get_active_objective())

    def test_flux_bounds_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        flux_bounds = list(reader.flux_bounds)
        self.assertEqual(len(flux_bounds), 0)

    def test_create_and_convert_model(self):
        reader = sbml.SBMLReader(self.doc)
        model = reader.create_model()
        sbml.convert_sbml_model(model)

        self.assertEqual(
            {entry.id for entry in model.compounds},
            {'Glucose', 'Glucose_6_P', 'H2O', 'Phosphate'})
        self.assertEqual(
            {entry.id for entry in model.reactions},
            {'G6Pase', 'Biomass'})
        self.assertEqual(
            {entry.id for entry in model.compartments},
            {'cell'})

        self.assertEqual(set(model.model), {'Biomass', 'G6Pase'})


class TestSBMLDatabaseL2V5(unittest.TestCase):
    """Test parsing of a simple level 2 version 5 SBML file"""

    def setUp(self):
        self.doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version5"
      level="2" version="5">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell"/>
   <compartment id="C_b" name="boundary"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose_LPAREN_c_RPAREN_" name="Glucose" compartment="C_c"
            charge="0">
    <notes>
     <body xmlns="http://www.w3.org/1999/xhtml">
      <p>FORMULA: C6H12O6</p>
      <p>Charge: "0"</p>
      <p>Additional notes..</p>
      <p>KEGG ID: C00031</p>
      <p>Custom: 123</p>
     </body>
    </notes>
   </species>
   <species id="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_" name="Glucose-6-P"
            compartment="C_c" charge="-2"/>
   <species id="M_H2O_LPAREN_c_RPAREN_" name="H2O" compartment="C_c"
            charge="0"/>
   <species id="M_Phosphate_LPAREN_c_RPAREN_" name="Phosphate"
            compartment="C_c" boundaryCondition="false"/>
   <species id="M_Biomass" name="Biomass" compartment="C_b"
            boundaryCondition="true"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true">
    <listOfReactants>
     <speciesReference species="M_Glucose_LPAREN_c_RPAREN_" stoichiometry="2"/>
     <speciesReference species="M_Phosphate_LPAREN_c_RPAREN_"
                       stoichiometry="2"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_H2O_LPAREN_c_RPAREN_" stoichiometry="2"/>
     <speciesReference species="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_"
                       stoichiometry="2"/>
    </listOfProducts>
    <notes>
     <body xmlns="http://www.w3.org/1999/xhtml">
      <p> Authors: Jane Doe ; John Doe </p>
      <p>CONFIDENCE: 3</p>
      <p>EC NUMBER: 3.1.3.9</p>
      <p>gene association :   b0822 </p>
      <p>SUBSYSTEM: "Glycolysis / Gluconeogenesis"</p>
      <p>Additional notes...</p>
     </body>
    </notes>
   </reaction>
   <reaction id="R_Biomass" reversible="false">
    <listOfReactants>
     <speciesReference species="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_"
                       stoichiometry="0.56"/>
     <speciesReference species="M_Glucose_LPAREN_c_RPAREN_"
                       stoichiometry="0.88"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_Biomass"/>
    </listOfProducts>
    <kineticLaw>
     <listOfParameters>
      <parameter id="LOWER_BOUND" value="0"/>
      <parameter id="UPPER_BOUND" value="1000"/>
      <parameter id="SOME_CUSTOM_PARAMETER" value="123.4"/>
      <parameter id="OBJECTIVE_COEFFICIENT" value="1"/>
     </listOfParameters>
    </kineticLaw>
   </reaction>
  </listOfReactions>
 </model>
</sbml>'''.encode('utf-8'))

    def test_model_name(self):
        reader = sbml.SBMLReader(self.doc)
        self.assertEqual(reader.id, 'test_model')
        self.assertEqual(reader.name, 'Test model')

    def test_compartment_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 2)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')
        self.assertEqual(compartments['C_b'].id, 'C_b')
        self.assertEqual(compartments['C_b'].name, 'boundary')

    def test_compartment_exists_with_ignore_boundary(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=True)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 1)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')

    def test_compounds_exist(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        species = {entry.id: entry for entry in reader.species}
        self.assertEqual(len(species), 5)

        gluc_species = species['M_Glucose_LPAREN_c_RPAREN_']
        self.assertEqual(gluc_species.id, 'M_Glucose_LPAREN_c_RPAREN_')
        self.assertEqual(gluc_species.name, 'Glucose')
        self.assertEqual(gluc_species.compartment, 'C_c')
        self.assertFalse(gluc_species.boundary)
        self.assertEqual(gluc_species.charge, 0)

        g6p_species = species['M_Glucose_6_DASH_P_LPAREN_c_RPAREN_']
        self.assertEqual(g6p_species.id, 'M_Glucose_6_DASH_P_LPAREN_c_RPAREN_')
        self.assertEqual(g6p_species.name, 'Glucose-6-P')
        self.assertEqual(g6p_species.compartment, 'C_c')
        self.assertFalse(g6p_species.boundary)
        self.assertEqual(g6p_species.charge, -2)

        h2o_species = species['M_H2O_LPAREN_c_RPAREN_']
        self.assertEqual(h2o_species.id, 'M_H2O_LPAREN_c_RPAREN_')
        self.assertEqual(h2o_species.name, 'H2O')
        self.assertEqual(h2o_species.compartment, 'C_c')
        self.assertFalse(h2o_species.boundary)
        self.assertEqual(h2o_species.charge, 0)

        self.assertFalse(species['M_Phosphate_LPAREN_c_RPAREN_'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

    def test_glucose_parse_notes(self):
        reader = sbml.SBMLReader(self.doc)
        species = reader.get_species('M_Glucose_LPAREN_c_RPAREN_')
        notes_dict = sbml.parse_xhtml_species_notes(species)
        self.assertEqual(notes_dict, {
            'formula': 'C6H12O6',
            'kegg': 'C00031',
            'charge': 0
        })

    def test_g6pase_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Both, [
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Phosphate_LPAREN_c_RPAREN_', 'C_c'), 2)
        ], [
            (Compound('M_H2O_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'), 2)
        ])
        self.assertEqual(reaction.equation, actual_equation)

    def test_g6pase_parse_notes(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        notes = sbml.parse_xhtml_reaction_notes(reaction)
        self.assertEqual(notes, {
            'subsystem': 'Glycolysis / Gluconeogenesis',
            'genes': 'b0822',
            'ec': '3.1.3.9',
            'confidence': 3,
            'authors': ['Jane Doe', 'John Doe']
        })

    def test_parse_reaction_cobra_flux_bounds(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        lower, upper = sbml.parse_flux_bounds(reaction)
        self.assertIsNone(lower)
        self.assertIsNone(upper)

        reaction = reader.get_reaction('R_Biomass')
        lower, upper = sbml.parse_flux_bounds(reaction)
        self.assertEqual(lower, 0)
        self.assertEqual(upper, 1000)

    def test_parse_reaction_cobra_objective(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        coeff = sbml.parse_objective_coefficient(reaction)
        self.assertIsNone(coeff)

        reaction = reader.get_reaction('R_Biomass')
        coeff = sbml.parse_objective_coefficient(reaction)
        self.assertEqual(coeff, 1)

    def test_biomass_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Forward, [
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'),
             Decimal('0.56')),
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), Decimal('0.88'))
        ], [
            (Compound('M_Biomass', 'C_b'), 1)
        ])
        self.assertEqual(reaction.equation, actual_equation)

    def test_objective_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = list(reader.objectives)
        self.assertEqual(len(objectives), 0)
        self.assertIsNone(reader.get_active_objective())

    def test_flux_bounds_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        flux_bounds = list(reader.flux_bounds)
        self.assertEqual(len(flux_bounds), 0)

    def test_create_and_convert_model(self):
        reader = sbml.SBMLReader(self.doc)
        model = reader.create_model()
        sbml.convert_sbml_model(model)

        self.assertEqual(
            {entry.id for entry in model.compounds},
            {'Glucose(c)', 'Glucose_6-P(c)', 'H2O(c)', 'Phosphate(c)'})
        self.assertEqual(
            {entry.id for entry in model.reactions},
            {'G6Pase', 'Biomass'})
        self.assertEqual(
            {entry.id for entry in model.compartments},
            {'c'})

        self.assertEqual(model.limits['Biomass'], ('Biomass', 0, 1000))
        self.assertEqual(model.biomass_reaction, 'Biomass')
        self.assertEqual(set(model.model), {'Biomass', 'G6Pase'})


class TestSBMLDatabaseL3V1(unittest.TestCase):
    """Test parsing of a simple level 3 version 1 SBML file"""

    def setUp(self):
        self.doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="3" version="1">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell" constant="true"/>
   <compartment id="C_b" name="boundary" constant="true"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose_LPAREN_c_RPAREN_" name="Glucose" compartment="C_c"
            constant="false" boundaryCondition="false"
            hasOnlySubstanceUnits="false"/>
   <species id="M_Glucose_6_DASH_P_LPAREN_c_RPAREN_" name="Glucose-6-P"
            compartment="C_c" constant="false" boundaryCondition="false"
            hasOnlySubstanceUnits="false"/>
   <species id="M_H2O_LPAREN_c_RPAREN_" name="H2O" compartment="C_c"
            constant="false" boundaryCondition="false"
            hasOnlySubstanceUnits="false"/>
   <species id="M_Phosphate_LPAREN_c_RPAREN_" name="Phosphate"
            compartment="C_c" constant="false" boundaryCondition="false"
            hasOnlySubstanceUnits="false"/>
   <species id="M_Biomass" name="Biomass" compartment="C_b" constant="false"
            boundaryCondition="true" hasOnlySubstanceUnits="false"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true" fast="false">
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
     <html:p>Glucose 6-phosphatase</html:p>
    </notes>
   </reaction>
   <reaction id="R_Biomass" reversible="false" fast="false">
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
 </model>
</sbml>'''.encode('utf-8'))

    def test_model_name(self):
        reader = sbml.SBMLReader(self.doc)
        self.assertEqual(reader.id, 'test_model')
        self.assertEqual(reader.name, 'Test model')

    def test_compartment_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 2)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')
        self.assertEqual(compartments['C_b'].id, 'C_b')
        self.assertEqual(compartments['C_b'].name, 'boundary')

    def test_compartment_exists_with_ignore_boundary(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=True)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 1)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')

    def test_compounds_exist(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        species = {entry.id: entry for entry in reader.species}
        self.assertEqual(len(species), 5)

        gluc_species = species['M_Glucose_LPAREN_c_RPAREN_']
        self.assertEqual(gluc_species.id, 'M_Glucose_LPAREN_c_RPAREN_')
        self.assertEqual(gluc_species.name, 'Glucose')
        self.assertEqual(gluc_species.compartment, 'C_c')
        self.assertFalse(gluc_species.boundary)

        g6p_species = species['M_Glucose_6_DASH_P_LPAREN_c_RPAREN_']
        self.assertEqual(g6p_species.id, 'M_Glucose_6_DASH_P_LPAREN_c_RPAREN_')
        self.assertEqual(g6p_species.name, 'Glucose-6-P')
        self.assertEqual(g6p_species.compartment, 'C_c')
        self.assertFalse(g6p_species.boundary)

        h2o_species = species['M_H2O_LPAREN_c_RPAREN_']
        self.assertEqual(h2o_species.id, 'M_H2O_LPAREN_c_RPAREN_')
        self.assertEqual(h2o_species.name, 'H2O')
        self.assertEqual(h2o_species.compartment, 'C_c')
        self.assertFalse(h2o_species.boundary)

        self.assertFalse(species['M_Phosphate_LPAREN_c_RPAREN_'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

    def test_g6pase_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Both, [
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Phosphate_LPAREN_c_RPAREN_', 'C_c'), 2)
        ], [
            (Compound('M_H2O_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'), 2)
        ])
        self.assertEqual(reaction.equation, actual_equation)

    def test_biomass_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Forward, [
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'),
             Decimal('0.56')),
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), Decimal('0.88'))
        ], [
            (Compound('M_Biomass', 'C_b'), 1)
        ])
        self.assertEqual(reaction.equation, actual_equation)

    def test_reaction_xml_notes(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        notes = reaction.xml_notes

        notes_tags = list(notes)
        self.assertEqual(len(notes_tags), 1)
        self.assertEqual(notes_tags[0].tag, '{http://www.w3.org/1999/xhtml}p')
        self.assertEqual(notes_tags[0].text, 'Glucose 6-phosphatase')

    def test_objective_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = list(reader.objectives)
        self.assertEqual(len(objectives), 0)
        self.assertIsNone(reader.get_active_objective())

    def test_flux_bounds_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        flux_bounds = list(reader.flux_bounds)
        self.assertEqual(len(flux_bounds), 0)

    def test_create_and_convert_model(self):
        reader = sbml.SBMLReader(self.doc)
        model = reader.create_model()
        sbml.convert_sbml_model(model)

        self.assertEqual(
            {entry.id for entry in model.compounds},
            {'Glucose(c)', 'Glucose_6-P(c)', 'H2O(c)', 'Phosphate(c)'})
        self.assertEqual(
            {entry.id for entry in model.reactions},
            {'G6Pase', 'Biomass'})
        self.assertEqual(
            {entry.id for entry in model.compartments},
            {'c'})

        self.assertEqual(set(model.model), {'Biomass', 'G6Pase'})


class TestSBMLDatabaseL3V1WithFBCV1(unittest.TestCase):
    """Test parsing of a level 3 version 1 SBML file with FBC version 1"""

    def setUp(self):
        self.doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version1"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="3" version="1"
      fbc:required="false">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell" constant="true"/>
   <compartment id="C_b" name="boundary" constant="true"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose_LPAREN_c_RPAREN_" name="Glucose" compartment="C_c"
            constant="false" boundaryCondition="false"
            hasOnlySubstanceUnits="false" fbc:charge="0"
            fbc:chemicalFormula="C6H12O6"/>
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
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true" fast="false">
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
   </reaction>
   <reaction id="R_Biomass" reversible="false" fast="false">
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
  <fbc:listOfFluxBounds>
   <fbc:fluxBound fbc:reaction="R_G6Pase" fbc:operation="greaterEqual"
                  fbc:value="-10"/>
   <fbc:fluxBound fbc:reaction="R_G6Pase" fbc:operation="lessEqual"
                  fbc:value="1000"/>
   <fbc:fluxBound fbc:reaction="R_Biomass" fbc:operation="greaterEqual"
                  fbc:value="0"/>
   <fbc:fluxBound fbc:reaction="R_Biomass" fbc:operation="lessEqual"
                  fbc:value="1000"/>
  </fbc:listOfFluxBounds>
 </model>
</sbml>'''.encode('utf-8'))

    def test_model_name(self):
        reader = sbml.SBMLReader(self.doc)
        self.assertEqual(reader.id, 'test_model')
        self.assertEqual(reader.name, 'Test model')

    def test_compartment_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 2)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')
        self.assertEqual(compartments['C_b'].id, 'C_b')
        self.assertEqual(compartments['C_b'].name, 'boundary')

    def test_compounds_exist(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        species = {entry.id: entry for entry in reader.species}
        self.assertEqual(len(species), 5)

        gluc_species = species['M_Glucose_LPAREN_c_RPAREN_']
        self.assertEqual(gluc_species.id, 'M_Glucose_LPAREN_c_RPAREN_')
        self.assertEqual(gluc_species.name, 'Glucose')
        self.assertEqual(gluc_species.compartment, 'C_c')
        self.assertFalse(gluc_species.boundary)
        self.assertEqual(gluc_species.formula, 'C6H12O6')
        self.assertEqual(gluc_species.charge, 0)

        g6p_species = species['M_Glucose_6_DASH_P_LPAREN_c_RPAREN_']
        self.assertEqual(g6p_species.id, 'M_Glucose_6_DASH_P_LPAREN_c_RPAREN_')
        self.assertEqual(g6p_species.name, 'Glucose-6-P')
        self.assertEqual(g6p_species.compartment, 'C_c')
        self.assertFalse(g6p_species.boundary)
        self.assertEqual(g6p_species.formula, 'C6H11O9P')
        self.assertEqual(g6p_species.charge, -2)

        h2o_species = species['M_H2O_LPAREN_c_RPAREN_']
        self.assertEqual(h2o_species.id, 'M_H2O_LPAREN_c_RPAREN_')
        self.assertEqual(h2o_species.name, 'H2O')
        self.assertEqual(h2o_species.compartment, 'C_c')
        self.assertFalse(h2o_species.boundary)
        self.assertEqual(h2o_species.formula, 'H2O')
        self.assertEqual(h2o_species.charge, 0)

        self.assertFalse(species['M_Phosphate_LPAREN_c_RPAREN_'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

        self.assertIsNone(species['M_Biomass'].formula)
        self.assertIsNone(species['M_Biomass'].charge)

    def test_g6pase_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Both, [
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Phosphate_LPAREN_c_RPAREN_', 'C_c'), 2)
        ], [
            (Compound('M_H2O_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'), 2)
        ])
        self.assertEqual(reaction.equation, actual_equation)

        self.assertEqual(reaction.properties['lower_flux'], -10)
        self.assertEqual(reaction.properties['upper_flux'], 1000)

    def test_biomass_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Forward, [
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'),
             Decimal('0.56')),
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), Decimal('0.88'))
        ], [
            (Compound('M_Biomass', 'C_b'), 1)
        ])
        self.assertEqual(reaction.equation, actual_equation)

        self.assertEqual(reaction.properties['lower_flux'], 0)
        self.assertEqual(reaction.properties['upper_flux'], 1000)

    def test_objective_exists(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = {entry.id: entry for entry in reader.objectives}
        self.assertEqual(len(objectives), 1)

        objective = objectives['obj1']
        self.assertEqual(objective.name, 'Objective 1')
        self.assertEqual(objective.type, 'maximize')
        self.assertEqual(dict(objective.reactions), {'R_Biomass': 1})

    def test_active_objective(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = {entry.id: entry for entry in reader.objectives}
        self.assertEqual(reader.get_active_objective(), objectives['obj1'])

    def test_flux_bounds_exists(self):
        reader = sbml.SBMLReader(self.doc)
        flux_bounds = list(reader.flux_bounds)
        self.assertEqual(len(flux_bounds), 4)

        biomass_bounds = set(
            (b.operation, b.value) for b in flux_bounds
            if b.reaction == 'R_Biomass')
        self.assertEqual(biomass_bounds, {
            ('greaterEqual', 0),
            ('lessEqual', 1000)})

        g6pase_bounds = set(
            (b.operation, b.value) for b in flux_bounds
            if b.reaction == 'R_G6Pase')
        self.assertEqual(g6pase_bounds, {
            ('greaterEqual', -10),
            ('lessEqual', 1000)})

    def test_create_and_convert_model(self):
        reader = sbml.SBMLReader(self.doc)
        model = reader.create_model()
        sbml.convert_sbml_model(model)

        self.assertEqual(
            {entry.id for entry in model.compounds},
            {'Glucose(c)', 'Glucose_6-P(c)', 'H2O(c)', 'Phosphate(c)'})
        self.assertEqual(
            {entry.id for entry in model.reactions},
            {'G6Pase', 'Biomass'})
        self.assertEqual(
            {entry.id for entry in model.compartments},
            {'c'})

        self.assertEqual(model.limits['Biomass'], ('Biomass', 0, 1000))
        self.assertEqual(model.limits['G6Pase'], ('G6Pase', -10, 1000))
        self.assertEqual(set(model.model), {'Biomass', 'G6Pase'})
        self.assertEqual(model.biomass_reaction, 'Biomass')


class TestSBMLDatabaseL3V1WithFBCV2(unittest.TestCase):
    """Test parsing of a level 3 version 1 SBML file with FBC version 2"""

    def setUp(self):
        self.doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
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

    def test_model_name(self):
        reader = sbml.SBMLReader(self.doc)
        self.assertEqual(reader.id, 'test_model')
        self.assertEqual(reader.name, 'Test model')

    def test_compartment_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        compartments = {entry.id: entry for entry in reader.compartments}
        self.assertEqual(len(compartments), 2)
        self.assertEqual(compartments['C_c'].id, 'C_c')
        self.assertEqual(compartments['C_c'].name, 'cell')
        self.assertEqual(compartments['C_b'].id, 'C_b')
        self.assertEqual(compartments['C_b'].name, 'boundary')

    def test_compounds_exist(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        species = {entry.id: entry for entry in reader.species}
        self.assertEqual(len(species), 5)

        gluc_species = species['M_Glucose_LPAREN_c_RPAREN_']
        self.assertEqual(gluc_species.id, 'M_Glucose_LPAREN_c_RPAREN_')
        self.assertEqual(gluc_species.name, 'Glucose')
        self.assertEqual(gluc_species.compartment, 'C_c')
        self.assertFalse(gluc_species.boundary)
        self.assertEqual(gluc_species.formula, 'C6H12O6')
        self.assertEqual(gluc_species.charge, 0)
        self.assertIsNotNone(gluc_species.xml_notes)

        g6p_species = species['M_Glucose_6_DASH_P_LPAREN_c_RPAREN_']
        self.assertEqual(g6p_species.id, 'M_Glucose_6_DASH_P_LPAREN_c_RPAREN_')
        self.assertEqual(g6p_species.name, 'Glucose-6-P')
        self.assertEqual(g6p_species.compartment, 'C_c')
        self.assertFalse(g6p_species.boundary)
        self.assertEqual(g6p_species.formula, 'C6H11O9P')
        self.assertEqual(g6p_species.charge, -2)
        self.assertIsNone(g6p_species.xml_notes)

        h2o_species = species['M_H2O_LPAREN_c_RPAREN_']
        self.assertEqual(h2o_species.id, 'M_H2O_LPAREN_c_RPAREN_')
        self.assertEqual(h2o_species.name, 'H2O')
        self.assertEqual(h2o_species.compartment, 'C_c')
        self.assertFalse(h2o_species.boundary)
        self.assertEqual(h2o_species.formula, 'H2O')
        self.assertEqual(h2o_species.charge, 0)

        self.assertFalse(species['M_Phosphate_LPAREN_c_RPAREN_'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

        self.assertIsNone(species['M_Biomass'].formula)
        self.assertIsNone(species['M_Biomass'].charge)

    def test_g6pase_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc)
        reaction = reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Both, [
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Phosphate_LPAREN_c_RPAREN_', 'C_c'), 2)
        ], [
            (Compound('M_H2O_LPAREN_c_RPAREN_', 'C_c'), 2),
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'), 2)
        ])
        self.assertEqual(reaction.equation, actual_equation)

        self.assertEqual(reaction.properties['lower_flux'], -10)
        self.assertEqual(reaction.properties['upper_flux'], 1000)

        self.assertIsNotNone(reaction.xml_notes)
        self.assertIsNotNone(reaction.xml_annotation)

    def test_biomass_reaction_exists(self):
        reader = sbml.SBMLReader(self.doc, ignore_boundary=False)
        reaction = reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Direction.Forward, [
            (Compound('M_Glucose_6_DASH_P_LPAREN_c_RPAREN_', 'C_c'),
             Decimal('0.56')),
            (Compound('M_Glucose_LPAREN_c_RPAREN_', 'C_c'), Decimal('0.88'))
        ], [
            (Compound('M_Biomass', 'C_b'), 1)
        ])
        self.assertEqual(reaction.equation, actual_equation)

        self.assertEqual(reaction.properties['lower_flux'], 0)
        self.assertEqual(reaction.properties['upper_flux'], 1000)

        self.assertIsNone(reaction.xml_notes)
        self.assertIsNone(reaction.xml_annotation)

    def test_objective_exists(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = {entry.id: entry for entry in reader.objectives}
        self.assertEqual(len(objectives), 1)

        objective = objectives['obj1']
        self.assertEqual(objective.name, 'Objective 1')
        self.assertEqual(objective.type, 'maximize')
        self.assertEqual(dict(objective.reactions), {'R_Biomass': 1})

    def test_active_objective(self):
        reader = sbml.SBMLReader(self.doc)
        objectives = {entry.id: entry for entry in reader.objectives}
        self.assertEqual(reader.get_active_objective(),
                         objectives['obj1'])

    def test_flux_bounds_not_present(self):
        reader = sbml.SBMLReader(self.doc)
        flux_bounds = list(reader.flux_bounds)
        self.assertEqual(len(flux_bounds), 0)

    def test_create_and_convert_model(self):
        reader = sbml.SBMLReader(self.doc)
        model = reader.create_model()
        sbml.convert_sbml_model(model)

        self.assertEqual(
            {entry.id for entry in model.compounds},
            {'Glucose(c)', 'Glucose_6-P(c)', 'H2O(c)', 'Phosphate(c)'})
        self.assertEqual(
            {entry.id for entry in model.reactions},
            {'G6Pase', 'Biomass'})
        self.assertEqual(
            {entry.id for entry in model.compartments},
            {'c'})

        self.assertEqual(model.limits['Biomass'], ('Biomass', 0, 1000))
        self.assertEqual(model.limits['G6Pase'], ('G6Pase', -10, 1000))
        self.assertEqual(set(model.model), {'Biomass', 'G6Pase'})
        self.assertEqual(model.biomass_reaction, 'Biomass')


class TestModelExtracellularCompartment(unittest.TestCase):
    def setUp(self):
        self.model = native.NativeModel()
        self.model.reactions.update([
            entry.DictReactionEntry({
                'id': 'EX_g6p',
                'equation': parse_reaction('g6p[e] <=>')
            }),
            entry.DictReactionEntry({
                'id': 'EX_glc-D',
                'equation': parse_reaction('glc-D[e] <=>')
            }),
            entry.DictReactionEntry({
                'id': 'SINK_glc-D',
                'equation': parse_reaction('glc-D[c] <=>')
            }),
            entry.DictReactionEntry({
                'id': 'rxn_1',
                'equation': parse_reaction('g6p[c] <=> glc-D[c]')
            }),
            entry.DictReactionEntry({
                'id': 'TP_glc-D',
                'equation': parse_reaction('glc-D[c] <=> glc-D[e]')
            })
        ])

    def test_detect_model(self):
        """Test that the extracellular compartment is detected.

        Since the 'e' compartment occurs more often in exchange reactions,
        this compartment should be detected.
        """
        extracellular = sbml.detect_extracellular_compartment(self.model)
        self.assertEqual(extracellular, 'e')

    def test_convert_model(self):
        self.model.extracellular_compartment = 'e'
        sbml.convert_exchange_to_compounds(self.model)

        self.assertEqual(
            {entry.id for entry in self.model.reactions},
            {'rxn_1', 'TP_glc-D', 'SINK_glc-D'})
        self.assertEqual(
            set(itervalues(self.model.exchange)),
            {
                (Compound('g6p', 'e'), 'EX_g6p', None, None),
                (Compound('glc-D', 'e'), 'EX_glc-D', None, None)
            })


class TestMergeEquivalentCompounds(unittest.TestCase):
    def test_merge(self):
        model = native.NativeModel()
        model.compounds.update([
            entry.DictCompoundEntry({
                'id': 'g6p_c',
                'name': 'G6P'
            }),
            entry.DictCompoundEntry({
                'id': 'g6p_e',
                'name': 'G6P'
            })
        ])
        model.reactions.update([
            entry.DictReactionEntry({
                'id': 'TP_g6p',
                'equation': parse_reaction('g6p_c[c] <=> g6p_e[e]')
            })
        ])
        exchange_compound = Compound('g6p_e', 'e')
        model.exchange[exchange_compound] = (
            exchange_compound, 'EX_g6p_e', -10, 0)

        sbml.merge_equivalent_compounds(model)

        self.assertEqual({entry.id for entry in model.compounds}, {'g6p'})
        self.assertEqual(
            model.reactions['TP_g6p'].equation,
            parse_reaction('g6p[c] <=> g6p[e]'))

        new_exchange_compound = Compound('g6p', 'e')
        self.assertEqual(
            model.exchange[new_exchange_compound],
            (new_exchange_compound, 'EX_g6p_e', -10, 0))


if __name__ == '__main__':
    unittest.main()
