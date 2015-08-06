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

from psamm.datasource import sbml
from psamm.reaction import Reaction, Compound

from decimal import Decimal
from fractions import Fraction
from six import StringIO


class TestSBMLDatabaseL1V2(unittest.TestCase):
    """Test parsing of a simple level 1 version 2 SBML file"""

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level1"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="1" version="2">
 <model name="Test model">
  <listOfCompartments>
   <compartment name="cell"/>
  </listOfCompartments>
  <listOfSpecies>
   <species name="Glucose" compartment="cell" initialAmount="1"/>
   <species name="Glucose_6_P" compartment="cell" initialAmount="1"/>
   <species name="H2O" compartment="cell" initialAmount="1"/>
   <species name="Phosphate" compartment="cell" initialAmount="1" boundaryCondition="false"/>
   <species name="Biomass" compartment="cell" initialAmount="1" boundaryCondition="true"/>
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
     <speciesReference species="Glucose_6_P" stoichiometry="56" denominator="100"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="Biomass"/>
    </listOfProducts>
   </reaction>
  </listOfReactions>
 </model>
</sbml>''')
        self.reader = sbml.SBMLReader(s)

    def test_model_name(self):
        self.assertEqual(self.reader.name, 'Test model')

    def test_compounds_exist(self):
        species = {entry.id: entry for entry in self.reader.species}
        self.assertEqual(len(species), 5)

        self.assertEqual(species['Glucose'].id, 'Glucose')
        self.assertEqual(species['Glucose'].name, 'Glucose')
        self.assertEqual(species['Glucose'].compartment, 'cell')
        self.assertFalse(species['Glucose'].boundary)

        self.assertEqual(species['Glucose_6_P'].id, 'Glucose_6_P')
        self.assertEqual(species['Glucose_6_P'].name, 'Glucose_6_P')
        self.assertEqual(species['Glucose_6_P'].compartment, 'cell')
        self.assertFalse(species['Glucose_6_P'].boundary)

        self.assertEqual(species['H2O'].id, 'H2O')
        self.assertEqual(species['H2O'].name, 'H2O')
        self.assertEqual(species['H2O'].compartment, 'cell')
        self.assertFalse(species['H2O'].boundary)

        self.assertFalse(species['Phosphate'].boundary)
        self.assertTrue(species['Biomass'].boundary)

    def test_g6pase_reaction_exists(self):
        reaction = self.reader.get_reaction('G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Bidir,
                                   [(Compound('Glucose', 'cell'), 2),
                                    (Compound('Phosphate', 'cell'), 2)],
                                   [(Compound('H2O', 'cell'), 2),
                                    (Compound('Glucose_6_P', 'cell'), 2)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_biomass_reaction_exists(self):
        reaction = self.reader.get_reaction('Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Right,
                                   [(Compound('Glucose_6_P', 'cell'),
                                     Fraction(56, 100))],
                                   [(Compound('Biomass', 'cell'), 1)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_reaction_xml_notes(self):
        reaction = self.reader.get_reaction('G6Pase')
        notes = reaction.xml_notes

        notes_tags = list(notes)
        self.assertEqual(len(notes_tags), 1)
        self.assertEqual(notes_tags[0].tag, '{http://www.w3.org/1999/xhtml}p')
        self.assertEqual(notes_tags[0].text, 'Glucose 6-phosphatase')


class TestSBMLDatabaseL2V5(unittest.TestCase):
    """Test parsing of a simple level 2 version 5 SBML file"""

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version5"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="2" version="5">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose" name="Glucose" compartment="C_c"/>
   <species id="M_Glucose_6_P" name="Glucose-6-P" compartment="C_c"/>
   <species id="M_H2O" name="H2O" compartment="C_c"/>
   <species id="M_Phosphate" name="Phosphate" compartment="C_c" boundaryCondition="false"/>
   <species id="M_Biomass" name="Biomass" compartment="C_c" boundaryCondition="true"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true">
    <listOfReactants>
     <speciesReference species="M_Glucose" stoichiometry="2"/>
     <speciesReference species="M_Phosphate" stoichiometry="2"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_H2O" stoichiometry="2"/>
     <speciesReference species="M_Glucose_6_P" stoichiometry="2"/>
    </listOfProducts>
    <notes>
     <html:p>Glucose 6-phosphatase</html:p>
    </notes>
   </reaction>
   <reaction id="R_Biomass" reversible="false">
    <listOfReactants>
     <speciesReference species="M_Glucose_6_P" stoichiometry="0.56"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_Biomass"/>
    </listOfProducts>
   </reaction>
  </listOfReactions>
 </model>
</sbml>''')
        self.reader = sbml.SBMLReader(s)

    def test_model_name(self):
        self.assertEqual(self.reader.id, 'test_model')
        self.assertEqual(self.reader.name, 'Test model')

    def test_compounds_exist(self):
        species = {entry.id: entry for entry in self.reader.species}
        self.assertEqual(len(species), 5)

        self.assertEqual(species['M_Glucose'].id, 'M_Glucose')
        self.assertEqual(species['M_Glucose'].name, 'Glucose')
        self.assertEqual(species['M_Glucose'].compartment, 'C_c')
        self.assertFalse(species['M_Glucose'].boundary)

        self.assertEqual(species['M_Glucose_6_P'].id, 'M_Glucose_6_P')
        self.assertEqual(species['M_Glucose_6_P'].name, 'Glucose-6-P')
        self.assertEqual(species['M_Glucose_6_P'].compartment, 'C_c')
        self.assertFalse(species['M_Glucose_6_P'].boundary)

        self.assertEqual(species['M_H2O'].id, 'M_H2O')
        self.assertEqual(species['M_H2O'].name, 'H2O')
        self.assertEqual(species['M_H2O'].compartment, 'C_c')
        self.assertFalse(species['M_H2O'].boundary)

        self.assertFalse(species['M_Phosphate'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

    def test_g6pase_reaction_exists(self):
        reaction = self.reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Bidir,
                                   [(Compound('M_Glucose', 'C_c'), 2),
                                    (Compound('M_Phosphate', 'C_c'), 2)],
                                   [(Compound('M_H2O', 'C_c'), 2),
                                    (Compound('M_Glucose_6_P', 'C_c'), 2)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_biomass_reaction_exists(self):
        reaction = self.reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Right,
                                   [(Compound('M_Glucose_6_P', 'C_c'),
                                     Decimal('0.56'))],
                                   [(Compound('M_Biomass', 'C_c'), 1)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_reaction_xml_notes(self):
        reaction = self.reader.get_reaction('R_G6Pase')
        notes = reaction.xml_notes

        notes_tags = list(notes)
        self.assertEqual(len(notes_tags), 1)
        self.assertEqual(notes_tags[0].tag, '{http://www.w3.org/1999/xhtml}p')
        self.assertEqual(notes_tags[0].text, 'Glucose 6-phosphatase')


class TestSBMLDatabaseL3V1(unittest.TestCase):
    """Test parsing of a simple level 3 version 1 SBML file"""

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="3" version="1">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell" constant="true"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose" name="Glucose" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false"/>
   <species id="M_Glucose_6_P" name="Glucose-6-P" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false"/>
   <species id="M_H2O" name="H2O" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false"/>
   <species id="M_Phosphate" name="Phosphate" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false"/>
   <species id="M_Biomass" name="Biomass" compartment="C_c" constant="false" boundaryCondition="true" hasOnlySubstanceUnits="false"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true" fast="false">
    <listOfReactants>
     <speciesReference species="M_Glucose" stoichiometry="2" constant="true"/>
     <speciesReference species="M_Phosphate" stoichiometry="2" constant="true"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_H2O" stoichiometry="2" constant="true"/>
     <speciesReference species="M_Glucose_6_P" stoichiometry="2" constant="true"/>
    </listOfProducts>
    <notes>
     <html:p>Glucose 6-phosphatase</html:p>
    </notes>
   </reaction>
   <reaction id="R_Biomass" reversible="false" fast="false">
    <listOfReactants>
     <speciesReference species="M_Glucose_6_P" stoichiometry="0.56" constant="true"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_Biomass" stoichiometry="1" constant="true"/>
    </listOfProducts>
   </reaction>
  </listOfReactions>
 </model>
</sbml>''')
        self.reader = sbml.SBMLReader(s)

    def test_model_name(self):
        self.assertEqual(self.reader.id, 'test_model')
        self.assertEqual(self.reader.name, 'Test model')

    def test_compounds_exist(self):
        species = {entry.id: entry for entry in self.reader.species}
        self.assertEqual(len(species), 5)

        self.assertEqual(species['M_Glucose'].id, 'M_Glucose')
        self.assertEqual(species['M_Glucose'].name, 'Glucose')
        self.assertEqual(species['M_Glucose'].compartment, 'C_c')
        self.assertFalse(species['M_Glucose'].boundary)

        self.assertEqual(species['M_Glucose_6_P'].id, 'M_Glucose_6_P')
        self.assertEqual(species['M_Glucose_6_P'].name, 'Glucose-6-P')
        self.assertEqual(species['M_Glucose_6_P'].compartment, 'C_c')
        self.assertFalse(species['M_Glucose_6_P'].boundary)

        self.assertEqual(species['M_H2O'].id, 'M_H2O')
        self.assertEqual(species['M_H2O'].name, 'H2O')
        self.assertEqual(species['M_H2O'].compartment, 'C_c')
        self.assertFalse(species['M_H2O'].boundary)

        self.assertFalse(species['M_Phosphate'].boundary)
        self.assertTrue(species['M_Biomass'].boundary)

    def test_g6pase_reaction_exists(self):
        reaction = self.reader.get_reaction('R_G6Pase')
        self.assertTrue(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Bidir,
                                   [(Compound('M_Glucose', 'C_c'), 2),
                                    (Compound('M_Phosphate', 'C_c'), 2)],
                                   [(Compound('M_H2O', 'C_c'), 2),
                                    (Compound('M_Glucose_6_P', 'C_c'), 2)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_biomass_reaction_exists(self):
        reaction = self.reader.get_reaction('R_Biomass')
        self.assertFalse(reaction.reversible)

        # Compare equation of reaction
        actual_equation = Reaction(Reaction.Right,
                                   [(Compound('M_Glucose_6_P', 'C_c'),
                                     Decimal('0.56'))],
                                   [(Compound('M_Biomass', 'C_c'), 1)])
        self.assertEqual(reaction.equation, actual_equation)

    def test_reaction_xml_notes(self):
        reaction = self.reader.get_reaction('R_G6Pase')
        notes = reaction.xml_notes

        notes_tags = list(notes)
        self.assertEqual(len(notes_tags), 1)
        self.assertEqual(notes_tags[0].tag, '{http://www.w3.org/1999/xhtml}p')
        self.assertEqual(notes_tags[0].text, 'Glucose 6-phosphatase')


if __name__ == '__main__':
    unittest.main()
