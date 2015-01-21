#!/usr/bin/env python

import unittest

from metnet.datasource import sbml
from metnet.reaction import Compound

from decimal import Decimal
from fractions import Fraction
from StringIO import StringIO

class TestSBMLDatabaseL1V2(unittest.TestCase):
    '''Test parsing of a simple level 1 version 2 SBML file'''

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level1" level="1" version="2">
 <model>
  <listOfCompartments>
   <compartment name="cell"/>
  </listOfCompartments>
  <listOfSpecies>
   <species name="Glucose" compartment="cell" initialAmount="1"/>
   <species name="Glucose_6_P" compartment="cell" initialAmount="1"/>
   <species name="H2O" compartment="cell" initialAmount="1"/>
   <species name="Phosphate" compartment="cell" initialAmount="1"/>
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
        self.database = sbml.SBMLDatabase(s)

    def test_compounds_exist(self):
        self.assertEquals(len(list(self.database.compounds)), 5)
        self.assertIn(Compound('Glucose', 'cell'), self.database.compounds)
        self.assertIn(Compound('Glucose_6_P', 'cell'), self.database.compounds)
        self.assertIn(Compound('H2O', 'cell'), self.database.compounds)
        self.assertIn(Compound('Phosphate', 'cell'), self.database.compounds)
        self.assertIn(Compound('Biomass', 'cell'), self.database.compounds)

    def test_g6pase_reaction_exists(self):
        self.assertTrue(self.database.is_reversible('G6Pase'))

        # Read stoichiometric valus of reaction
        values = dict(self.database.get_reaction_values('G6Pase'))
        self.assertEquals(values, { Compound('Glucose', 'cell'): -2,
                                    Compound('Phosphate', 'cell'): -2,
                                    Compound('Glucose_6_P', 'cell'): 2,
                                    Compound('H2O', 'cell'): 2 })

    def test_biomass_reaction_exists(self):
        self.assertFalse(self.database.is_reversible('Biomass'))

        # Read stoichiometric values of reaction
        values = dict(self.database.get_reaction_values('Biomass'))
        self.assertEquals(values, { Compound('Glucose_6_P', 'cell'): -Fraction(56, 100),
                                    Compound('Biomass', 'cell'): 1 })

class TestSBMLDatabaseL2V5(unittest.TestCase):
    '''Test parsing of a simple level 2 version 5 SBML file'''

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level2/version5" level="2" version="5">
 <model>
  <listOfCompartments>
   <compartment id="C_c" name="cell"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose" name="Glucose" compartment="C_c"/>
   <species id="M_Glucose_6_P" name="Glucose-6-P" compartment="C_c"/>
   <species id="M_H2O" name="H2O" compartment="C_c"/>
   <species id="M_Phosphate" name="Phosphate" compartment="C_c"/>
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
        self.database = sbml.SBMLDatabase(s)

    def test_compounds_exist(self):
        self.assertEquals(len(list(self.database.compounds)), 5)
        self.assertIn(Compound('M_Glucose', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Glucose_6_P', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_H2O', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Phosphate', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Biomass', 'C_c'), self.database.compounds)

    def test_g6pase_reaction_exists(self):
        self.assertTrue(self.database.is_reversible('R_G6Pase'))

        # Read stoichiometric valus of reaction
        values = dict(self.database.get_reaction_values('R_G6Pase'))
        self.assertEquals(values, { Compound('M_Glucose', 'C_c'): -2,
                                    Compound('M_Phosphate', 'C_c'): -2,
                                    Compound('M_Glucose_6_P', 'C_c'): 2,
                                    Compound('M_H2O', 'C_c'): 2 })

    def test_biomass_reaction_exists(self):
        self.assertFalse(self.database.is_reversible('R_Biomass'))

        # Read stoichiometric values of reaction
        values = dict(self.database.get_reaction_values('R_Biomass'))
        self.assertEquals(values, { Compound('M_Glucose_6_P', 'C_c'): -Decimal('0.56'),
                                    Compound('M_Biomass', 'C_c'): 1 })

class TestSBMLDatabaseL3V1(unittest.TestCase):
    '''Test parsing of a simple level 3 version 1 SBML file'''

    def setUp(self):
        s = StringIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
 <model>
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
        self.database = sbml.SBMLDatabase(s)

    def test_compounds_exist(self):
        self.assertEquals(len(list(self.database.compounds)), 5)
        self.assertIn(Compound('M_Glucose', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Glucose_6_P', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_H2O', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Phosphate', 'C_c'), self.database.compounds)
        self.assertIn(Compound('M_Biomass', 'C_c'), self.database.compounds)

    def test_g6pase_reaction_exists(self):
        self.assertTrue(self.database.is_reversible('R_G6Pase'))

        # Read stoichiometric valus of reaction
        values = dict(self.database.get_reaction_values('R_G6Pase'))
        self.assertEquals(values, { Compound('M_Glucose', 'C_c'): -2,
                                    Compound('M_Phosphate', 'C_c'): -2,
                                    Compound('M_Glucose_6_P', 'C_c'): 2,
                                    Compound('M_H2O', 'C_c'): 2 })

    def test_biomass_reaction_exists(self):
        self.assertFalse(self.database.is_reversible('R_Biomass'))

        # Read stoichiometric values of reaction
        values = dict(self.database.get_reaction_values('R_Biomass'))
        self.assertEquals(values, { Compound('M_Glucose_6_P', 'C_c'): -Decimal('0.56'),
                                    Compound('M_Biomass', 'C_c'): 1 })
if __name__ == '__main__':
    unittest.main()
