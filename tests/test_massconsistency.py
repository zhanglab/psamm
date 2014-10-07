#!/usr/bin/env python

import unittest

from metnet import metabolicmodel
from metnet import massconsistency
from metnet import lpsolver
from metnet.reaction import ModelSEED

class TestMassConsistency(unittest.TestCase):
    '''Test fastcore using a simple model'''

    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = metabolicmodel.DictDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.load_model_from_file(iter(self.database.reactions))

        self.masscons = massconsistency.MassConsistencyCheck(lpsolver.CplexSolver(None))

    def test_mass_consistent_is_consistent(self):
        exchange = { 'rxn_1', 'rxn_6' }
        self.assertTrue(self.masscons.is_consistent(self.model, exchange, set()))

    def test_mass_inconsistent_is_consistent(self):
        exchange = { 'rxn_1', 'rxn_6' }
        self.database.set_reaction('rxn_7', ModelSEED.parse('|D| => (2) |C|'))
        self.model.add_reaction('rxn_7')
        self.assertFalse(self.masscons.is_consistent(self.model, exchange, set()))

if __name__ == '__main__':
    unittest.main()
