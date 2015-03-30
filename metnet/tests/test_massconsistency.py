#!/usr/bin/env python

import unittest

from metnet.metabolicmodel import MetabolicModel
from metnet.database import DictDatabase
from metnet import massconsistency
from metnet.lpsolver import cplex
from metnet.datasource.modelseed import parse_reaction

class TestMassConsistency(unittest.TestCase):
    '''Test fastcore using a simple model'''

    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_5', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_6', parse_reaction('|D| =>'))
        self.model = MetabolicModel.load_model(self.database, self.database.reactions)

        self.solver = cplex.Solver()

    def test_mass_consistent_is_consistent(self):
        exchange = { 'rxn_1', 'rxn_6' }
        self.assertTrue(massconsistency.is_consistent(
            self.model, self.solver, exchange, set()))

    def test_mass_inconsistent_is_consistent(self):
        exchange = { 'rxn_1', 'rxn_6' }
        self.database.set_reaction('rxn_7', parse_reaction('|D| => (2) |C|'))
        self.model.add_reaction('rxn_7')
        self.assertFalse(massconsistency.is_consistent(
            self.model, self.solver, exchange, set()))

if __name__ == '__main__':
    unittest.main()
