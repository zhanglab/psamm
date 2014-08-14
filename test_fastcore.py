#!/usr/bin/env python

import unittest

import metabolicmodel
import fastcore
from reaction import ModelSEED

class TestFastcore(unittest.TestCase):
    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.load_model_from_file(iter(self.database.reactions))

    def test_lp10(self):
        result = fastcore.fastcore_lp10_cplex(self.model, { 'rxn_6' }, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5' }, 0.001)
        supp = fastcore.support_set(result, 0.99*0.001)
        self.assertEqual(supp, { 'rxn_1', 'rxn_3', 'rxn_6' })

    def test_lp7(self):
        result = fastcore.fastcore_lp7_cplex(self.model, self.model.reaction_set, 0.001)
        supp = fastcore.support_set_positive(result, 0.001*0.99)
        self.assertEqual(supp, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

        result = fastcore.fastcore_lp7_cplex(self.model, {'rxn_5'}, 0.001)
        supp = fastcore.support_set_positive(result, 0.001*0.99)
        self.assertEqual(supp, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_find_sparse_mode(self):
        core = { 'rxn_1' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

        core = { 'rxn_2' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, set())

        core = { 'rxn_3' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

        core = { 'rxn_4' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

        core = { 'rxn_5' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

        core = { 'rxn_6' }
        mode = fastcore.find_sparse_mode(self.model, core, self.model.reaction_set - core, False, 0.001)
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

    def test_fastcc_inconsistent(self):
        self.assertEqual(list(fastcore.fastcc(self.model, 0.001)), ['rxn_2'])

    def test_fastcc_is_consistent_on_inconsistent(self):
        self.assertEqual(fastcore.fastcc_is_consistent(self.model, 0.001), False)

    def test_fastcc_is_consistent_on_consistent(self):
        self.model.remove_reaction('rxn_2')
        self.assertEqual(fastcore.fastcc_is_consistent(self.model, 0.001), True)

    def test_fastcc_consistent_subset(self):
        self.assertEqual(fastcore.fastcc_consistent_subset(self.model, 0.001), set(['rxn_1', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6']))

    def test_fastcore_global_inconsistent(self):
        self.database.set_reaction('rxn_7', ModelSEED.parse('|E| <=>'))
        self.assertRaises(Exception, fastcore.fastcore, (self.model, { 'rxn_7' }, 0.001))

if __name__ == '__main__':
    unittest.main()