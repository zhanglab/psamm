#!/usr/bin/env python

import unittest

from metnet import metabolicmodel
from metnet import fastcore
from metnet import lpsolver
from metnet.reaction import ModelSEED

class TestFastcoreSimpleVlassisModel(unittest.TestCase):
    '''Test fastcore using the simple model in Vlassis et al. 2014.'''

    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.get_model(self.database.reactions)
        self.fastcore = fastcore.Fastcore(lpsolver.CplexSolver(None))

    def test_lp10(self):
        result = self.fastcore.lp10(self.model, { 'rxn_6' }, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5' },
                                    epsilon=0.001, scaling=1e3)
        supp = set(fastcore.support(result, 0.999*0.001))
        self.assertEqual(supp, { 'rxn_1', 'rxn_3', 'rxn_6' })

    def test_lp10_weighted(self):
        weights = { 'rxn_3': 1 }
        result = self.fastcore.lp10(self.model, { 'rxn_6' }, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5' },
                                    epsilon=0.001, scaling=1e3, weights=weights)
        supp = set(fastcore.support(result, 0.999*0.001))
        self.assertEqual(supp, { 'rxn_1', 'rxn_3', 'rxn_6' })

        weights = { 'rxn_3': 3 }
        result = self.fastcore.lp10(self.model, { 'rxn_6' }, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5' },
                                    epsilon=0.001, scaling=1e3, weights=weights)
        supp = set(fastcore.support(result, 0.999*0.001))
        self.assertEqual(supp, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_lp7(self):
        result = self.fastcore.lp7(self.model, set(self.model.reaction_set), 0.001)
        supp = set(fastcore.support_positive(result, 0.001*0.999))
        self.assertEqual(supp, { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

        result = self.fastcore.lp7(self.model, {'rxn_5'}, 0.001)
        supp = set(fastcore.support_positive(result, 0.001*0.999))
        self.assertEqual(supp, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_find_sparse_mode_singleton(self):
        core = { 'rxn_1' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

        core = { 'rxn_2' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, set())

        core = { 'rxn_3' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

        core = { 'rxn_4' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

        core = { 'rxn_5' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

        core = { 'rxn_6' }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3))
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

    def test_find_sparse_mode_weighted(self):
        core = { 'rxn_1' }
        weights = { 'rxn_3': 1 }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3, weights=weights))
        self.assertEqual(mode, { 'rxn_1', 'rxn_3', 'rxn_6' })

        weights = { 'rxn_3': 3 }
        mode = set(self.fastcore.find_sparse_mode(self.model, core, set(self.model.reaction_set) - core,
                                                    epsilon=0.001, scaling=1e3, weights=weights))
        self.assertEqual(mode, { 'rxn_1', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_fastcc_inconsistent(self):
        self.assertEqual(set(self.fastcore.fastcc(self.model, 0.001)), { 'rxn_2' })

    def test_fastcc_is_consistent_on_inconsistent(self):
        self.assertFalse(self.fastcore.fastcc_is_consistent(self.model, 0.001))

    def test_fastcc_is_consistent_on_consistent(self):
        self.model.remove_reaction('rxn_2')
        self.assertTrue(self.fastcore.fastcc_is_consistent(self.model, 0.001))

    def test_fastcc_consistent_subset(self):
        self.assertEqual(self.fastcore.fastcc_consistent_subset(self.model, 0.001), set(['rxn_1', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6']))

    def test_fastcore_global_inconsistent(self):
        self.database.set_reaction('rxn_7', ModelSEED.parse('|E| <=>'))
        with self.assertRaises(Exception):
            self.fastcore.fastcore(self.model, { 'rxn_7' }, 0.001)

class TestFastcoreTinyBiomassModel(unittest.TestCase):
    '''Test fastcore using a model with tiny values in biomass reaction'''

    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('(0.000001) |A| =>'))
        self.model = self.database.get_model(self.database.reactions)
        self.fastcore = fastcore.Fastcore(lpsolver.CplexSolver(None))

    def test_fastcc_is_consistent(self):
        self.assertTrue(self.fastcore.fastcc_is_consistent(self.model, 0.001))

    def test_fastcore_induced_model(self):
        core = { 'rxn_2' }
        self.assertEquals(set(self.fastcore.fastcore(self.model, core, 0.001)), { 'rxn_1', 'rxn_2' })

    def test_fastcore_induced_model_high_epsilon(self):
        core = { 'rxn_2' }
        self.assertEquals(set(self.fastcore.fastcore(self.model, core, 0.1)), { 'rxn_1', 'rxn_2' })

class TestFlippingModel(unittest.TestCase):
    '''Test fastcore on a model that has to flip'''

    def setUp(self):
        # TODO use mock model instead of actual model
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('|A| <=>'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|C| <=> |B|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|C| <=>'))
        self.model = self.database.get_model(self.database.reactions)
        self.fastcore = fastcore.Fastcore(lpsolver.CplexSolver(None))

    def test_fastcore_induced_model(self):
        core = { 'rxn_2', 'rxn_3' }
        self.assertEquals(set(self.fastcore.fastcore(self.model, core, 0.001)),
                            { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4' })

if __name__ == '__main__':
    unittest.main()
