#!/usr/bin/env python

import unittest

from metnet import metabolicmodel
from metnet import fluxanalysis
from metnet import lpsolver
from metnet.reaction import ModelSEED

class TestFluxBalance(unittest.TestCase):
    def setUp(self):
        self.database = metabolicmodel.DictDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.get_model(self.database.reactions)
        self.solver = lpsolver.CplexSolver(None)

    def test_flux_balance_rxn_1(self):
        fluxes = dict(fluxanalysis.flux_balance(self.model, 'rxn_1', solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_6'], 1000)

    def test_flux_balance_rxn_2(self):
        fluxes = dict(fluxanalysis.flux_balance(self.model, 'rxn_2', solver=self.solver))
        self.assertEqual(fluxes['rxn_2'], 0)

    def test_flux_balance_rxn_3(self):
        fluxes = dict(fluxanalysis.flux_balance(self.model, 'rxn_3', solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_3'], 1000)
        self.assertEqual(fluxes['rxn_6'], 1000)

    def test_flux_balance_rxn_6(self):
        fluxes = dict(fluxanalysis.flux_balance(self.model, 'rxn_6', solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_6'], 1000)

class TestNaiveConsistency(unittest.TestCase):
    def setUp(self):
        self.database = metabolicmodel.DictDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.get_model(self.database.reactions)
        self.solver = lpsolver.CplexSolver(None)

    def test_check_inconsistent(self):
        self.model.remove_reaction('rxn_2')
        core = self.model.reaction_set
        inconsistent = set(fluxanalysis.naive_consistency_check(self.model, core, 0.001, solver=self.solver))
        self.assertEqual(inconsistent, {})

    def test_check_inconsistent(self):
        core = set(self.model.reactions)
        inconsistent = set(fluxanalysis.naive_consistency_check(self.model, core, 0.001, solver=self.solver))
        self.assertEqual(inconsistent, { 'rxn_2' })

if __name__ == '__main__':
    unittest.main()
