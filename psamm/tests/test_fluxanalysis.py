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

from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm import fluxanalysis
from psamm.datasource.modelseed import parse_reaction
from psamm.lpsolver import generic


class TestFluxBalance(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_5', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_6', parse_reaction('|D| =>'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_flux_balance_rxn_1(self):
        fluxes = dict(fluxanalysis.flux_balance(
            self.model, 'rxn_1', tfba=False, solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_6'], 1000)

    def test_flux_balance_rxn_2(self):
        fluxes = dict(fluxanalysis.flux_balance(
            self.model, 'rxn_2', tfba=False, solver=self.solver))
        self.assertEqual(fluxes['rxn_2'], 0)

    def test_flux_balance_rxn_3(self):
        fluxes = dict(fluxanalysis.flux_balance(
            self.model, 'rxn_3', tfba=False, solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_3'], 1000)
        self.assertEqual(fluxes['rxn_6'], 1000)

    def test_flux_balance_rxn_6(self):
        fluxes = dict(fluxanalysis.flux_balance(
            self.model, 'rxn_6', tfba=False, solver=self.solver))
        self.assertEqual(fluxes['rxn_1'], 500)
        self.assertEqual(fluxes['rxn_2'], 0)
        self.assertEqual(fluxes['rxn_6'], 1000)


class TestFluxBalanceThermodynamic(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('ex_A', parse_reaction('|A| <=>'))
        self.database.set_reaction('ex_D', parse_reaction('|D| <=>'))
        self.database.set_reaction('rxn_1', parse_reaction('|A| => |B|'))
        self.database.set_reaction('rxn_2', parse_reaction('|B| <=> |C|'))
        self.database.set_reaction('rxn_3', parse_reaction('|C| <=> |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|D| <=> |E|'))
        self.database.set_reaction('rxn_5', parse_reaction('|E| => |B|'))

        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)
        self.model.limits['ex_A'].lower = -10 # Low uptake
        self.model.limits['ex_D'].lower = 0 # No uptake

        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MILP solver for tests')

    def test_flux_balance_tfba_exchange_d(self):
        fluxes = dict(fluxanalysis.flux_balance(
            self.model, 'ex_D', tfba=True, solver=self.solver))
        self.assertEquals(fluxes['ex_A'], -10)
        self.assertEquals(fluxes['ex_D'], 10)
        self.assertEquals(fluxes['rxn_2'], 10)
        self.assertEquals(fluxes['rxn_4'], 0)
        self.assertEquals(fluxes['rxn_5'], 0)


class TestFluxVariability(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_5', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_6', parse_reaction('|D| =>'))
        self.database.set_reaction('rxn_7', parse_reaction('|E| => |F|'))
        self.database.set_reaction('rxn_8', parse_reaction('|F| => |E|'))

        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)
        self.model.limits['rxn_5'].upper = 100

        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MILP solver for tests')

    def test_flux_variability(self):
        fluxes = dict(fluxanalysis.flux_variability(
            self.model, self.model.reactions, {'rxn_6': 200},
            tfba=False, solver=self.solver))

        self.assertEqual(fluxes['rxn_1'][0], 100)

        self.assertEqual(fluxes['rxn_2'][0], 0)
        self.assertEqual(fluxes['rxn_2'][1], 0)

        self.assertEqual(fluxes['rxn_5'][0], 0)
        self.assertEqual(fluxes['rxn_5'][1], 100)

        self.assertEqual(fluxes['rxn_6'][0], 200)

        self.assertGreater(fluxes['rxn_7'][1], 0)
        self.assertGreater(fluxes['rxn_8'][1], 0)

    def test_flux_variability_thermodynamic(self):
        fluxes = dict(fluxanalysis.flux_variability(
            self.model, self.model.reactions, {'rxn_6': 200},
            tfba=True, solver=self.solver))

        self.assertEqual(fluxes['rxn_1'][0], 100)

        self.assertEqual(fluxes['rxn_2'][0], 0)
        self.assertEqual(fluxes['rxn_2'][1], 0)

        self.assertEqual(fluxes['rxn_5'][0], 0)
        self.assertEqual(fluxes['rxn_5'][1], 100)

        self.assertEqual(fluxes['rxn_6'][0], 200)

        self.assertEqual(fluxes['rxn_7'][1], 0)
        self.assertEqual(fluxes['rxn_8'][1], 0)


class TestFluxConsistency(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_5', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_6', parse_reaction('|D| =>'))
        self.database.set_reaction('rxn_7', parse_reaction('|E| => |F|'))
        self.database.set_reaction('rxn_8', parse_reaction('|F| => |E|'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MILP solver for tests')

    def test_check_on_consistent(self):
        self.model.remove_reaction('rxn_2')
        core = self.model.reactions
        inconsistent = set(fluxanalysis.consistency_check(
            self.model, core, epsilon=0.001, tfba=False, solver=self.solver))
        self.assertEqual(inconsistent, set())

    def test_check_on_inconsistent(self):
        core = set(self.model.reactions)
        inconsistent = set(fluxanalysis.consistency_check(
            self.model, core, epsilon=0.001, tfba=False, solver=self.solver))
        self.assertEqual(inconsistent, {'rxn_2'})

    def test_check_on_inconsistent_with_thermodynamic_constraints(self):
        self.model.remove_reaction('rxn_2')
        core = self.model.reactions
        inconsistent = set(fluxanalysis.consistency_check(
            self.model, core, epsilon=0.001, tfba=True, solver=self.solver))
        self.assertEqual(inconsistent, {'rxn_7', 'rxn_8'})


if __name__ == '__main__':
    unittest.main()
