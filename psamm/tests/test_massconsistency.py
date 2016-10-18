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
from psamm import massconsistency
from psamm.datasource.reaction import parse_reaction
from psamm.reaction import Compound
from psamm.lpsolver import generic

from six import iteritems


class TestMassConsistency(unittest.TestCase):
    """Test mass consistency using a simple model"""

    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |D|'))
        self.database.set_reaction('rxn_4', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_5', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_6', parse_reaction('|D| =>'))
        self.model = MetabolicModel.load_model(self.database, self.database.reactions)

        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

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

    def test_mass_consistent_reactions_returns_compounds(self):
        exchange = { 'rxn_1', 'rxn_6' }
        _, compounds = massconsistency.check_reaction_consistency(
            self.model, exchange=exchange, solver=self.solver)
        for c, value in compounds:
            self.assertIn(c, self.model.compounds)
            self.assertGreaterEqual(value, 1.0)

    def test_mass_consistent_reactions_returns_reactions(self):
        exchange = { 'rxn_1', 'rxn_6' }
        reactions, _ = massconsistency.check_reaction_consistency(
            self.model, exchange=exchange, solver=self.solver)
        for r, residual in reactions:
            self.assertIn(r, self.model.reactions)


class TestMassConsistencyZeroMass(unittest.TestCase):
    """Test mass consistency using a model with zero-mass compound"""

    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction(
            '|A| + |B| => |C|'))
        self.database.set_reaction('rxn_2', parse_reaction(
            '|C| + |Z| => |A| + |B|'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_is_consistent_with_zeromass(self):
        consistent = massconsistency.is_consistent(
            self.model, solver=self.solver, zeromass={Compound('Z')})
        self.assertTrue(consistent)

    def test_compound_consistency_with_zeromass(self):
        compounds = dict(massconsistency.check_compound_consistency(
            self.model, solver=self.solver, zeromass={Compound('Z')}))
        for c, value in iteritems(compounds):
            self.assertGreaterEqual(value, 1)

    def test_reaction_consistency_with_zeromass(self):
        reactions, _ = massconsistency.check_reaction_consistency(
            self.model, solver=self.solver, zeromass={Compound('Z')})
        reactions = dict(reactions)

        for r, value in iteritems(reactions):
            self.assertAlmostEqual(value, 0)


if __name__ == '__main__':
    unittest.main()
