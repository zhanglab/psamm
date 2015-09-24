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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import unittest

from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm import fluxcoupling
from psamm.datasource.modelseed import parse_reaction
from psamm.lpsolver import generic


class TestFluxCouplingBurgardModel(unittest.TestCase):
    """Test case based on the simple model in [Burgard04]_."""

    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('|A| => |B|'))
        self.database.set_reaction('rxn_2', parse_reaction('|B| => |G|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |G|'))
        self.database.set_reaction('rxn_4', parse_reaction('|B| => |H|'))
        self.database.set_reaction('rxn_5', parse_reaction('|B| => |C| + |F|'))
        self.database.set_reaction('rxn_6', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_7', parse_reaction('|E| <=> |D|'))
        self.database.set_reaction('rxn_9', parse_reaction('|I| => |J|'))
        self.database.set_reaction('rxn_10', parse_reaction('|J| => |K|'))
        self.database.set_reaction('rxn_A', parse_reaction('=> |A|'))
        self.database.set_reaction('rxn_G', parse_reaction('|G| =>'))
        self.database.set_reaction('rxn_E', parse_reaction('|E| =>'))
        self.database.set_reaction(
            'rxn_bio', parse_reaction('|D| + (2.5) |F| =>'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_flux_coupling(self):
        fcp = fluxcoupling.FluxCouplingProblem(self.model, {}, self.solver)
        reactions = sorted(self.model.reactions)
        couplings = {}
        for r1 in reactions:
            for r2 in reactions:
                couplings[r1, r2] = fcp.solve(r1, r2)

        inconsistent = {'rxn_4', 'rxn_9', 'rxn_10'}
        coupled_set = {'rxn_5', 'rxn_6', 'rxn_7', 'rxn_bio', 'rxn_E'}

        for i, r1 in enumerate(reactions):
            for r2 in reactions[i:]:
                if r1 == r2:
                    # A reaction is always fully coupled to itself
                    # unless it is inconsistent.
                    lower, upper = couplings[r1, r2]
                    if r1 in inconsistent:
                        self.assertEqual((lower, upper), (None, None))
                    else:
                        print('{}, {}: {}, {}'.format(r1, r2, lower, upper))
                        self.assertAlmostEqual(lower, upper)
                elif r1 in inconsistent or r2 in inconsistent:
                    # Couplings with inconsistent reactions are always
                    # fully coupled at 0.0 or uncoupled.
                    coupling = couplings[r1, r2]
                    self.assertIn(coupling, ((0.0, 0.0), (None, None)))
                elif r1 in coupled_set and r2 in coupled_set:
                    lower, upper = couplings[r1, r2]
                    self.assertAlmostEqual(lower, upper)
                else:
                    # At least one of the values must be 0 or None if the
                    # reactions are not in the coupled set.
                    coupling = couplings[r1, r2]
                    self.assertTrue(any(v in (0.0, None) for v in coupling))

        # Couplings with rxn_7 are negative
        for r in coupled_set:
            if r != 'rxn_7':
                self.assertLess(couplings['rxn_7', r][0], 0.0)


class TestFluxCouplingClass(unittest.TestCase):
    def test_uncoupled(self):
        cc = fluxcoupling.classify_coupling((None, None))
        self.assertEqual(cc, fluxcoupling.CouplingClass.Uncoupled)

    def test_directional_reverse(self):
        cc = fluxcoupling.classify_coupling((None, 1000))
        self.assertEqual(cc, fluxcoupling.CouplingClass.DirectionalReverse)

        cc = fluxcoupling.classify_coupling((-1000, None))
        self.assertEqual(cc, fluxcoupling.CouplingClass.DirectionalReverse)

    def test_inconsistent(self):
        cc = fluxcoupling.classify_coupling((0.0, 0.0))
        self.assertEqual(cc, fluxcoupling.CouplingClass.Inconsistent)

    def test_directional_forward(self):
        cc = fluxcoupling.classify_coupling((-1000, 100.0))
        self.assertEqual(cc, fluxcoupling.CouplingClass.DirectionalForward)

    def test_full(self):
        cc = fluxcoupling.classify_coupling((-1000, -1000.0))
        self.assertEqual(cc, fluxcoupling.CouplingClass.Full)


if __name__ == '__main__':
    unittest.main()
