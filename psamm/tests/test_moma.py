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
# Copyright 2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest

from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm import moma
from psamm.datasource.reaction import parse_reaction
from psamm.lpsolver import generic


class TestLinearMOMA(unittest.TestCase):
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

    def test_moma_fba(self):
        p = moma.MOMAProblem(self.model, self.solver)
        fluxes = p.get_fba_flux('rxn_6')
        self.assertAlmostEqual(fluxes['rxn_1'], 500)
        self.assertAlmostEqual(fluxes['rxn_2'], 0)
        self.assertAlmostEqual(fluxes['rxn_6'], 1000)

    def test_moma_minimal_fba(self):
        p = moma.MOMAProblem(self.model, self.solver)
        fluxes = p.get_minimal_fba_flux('rxn_6')
        self.assertAlmostEqual(fluxes['rxn_1'], 500)
        self.assertAlmostEqual(fluxes['rxn_2'], 0)
        self.assertAlmostEqual(fluxes['rxn_3'], 1000)
        self.assertAlmostEqual(fluxes['rxn_6'], 1000)

    def test_linear_moma(self):
        p = moma.MOMAProblem(self.model, self.solver)
        with p.constraints(p.get_flux_var('rxn_3') == 0):
            p.lin_moma({
                'rxn_3': 1000,
                'rxn_4': 0,
                'rxn_5': 0,
            })

        # The closest solution when these are constrained is for
        # rxn_6 to take on a flux of zero.
        self.assertAlmostEqual(p.get_flux('rxn_6'), 0)

    def test_linear_moma2(self):
        p = moma.MOMAProblem(self.model, self.solver)
        with p.constraints(p.get_flux_var('rxn_3') == 0):
            p.lin_moma2('rxn_6', 1000)

        self.assertAlmostEqual(p.get_flux('rxn_1'), 500)
        self.assertAlmostEqual(p.get_flux('rxn_2'), 0)
        self.assertAlmostEqual(p.get_flux('rxn_3'), 0)
        self.assertAlmostEqual(p.get_flux('rxn_4'), 1000)
        self.assertAlmostEqual(p.get_flux('rxn_5'), 1000)
        self.assertAlmostEqual(p.get_flux('rxn_6'), 1000)


class TestQuadraticMOMA(unittest.TestCase):
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
            self.solver = generic.Solver(quadratic=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MIQP solver for tests')

    def test_quadratic_moma(self):
        p = moma.MOMAProblem(self.model, self.solver)
        with p.constraints(p.get_flux_var('rxn_3') == 0):
            p.moma({
                'rxn_3': 1000,
                'rxn_4': 0,
                'rxn_5': 0,
            })

        # The closest solution when these are constrained is for
        # rxn_6 to take on a flux of zero.
        self.assertAlmostEqual(p.get_flux('rxn_6'), 0)

    def test_quadratic_moma2(self):
        p = moma.MOMAProblem(self.model, self.solver)
        with p.constraints(p.get_flux_var('rxn_3') == 0):
            p.moma2('rxn_6', 1000)

        self.assertAlmostEqual(p.get_flux('rxn_1'), 500)
        self.assertAlmostEqual(p.get_flux('rxn_2'), 0)
        self.assertAlmostEqual(p.get_flux('rxn_3'), 0)
        self.assertAlmostEqual(p.get_flux('rxn_4'), 1000)
        self.assertAlmostEqual(p.get_flux('rxn_5'), 1000)
        self.assertAlmostEqual(p.get_flux('rxn_6'), 1000)


if __name__ == '__main__':
    unittest.main()
