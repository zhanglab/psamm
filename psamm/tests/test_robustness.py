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
# Copyright 2020  Elysha Sameth <esameth1@my.uri.edu>

from __future__ import unicode_literals

import unittest

from psamm import fluxanalysis
from psamm.lpsolver import generic
from psamm.datasource import native
from psamm.database import DictDatabase
from psamm.metabolicmodel import MetabolicModel
from psamm.datasource.reaction import parse_reaction
from psamm.commands.robustness import RobustnessTaskHandler, \
    RobustnessTaskHandlerFva

from six import itervalues


class TestRobustnessTaskHandler(unittest.TestCase):
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

        self._problem = fluxanalysis.FluxBalanceProblem(
            self.model, self.solver)
        self._reactions = list(self.model.reactions)

    def test1_loop_removal_none(self):
        fluxes = RobustnessTaskHandler.handle_task(RobustnessTaskHandler(
            self.model, self.solver, 'none', self._reactions),
            ['rxn_6', 1000], 'rxn_6')

        self.assertAlmostEqual(fluxes['rxn_1'], 500)
        self.assertAlmostEqual(fluxes['rxn_2'], 0)
        self.assertAlmostEqual(fluxes['rxn_6'], 1000)

    def test2_loop_removal_l1min(self):
        fluxes = RobustnessTaskHandler.handle_task(RobustnessTaskHandler(
            self.model, self.solver, 'l1min', self._reactions),
            ['rxn_6', 1000], 'rxn_6')

        self.assertAlmostEqual(fluxes['rxn_1'], 500)
        self.assertAlmostEqual(fluxes['rxn_2'], 0)
        self.assertAlmostEqual(fluxes['rxn_3'], 1000)
        self.assertAlmostEqual(fluxes['rxn_4'], 0)
        self.assertAlmostEqual(fluxes['rxn_5'], 0)
        self.assertAlmostEqual(fluxes['rxn_6'], 1000)

    def test3_loop_removal_tfba(self):
        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MIQP solver for tests')

        fluxes = RobustnessTaskHandler.handle_task(RobustnessTaskHandler(
            self.model, self.solver, 'tfba', self._reactions),
            ['rxn_3', 100], 'rxn_1')

        self.assertAlmostEqual(fluxes['rxn_1'], 500)
        self.assertAlmostEqual(fluxes['rxn_2'], 0)
        self.assertAlmostEqual(fluxes['rxn_3'], 100)
        self.assertAlmostEqual(fluxes['rxn_4'], 900)
        self.assertAlmostEqual(fluxes['rxn_5'], 900)
        self.assertAlmostEqual(fluxes['rxn_6'], 1000)

    def test4_no_reactions(self):
        self._reactions2 = None
        fluxes = RobustnessTaskHandler.handle_task(RobustnessTaskHandler(
            self.model, self.solver, 'none', self._reactions2),
            ['rxn_6', 10], 'rxn_6')
        self.assertEqual(fluxes, 10)


class TestRobustnessTaskHandlerFva(unittest.TestCase):
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
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

        self._problem = fluxanalysis.FluxBalanceProblem(
            self.model, self.solver)
        self._reactions = list(self.model.reactions)

    def test1_loop_removal_none(self):
        fluxes = RobustnessTaskHandlerFva.handle_task(RobustnessTaskHandlerFva(
            self.model, self.solver, 'none', self._reactions),
            ['rxn_6', 200], 'rxn_6')

        for bounds in itervalues(fluxes):
            self.assertEqual(len(bounds), 2)

        self.assertAlmostEqual(fluxes['rxn_1'][0], 100)

        self.assertAlmostEqual(fluxes['rxn_2'][0], 0)
        self.assertAlmostEqual(fluxes['rxn_2'][1], 0)

        self.assertAlmostEqual(fluxes['rxn_5'][0], 0)
        self.assertAlmostEqual(fluxes['rxn_5'][1], 100)

        self.assertAlmostEqual(fluxes['rxn_6'][0], 200)

        self.assertGreater(fluxes['rxn_7'][1], 0)
        self.assertGreater(fluxes['rxn_8'][1], 0)

    def test2_loop_removal_tfba(self):
        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MIQP solver for tests')

        fluxes = RobustnessTaskHandlerFva.handle_task(RobustnessTaskHandlerFva(
            self.model, self.solver, 'tfba', self._reactions),
            ['rxn_6', 200], 'rxn_3')

        for bounds in itervalues(fluxes):
            self.assertEqual(len(bounds), 2)

        self.assertAlmostEqual(fluxes['rxn_1'][0], 100)

        self.assertAlmostEqual(fluxes['rxn_2'][0], 0)
        self.assertAlmostEqual(fluxes['rxn_2'][1], 0)

        self.assertAlmostEqual(fluxes['rxn_3'][0], 200)

        self.assertAlmostEqual(fluxes['rxn_5'][0], 0)
        self.assertAlmostEqual(fluxes['rxn_5'][1], 0)

        self.assertAlmostEqual(fluxes['rxn_6'][0], 200)

        self.assertAlmostEqual(fluxes['rxn_7'][1], 0)
        self.assertAlmostEqual(fluxes['rxn_8'][1], 0)

    def test3_no_reactions(self):
        self._reactions2 = None
        fluxes = RobustnessTaskHandlerFva.handle_task(RobustnessTaskHandlerFva(
            self.model, self.solver, 'none', self._reactions2),
            ['rxn_6', 200], 'rxn_6')

        for bounds in itervalues(fluxes):
            self.assertEqual(len(bounds), 2)

        self.assertAlmostEqual(fluxes['rxn_6'][0], 200)
        self.assertAlmostEqual(fluxes['rxn_6'][1], 200)


if __name__ == '__main__':
    unittest.main()
