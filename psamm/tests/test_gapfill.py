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

from psamm.database import DictDatabase
from psamm.metabolicmodel import MetabolicModel
from psamm.datasource.reaction import parse_reaction
from psamm.lpsolver import generic
from psamm.gapfill import gapfind, gapfill, GapFillError
from psamm.reaction import Compound


class TestGapfind(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> |A|'))
        self.database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', parse_reaction('|A| => |C|'))
        self.database.set_reaction('rxn_4', parse_reaction('|C| => |D|'))
        self.database.set_reaction('rxn_5', parse_reaction('|D| <=> |E|'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MILP solver for tests')

    def test_gapfind(self):
        self.model.remove_reaction('rxn_4')
        compounds = set(gapfind(self.model, self.solver, epsilon=0.1))
        self.assertEqual(compounds, {Compound('D'), Compound('E')})

    def test_gapfill_add_reaction(self):
        core = set(self.model.reactions) - {'rxn_4'}
        blocked = {Compound('D'), Compound('E')}

        add, rev = gapfill(
            self.model, core, blocked, {}, self.solver, epsilon=0.1)
        self.assertEqual(set(rev), set())
        self.assertEqual(set(add), {'rxn_4'})

    def test_gapfill_exclude_addition(self):
        core = set(self.model.reactions) - {'rxn_4'}
        blocked = {Compound('D'), Compound('E')}
        exclude = {'rxn_4'}

        with self.assertRaises(GapFillError):
            gapfill(self.model, core, blocked, exclude, self.solver,
                    epsilon=0.1)

    def test_gapfill_reverse_reaction(self):
        self.model.database.set_reaction('rxn_4', parse_reaction('|D| => |C|'))
        core = set(self.model.reactions)
        blocked = {Compound('D'), Compound('E')}

        add, rev = gapfill(
            self.model, core, blocked, {}, self.solver, epsilon=0.1)
        self.assertEqual(set(rev), {'rxn_4'})
        self.assertEqual(set(add), set())

    def test_gapfill_exclude_reversal(self):
        self.model.database.set_reaction('rxn_4', parse_reaction('D => C'))
        core = set(self.model.reactions)
        blocked = {Compound('D'), Compound('E')}
        exclude = {'rxn_4'}

        with self.assertRaises(GapFillError):
            gapfill(self.model, core, blocked, exclude, self.solver,
                    epsilon=0.1)
