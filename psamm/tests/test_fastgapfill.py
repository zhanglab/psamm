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
# Copyright 2016-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016  Chao liu <lcddzyx@gmail.com>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest

from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm import fastgapfill
from psamm.lpsolver import generic
from psamm.datasource.reaction import parse_reaction


class TestFastGapFill(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('|A| <=>'))
        self._database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self._database.set_reaction('rxn_3', parse_reaction('|C| <=> |B|'))
        self._database.set_reaction('rxn_4', parse_reaction('|C| <=>'))
        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)

        try:
            self._solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_fastgapfill(self):
        core = {'rxn_2', 'rxn_3'}
        induced = fastgapfill.fastgapfill(
            self._mm,
            core,
            epsilon=0.001,
            solver=self._solver)

        self.assertEqual(
            set(induced),
            {'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4'})


if __name__ == '__main__':
    unittest.main()
