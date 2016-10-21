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
# Copyright 2016  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016  Chao liu <lcddzyx@gmail.com>

import os
import unittest
import tempfile
import shutil

from psamm.datasource.native import NativeModel
from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm.reaction import Compound
from psamm import fastgapfill
from psamm.lpsolver import generic
from psamm.datasource.reaction import parse_reaction


class TestCreateExtendedModel(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'reactions:',
                '  - id: rxn_1',
                '    equation: A[e] <=> B[c]',
                '  - id: rxn_2',
                '    equation: A[e] => C[c]',
                '  - id: rxn_3',
                '    equation: A[e] => D[e]',
                '  - id: rxn_4',
                '    equation: C[c] => D[e]',
                'compounds:',
                '  - id: A',
                '  - id: B',
                '  - id: C',
                '  - id: D',
                'media:',
                '  - compartment: e',
                '    compounds:',
                '      - id: A',
                '        reaction: rxn_5',
                '      - id: D',
                '        reaction: rxn_6',
                'model:',
                '  - reactions:',
                '     - rxn_1',
                '     - rxn_2',
            ]))
        self._model = NativeModel.load_model_from_path(self._model_dir)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_create_model_extended(self):
        expected_reactions = set([
            'rxn_1',
            'rxn_2',
            'rxn_3',
            'rxn_4',
            'rxn_5',
            'rxn_6',
            ('rxntp', Compound('B', 'c')),
            ('rxntp', Compound('C', 'c')),
            ('rxnex', Compound('A', 'e')),
            ('rxnex', Compound('B', 'c')),
            ('rxnex', Compound('C', 'c')),
            ('rxnex', Compound('D', 'e'))
        ])

        expected_weights = {
            'rxn_3': 5.6,
            'rxn_4': 1.0,
            ('rxntp', Compound('B', 'c')): 3.0,
            ('rxntp', Compound('C', 'c')): 3.0,
            ('rxnex', Compound('A', 'e')): 2.0,
            ('rxnex', Compound('B', 'c')): 2.0,
            ('rxnex', Compound('C', 'c')): 2.0,
            ('rxnex', Compound('D', 'e')): 2.0
        }
        penalties = {'rxn_3': 5.6}

        model_extended, weights = fastgapfill.create_extended_model(
            self._model,
            db_penalty=1.0,
            ex_penalty=2.0,
            tp_penalty=3.0,
            penalties=penalties)
        self.assertEqual(set(model_extended.reactions), expected_reactions)
        self.assertEqual(weights, expected_weights)


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
