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
import math

from psamm.formula import Formula
from psamm.datasource.native import NativeModel
from psamm import balancecheck
from psamm.datasource.reaction import parse_reaction


class TestBalanceCheckWithModel(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'reactions:',
                '  - id: rxn_1',
                '    equation: A[e] => B[c]',
                '  - id: rxn_2',
                '    equation: B[c] => C[e]',
                '  - id: rxn_3',
                '    equation: A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]',
                'compounds:',
                '  - id: A',
                '    formula: C6H12O6',
                '    charge: 1',
                '  - id: B',
                '    formula: O2',
                '    charge: 1',
                '  - id: C',
                '    formula: CO2',
                '    charge: -1',
                '  - id: D',
                '    formula: H2O',
            ]))
        self._model = NativeModel.load_model_from_path(self._model_dir)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_charge_balance(self):
        d = {reaction.id: value for reaction, value in
             balancecheck.charge_balance(self._model)}
        self.assertEqual(d['rxn_1'], 0)
        self.assertEqual(d['rxn_2'], -2)
        self.assertTrue(math.isnan(d['rxn_3']))

    def test_formula_balance(self):
        d = {reaction.id: value for reaction, value in
             balancecheck.formula_balance(self._model)}
        self.assertEqual(d['rxn_1'][0], Formula.parse('C6H12O6'))
        self.assertEqual(d['rxn_1'][1], Formula.parse('O2'))

        self.assertEqual(d['rxn_2'][0], Formula.parse('O2'))
        self.assertEqual(d['rxn_2'][1], Formula.parse('CO2'))

        self.assertEqual(d['rxn_3'][0], Formula.parse('C6H12O18'))
        self.assertEqual(d['rxn_3'][1], Formula.parse('C6H12O18'))


class TestBalanceCheckWithReaction(unittest.TestCase):
    def test_reaction_charge_zero_sum(self):
        reaction = parse_reaction('A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]')
        compound_charge = {'A': 6, 'B': -1, 'C': 1, 'D': -1}
        charge_sum = balancecheck.reaction_charge(reaction, compound_charge)
        self.assertEqual(charge_sum, 0)

    def test_reaction_charge_non_zero_sum(self):
        reaction = parse_reaction('A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]')
        compound_charge = {'A': 1, 'B': -1, 'C': 1, 'D': -1}
        charge_sum = balancecheck.reaction_charge(reaction, compound_charge)
        self.assertEqual(charge_sum, 5)

    def test_reaction_charge_nan_sum(self):
        reaction = parse_reaction('A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]')
        compound_charge = {'A': 1, 'B': -1, 'C': 1}
        charge_sum = balancecheck.reaction_charge(reaction, compound_charge)
        self.assertTrue(math.isnan(charge_sum))

    def test_reaction_formula_normal_return(self):
        reaction = parse_reaction('A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]')
        compound_formula = {
            'A': Formula.parse('C6H12O6'),
            'B': Formula.parse('O2'),
            'C': Formula.parse('CO2'),
            'D': Formula.parse('H2O')
        }
        left_form, right_form = balancecheck.reaction_formula(
            reaction, compound_formula)
        self.assertEqual(left_form, Formula.parse('C6H12O18'))
        self.assertEqual(right_form, Formula.parse('C6H12O18'))

    def test_reaction_formula_none_return(self):
        reaction = parse_reaction('A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]')
        compound_formula = {
            'A': Formula.parse('C6H12O6'),
            'B': Formula.parse('O2'),
            'C': Formula.parse('CO2'),
        }
        result = balancecheck.reaction_formula(reaction, compound_formula)
        self.assertIsNone(result)
