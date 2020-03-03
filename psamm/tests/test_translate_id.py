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
# Copyright 2020  Jing Wang <wjingsjtu@gmail.com>

import unittest
import tempfile
import shutil
import os

from psamm import translate_id as tr_id
from psamm.datasource.native import ModelReader, NativeModel
from psamm.datasource.reaction import parse_reaction


class TestTranslateId(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'biomass: rxn_1',
                'compartments:',
                '  - id: c',
                '    adjacent_to: e',
                '  - id: e',
                'reactions:',
                '  - id: rxn_1',
                '    name: rxn_1',
                '    equation: A[e] => B[c]',
                '    genes: gene1 and gene2 or gene3',
                '  - id: rxn_2',
                '    name: rxn_2',
                '    equation: B[c] => C[e]',
                '    genes: gene5 or gene6',
                '  - id: rxn_3',
                '    name: rxn_3',
                '    equation: A[e] + (6) B[c] <=> (6) C[e] + (6) D[c]',
                '    genes: gene7',
                'compounds:',
                '  - id: A',
                '    name: compound_A',
                '    formula: C6H12O6',
                '    charge: 1',
                '  - id: B',
                '    name: compound_B',
                '    formula: O2',
                '    charge: 1',
                '    kegg: C00010',
                '  - id: C',
                '    name: compound_C',
                '    formula: CO2',
                '    charge: -1',
                '  - id: D',
                '    name: compound_D',
                '    formula: H2O',
                'exchange:',
                '  - compartment: e',
                '  - compounds:',
                '    - id: A',
                '      upper: 100',
                '    - id: C',
                'limits:',
                '  - reaction: rxn_1',
                '    lower: -50',
                '  - reaction: rxn_2',
                '    upper: 100',
            ]))
        self._model = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'model.yaml')).create_model()

        with open(os.path.join(self._model_dir, 'newmodel.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'biomass: rxn_1a',
                'extracellular: s',
                'compartments:',
                '  - id: c',
                '    adjacent_to: s',
                '  - id: s',
                'reactions:',
                '  - id: rxn_1a',
                '    name: rxn_1',
                '    equation: A1[s] => B2[c]',
                '    genes: gene1 and gene2 or gene3',
                '  - id: rxn_2b',
                '    name: rxn_2',
                '    equation: B2[c] => A1[s]',
                '    genes: gene5 or gene6',
                '  - id: rxn_3c',
                '    name: rxn_3',
                '    equation: A1[s] + (6) B2[c] <=> (6) A1[s] + (6) D4[c]',
                '    genes: gene7',
                'compounds:',
                '  - id: A1',
                '    name: compound_A',
                '    formula: C6H12O6',
                '    charge: 1',
                '  - id: B2',
                '    name: compound_B',
                '    formula: O2',
                '    charge: 1',
                '    kegg: C00010',
                '  - id: C3',
                '    name: compound_C',
                '    formula: CO2',
                '    charge: -1',
                '  - id: D4',
                '    name: compound_D',
                '    formula: H2O',
                'exchange:',
                '  - compartment: s',
                '  - compounds:',
                '    - id: A1',
                '      upper: 100',
                '      reaction: EX_A1(s)',
                'limits:',
                '  - reaction: rxn_1a',
                '    lower: -50',
            ]))
        self._newmodel = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'newmodel.yaml')).create_model()

        with open(os.path.join(self._model_dir, 'compount_map.tsv'), 'w') as f:
            f.write('\n'.join([
                'id1\tid2\tp',
                'A\tA1\t1.0',
                'B\tB2\t0.9',
                'C\tA1\t0.001',
                'D\tD4\t0.32',
                'F\t\t',
            ]))

        with open(os.path.join(self._model_dir, 'reaction_map.tsv'), 'w') as f:
            f.write('\n'.join([
                'id1\tid2\tp',
                'rxn_1\trxn_1a\t1.0',
                'rxn_2\trxn_1a\t0.9',
                'rxn_4\t\t',
                'rxn_3\trxn_3c\t0.32',
            ]))

        with open(
                os.path.join(self._model_dir, 'compartment_map.tsv'),
                'w') as f:
            f.write('\n'.join([
                'e\ts',
            ]))

        self._compound_map = tr_id.read_mapping(
            os.path.join(self._model_dir, 'compount_map.tsv'))
        self._reaction_map = tr_id.read_mapping(
            os.path.join(self._model_dir, 'reaction_map.tsv'))
        self._compartment_map = tr_id.read_mapping(
            os.path.join(self._model_dir, 'compartment_map.tsv'))
        self._translated = tr_id.TranslatedModel(
            self._model, self._compound_map,
            self._reaction_map, self._compartment_map)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_instance(self):
        self.assertIsInstance(self._translated, NativeModel)
        self.assertIsInstance(self._translated, tr_id.TranslatedModel)

    def test_translate_compound(self):
        for cpd in ['A', 'B', 'C', 'D', 'C3']:
            self.assertNotIn(cpd, self._translated.compounds)
        for cpd in ['A1', 'B2', 'D4']:
            self.assertIn(cpd, self._translated.compounds)

    def test_translate_compartment(self):
        self.assertIn('c', self._translated.compartments)
        self.assertIn('s', self._translated.compartments)
        self.assertNotIn('e', self._translated.compartments)
        self.assertSetEqual(
            self._translated._compartment_boundaries, set({('c', 's')}))
        self.assertEqual(self._translated.extracellular_compartment, 's')

    def test_translate_reaction(self):
        for rxn in ['rxn_1', 'rxn_2', 'rxn_3', 'rxn_2b']:
            self.assertNotIn(rxn, self._translated.reactions)
        for rxn in ['rxn_1a', 'rxn_3c']:
            self.assertIn(rxn, self._translated.reactions)
            self.assertEqual(self._translated.reactions[rxn].equation,
                             self._newmodel.reactions[rxn].equation)

    def test_translate_exchange(self):
        self.assertDictEqual(self._translated.exchange,
                             self._newmodel.exchange)

    def test_translate_limits(self):
        self.assertDictEqual(self._translated.limits,
                             self._newmodel.limits)

    def test_write_model(self):
        self._translated.write_model(
            os.path.join(self._model_dir, 'translated'))

    def test_no_compartment_change(self):
        model = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'model.yaml')).create_model()
        translated = tr_id.TranslatedModel(
            model,
            self._translated.cpd_mapping_id,
            self._translated.rxn_mapping_id)
        self.assertEqual(translated.reactions['rxn_1a'].equation,
                         parse_reaction('A1[e] => B2[c]'))
        self.assertEqual(translated.reactions['rxn_3c'].equation,
                         parse_reaction(
                             'A1[e] + (6) B2[c] <=> (6) A1[e] + (6) D4[c]'))
        self.assertIn('c', translated.compartments)
        self.assertNotIn('s', translated.compartments)
        self.assertIn('e', translated.compartments)
        self.assertSetEqual(
            translated._compartment_boundaries, set({('c', 'e')}))
        self.assertEqual(translated.extracellular_compartment, 'e')
