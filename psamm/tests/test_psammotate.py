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
# Copyright 2020 Keith Dufault-Thompson <keitht547@uri.edu>

import unittest
import tempfile
import shutil
import os

from psamm.commands import psammotate
from psamm.datasource.native import ModelReader


class TestPsammotate(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'reactions:',
                '  - id: rxn_1',
                '    equation: A[e] => B[c]',
                '    genes: gene_1',
                '  - id: rxn_2',
                '    equation: B[c] => C[c]',
                '    genes: gene_2 and gene_3',
                '  - id: rxn_3',
                '    equation: C[c] <=> D[c]',
                '    genes: gene_3 and gene_4',
                '  - id: rxn_4',
                '    equation: D[c] => E[e]',
                '    genes: gene_3 or gene_4',
                '  - id: rxn_5',
                '    equation: C[c] => D[e]',
                '    genes: (gene_3 or gene_4) and (gene_1 or gene_5)',
                '  - id: rxn_6',
                '    equation: D[c] => A[e]',
                '    genes: Biomass',
                '  - id: rxn_7',
                '    equation: E[e] => A[e]',
                '  - id: rxn_8',
                '    equation: E[e] =>'
            ]))
        self._model = ModelReader.reader_from_path(
            self._model_dir).create_model()
        self.rxn_entry_dict = {}
        for rx in self._model.reactions:
            self.rxn_entry_dict[rx.id] = rx

        with open(os.path.join(self._model_dir, 'app.tsv'), 'w') as f:
            f.write('\n'.join([
                'gene_1\tgene_a',
                'gene_2\tgene_b',
                'gene_3\tgene_c',
                'gene_4\t-',
                'gene_5\tgene_e, gene_f',
                '-\tgene_g'
            ]))

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_app_reader(self):
        translation_dict = psammotate.app_reader(
            open(os.path.join(self._model_dir, 'app.tsv')), 1, 0)
        self.assertTrue('-' not in translation_dict)
        self.assertTrue(translation_dict['gene_1'] == ['gene_a'])
        self.assertTrue(translation_dict['gene_2'] == ['gene_b'])
        self.assertTrue(translation_dict['gene_3'] == ['gene_c'])
        self.assertTrue(translation_dict['gene_4'] == ['-'])
        self.assertTrue(translation_dict['gene_5'] == ['gene_e', 'gene_f'])

    def test_model_loader_ignorena(self):
        translation_dict = psammotate.app_reader(
            open(os.path.join(self._model_dir, 'app.tsv')), 1, 0)
        translated_genes = psammotate.model_loader(
            self._model, True, translation_dict)
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_1']],
                         ['gene_1', 'gene_a', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_2']],
                         ['gene_2 and gene_3', 'gene_b and gene_c', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_3']],
                         ['gene_3 and gene_4', 'gene_c and -', False])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_4']],
                         ['gene_3 or gene_4', 'gene_c or -', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_5']],
                         ['(gene_3 or gene_4) and (gene_1 or gene_5)',
                          '(gene_c or -) and (gene_a or (gene_e or gene_f))',
                          True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_7']],
                         [None, None, False])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_8']],
                         [None, None, False])

    def test_model_loader(self):
        translation_dict = psammotate.app_reader(
            open(os.path.join(self._model_dir, 'app.tsv')), 1, 0)
        translated_genes = psammotate.model_loader(
            self._model, False, translation_dict)
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_1']],
                         ['gene_1', 'gene_a', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_2']],
                         ['gene_2 and gene_3', 'gene_b and gene_c', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_3']],
                         ['gene_3 and gene_4', 'gene_c and -', False])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_4']],
                         ['gene_3 or gene_4', 'gene_c or -', True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_5']],
                         ['(gene_3 or gene_4) and (gene_1 or gene_5)',
                          '(gene_c or -) and (gene_a or (gene_e or gene_f))',
                          True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_7']],
                         [None, None, True])
        self.assertEqual(translated_genes[self.rxn_entry_dict['rxn_8']],
                         [None, None, True])

    def test_get_gene_list_and(self):
        genes = 'gene_a and gene_b'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b']))

    def test_get_gene_list_or(self):
        genes = 'gene_a or gene_b'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b']))

    def test_get_gene_list_nested_or(self):
        genes = 'gene_a or (gene_b or gene_c)'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b', 'gene_c']))

    def test_get_gene_list_nested_and(self):
        genes = 'gene_a or (gene_b and gene_c)'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b', 'gene_c']))

    def test_get_gene_list_nested_and_or(self):
        genes = '(gene_a or gene_d) or (gene_b and gene_c)'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b', 'gene_c', 'gene_d']))

    def test_get_gene_list_multiple_nests(self):
        genes = '(gene_a or (gene_b and gene_c)) or gene_d'
        self.assertEqual(psammotate.get_gene_list(genes),
                         frozenset(['gene_a', 'gene_b', 'gene_c', 'gene_d']))

    def test_remove_gaps_single(self):
        genes = 'gene_a'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a')

    def test_remove_gaps_or(self):
        genes = 'gene_a or gene_b'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a or gene_b')

    def test_remove_gaps_and(self):
        genes = 'gene_a and gene_b'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a and gene_b')

    def test_remove_gaps_and_missing(self):
        genes = 'gene_a and -'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a and -')

    def test_remove_gaps_or_missing(self):
        genes = 'gene_a or -'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a')

    def test_remove_gaps_or_missing_nested(self):
        genes = 'gene_a or (- or -)'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a')

    def test_remove_gaps_and_missing_nested(self):
        genes = 'gene_a and (- or -)'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a and -')

    def test_remove_gaps_double_nested(self):
        genes = 'gene_a or (- or (- or gene_b))'
        self.assertEqual(psammotate.remove_gap(genes), 'gene_a or ((gene_b))')


if __name__ == '__main__':
    unittest.main()
