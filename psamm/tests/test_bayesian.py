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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest
import tempfile
import os
import shutil
from collections import namedtuple
import pandas as pd

from psamm import bayesian, bayesian_util as util
from psamm.datasource.native import (
    ModelReader, CompoundEntry, Formula, ReactionEntry)
from psamm.datasource.reaction import parse_reaction


class TestBayesianPredictor(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model1.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
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
                '    name: A',
                '    formula: C6H12O6',
                '    charge: 1',
                '  - id: B',
                '    name: B',
                '    formula: O2',
                '    charge: 1',
                '    kegg: C00010',
                '  - id: C',
                '    name: C',
                '    formula: CO2',
                '    charge: -1',
                '  - id: D',
                '    name: D',
                '    formula: H2O',
            ]))
        self._model1 = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'model1.yaml')).create_model()
        self._model1 = bayesian.MappingModel(self._model1)

        with open(os.path.join(self._model_dir, 'model2.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'reactions:',
                '  - id: rxn_1',
                '    name: rxn_1',
                '    equation: A[e] => B2[c]',
                '    genes: gene1 and gene2 or gene3',
                '  - id: rnx_2',
                '    equation: B2[c] => C[e]',
                '    genes: gene5 and gene6',
                '  - id: rxn_3',
                '    name: rxn_3',
                '    equation: A[e] <=> (6) C[e] + (6) D4[c]',
                '  - id: rxn_4',
                '    equation: E[c] => F[c]',
                'compounds:',
                '  - id: A',
                '    name: a',
                '    formula: C6H11O6',
                '    charge: 0',
                '  - id: B2',
                '    name: -B ()',
                '    formula: O2',
                '    charge: 1',
                '    kegg: C00010',
                '  - id: C',
                '    name: C',
                '    charge: -1',
                '  - id: D4',
                '    name: D',
                '    formula: H2O',
                '    charge: -1',
            ]))
        self._model2 = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'model2.yaml')).create_model()
        self._model2 = bayesian.MappingModel(self._model2)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_summary(self):
        self._model1.print_summary()
        self._model2.print_summary()

    def test_consistency_check(self):
        self.assertTrue(self._model1.check_reaction_compounds())
        self.assertFalse(self._model2.check_reaction_compounds())

    def test_compound_map(self):
        p = bayesian.BayesianCompoundPredictor(
            self._model1, self._model2, 1, self._model_dir, True, True
        )
        self.assertEqual(p.model1, self._model1)
        self.assertEqual(p.model2, self._model2)
        self.assertEqual(p.map('A', 'A'), 1)
        self.assertIsInstance(p, bayesian.BayesianCompoundPredictor)
        best_map = p.get_best_map()
        self.assertIsInstance(best_map, pd.DataFrame)
        self.assertIsInstance(p.get_raw_map(), pd.DataFrame)
        self.assertIsInstance(p.get_cpd_pred(), pd.Series)
        self.assertLess(best_map.loc[('B', 'B2'), 'p_id'],
                        best_map.loc[('A', 'A'), 'p_id'])
        self.assertIn(('A', 'A'), best_map.index)
        self.assertIn(('B', 'B2'), best_map.index)
        self.assertIn(('C', 'C'), best_map.index)
        self.assertIn(('D', 'D4'), best_map.index)

    def test_reaction_map(self):
        cpd_pred = bayesian.BayesianCompoundPredictor(
            self._model1, self._model2, 1, self._model_dir, False, False
        ).get_cpd_pred()
        p = bayesian.BayesianReactionPredictor(
            self._model1, self._model2, cpd_pred,
            1, self._model_dir, True, True
        )
        self.assertIsInstance(p, bayesian.BayesianReactionPredictor)
        best_map = p.get_best_map()
        self.assertIsInstance(best_map, pd.DataFrame)
        self.assertIsInstance(p.get_raw_map(), pd.DataFrame)
        self.assertIn(('rxn_1', 'rxn_1'), best_map.index)
        self.assertIn(('rxn_2', 'rnx_2'), best_map.index)
        self.assertIn(('rxn_3', 'rxn_3'), best_map.index)

        p = bayesian.BayesianReactionPredictor(
            self._model1, self._model2, cpd_pred,
            1, self._model_dir, True, False
        )
        self.assertEqual(p.model1, self._model1)
        self.assertEqual(p.model2, self._model2)
        self.assertAlmostEqual(p.map('rxn_1', 'rxn_1'), 1, 1)
        self.assertIsInstance(p, bayesian.BayesianReactionPredictor)
        best_map = p.get_best_map()
        self.assertIsInstance(best_map, pd.DataFrame)
        self.assertIsInstance(p.get_raw_map(), pd.DataFrame)

    def test_compound_id_likelihood(self):
        match, unmatch = bayesian.compound_id_likelihood(
            CompoundEntry({'id': 'A'}), CompoundEntry({'id': 'A'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.compound_id_likelihood(
            CompoundEntry({'id': 'A'}), CompoundEntry({'id': 'B'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertLess(match, unmatch)

    def test_compound_name_likelihood(self):
        match, unmatch = bayesian.compound_name_likelihood(
            CompoundEntry({'id': 'A', 'name': 'Test'}),
            CompoundEntry({'id': 'B', 'name': 'te-s(t)'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.compound_name_likelihood(
            CompoundEntry({'id': 'A', 'name': 'Test'}),
            CompoundEntry({'id': 'B', 'name': 'Test1'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertLess(match, unmatch)

    def test_compound_charge_likelihood(self):
        match, unmatch = bayesian.compound_charge_likelihood(
            CompoundEntry({'id': 'A', 'charge': 1}),
            CompoundEntry({'id': 'B', 'charge': 1}),
            1.0 / 500, 0.5, 0.5
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.compound_charge_likelihood(
            CompoundEntry({'id': 'A', 'charge': 1}),
            CompoundEntry({'id': 'B', 'charge': -1}),
            1.0 / 500, 0.5, 0.5
        )
        self.assertLess(match, unmatch)
        match, unmatch = bayesian.compound_charge_likelihood(
            CompoundEntry({'id': 'A'}),
            CompoundEntry({'id': 'B', 'charge': -1}),
            1.0 / 500, 0.5, 0.5
        )
        self.assertEqual(match, unmatch)

    def test_compound_formula_likelihood(self):
        match, unmatch = bayesian.compound_formula_likelihood(
            CompoundEntry({'id': 'A', 'formula': Formula.parse('C3H4')}),
            CompoundEntry(
                {'id': 'B', 'formula': Formula.parse('C3H5'), 'charge': 1}),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.compound_formula_likelihood(
            CompoundEntry({'id': 'A', 'formula': Formula.parse('C3H4')}),
            CompoundEntry({'id': 'B', 'formula': Formula.parse('C3H5')}),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertLess(match, unmatch)
        match, unmatch = bayesian.compound_formula_likelihood(
            CompoundEntry({'id': 'A', 'formula': Formula.parse('C3H4')}),
            CompoundEntry({'id': 'B', 'formula': None}),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertEqual(match, unmatch)

    def test_compound_kegg_likelihood(self):
        CompoundEntry = namedtuple('CompoundEntry', ['id', 'kegg'])
        match, unmatch = bayesian.compound_kegg_likelihood(
            CompoundEntry(id='A', kegg='C00001'),
            CompoundEntry(id='B', kegg='C00001'),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.compound_kegg_likelihood(
            CompoundEntry(id='A', kegg='C00001'),
            CompoundEntry(id='B', kegg='C00003'),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertLess(match, unmatch)
        match, unmatch = bayesian.compound_kegg_likelihood(
            CompoundEntry(id='A', kegg=None),
            CompoundEntry(id='B', kegg='C00001'),
            1.0 / 500, 1.0 / 1000, 0.5
        )
        self.assertEqual(match, unmatch)

    def test_reaction_id_likelihood(self):
        match, unmatch = bayesian.reaction_id_likelihood(
            ReactionEntry({'id': 'r1'}),
            ReactionEntry({'id': 'r1'}),
            1.0 / 500, 1.0 / 500, 499.0 / 500
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.reaction_id_likelihood(
            ReactionEntry({'id': 'r1'}),
            ReactionEntry({'id': 'r2'}),
            1.0 / 500, 1.0 / 500, 499.0 / 500
        )
        self.assertLess(match, unmatch)

    def test_reaction_name_likelihood(self):
        match, unmatch = bayesian.reaction_name_likelihood(
            ReactionEntry({'id': 'r1', 'name': 'Test'}),
            ReactionEntry({'id': 'r1', 'name': 'te-s(t)'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.reaction_name_likelihood(
            ReactionEntry({'id': 'r1', 'name': 'Test'}),
            ReactionEntry({'id': 'r2', 'name': 'Test2'}),
            1.0 / 500, 1.0 / 1000
        )
        self.assertLess(match, unmatch)

    def test_reaction_genes_likelihood(self):
        match, unmatch = bayesian.reaction_genes_likelihood(
            ReactionEntry({'id': 'r1', 'genes': 'g1 and g2 or g3'}),
            ReactionEntry({'id': 'r1', 'genes': 'g2 and g1'}),
            1.0 / 500, 1.0 / 1000, 0.8
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.reaction_genes_likelihood(
            ReactionEntry({'id': 'r1', 'genes': 'g1 and g2 or g3'}),
            ReactionEntry({'id': 'r1', 'genes': 'g2 and g4'}),
            1.0 / 500, 1.0 / 1000, 0.8
        )
        self.assertLessEqual(match, unmatch)
        match, unmatch = bayesian.reaction_genes_likelihood(
            ReactionEntry({'id': 'r1', 'genes': 'g1 and g2 or g3'}),
            ReactionEntry({'id': 'r1', 'genes': 'g2 and g4'}),
            1.0 / 500, 1.0 / 1000, 0.8, {'g1': 'g4', 'g4': 'g1'}
        )
        self.assertGreater(match, unmatch)
        match, unmatch = bayesian.reaction_genes_likelihood(
            ReactionEntry({'id': 'r1', 'genes': None}),
            ReactionEntry({'id': 'r1', 'genes': 'g4 and g5'}),
            1.0 / 500, 1.0 / 1000, 0.8
        )
        self.assertEqual(match, unmatch)

    def test_reaction_equation_likelihood(self):
        cpd_map = {'A': set(['A1']),
                   'B': set(['B1']),
                   'C': set(['C1']),
                   'D': set(['D1'])}
        cpd_score = {'A': 1, 'B': 1, 'C': 0.9, 'D': 0.9}
        eq1 = parse_reaction('A[c] + C[c] + B[c] => D[c]')
        eq2 = parse_reaction('A1[c] + C1[c] + B1[c] => D1[c]')
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': eq2}),
            cpd_map, cpd_score
        )
        self.assertGreater(m, un)
        eq2 = parse_reaction('D1[c] + F1[c] <=> A1[c] + C1[c] + B1[c]')
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': eq2}),
            cpd_map, cpd_score
        )
        self.assertGreater(m, un)
        eq2 = parse_reaction('A1[c] => D1[c] + C1[c]')
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': eq2}),
            cpd_map, cpd_score
        )
        self.assertLess(m, un)
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': None}),
            cpd_map, cpd_score
        )
        self.assertEqual(m, un)
        eq2 = parse_reaction('A1[s] + C1[s] + B1[s] => D1[s]')
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': eq2}),
            cpd_map, cpd_score
        )
        self.assertLess(m, un)
        m, un = bayesian.reaction_equation_compound_mapping_likelihood(
            ReactionEntry({'id': 'r1', 'equation': eq1}),
            ReactionEntry({'id': 'r2', 'equation': eq2}),
            cpd_map, cpd_score,
            {'c': 's'}
        )
        self.assertGreater(m, un)


class Test_bayesian_util(unittest.TestCase):
    def test_name_equals(self):
        self.assertFalse(util.name_equals('Test', None))

    def test_name_similar(self):
        score = util.name_similar('kit-ten', 'sitt(i)ng')
        self.assertEqual(score, 1 - 3.0 / 7)
        score = util.name_similar('kitten', None)
        self.assertEqual(score, 0.01)
        score = util.name_similar('kitten', '')
        self.assertEqual(score, 0.01)

    def test_fomula_equals(self):
        self.assertFalse(util.formula_equals(None, None, 0, 1))

    def test_formula_exact(self):
        self.assertTrue(
            util.formula_exact(Formula.parse('C3H4O5'),
                               Formula.parse('C3H4O5'))
        )
        self.assertFalse(
            util.formula_exact(Formula.parse('C2H3O6'), None)
        )
        self.assertFalse(
            util.formula_exact(Formula.parse('C2H3O6'),
                               Formula.parse('C5H11'))
        )

    def test_pairwise_distance(self):
        for pair, score in util.pairwise_distance(['kitten'], ['sitting'],
                                                  util.jaccard):
            self.assertTupleEqual(
                (pair, score),
                (('kitten', 'sitting'), 3.0 / 7)
            )

    def test_gene_equals(self):
        g1 = '(A1 and B21) or (A12 and B2)'
        g2 = 'A12 and B2'
        g3 = 'A1 and B2'
        g4 = '(C and D) or (E and F)'
        gmap = {'A12': 'D', 'B2': 'C', 'C': 'B2', 'D': 'A12'}
        self.assertTrue(util.genes_equals(g1, g2))
        self.assertFalse(util.genes_equals(g1, g3))
        self.assertFalse(util.genes_equals(g1, None))
        self.assertTrue(util.genes_equals(g1, g4, gmap))
