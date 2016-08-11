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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2016  Chao Liu <lcddzyx@gmail.com>

import os
import unittest
import tempfile
import shutil

from psamm.datasource.native import NativeModel
from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm.datasource.reaction import parse_reaction
from psamm.expression import boolean
from psamm.lpsolver import generic
from psamm import fluxanalysis, randomsparse


class TestGetGeneAssociation(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('\n'.join([
                '---',
                'reactions:',
                '  - id: rxn_1',
                '    equation: A[e] <=> B[c]',
                '    genes: gene_1 and gene_2',
                '  - id: rxn_2',
                '    equation: A[e] => C[c]',
                '    genes: gene_3',
                '  - id: rxn_3',
                '    equation: A[e] => D[e]',
                '  - id: rxn_4',
                '    equation: C[c] => D[e]',
                '    genes: [gene_4, gene_5, gene_6]',
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
                '     - rxn_3',
                '     - rxn_4',
            ]))
        self._model = NativeModel.load_model_from_path(self._model_dir)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_get_gene_association(self):
        expected_association = {
            'rxn_1': boolean.Expression('gene_1 and gene_2'),
            'rxn_2': boolean.Expression('gene_3'),
            'rxn_4': boolean.Expression('gene_4 and gene_6 and gene_5')
        }

        self.assertEqual(
            dict(randomsparse.get_gene_associations(self._model)),
            expected_association)


class TestGetExchangeReactions(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('|A| <=>'))
        self._database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self._database.set_reaction('rxn_3', parse_reaction('|C| <=> |B|'))
        self._database.set_reaction('rxn_4', parse_reaction('|C| <=>'))
        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)

    def test_get_exchange_reactions(self):
        expected_reactions = {'rxn_1', 'rxn_4'}
        self.assertEqual(
            set(randomsparse.get_exchange_reactions(self._mm)),
            expected_reactions)


class TestReactionDeletionStrategy(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('|A| <=>'))
        self._database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self._database.set_reaction('rxn_3', parse_reaction('|C| <=> |B|'))
        self._database.set_reaction('rxn_4', parse_reaction('|C| <=>'))
        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)

        self._strategy = randomsparse.ReactionDeletionStrategy(self._mm)

    def test_method_get_all(self):
        expected_total = {'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4'}
        self.assertEqual(set(self._strategy.entities), expected_total)

    def test_method_tests(self):
        expected_reactions = {
            'rxn_1': {'rxn_1'},
            'rxn_2': {'rxn_2'},
            'rxn_3': {'rxn_3'},
            'rxn_4': {'rxn_4'}
        }
        self.assertEqual(dict(self._strategy.iter_tests()), expected_reactions)


class TestGeneDeletionStrategy(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('|A| <=>'))
        self._database.set_reaction('rxn_2', parse_reaction('|A| <=> |B|'))
        self._database.set_reaction('rxn_3', parse_reaction('|C| <=> |B|'))
        self._database.set_reaction('rxn_4', parse_reaction('|C| <=>'))
        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)
        self._assoc = {'rxn_2': boolean.Expression('gene_1 or gene_2')}

        self._strategy = randomsparse.GeneDeletionStrategy(
            self._mm, self._assoc)

    def test_method_get_all(self):
        expected_total = {'gene_1', 'gene_2'}
        self.assertEqual(set(self._strategy.entities), expected_total)

    def test_method_tests_and_delete(self):
        expected_genes = {'gene_1', 'gene_2'}
        expected_reaction_set = {'rxn_2'}
        test_dict = {}

        for i, (gene, deleted_reac) in enumerate(self._strategy.iter_tests()):
            self._strategy.delete(gene, deleted_reac)
            if i == 0:
                self.assertEqual(deleted_reac, set())
            else:
                self.assertEqual(deleted_reac, {'rxn_2'})
            test_dict[gene] = deleted_reac

        self.assertTrue(all(x in test_dict for x in expected_genes))
        self.assertTrue(expected_reaction_set in test_dict.values())


class TestRandomSparse(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('A[e] => B[e]'))
        self._database.set_reaction('rxn_2', parse_reaction('B[e] => C[e]'))
        self._database.set_reaction('rxn_3', parse_reaction('B[e] => D[e]'))
        self._database.set_reaction('rxn_4', parse_reaction('C[e] => E[e]'))
        self._database.set_reaction('rxn_5', parse_reaction('D[e] => E[e]'))
        self._database.set_reaction('rxn_6', parse_reaction('E[e] =>'))
        self._database.set_reaction('ex_A', parse_reaction('A[e] <=>'))

        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)
        self._assoc = {
            'rxn_1': boolean.Expression('gene_1'),
            'rxn_2': boolean.Expression('gene_2 or gene_3'),
            'rxn_5': boolean.Expression('gene_3 and gene_4')
        }
        self._obj_reaction = 'rxn_6'

        try:
            self._solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

        self._prob = fluxanalysis.FluxBalanceProblem(self._mm, self._solver)

    def test_random_sparse_reaction_strategy(self):
        expected = [
            ({'rxn_1', 'rxn_6', 'ex_A', 'rxn_2', 'rxn_4'}, {'rxn_3', 'rxn_5'}),
            ({'rxn_1', 'rxn_6', 'ex_A', 'rxn_3', 'rxn_5'}, {'rxn_2', 'rxn_4'})
        ]

        strategy = randomsparse.ReactionDeletionStrategy(self._mm)
        essential, deleted = randomsparse.random_sparse(
            strategy,
            self._prob,
            self._obj_reaction,
            flux_threshold=100)

        self.assertTrue((essential, deleted) in expected)

    def test_random_sparse_gene_strategy(self):
        expected = [
            ({'gene_1', 'gene_2'}, {'gene_3', 'gene_4'}),
            ({'gene_1', 'gene_3'}, {'gene_2', 'gene_4'})
        ]
        strategy = randomsparse.GeneDeletionStrategy(
            self._mm, self._assoc)
        essential, deleted = randomsparse.random_sparse(
            strategy,
            self._prob,
            self._obj_reaction,
            flux_threshold=100)

        self.assertTrue((essential, deleted) in expected)
