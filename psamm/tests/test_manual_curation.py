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

from psamm import manual_curation as curation
from psamm.datasource.native import ModelReader
from psamm.datasource.reaction import Compound


class TestCurator(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'compound.tsv'), 'w') as f:
            f.write('id1\tid2\tp\n')
            f.write('%s\t%s\t%f\n' % ('A', 'A1', 1.0))
            f.write('%s\t%s\t%f\n' % ('B', 'B1', 1.0))
            f.write('%s\t%s\t%f\n' % ('C', 'C1', 0.03))
            f.write('%s\t%s\t%f\n' % ('D', 'D1', 0.00005))
        with open(os.path.join(self._model_dir, 'reaction.tsv'), 'w') as f:
            f.write('id1\tid2\tp\n')
            f.write('%s\t%s\t%f\n' % ('rxn1', 'rxn1', 1.0))
            f.write('%s\t%s\t%f\n' % ('rxn2', 'rxn2', 1.0))
            f.write('%s\t%s\t%f\n' % ('rxn3', 'rxn3', 0.03))
            f.write('%s\t%s\t%f\n' % ('rxn4', 'rxn4', 0.00005))
        self.curator = curation.Curator(
            os.path.join(self._model_dir, 'compound.tsv'),
            os.path.join(self._model_dir, 'reaction.tsv'),
            os.path.join(self._model_dir, 'curated_compound.tsv'),
            os.path.join(self._model_dir, 'curated_reaction.tsv')
        )

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_add_mapping(self):
        self.curator.add_mapping(('A', 'A1'), 'c', True)
        self.curator.add_mapping(('D', 'D1'), 'c', False)
        self.curator.add_mapping(('rxn1', 'rxn1'), 'r', True)
        self.curator.add_mapping(('rxn4', 'rxn4'), 'r', False)
        self.curator.save()

    def test_compound_checked(self):
        self.assertFalse(self.curator.compound_checked(('B', 'B1')))
        self.curator.add_mapping(('B', 'B1'), 'c', True)
        self.assertFalse(self.curator.compound_checked(('C', 'C1')))
        self.curator.save()
        curator = self.curator = curation.Curator(
            self.curator.compound_map_file,
            self.curator.reaction_map_file,
            self.curator.curated_compound_map_file,
            self.curator.curated_reaction_map_file
        )
        self.assertTrue(curator.compound_checked(('B', 'B1')))
        self.assertTrue(curator.compound_checked('B'))

    def test_compound_ignored(self):
        self.assertFalse(self.curator.compound_ignored('C'))
        self.curator.add_ignore('C', 'c')
        self.assertTrue(self.curator.compound_ignored('C'))
        self.curator.save()
        curator = self.curator = curation.Curator(
            self.curator.compound_map_file,
            self.curator.reaction_map_file,
            self.curator.curated_compound_map_file,
            self.curator.curated_reaction_map_file
        )
        self.assertTrue(curator.compound_ignored('C'))

    def test_reaction_checked(self):
        self.assertFalse(self.curator.reaction_checked(('rxn2', 'rxn2')))
        self.curator.add_mapping(('rxn2', 'rxn2'), 'r', False)
        self.assertFalse(self.curator.reaction_checked(('rxn3', 'rxn3')))
        self.curator.save()
        curator = self.curator = curation.Curator(
            self.curator.compound_map_file,
            self.curator.reaction_map_file,
            self.curator.curated_compound_map_file,
            self.curator.curated_reaction_map_file
        )
        self.assertTrue(curator.reaction_checked(('rxn2', 'rxn2')))

    def test_reaction_ignored(self):
        self.assertFalse(self.curator.reaction_ignored('rxn3'))
        self.curator.add_ignore('rxn3', 'r')
        self.assertTrue(self.curator.reaction_ignored('rxn3'))
        self.curator.save()
        curator = self.curator = curation.Curator(
            self.curator.compound_map_file,
            self.curator.reaction_map_file,
            self.curator.curated_compound_map_file,
            self.curator.curated_reaction_map_file
        )
        self.assertTrue(curator.reaction_ignored('rxn3'))


class TestSearch(unittest.TestCase):
    def setUp(self):
        self._model_dir = tempfile.mkdtemp()
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
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
            ]))
        self._model = ModelReader.reader_from_path(
            os.path.join(self._model_dir, 'model.yaml')).create_model()

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_search_compound(self):
        curation.search_compound(self._model, ['C', 'Z'])

    def test_search_reaction(self):
        cpds = curation.search_reaction(self._model, ['rxn_1', 'test'])
        self.assertListEqual(
            next(cpds),
            [Compound('A', 'e'), Compound('B', 'c')])

    def test_filter_search_term(self):
        self.assertEqual(curation.filter_search_term('Test-(tes)_t'),
                         'testtest')
