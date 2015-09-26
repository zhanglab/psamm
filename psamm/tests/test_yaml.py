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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import os
import shutil
import tempfile
import unittest

from psamm.datasource import native
from psamm.reaction import Reaction, Compound

from six import StringIO


class TestYAMLDataSource(unittest.TestCase):
    def test_parse_reaction_list(self):
        reactions = list(native.parse_reaction_list('./test.yaml', [
            {
                'id': 'rxn1',
                'equation': {
                    'reversible': True,
                    'left': [
                        { 'id': 'A', 'value': 1 },
                        { 'id': 'B', 'value': 2 } ],
                    'right': [
                        { 'id': 'C', 'value': 1 }
                    ]
                }
            }
        ]))

        self.assertEqual(len(reactions), 1)

        reaction = Reaction(Reaction.Bidir,
                            [(Compound('A'), 1), (Compound('B'), 2)],
                            [(Compound('C'), 1)])
        self.assertEqual(reactions[0].equation, reaction)

    def test_parse_reaction_list_missing_value(self):
        with self.assertRaises(native.ParseError):
            reactions = list(native.parse_reaction_list('./test.yaml', [
                {
                    'id': 'rxn1',
                    'equation': {
                        'left': [
                            { 'id': 'A' }
                        ]
                    }
                }
            ]))

    def test_parse_medium_table(self):
        table = '''
ac      e
glcD    e       -10
co2     e       -       50
'''

        medium = list(native.parse_medium_table_file(StringIO(table.strip())))
        self.assertEqual(len(medium), 3)
        self.assertEqual(medium[0], (Compound('ac', 'e'), None, None, None))
        self.assertEqual(medium[1], (Compound('glcD', 'e'), None, -10, None))
        self.assertEqual(medium[2], (Compound('co2', 'e'), None, None, 50))

    def test_parse_medium(self):
        medium = list(native.parse_medium({
            'compartment': 'e',
            'compounds': [
                {'id': 'ac'},
                {'id': 'glcD', 'lower': -10},
                {'id': 'co2', 'upper': 50},
                {'id': 'compound_x', 'compartment': 'c'},
                {'id': 'compound_y', 'reaction': 'EX_cpdy'}
            ]
        }))

        self.assertEqual(len(medium), 5)
        self.assertEqual(medium[0], (Compound('ac', 'e'), None, None, None))
        self.assertEqual(medium[1], (Compound('glcD', 'e'), None, -10, None))
        self.assertEqual(medium[2], (Compound('co2', 'e'), None, None, 50))
        self.assertEqual(
            medium[3], (Compound('compound_x', 'c'), None, None, None))
        self.assertEqual(
            medium[4], (Compound('compound_y', 'e'), 'EX_cpdy', None, None))


class TestYAMLFileSystemData(unittest.TestCase):
    """Test loading files from file system."""

    def setUp(self):
        self._model_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def write_model_file(self, filename, contents):
        path = os.path.join(self._model_dir, filename)
        with open(path, 'w') as f:
            f.write(contents)
        return path

    def test_parse_compound_tsv_file(self):
        path = self.write_model_file('compounds.tsv', '\n'.join([
            'id\tname\tformula\tcharge\tkegg\tcas',
            'h2o\tWater\tH2O\t0\tC00001\t7732-18-5'
        ]))

        compounds = list(native.parse_compound_file(path, 'tsv'))
        self.assertEqual(len(compounds), 1)

        self.assertEqual(compounds[0].id, 'h2o')
        self.assertEqual(compounds[0].name, 'Water')
        self.assertEqual(compounds[0].properties['name'], 'Water')
        self.assertEqual(compounds[0].formula, 'H2O')
        self.assertEqual(compounds[0].charge, 0)
        self.assertEqual(compounds[0].kegg, 'C00001')
        self.assertEqual(compounds[0].cas, '7732-18-5')

    def test_parse_compound_yaml_file(self):
        path = self.write_model_file('compounds.yaml', '\n'.join([
            '- id: h2o',
            '  name: Water',
            '  formula: H2O',
            '  user_annotation: ABC'
        ]))

        compounds = list(native.parse_compound_file(path, 'yaml'))
        self.assertEqual(len(compounds), 1)

        self.assertEqual(compounds[0].id, 'h2o')
        self.assertEqual(compounds[0].name, 'Water')
        self.assertEqual(compounds[0].formula, 'H2O')
        self.assertEqual(compounds[0].properties['user_annotation'], 'ABC')

    def test_parse_reaction_table_file(self):
        path = self.write_model_file('reactions.tsv', '\n'.join([
            'id\tname\tequation',
            'rxn_1\tReaction 1\t|A| => (2) |B|',
            'rxn_2\tSecond reaction\t<=> |C|'
        ]))

        reactions = list(native.parse_reaction_file(path))
        self.assertEqual(len(reactions), 2)

        self.assertEqual(reactions[0].id, 'rxn_1')
        self.assertEqual(reactions[0].name, 'Reaction 1')
        self.assertEqual(reactions[0].properties['name'], 'Reaction 1')
        self.assertEqual(reactions[0].equation, Reaction(
            Reaction.Right, [(Compound('A'), 1)], [(Compound('B'), 2)]))
        self.assertEqual(reactions[0].ec, None)
        self.assertEqual(reactions[0].genes, None)
        self.assertEqual(reactions[0].filemark.filecontext.filepath, path)
        self.assertEqual(reactions[0].filemark.line, 2)

        self.assertEqual(reactions[1].id, 'rxn_2')
        self.assertEqual(reactions[1].name, 'Second reaction')
        self.assertEqual(reactions[1].equation, Reaction(
            Reaction.Bidir, [], [(Compound('C'), 1)]))
        self.assertEqual(reactions[1].filemark.filecontext.filepath, path)
        self.assertEqual(reactions[1].filemark.line, 3)

    def test_parse_reaction_yaml_file(self):
        path = self.write_model_file('reaction.yaml', '\n'.join([
            '- id: rxn_1',
            '  equation: "|A| => |B|"',
            '- id: rxn_2',
            '  equation:',
            '    reversible: no',
            '    left:',
            '      - id: B',
            '        value: 1',
            '    right:',
            '      - id: A',
            '        value: 1'
        ]))

        reactions = list(native.parse_reaction_file(path))
        self.assertEqual(len(reactions), 2)

        self.assertEqual(reactions[0].id, 'rxn_1')
        self.assertEqual(reactions[1].id, 'rxn_2')

    def test_parse_medium_table_file(self):
        path = self.write_model_file('medium.tsv', '\n'.join([
            '',
            '# comment',
            'cpd_A\tc',
            'cpd_B\te\t-1000',
            'cpd_C\te\t-\t20  # line comment',
            'cpd_D\te\t-100\t-10'
        ]))

        medium = list(native.parse_medium_file(path))
        self.assertEqual(medium, [
            (Compound('cpd_A', 'c'), None, None, None),
            (Compound('cpd_B', 'e'), None, -1000, None),
            (Compound('cpd_C', 'e'), None, None, 20),
            (Compound('cpd_D', 'e'), None, -100, -10)
        ])

    def test_parse_medium_yaml_file(self):
        path = self.write_model_file('medium.yaml', '\n'.join([
            'compartment: e',
            'compounds:',
            '  - id: cpd_A',
            '    reaction: EX_A',
            '    lower: -40',
            '  - id: cpd_B',
            '    upper: 100',
            '  - id: cpd_C',
            '    lower: -100.0',
            '    upper: 500.0',
            '  - id: cpd_D',
            '    compartment: c'
        ]))

        medium = list(native.parse_medium_file(path))
        self.assertEqual(medium, [
            (Compound('cpd_A', 'e'), 'EX_A', -40, None),
            (Compound('cpd_B', 'e'), None, None, 100),
            (Compound('cpd_C', 'e'), None, -100, 500),
            (Compound('cpd_D', 'c'), None, None, None)
        ])

    def test_parse_medium_yaml_list(self):
        self.write_model_file('medium.yaml', '\n'.join([
            'compartment: e',
            'compounds:',
            '  - id: cpd_A',
            '    lower: -42'
        ]))

        path = os.path.join(self._model_dir, 'fake.yaml')
        medium = list(native.parse_medium_list(path, [
            {'include': 'medium.yaml'},
            {
                'compartment': 'e',
                'compounds': [
                    {'id': 'cpd_B', 'upper': 767}
                ]
            }
        ]))

        self.assertEqual(medium, [
            (Compound('cpd_A', 'e'), None, -42, None),
            (Compound('cpd_B', 'e'), None, None, 767)
        ])

    def test_parse_limits_table_file(self):
        path = self.write_model_file('limits_1.tsv', '\n'.join([
            '# comment',
            '',
            'rxn_1',
            'rxn_2\t-100',
            'rxn_3\t-\t3e-1',
            'rxn_4\t-1000\t-100  # line comment'
        ]))

        limits = list(native.parse_limits_file(path))
        self.assertEqual(limits[0], ('rxn_1', None, None))
        self.assertEqual(limits[1], ('rxn_2', -100, None))
        self.assertEqual(limits[2][0], 'rxn_3')
        self.assertEqual(limits[2][1], None)
        self.assertAlmostEqual(limits[2][2], 3e-1)
        self.assertEqual(limits[3], ('rxn_4', -1000, -100))

    def test_parse_limits_table_file_too_many_fields(self):
        path = self.write_model_file('limits.tsv', '\n'.join([
            'rxn_1\t-\t100\ttext']))
        with self.assertRaises(native.ParseError):
            list(native.parse_limits_file(path))

    def test_parse_limits_yaml_file(self):
        path = self.write_model_file('limits_1.yaml', '\n'.join([
            '- include: limits_2.yml',
            '- reaction: rxn_3',
            '  lower: -1000',
            '- reaction: rxn_4',
            '  upper: 0'
        ]))

        self.write_model_file('limits_2.yml', '\n'.join([
            '- reaction: rxn_1',
            '  lower: -1',
            '  upper: 25.5',
            '- reaction: rxn_2'
        ]))

        limits = list(native.parse_limits_file(path))
        self.assertEqual(limits, [
            ('rxn_1', -1, 25.5),
            ('rxn_2', None, None),
            ('rxn_3', -1000, None),
            ('rxn_4', None, 0)
        ])

    def test_parse_model_table_file(self):
        path = self.write_model_file('model_1.tsv', '\n'.join([
            '# comment',
            'rxn_1',
            'rxn_2',
            'rxn_3',
            'rxn_4  # line comment']))

        reactions = list(native.parse_model_file(path))
        self.assertEqual(reactions, ['rxn_1', 'rxn_2', 'rxn_3', 'rxn_4'])

    def test_parse_model_yaml_file(self):
        path = self.write_model_file('model_1.yaml', '''---
            - include: model_2.yaml
            - reactions:
              - rxn_3
              - rxn_4''')

        self.write_model_file('model_2.yaml', '''---
            - groups:
              - name: First group
                reactions: [rxn_1]
              - name: Second group
                reactions: [rxn_2]''')

        reactions = list(native.parse_model_file(path))
        self.assertEqual(reactions, ['rxn_1', 'rxn_2', 'rxn_3', 'rxn_4'])


if __name__ == '__main__':
    unittest.main()
