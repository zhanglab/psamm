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
import math
from decimal import Decimal
from collections import OrderedDict

from psamm.datasource import native, context, entry
from psamm.reaction import Reaction, Compound, Direction
from psamm.formula import Formula

from six import StringIO
import yaml


class TestYAMLDataSource(unittest.TestCase):
    def test_parse_reaction(self):
        reaction = native.parse_reaction({
            'id': 'reaction_123',
            'equation': 'A + 2 B => C[e]'
        }, default_compartment='p')

        self.assertEqual(reaction.id, 'reaction_123')
        self.assertEqual(reaction.equation, Reaction(
            Direction.Forward,
            [(Compound('A', 'p'), 1), (Compound('B', 'p'), 2)],
            [(Compound('C', 'e'), 1)]))

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

        reaction = Reaction(Direction.Both,
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
        }, 'e'))

        self.assertEqual(len(medium), 5)
        self.assertEqual(medium[0], (Compound('ac', 'e'), None, None, None))
        self.assertEqual(medium[1], (Compound('glcD', 'e'), None, -10, None))
        self.assertEqual(medium[2], (Compound('co2', 'e'), None, None, 50))
        self.assertEqual(
            medium[3], (Compound('compound_x', 'c'), None, None, None))
        self.assertEqual(
            medium[4], (Compound('compound_y', 'e'), 'EX_cpdy', None, None))

    def test_parse_normal_float(self):
        v = native.yaml_load('-23.456')
        self.assertEqual(v, Decimal('-23.456'))
        self.assertIsInstance(v, Decimal)

    def test_parse_special_float(self):
        self.assertEqual(native.yaml_load('.inf'), Decimal('Infinity'))
        self.assertEqual(native.yaml_load('-.inf'), -Decimal('Infinity'))
        self.assertTrue(math.isnan(native.yaml_load('.nan')))


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

    def test_long_string_model(self):
        long_string = '''---
            name: Test model
            biomass: rxn_1
            reactions:
              - id: rxn_1
                equation: '|A_\u2206[e]| => |B[c]|'
                genes:
                  - gene_1
                  - gene_2
              - id: rxn_2_\u03c0
                equation: '|B[c]| => |C[e]|'
                genes: 'gene_3 or (gene_4 and gene_5)'
            compounds:
              - id: A_\u2206
              - id: B
              - id: C
            media:
              - compartment: e
                compounds:
                  - id: A_\u2206
                  - id: C
            limits:
              - reaction: rxn_2_\u03c0
                upper: 100
            '''
        m = native.NativeModel(long_string)
        self.assertEqual(m.name, 'Test model')
        self.assertEqual(m.biomass_reaction, 'rxn_1')
        self.assertEqual(m.extracellular_compartment, 'e')
        reactions = list(m.parse_reactions())
        self.assertEqual(reactions[0].id, 'rxn_1')

    def test_long_string_model_include(self):
        longString = '''---
            name: Test model
            biomass: rxn_1
            reactions:
              - include: medium.yaml
            '''
        m = native.NativeModel(longString)
        with self.assertRaises(context.ContextError):
            list(m.parse_reactions())

    def test_dict_model(self):
        dict_model = {
            'name': 'Test model',
            'biomass': 'rxn_1',
            'reactions': [
                {'id': 'rxn_1',
                 'equation': '|A_\u2206[e]| => |B[c]|',
                 'genes': [
                    'gene_1',
                    'gene_2']
                },
                {'id': 'rxn_2_\u03c0',
                 'equation': '|B[c]| => |C[e]|',
                 'genes': 'gene_3 or (gene_4 and gene_5)'
                }
                ],
            'compounds': [
                {'id': 'A_\u2206'},
                {'id': 'B'},
                {'id': 'C'}
              ],
            'media': [
                {'compartment': 'e',
                 'compounds':[
                    {'id': 'A_\u2206'},
                    {'id': 'C'}]
                }],
            'limits':[
                {'reaction': 'rxn_2_\u03c0',
                 'upper': 100
                }]
            }

        dmodel = native.NativeModel(dict_model)
        self.assertEqual(dmodel.name, 'Test model')
        self.assertEqual(dmodel.biomass_reaction, 'rxn_1')
        self.assertEqual(dmodel.extracellular_compartment, 'e')

    def test_parse_model_file(self):
        path = self.write_model_file('model.yaml', '\n'.join([
            'name: Test model',
            'biomass: biomass_reaction_id',
            'extracellular: Extra',
            'default_flux_limit: 500'
        ]))

        model = native.NativeModel.load_model_from_path(path)

        self.assertEqual(model.name, 'Test model')
        self.assertEqual(model.biomass_reaction, 'biomass_reaction_id')
        self.assertEqual(model.extracellular_compartment, 'Extra')
        self.assertEqual(model.default_flux_limit, 500)

    def test_bad_path(self):
        with self.assertRaises(native.ParseError):
            native.NativeModel.load_model_from_path('/nope/nreal/path')

    def test_invalid_model_type(self):
        with self.assertRaises(ValueError):
            native.NativeModel(42.2)

    def test_parse_model_file_with_media(self):
        path = self.write_model_file('model.yaml', '\n'.join([
            'extracellular: Ex',
            'media:',
            ' - compounds:',
            '    - id: A',
            '    - id: B',
            '      compartment: c'
        ]))

        model = native.NativeModel.load_model_from_path(path)

        medium = list(model.parse_medium())
        self.assertEqual(medium[0][0], Compound('A', 'Ex'))
        self.assertEqual(medium[1][0], Compound('B', 'c'))

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
        self.assertEqual(compounds[0].properties['kegg'], 'C00001')
        self.assertEqual(compounds[0].properties['cas'], '7732-18-5')

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
            Direction.Forward, [(Compound('A'), 1)], [(Compound('B'), 2)]))
        self.assertEqual(reactions[0].properties.get('ec'), None)
        self.assertEqual(reactions[0].genes, None)
        self.assertEqual(reactions[0].filemark.filecontext.filepath, path)
        self.assertEqual(reactions[0].filemark.line, 2)

        self.assertEqual(reactions[1].id, 'rxn_2')
        self.assertEqual(reactions[1].name, 'Second reaction')
        self.assertEqual(reactions[1].equation, Reaction(
            Direction.Both, [], [(Compound('C'), 1)]))
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

        medium = list(native.parse_medium_file(path, 'e'))
        self.assertEqual(medium, [
            (Compound('cpd_A', 'c'), None, None, None),
            (Compound('cpd_B', 'e'), None, -1000, None),
            (Compound('cpd_C', 'e'), None, None, 20),
            (Compound('cpd_D', 'e'), None, -100, -10)
        ])

    def test_get_limits_invalid_fixed(self):
        d = {
            'fixed' : 10,
            'upper' : 20
            }
        with self.assertRaises(native.ParseError):
            native.get_limits(d)

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
            '    compartment: c',
            '  - id: cpd_E',
            '    fixed: 100.0',
        ]))

        medium = list(native.parse_medium_file(path, 'e'))
        self.assertEqual(medium, [
            (Compound('cpd_A', 'e'), 'EX_A', -40, None),
            (Compound('cpd_B', 'e'), None, None, 100),
            (Compound('cpd_C', 'e'), None, -100, 500),
            (Compound('cpd_D', 'c'), None, None, None),
            (Compound('cpd_E', 'e'), None, 100, 100)
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
        ], 'e'))

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
            '  fixed: 10'
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
            ('rxn_4', 10, 10)
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


class TestCheckId(unittest.TestCase):
    def test_check_id_none(self):
        with self.assertRaises(native.ParseError):
            native._check_id(None, 'Compound')

    def test_check_id_not_string(self):
        with self.assertRaises(native.ParseError):
            native._check_id(False, 'Compound')

    def test_check_id_empty_string(self):
        with self.assertRaises(native.ParseError):
            native._check_id('', 'Reaction')

    def test_check_id_success(self):
        native._check_id(u'\u222b', 'Compound')


class TestNativeModelWriter(unittest.TestCase):
    def setUp(self):
        self.writer = native.ModelWriter()

    def test_convert_compound_entry(self):
        compound = entry.DictCompoundEntry({
            'id': 'c1',
            'name': 'Compound 1',
            'charge': 2,
            'formula': None,
            'custom': 'ABC'
        })
        d = self.writer.convert_compound_entry(compound)
        self.assertIsInstance(d, OrderedDict)
        self.assertEqual(d['id'], 'c1')
        self.assertEqual(d['name'], 'Compound 1')
        self.assertEqual(d['charge'], 2)
        self.assertEqual(d['custom'], 'ABC')
        self.assertNotIn('formula', d)

    def test_convert_reaction_entry(self):
        equation = Reaction(
            Direction.Both, {Compound('c1', 'c'): -1, Compound('c2', 'c'): 1})
        reaction = entry.DictReactionEntry({
            'id': 'rxn01234',
            'name': 'Test reaction',
            'equation': equation,
            'custom_property': -34,
            'another_custom_property': None
        })
        d = self.writer.convert_reaction_entry(reaction)
        self.assertIsInstance(d, OrderedDict)
        self.assertEqual(d['id'], 'rxn01234')
        self.assertEqual(d['name'], 'Test reaction')
        self.assertEqual(d['equation'], equation)
        self.assertEqual(d['custom_property'], -34)
        self.assertNotIn('another_custom_property', d)

    def test_write_compounds(self):
        stream = StringIO()
        compounds = [
            entry.DictCompoundEntry({
                'id': 'c1',
                'name': 'Compound 1',
                'charge': 2,
                'custom': 34.5
            }),
            entry.DictCompoundEntry({
                'id': 'c2',
                'name': 'Compound 2',
                'formula': Formula.parse('H2O')
            })
        ]
        self.writer.write_compounds(stream, compounds)

        self.assertEqual(yaml.safe_load(stream.getvalue()), [
            {
                'id': 'c1',
                'name': 'Compound 1',
                'charge': 2,
                'custom': 34.5
            },
            {
                'id': 'c2',
                'name': 'Compound 2',
                'formula': 'H2O'
            }
        ])

    def test_write_compounds_with_properties(self):
        stream = StringIO()
        compounds = [
            entry.DictCompoundEntry({
                'id': 'c1',
                'name': 'Compound 1',
                'charge': 2
            }),
            entry.DictCompoundEntry({
                'id': 'c2',
                'name': 'Compound 2',
                'formula': 'H2O'
            })
        ]
        self.writer.write_compounds(stream, compounds, properties={'name'})

        self.assertEqual(yaml.safe_load(stream.getvalue()), [
            {
                'id': 'c1',
                'name': 'Compound 1'
            },
            {
                'id': 'c2',
                'name': 'Compound 2'
            }
        ])

    def test_write_reactions(self):
        stream = StringIO()
        reactions = [
            entry.DictReactionEntry({
                'id': 'r1',
                'name': 'Reaction 1',
                'equation': Reaction(Direction.Both, {})
            }),
            entry.DictReactionEntry({
                'id': 'r2',
                'name': 'Reaction 2',
                'equation': Reaction(Direction.Forward, {
                    Compound('c1', 'c'): -1,
                    Compound('c2', 'c'): 2
                })
            })
        ]
        self.writer.write_reactions(stream, reactions)

        # The reaction equation for r1 is invalid (no compounds) and is
        # therefore skipped.
        self.assertEqual(yaml.safe_load(stream.getvalue()), [
            {
                'id': 'r1',
                'name': 'Reaction 1'
            },
            {
                'id': 'r2',
                'name': 'Reaction 2',
                'equation': '|c1[c]| => (2) |c2[c]|'
            }
        ])
