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


if __name__ == '__main__':
    unittest.main()
