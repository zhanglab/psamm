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


if __name__ == '__main__':
    unittest.main()
