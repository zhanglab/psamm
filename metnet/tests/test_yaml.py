#!/usr/bin/env python

import unittest

from metnet.datasource import native
from metnet.reaction import Reaction, Compound

class TestYAMLDataSource(unittest.TestCase):
    def test_parse_reaction_list(self):
        reactions = list(native.parse_reaction_list('./test.yaml', [
            {
                'id': 'rxn1',
                'reversible': True,
                'left': [
                    { 'id': 'A', 'value': 1 },
                    { 'id': 'B', 'value': 2 } ],
                'right': [
                    { 'id': 'C', 'value': 1 }
                ]
            }
        ]))

        self.assertEquals(len(reactions), 1)

        reaction = Reaction(Reaction.Bidir,
                            [(Compound('A'), 1), (Compound('B'), 2)],
                            [(Compound('C'), 1)])
        self.assertEquals(reactions[0].equation, reaction)

    def test_parse_reaction_list_missing_value(self):
        with self.assertRaises(native.ParseError):
            reactions = list(native.parse_reaction_list('./test.yaml', [
                {
                    'id': 'rxn1',
                    'left': [
                        { 'id': 'A' }
                    ]
                }
            ]))


if __name__ == '__main__':
    unittest.main()
