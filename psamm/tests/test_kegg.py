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

from psamm.datasource import kegg
from psamm.reaction import Reaction, Compound, Direction
from psamm.expression.affine import Expression

from six import StringIO


class TestKEGGEntryParser(unittest.TestCase):
    def setUp(self):
        self.f = StringIO('\n'.join([
            'ENTRY     ID001',
            'NAME      Test entry',
            'PROPERTY  This is a multi-line',
            '          property!',
            '///',
            'ENTRY     ID002',
            'NAME      Another entry',
            'PROPERTY  Single line property',
            'REFS      ref1: abcdef',
            '          ref2: defghi',
            '///'
        ]))

    def test_parse_entries(self):
        entries = list(kegg.parse_kegg_entries(self.f))
        self.assertEqual(len(entries), 2)

        self.assertEqual(entries[0].properties, {
            'entry': ['ID001'],
            'name': ['Test entry'],
            'property': ['This is a multi-line', 'property!']
        })

        self.assertEqual(entries[1].properties, {
            'entry': ['ID002'],
            'name': ['Another entry'],
            'property': ['Single line property'],
            'refs': [
                'ref1: abcdef',
                'ref2: defghi'
            ]
        })


class TestKEGGCompoundEntry(unittest.TestCase):
    def test_minimal_compound_entry(self):
        c = kegg.CompoundEntry(kegg.KEGGEntry({
            'entry': ['C00001        Compound']
        }))

        self.assertEqual(c.id, 'C00001')
        self.assertIsNone(c.properties['name'])
        self.assertIsNone(c.filemark)

    def test_complete_compound_entry(self):
        c = kegg.CompoundEntry(kegg.KEGGEntry({
            'entry': ['C00001        Compound'],
            'name': ['H2O;', 'Water'],
            'reaction': ['R00001 R00002', 'R00003'],
            'enzyme': ['1.2.3.4 2.3.4.5', '7.6.50.4 2.1.-,-'],
            'formula': ['H2O'],
            'exact_mass': ['18.01'],
            'mol_weight': ['18.01'],
            'pathway': [
                'map00001  First pathway',
                'map00002  Second pathway'
            ],
            'dblinks': [
                'CAS: 12345',
                'ChEBI: B2345'
            ],
            'comment': ['This information is purely for testing!']
        }))

        self.assertEqual(c.id, 'C00001')
        self.assertEqual(c.properties['name'], 'H2O')
        self.assertEqual(c.properties['names'], ['H2O', 'Water'])
        self.assertEqual(
            c.properties['reactions'], ['R00001', 'R00002', 'R00003'])
        self.assertEqual(c.properties['enzymes'], [
            '1.2.3.4', '2.3.4.5', '7.6.50.4', '2.1.-,-'])
        self.assertEqual(c.properties['formula'], 'H2O')
        self.assertAlmostEqual(c.properties['exact_mass'], 18.01)
        self.assertAlmostEqual(c.properties['mol_weight'], 18.01)
        self.assertEqual(c.properties['pathways'], [
            ('map00001', 'First pathway'),
            ('map00002', 'Second pathway')
        ])
        self.assertEqual(c.properties['dblinks'], [
            ('CAS', '12345'), ('ChEBI', 'B2345')
        ])
        self.assertEqual(
            c.properties['comment'], 'This information is purely for testing!')


class TestKEGGReactionParser(unittest.TestCase):
    def test_kegg_parse(self):
        r = kegg.parse_reaction('C00013 + C00001 <=> 2 C00009')
        self.assertEqual(
            r, Reaction(Direction.Both,
                        [(Compound('C00013'), 1), (Compound('C00001'), 1)],
                        [(Compound('C00009'), 2)]))

    def test_kegg_parse_with_count_expression(self):
        r = kegg.parse_reaction('2n C00404 + n C00001 <=> (n+1) C02174')
        self.assertEqual(
            r, Reaction(Direction.Both,
                        [(Compound('C00404'), Expression('2n')),
                         (Compound('C00001'), Expression('n'))],
                        [(Compound('C02174'), Expression('n+1'))]))

    def test_kegg_parse_with_compound_argument(self):
        r = kegg.parse_reaction('C00039(n) <=> C00013 + C00039(n+1)')
        self.assertEqual(
            r, Reaction(Direction.Both,
                        [(Compound('C00039', arguments=[Expression('n')]), 1)],
                        [(Compound('C00013'), 1),
                         (Compound('C00039',
                                   arguments=[Expression('n+1')]), 1)]))


if __name__ == '__main__':
    unittest.main()
