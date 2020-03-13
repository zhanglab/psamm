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
# Copyright 2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest

from psamm.datasource import entry


class TestDictEntries(unittest.TestCase):
    def test_create_compound_dict_entry(self):
        props = {
            'id': 'test_compound',
            'name': 'Test compound',
            'formula': 'H2O',
            'charge': 1,
            'custom': 'ABC'
        }
        e = entry.DictCompoundEntry(props)
        self.assertIsInstance(e, entry.CompoundEntry)
        self.assertEqual(e.id, 'test_compound')
        self.assertEqual(e.name, 'Test compound')
        self.assertEqual(e.formula, 'H2O')
        self.assertEqual(e.charge, 1)
        self.assertEqual(e.properties.get('custom'), 'ABC')

    def test_create_compound_dict_entry_copy(self):
        props = {
            'id': 'test_compound',
            'name': 'Test compound'
        }
        e = entry.DictCompoundEntry(props)
        e2 = entry.DictCompoundEntry(e)
        self.assertIsInstance(e2, entry.CompoundEntry)
        self.assertEqual(e2.id, 'test_compound')
        self.assertEqual(e2.name, 'Test compound')

    def test_create_compound_dict_entry_without_id(self):
        props = {
            'name': 'Compound 1',
            'formula': 'CO2'
        }
        with self.assertRaises(ValueError):
            e = entry.DictCompoundEntry(props)

    def test_create_compound_dict_entry_with_id_override(self):
        props = {
            'id': 'old_id',
            'name': 'Compound 1',
            'formula': 'CO2'
        }
        e = entry.DictCompoundEntry(props, id='new_id')
        self.assertEqual(e.id, 'new_id')
        self.assertEqual(e.properties['id'], 'new_id')

    def test_use_compound_dict_entry_setters(self):
        e = entry.DictCompoundEntry({}, id='new_id')
        e.formula = 'CO2'
        e.name = 'Compound 1'
        e.charge = 5
        self.assertEqual(e.formula, 'CO2')
        self.assertEqual(e.properties['formula'], 'CO2')
        self.assertEqual(e.name, 'Compound 1')
        self.assertEqual(e.charge, 5)

    def test_create_reaction_entry(self):
        props = {
            'id': 'reaction_1',
            'name': 'Reaction 1'
        }
        e = entry.DictReactionEntry(props)
        self.assertIsInstance(e, entry.ReactionEntry)
        self.assertEqual(e.id, 'reaction_1')
        self.assertEqual(e.name, 'Reaction 1')

    def test_create_reaction_entry_from_compound_entry(self):
        e = entry.DictCompoundEntry({
            'id': 'compound_1',
            'name': 'Compound 1'
        })
        with self.assertRaises(ValueError):
            e2 = entry.DictReactionEntry(e)

    def test_use_reaction_dict_entry_setters(self):
        e = entry.DictReactionEntry({}, id='reaction_1')
        e.name = 'Reaction 1'
        e.equation = 'A => B'
        e.genes = 'gene_1 and gene_2'
        self.assertEqual(e.name, 'Reaction 1')
        self.assertEqual(e.equation, 'A => B')
        self.assertEqual(e.genes, 'gene_1 and gene_2')
        self.assertEqual(e.properties['genes'], 'gene_1 and gene_2')

    def test_create_compartment_entry(self):
        props = {
            'id': 'c',
            'name': 'Cytosol'
        }
        e = entry.DictCompartmentEntry(props)
        self.assertIsInstance(e, entry.CompartmentEntry)
        self.assertEqual(e.id, 'c')
        self.assertEqual(e.name, 'Cytosol')


if __name__ == '__main__':
    unittest.main()
