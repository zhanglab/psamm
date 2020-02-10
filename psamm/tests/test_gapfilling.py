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
# Copyright 2016-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2016  Chao liu <lcddzyx@gmail.com>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest

from psamm.datasource import native
from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm.datasource.reaction import parse_reaction
from psamm import gapfilling


class TestAddReactions(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) A[c]'))
        self.database.set_reaction('rxn_2', parse_reaction('A[c] <=> B[c]'))
        self.database.set_reaction('rxn_3', parse_reaction('A[c] => D[e]'))
        self.database.set_reaction('rxn_4', parse_reaction('A[c] => C[c]'))
        self.database.set_reaction('rxn_5', parse_reaction('C[c] => D[e]'))
        self.database.set_reaction('rxn_6', parse_reaction('D[e] =>'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)

    def test_add_all_database_reactions(self):
        # Should get added
        self.database.set_reaction('rxn_7', parse_reaction('D[c] => E[c]'))
        # Not added because of compartment
        self.database.set_reaction('rxn_8', parse_reaction('D[c] => E[p]'))
        added = gapfilling.add_all_database_reactions(self.model, {'c', 'e'})
        self.assertEqual(added, {'rxn_7'})
        self.assertEqual(set(self.model.reactions), {
            'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6', 'rxn_7'
        })

    def test_add_all_database_reactions_none(self):
        added = gapfilling.add_all_database_reactions(self.model, {'c', 'e'})
        self.assertEqual(added, set())
        self.assertEqual(set(self.model.reactions), {
            'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6'})

    def test_add_all_transport_reactions(self):
        added = gapfilling.add_all_transport_reactions(
            self.model, {('e', 'c')})
        for reaction in added:
            compartments = tuple(c.compartment for c, _ in
                                 self.model.get_reaction_values(reaction))
            self.assertEqual(len(compartments), 2)
            self.assertTrue('e' in compartments)
            self.assertTrue(compartments[0] != compartments[1])


class TestCreateExtendedModel(unittest.TestCase):
    def setUp(self):
        reader = native.ModelReader({
            'compartments': [
                {
                    'id': 'c',
                    'adjacent_to': 'e'
                }, {
                    'id': 'e'
                }
            ],
            'reactions': [
                {
                    'id': 'rxn_1',
                    'equation': 'A[e] <=> B[c]'
                }, {
                    'id': 'rxn_2',
                    'equation': 'A[e] => C[c]'
                }, {
                    'id': 'rxn_3',
                    'equation': 'A[e] => D[e]'
                }, {
                    'id': 'rxn_4',
                    'equation': 'C[c] => D[e]'
                }
            ],
            'compounds': [
                {'id': 'A'},
                {'id': 'B'},
                {'id': 'C'},
                {'id': 'D'}
            ],
            'exchange': [{
                'compartment': 'e',
                'compounds': [
                    {
                        'id': 'A',
                        'reaction': 'rxn_5'
                    },
                    {
                        'id': 'D',
                        'reaction': 'rxn_6'
                    }
                ]
            }],
            'model': [{
                'reactions': [
                    'rxn_1',
                    'rxn_2'
                ]
            }]
        })
        self._model = reader.create_model()

    def test_create_model_extended(self):
        expected_reactions = set([
            'rxn_1',
            'rxn_2',
            'rxn_3',
            'rxn_4',
            'rxn_5',
            'rxn_6',
            'TP_A[c]_A[e]',
            'TP_B[c]_B[e]',
            'TP_C[c]_C[e]',
            'TP_D[c]_D[e]',
            'EX_B[e]',
            'EX_C[e]',
        ])

        expected_weights = {
            'rxn_3': 5.6,
            'rxn_4': 1.0,
            'TP_A[c]_A[e]': 3.0,
            'TP_B[c]_B[e]': 3.0,
            'TP_C[c]_C[e]': 3.0,
            'TP_D[c]_D[e]': 3.0,
            'EX_B[e]': 2.0,
            'EX_C[e]': 2.0
        }
        penalties = {'rxn_3': 5.6}

        model_extended, weights = gapfilling.create_extended_model(
            self._model,
            db_penalty=1.0,
            ex_penalty=2.0,
            tp_penalty=3.0,
            penalties=penalties)
        self.assertEqual(set(model_extended.reactions), expected_reactions)
        self.assertEqual(weights, expected_weights)


if __name__ == '__main__':
    unittest.main()
