#!/usr/bin/env python

import unittest

from metnet import metabolicmodel
from metnet.reaction import ModelSEED

class TestMetabolicDatabase(unittest.TestCase):
    def setUp(self):
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))

    def test_database_get_reaction(self):
        reaction = ModelSEED.parse('|A| => |D|')
        self.assertEqual(self.database.get_reaction('rxn_3'), reaction)

class TestMetabolicModel(unittest.TestCase):
    def setUp(self):
        # TODO use mock database instead of actual database
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.load_model_from_file(iter(self.database.reactions))

    def test_reaction_set(self):
        self.assertEqual(self.model.reaction_set, { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_compound_set(self):
        self.assertEqual(self.model.compound_set, { ('A', None), ('B', None), ('C', None), ('D', None) })

    def test_matrix_get_item(self):
        self.assertEqual(self.model.matrix[('A', None), 'rxn_1'], 2)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[('B', None), 'rxn_2'], 1)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_4'], -1)
        self.assertEqual(self.model.matrix[('C', None), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[('C', None), 'rxn_5'], -1)
        self.assertEqual(self.model.matrix[('D', None), 'rxn_5'], 1)

    def test_matrix_iteritems(self):
        matrix_set = { ((('A', None), 'rxn_1'), 2),
                         ((('A', None), 'rxn_2'), -1),
                         ((('B', None), 'rxn_2'), 1),
                         ((('A', None), 'rxn_3'), -1),
                         ((('D', None), 'rxn_3'), 1),
                         ((('A', None), 'rxn_4'), -1),
                         ((('C', None), 'rxn_4'), 1),
                         ((('C', None), 'rxn_5'), -1),
                         ((('D', None), 'rxn_5'), 1),
                         ((('D', None), 'rxn_6'), -1) }
        self.assertEqual(set(self.model.matrix.iteritems()), matrix_set)

    def test_limits_get_item(self):
        self.assertEqual(self.model.limits['rxn_1'], metabolicmodel.FluxBounds(0, 1000))
        self.assertEqual(self.model.limits['rxn_2'], metabolicmodel.FluxBounds(-1000, 1000))
        self.assertEqual(self.model.limits['rxn_3'], metabolicmodel.FluxBounds(0, 1000))

class TestMetabolicModel(unittest.TestCase):
    def setUp(self):
        # TODO use mock database instead of actual database
        self.database = metabolicmodel.MetabolicDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))

        model = self.database.load_model_from_file(iter(self.database.reactions))
        self.model = metabolicmodel.FlipableModelView(model)

    def test_flipable_model_view_matrix_get_item_after_flip(self):
        self.model.flip({ 'rxn_4' })
        self.assertEqual(self.model.matrix[('A', None), 'rxn_1'], 2)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[('C', None), 'rxn_4'], -1)

    def test_flipable_model_view_matrix_get_item_after_double_flip(self):
        self.model.flip({ 'rxn_4', 'rxn_5' })
        self.model.flip({ 'rxn_1', 'rxn_4', 'rxn_2' })
        self.assertEqual(self.model.matrix[('A', None), 'rxn_1'], -2)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_2'], 1)
        self.assertEqual(self.model.matrix[('B', None), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[('A', None), 'rxn_4'], -1)
        self.assertEqual(self.model.matrix[('C', None), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[('C', None), 'rxn_5'], 1)
        self.assertEqual(self.model.matrix[('D', None), 'rxn_5'], -1)

    def test_flipable_model_view_limits_get_item_after_flip(self):
        self.model.flip({ 'rxn_1', 'rxn_2' })
        self.assertEqual(self.model.limits['rxn_1'], metabolicmodel.FluxBounds(-1000, 0))
        self.assertEqual(self.model.limits['rxn_2'], metabolicmodel.FluxBounds(-1000, 1000))
        self.assertEqual(self.model.limits['rxn_3'], metabolicmodel.FluxBounds(0, 1000))


if __name__ == '__main__':
    unittest.main()
