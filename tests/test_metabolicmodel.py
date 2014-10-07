#!/usr/bin/env python

import unittest

from metnet import metabolicmodel
from metnet.reaction import ModelSEED, Compound

class TestMetabolicDatabase(unittest.TestCase):
    def setUp(self):
        self.database = metabolicmodel.DictDatabase()
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
        self.database = metabolicmodel.DictDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))
        self.model = self.database.get_model(self.database.reactions)

    def test_database_property(self):
        self.assertIs(self.model.database, self.database)

    def test_reaction_set(self):
        self.assertEqual(set(self.model.reactions), { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_compound_set(self):
        self.assertEqual(set(self.model.compounds), { Compound('A'), Compound('B'), Compound('C'), Compound('D') })

    def test_add_reaction_new(self):
        self.database.set_reaction('rxn_7', ModelSEED.parse('|D| => |E|'))
        self.model.add_reaction('rxn_7')
        self.assertIn('rxn_7', set(self.model.reactions))
        self.assertIn(Compound('E'), set(self.model.compounds))

    def test_add_reaction_existing(self):
        self.model.add_reaction('rxn_1')
        self.assertEqual(set(self.model.reactions), { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })
        self.assertEqual(set(self.model.compounds), { Compound('A'), Compound('B'), Compound('C'), Compound('D') })

    def test_add_reaction_invalid(self):
        with self.assertRaises(Exception):
            self.model.add_reaction('rxn_7')

    def test_remove_reaction_existing(self):
        self.model.remove_reaction('rxn_2')
        self.assertEqual(set(self.model.reactions), { 'rxn_1', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })
        self.assertEqual(set(self.model.compounds), { Compound('A'), Compound('C'), Compound('D') })

    def test_add_all_database_reactions(self):
        self.database.set_reaction('rxn_7', ModelSEED.parse('|D| => |E|'))
        added = self.model.add_all_database_reactions()
        self.assertEqual(added, { 'rxn_7' })
        self.assertEqual(set(self.model.reactions), { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6', 'rxn_7' })

    def test_add_all_database_reactions_none(self):
        added = self.model.add_all_database_reactions()
        self.assertEqual(added, set())
        self.assertEqual(set(self.model.reactions), { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_matrix_get_item(self):
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_1'], 2)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[Compound('B'), 'rxn_2'], 1)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_4'], -1)
        self.assertEqual(self.model.matrix[Compound('C'), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[Compound('C'), 'rxn_5'], -1)
        self.assertEqual(self.model.matrix[Compound('D'), 'rxn_5'], 1)

    def test_matrix_get_item_invalid_key(self):
        with self.assertRaises(KeyError):
            a = self.model.matrix[Compound('A'), 'rxn_5']
        with self.assertRaises(KeyError):
            b = self.model.matrix['rxn_1']

    def test_matrix_set_item_is_invalid(self):
        with self.assertRaises(TypeError):
            self.model.matrix[Compound('A'), 'rxn_1'] = 4

    def test_matrix_iter(self):
        matrix_keys = { (Compound('A'), 'rxn_1'),
                        (Compound('A'), 'rxn_2'),
                        (Compound('B'), 'rxn_2'),
                        (Compound('A'), 'rxn_3'),
                        (Compound('D'), 'rxn_3'),
                        (Compound('A'), 'rxn_4'),
                        (Compound('C'), 'rxn_4'),
                        (Compound('C'), 'rxn_5'),
                        (Compound('D'), 'rxn_5'),
                        (Compound('D'), 'rxn_6') }
        self.assertEqual(set(iter(self.model.matrix)), matrix_keys)

    def test_matrix_len(self):
        self.assertEqual(len(self.model.matrix), 10)

    def test_limits_get_item(self):
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 1000))
        self.assertEqual(self.model.limits['rxn_2'].bounds, (-1000, 1000))
        self.assertEqual(self.model.limits['rxn_3'].bounds, (0, 1000))

    def test_limits_get_item_invalid_key(self):
        with self.assertRaises(KeyError):
            a = self.model.limits['rxn_7']

    def test_limits_set_item_is_invalid(self):
        with self.assertRaises(TypeError):
            self.model.limits['rxn_1'] = None

    def test_limits_set_upper_flux_bounds(self):
        self.model.limits['rxn_1'].upper = 500
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 500))

    def test_limits_set_upper_flux_bounds_to_invalid(self):
        with self.assertRaises(ValueError):
            self.model.limits['rxn_1'].upper = -20

    def test_limits_delete_upper_flux_bounds(self):
        self.model.limits['rxn_1'].upper = 500
        del self.model.limits['rxn_1'].upper
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 1000))

    def test_limits_delete_upper_not_set(self):
        del self.model.limits['rxn_1'].upper
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 1000))

    def test_limits_set_lower_flux_bounds(self):
        self.model.limits['rxn_1'].lower = 500
        self.assertEqual(self.model.limits['rxn_1'].bounds, (500, 1000))

    def test_limits_set_lower_flux_bounds_to_invalid(self):
        with self.assertRaises(ValueError):
            self.model.limits['rxn_1'].lower = 1001

    def test_limits_delete_lower_flux_bounds(self):
        self.model.limits['rxn_1'].lower = 500
        del self.model.limits['rxn_1'].lower
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 1000))

    def test_limits_delete_lower_not_set(self):
        del self.model.limits['rxn_1'].lower
        self.assertEqual(self.model.limits['rxn_1'].bounds, (0, 1000))

    def test_limits_set_both_flux_bounds(self):
        self.model.limits['rxn_2'].bounds = 1001, 1002
        self.assertEqual(self.model.limits['rxn_2'].lower, 1001)
        self.assertEqual(self.model.limits['rxn_2'].upper, 1002)
        self.assertEqual(self.model.limits['rxn_2'].bounds, (1001, 1002))

    def test_limits_set_both_flux_bounds_to_invalid(self):
        with self.assertRaises(ValueError):
            self.model.limits['rxn_2'].bounds = 10, -1000

    def test_limits_delete_both_flux_bounds(self):
        self.model.limits['rxn_2'].bounds = -10, 800
        del self.model.limits['rxn_2'].bounds
        self.assertEqual(self.model.limits['rxn_2'].bounds, (-1000, 1000))

    def test_limits_set_both_flux_bounds_delete_lower(self):
        self.model.limits['rxn_2'].bounds = -10, 800
        del self.model.limits['rxn_2'].lower
        self.assertEqual(self.model.limits['rxn_2'].bounds, (-1000, 800))

    def test_limits_set_both_flux_bounds_delete_upper(self):
        self.model.limits['rxn_2'].bounds = -10, 800
        del self.model.limits['rxn_2'].upper
        self.assertEqual(self.model.limits['rxn_2'].bounds, (-10, 1000))

    def test_limits_iter(self):
        self.assertEqual(set(iter(self.model.limits)), { 'rxn_1', 'rxn_2', 'rxn_3', 'rxn_4', 'rxn_5', 'rxn_6' })

    def test_limits_len(self):
        self.assertEqual(len(self.model.limits), 6)

class TestMetabolicModelFlipableView(unittest.TestCase):
    def setUp(self):
        # TODO use mock database instead of actual database
        self.database = metabolicmodel.DictDatabase()
        self.database.set_reaction('rxn_1', ModelSEED.parse('=> (2) |A|'))
        self.database.set_reaction('rxn_2', ModelSEED.parse('|A| <=> |B|'))
        self.database.set_reaction('rxn_3', ModelSEED.parse('|A| => |D|'))
        self.database.set_reaction('rxn_4', ModelSEED.parse('|A| => |C|'))
        self.database.set_reaction('rxn_5', ModelSEED.parse('|C| => |D|'))
        self.database.set_reaction('rxn_6', ModelSEED.parse('|D| =>'))

        model = self.database.get_model(self.database.reactions)
        self.model = metabolicmodel.FlipableModelView(model)

    def test_flipable_model_view_matrix_get_item_after_flip(self):
        self.model.flip({ 'rxn_4' })
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_1'], 2)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[Compound('C'), 'rxn_4'], -1)

    def test_flipable_model_view_matrix_get_item_after_double_flip(self):
        self.model.flip({ 'rxn_4', 'rxn_5' })
        self.model.flip({ 'rxn_1', 'rxn_4', 'rxn_2' })
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_1'], -2)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_2'], 1)
        self.assertEqual(self.model.matrix[Compound('B'), 'rxn_2'], -1)
        self.assertEqual(self.model.matrix[Compound('A'), 'rxn_4'], -1)
        self.assertEqual(self.model.matrix[Compound('C'), 'rxn_4'], 1)
        self.assertEqual(self.model.matrix[Compound('C'), 'rxn_5'], 1)
        self.assertEqual(self.model.matrix[Compound('D'), 'rxn_5'], -1)

    def test_flipable_model_view_limits_get_item_after_flip(self):
        self.model.flip({ 'rxn_1', 'rxn_2' })
        self.assertEqual(self.model.limits['rxn_1'].bounds, (-1000, 0))
        self.assertEqual(self.model.limits['rxn_2'].bounds, (-1000, 1000))
        self.assertEqual(self.model.limits['rxn_3'].bounds, (0, 1000))

    def test_flipable_model_view_limits_set_item_after_flip(self):
        self.model.flip({ 'rxn_1' })
        self.model.limits['rxn_1'].bounds = -20, 500
        self.assertEqual(self.model.limits['rxn_1'].bounds, (-20, 500))

        self.model.flip({ 'rxn_1' })
        self.assertEqual(self.model.limits['rxn_1'].bounds, (-500, 20))

if __name__ == '__main__':
    unittest.main()
