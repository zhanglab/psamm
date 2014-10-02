#!/usr/bin/env python

import unittest

from metnet.reaction import Reaction, Compound
from metnet.reaction import ModelSEED
from metnet.reaction import KEGG
from metnet.reaction import SudenSimple
from metnet.reaction import MetNet

from metnet.expression.affine import Expression

from decimal import Decimal

class TestCompound(unittest.TestCase):
    def test_compound_init_no_arguments(self):
        with self.assertRaises(TypeError):
            c = Compound()

    def test_compound_init_with_name(self):
        c = Compound('Phosphate')

    def test_compound_init_with_compartment(self):
        c = Compound('Phosphate', compartment='e')

    def test_compound_init_with_argument(self):
        c = Compound('Menaquinone', arguments=[8])

    def test_compound_init_with_compartment_and_argument(self):
        c = Compound('Menaquinone', compartment='e', arguments=[8])

    def test_compound_name(self):
        c = Compound('Phosphate')
        self.assertEquals(c.name, 'Phosphate')

    def test_compound_compartment(self):
        c = Compound('Phosphate', compartment='e')
        self.assertEquals(c.compartment, 'e')

    def test_compound_compartment_none(self):
        c = Compound('Phosphate')
        self.assertIsNone(c.compartment)

    def test_compound_arguments(self):
        c = Compound('Menaquinone', arguments=[8])
        self.assertEquals(list(c.arguments), [8])

    def test_compound_arguments_none(self):
        c = Compound('Phosphate')
        self.assertEquals(list(c.arguments), [])

    def test_compound_id_basic(self):
        self.assertEquals(Compound('C1234').id, 'C1234')

    def test_compound_id_from_invalid_characters(self):
        self.assertEquals(Compound('D-Glucose 1-phosphate').id, 'D_Glucose_1_phosphate')

    def test_compound_id_from_digit_prefix(self):
        self.assertEquals(Compound('2-Oxoglutarate').id, '_2_Oxoglutarate')

    def test_compound_id_with_compartment(self):
        self.assertEquals(Compound('L-Glutamine', 'e').id, 'L_Glutamine_e')

    def test_compound_translate_name(self):
        c = Compound('Pb')
        self.assertEquals(c.translate(lambda x: x.lower()), Compound('pb'))
        self.assertIsNot(c.translate(lambda x: x.lower()), c)

    def test_compound_in_compartment_when_unassigned(self):
        c = Compound('H+')
        self.assertEquals(c.in_compartment('e'), Compound('H+', 'e'))
        self.assertIsNot(c.in_compartment('e'), c)

    def test_compound_in_compartment_when_assigned(self):
        c = Compound('H+', 'c')
        self.assertEquals(c.in_compartment('p'), Compound('H+', 'p'))
        self.assertIsNot(c.in_compartment('p'), c)

    def test_compound_equals_other_with_same_name(self):
        c = Compound('Phosphate')
        self.assertEquals(c, Compound('Phosphate'))

    def test_compound_not_equals_other_with_name(self):
        c = Compound('H+')
        self.assertNotEquals(c, Compound('Phosphate'))

    def test_compound_equals_other_with_compartment(self):
        c = Compound('Phosphate', 'e')
        self.assertEquals(c, Compound('Phosphate', 'e'))

    def test_compound_not_equals_other_with_compartment(self):
        c = Compound('Phosphate', 'e')
        self.assertNotEquals(c, Compound('Phosphate', None))

    def test_compound_equals_other_with_arguments(self):
        c = Compound('Polyphosphate', arguments=[5])
        self.assertEquals(c, Compound('Polyphosphate', arguments=[5]))

    def test_compound_not_equals_other_with_arguments(self):
        c = Compound('Polyphosphate', arguments=[4])
        self.assertNotEquals(c, Compound('Polyphosphate', arguments=[5]))

    def test_compounds_sorted(self):
        l = sorted([Compound('A', arguments=[10]), Compound('A'), Compound('B'),
                    Compound('A', arguments=[4]), Compound('A', 'e')])
        self.assertEquals(l, [Compound('A'), Compound('A', arguments=[4]),
                              Compound('A', arguments=[10]), Compound('A', 'e'), Compound('B')])

    def test_compound_str_basic(self):
        self.assertEquals(str(Compound('Phosphate')), 'Phosphate')

    def test_compound_str_with_compartment(self):
        self.assertEquals(str(Compound('Phosphate', 'e')), 'Phosphate[e]')

    def test_compound_str_with_argument(self):
        self.assertEquals(str(Compound('Polyphosphate', arguments=[3])), 'Polyphosphate(3)')

    def test_compound_str_with_compartment_and_argument(self):
        self.assertEquals(str(Compound('Polyphosphate', 'p', [3])), 'Polyphosphate(3)[p]')


class TestReaction(unittest.TestCase):
    def test_reaction_init_empty_bidir(self):
        r = Reaction(Reaction.Bidir, [], [])

    def test_reaction_init_empty_left(self):
        r = Reaction(Reaction.Left, [], [])

    def test_reaction_init_empty_right(self):
        r = Reaction(Reaction.Right, [], [])

    def test_reaction_init_empty_invalid_direction(self):
        with self.assertRaises(ValueError):
            r = Reaction(None, [], [])

    def test_reaction_init_left_empty(self):
        r = Reaction(Reaction.Bidir, [], [(Compound('A'), 1)])

    def test_reaction_init_right_empty(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [])

    def test_reaction_init_none_empty(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [(Compound('B'), 1)])

    def test_reaction_direction_property(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [(Compound('B'), 1)])
        self.assertEquals(r.direction, Reaction.Bidir)

    def test_reaction_left_property(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [(Compound('B'), 1)])
        self.assertEquals(list(r.left), [(Compound('A'), 1)])

    def test_reaction_right_property(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [(Compound('B'), 1)])
        self.assertEquals(list(r.right), [(Compound('B'), 1)])

    def test_reaction_compounds_property(self):
        r = Reaction(Reaction.Bidir, [(Compound('A'), 1)], [(Compound('B'), 1)])
        self.assertEquals(list(r.compounds), [(Compound('A'), 1), (Compound('B'), 1)])

    def test_reaction_normalized_of_bidir(self):
        r = Reaction(Reaction.Bidir, [(Compound('Au'), 1)], [(Compound('Pb'), 1)]).normalized()
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('Au'), 1)], [(Compound('Pb'), 1)]))

    def test_reaction_normalized_of_right(self):
        r = Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)]).normalized()
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)]))

    def test_reaction_normalized_of_left(self):
        r = Reaction(Reaction.Left, [(Compound('Au'), 1)], [(Compound('Pb'), 1)]).normalized()
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

    def test_reaction_translated_compounds(self):
        r = Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)])
        rt = r.translated_compounds(lambda name: name.lower())
        self.assertEquals(rt, Reaction(Reaction.Right, [(Compound('pb'), 1)], [(Compound('au'), 1)]))

    def test_reaction_format_simple(self):
        r = Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)])
        self.assertEquals(r.format(), '|Pb| => |Au|')

    def test_reaction_format_with_arguments(self):
        pp1 = Compound('Polyphosphate', arguments=[4])
        pp2 = Compound('Polyphosphate', arguments=[5])
        r = Reaction(Reaction.Right, [(Compound('ATP'), 1), (pp1, 1)], [(Compound('ADP'), 1), (pp2, 1)])
        self.assertEquals(r.format(), '|ATP| + |Polyphosphate(4)| => |ADP| + |Polyphosphate(5)|')

    def test_reaction_equals_other(self):
        r = Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)])
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

    def test_reaction_not_equals_other_with_different_compounds(self):
        r = Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)])
        self.assertNotEquals(r, Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

    def test_reaction_not_equals_other_with_different_direction(self):
        r = Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)])
        self.assertNotEquals(r, Reaction(Reaction.Left, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

class TestModelSEED(unittest.TestCase):
    def test_modelseed_parse(self):
        r = ModelSEED.parse('|H2O| + |PPi| => (2) |Phosphate| + (2) |H+|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2O'), 1), (Compound('PPi'), 1)],
                                      [(Compound('Phosphate'), 2), (Compound('H+'), 2)]))

    def test_modelseed_parse_with_decimal(self):
        r = ModelSEED.parse('|H2| + (0.5) |O2| => |H2O|')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('H2'), 1), (Compound('O2'), Decimal('0.5'))],
                                      [(Compound('H2O'), 1)]))

    def test_modelseed_format(self):
        r = Reaction(Reaction.Left, [(Compound('H2O'), 2)], [(Compound('H2'), 2), (Compound('O2'), 1)])
        self.assertEquals(r.format(), '(2) |H2O| <= (2) |H2| + |O2|')

class TestSudenSimple(unittest.TestCase):
    def test_sudensimple_parse(self):
        r = SudenSimple.parse('1 H2O + 1 PPi <=> 2 Phosphate + 2 proton')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('H2O'), 1), (Compound('PPi'), 1)],
                                      [(Compound('Phosphate'), 2), (Compound('proton'), 2)]))

    def test_sudensimple_parse_with_decimal(self):
        r = SudenSimple.parse('1 H2 + 0.5 O2 <=> 1 H2O')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('H2'), 1), (Compound('O2'), Decimal('0.5'))],
                                      [(Compound('H2O'), 1)]))

class TestMetNet(unittest.TestCase):
    def test_metnet_parse_with_global_compartment(self):
        r = MetNet.parse('[c] : akg + ala-L <==> glu-L + pyr')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('cpd_akg', 'c'), 1), (Compound('cpd_ala-L', 'c'), 1)],
                                      [(Compound('cpd_glu-L', 'c'), 1), (Compound('cpd_pyr', 'c'), 1)]))

    def test_metnet_parse_with_local_compartment(self):
        r =  MetNet.parse('(2) ficytcc553[c] + so3[c] + h2o[c] --> (2) focytcc553[c] + so4[c] + (2) h[e]')
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('cpd_ficytcc553', 'c'), 2), (Compound('cpd_so3', 'c'), 1),
                                                       (Compound('cpd_h2o', 'c'), 1)],
                                      [(Compound('cpd_focytcc553', 'c'), 2), (Compound('cpd_so4', 'c'), 1),
                                       (Compound('cpd_h', 'e'), 2)]))

class TestKEGG(unittest.TestCase):
    def test_kegg_parse(self):
        r = KEGG.parse('C00013 + C00001 <=> 2 C00009')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00013'), 1), (Compound('C00001'), 1)],
                                      [(Compound('C00009'), 2)]))

    def test_kegg_parse_with_count_expression(self):
        r = KEGG.parse('2n C00404 + n C00001 <=> (n+1) C02174')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00404'), Expression('2n')),
                                                       (Compound('C00001'), Expression('n'))],
                                      [(Compound('C02174'), Expression('n+1'))]))

    def test_kegg_parse_with_compound_argument(self):
        r = KEGG.parse('C00039(n) <=> C00013 + C00039(n+1)')
        self.assertEquals(r, Reaction(Reaction.Bidir, [(Compound('C00039', arguments=[Expression('n')]), 1)],
                                      [(Compound('C00013'), 1), (Compound('C00039', arguments=[Expression('n+1')]), 1)]))

if __name__ == '__main__':
    unittest.main()
