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

from psamm.reaction import Reaction, Compound


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
        self.assertEquals(list(r.compounds),
                          [(Compound('A'), -1), (Compound('B'), 1)])

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
        self.assertEquals(str(r), '|Pb| => |Au|')

    def test_reaction_format_with_arguments(self):
        pp1 = Compound('Polyphosphate', arguments=[4])
        pp2 = Compound('Polyphosphate', arguments=[5])
        r = Reaction(Reaction.Right, [(Compound('ATP'), 1), (pp1, 1)], [(Compound('ADP'), 1), (pp2, 1)])
        self.assertEquals(str(r), '|ATP| + |Polyphosphate(4)| => |ADP| + |Polyphosphate(5)|')

    def test_reaction_equals_other(self):
        r = Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)])
        self.assertEquals(r, Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

    def test_reaction_not_equals_other_with_different_compounds(self):
        r = Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)])
        self.assertNotEquals(r, Reaction(Reaction.Right, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))

    def test_reaction_not_equals_other_with_different_direction(self):
        r = Reaction(Reaction.Right, [(Compound('Au'), 1)], [(Compound('Pb'), 1)])
        self.assertNotEquals(r, Reaction(Reaction.Left, [(Compound('Pb'), 1)], [(Compound('Au'), 1)]))


if __name__ == '__main__':
    unittest.main()
