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

from psamm.formula import Formula, FormulaElement, Atom, Radical


class TestFormulaElement(unittest.TestCase):
    def test_add_formula_elements(self):
        e1 = FormulaElement()
        e2 = FormulaElement()
        self.assertEqual(e1 + e2, Formula({e1: 1, e2: 1}))

    def test_add_formula_element_to_self(self):
        e1 = FormulaElement()
        self.assertEqual(e1 + e1, Formula({e1: 2}))

    def test_add_formula_element_and_number(self):
        with self.assertRaises(TypeError):
            f = FormulaElement() + 42

    def test_merge_formula_elements(self):
        e1 = FormulaElement()
        e2 = FormulaElement()
        self.assertEqual(e1 | e2, Formula({e1: 1, e2: 1}))

    def test_merge_formula_element_to_self(self):
        e1 = FormulaElement()
        self.assertEqual(e1 | e1, Formula({e1: 2}))

    def test_substitute_into_formula_element(self):
        e1 = FormulaElement()
        self.assertEqual(
            e1.substitute(lambda v: {'x': 42}.get(v.symbol, v)), e1)


class TestAtom(unittest.TestCase):
    def test_atom_symbol(self):
        a = Atom('H')
        self.assertEqual(a.symbol, 'H')

    def test_atom_symbol_wide(self):
        a = Atom('Zn')
        self.assertEqual(a.symbol, 'Zn')

    def test_atom_symbol_non_standard(self):
        a = Atom('X')
        self.assertEqual(a.symbol, 'X')

    def test_atom_singleton_property(self):
        a = Atom.C
        self.assertEqual(a.symbol, 'C')

    def test_atom_to_string(self):
        a = Atom('C')
        self.assertEqual(str(a), 'C')

    def test_atom_repr(self):
        a = Atom('Si')
        self.assertEqual(repr(a), "Atom('Si')")

    def test_atom_equals(self):
        a1 = Atom('H')
        self.assertEqual(a1, Atom.H)
        self.assertNotEqual(a1, Atom.Zn)

    def test_atom_ordered(self):
        a1 = Atom('H')
        a2 = Atom('C')
        a3 = Atom('Zn')
        self.assertGreater(a1, a2)
        self.assertLess(a1, a3)


class TestFormula(unittest.TestCase):
    def test_formula_merge_same_formulas_with_same_atoms(self):
        f = Formula({Atom.H: 2, Atom.O: 1}) | Formula({Atom.N: 1, Atom.O: 2})
        self.assertEqual(f, Formula({Atom.H: 2, Atom.N: 1, Atom.O: 3}))

    def test_formula_merge_formulas_that_cancel_out(self):
        f = Formula({Atom.H: 3}) | Formula({Atom.H: -3})
        self.assertEqual(f, Formula())

    def test_formula_multiply_number(self):
        f = Formula({Atom.H: 2, Atom.O: 1}) * 4
        self.assertEqual(f, Formula({Atom.H: 8, Atom.O: 4}))

    def test_formula_multiply_one(self):
        f = Formula({Atom.H: 2, Atom.O: 1}) * 1
        self.assertEqual(f, f)

    def test_formula_multiply_zero(self):
        f = Formula({Atom.H: 2, Atom.O: 1}) * 0
        self.assertEqual(f, Formula())

    def test_formula_right_multiply_number(self):
        f = 2 * Formula({Atom.H: 2, Atom.O: 1})
        self.assertEqual(f, Formula({Atom.H: 4, Atom.O: 2}))

    def test_formula_repeat(self):
        f = Formula({Atom.H: 2, Atom.O: 1}).repeat(4)
        self.assertEqual(f, Formula({Formula({Atom.H: 2, Atom.O: 1}): 4}))

    def test_formula_equals_other_formula(self):
        f = Formula({Atom.H: 2, Atom.O: 1})
        self.assertEqual(f, Formula({Atom.O: 1, Atom.H: 2}))

    def test_formula_not_equals_other_with_distinct_elements(self):
        f = Formula({Atom.Au: 1})
        self.assertNotEqual(f, Formula({Atom.Ag: 1}))

    def test_formula_not_equals_other_with_different_number(self):
        f = Formula({Atom.Au: 1})
        self.assertNotEqual(f, Formula({Atom.Au: 2}))

    def test_formula_items(self):
        f = Formula({Atom.H: 12, Atom.C: 6, Atom.O: 6})
        self.assertEqual(dict(f.items()), {
            Atom.C: 6,
            Atom.H: 12,
            Atom.O: 6
        })

    def test_formula_contains(self):
        f = Formula({Atom.H: 12, Atom.C: 6, Atom.O: 6})
        self.assertIn(Atom.C, f)
        self.assertNotIn(Atom.Ag, f)

    def test_formula_to_string(self):
        f = Formula({Atom.H: 12, Atom.C: 6, Atom.O: 6})
        self.assertEqual(str(f), 'C6H12O6')

    def test_formula_to_string_with_group(self):
        f = Formula({
            Atom.C: 1,
            Atom.H: 3,
            Formula({Atom.C: 1, Atom.H: 2}): 14,
            Formula({
                Atom.C: 1, Atom.O: 1,
                Formula({Atom.H: 1, Atom.O: 1}): 1
            }): 1
        })
        # Ideally: self.assertEqual(str(f), 'CH3(CH2)14COOH')
        # The two subgroups are unordered so we cannot assert a specfic string
        # at this point.

    def test_formula_flattened(self):
        f = Formula({
            Atom.C: 1,
            Atom.H: 3,
            Formula({Atom.C: 1, Atom.H: 2}): 14,
            Formula({
                Atom.C: 1, Atom.O: 1,
                Formula({Atom.H: 1, Atom.O: 1}): 1
            }): 1
        })
        self.assertEqual(f.flattened(), Formula({
            Atom.C: 16,
            Atom.H: 32,
            Atom.O: 2
        }))

    def test_formula_substitute_non_positive(self):
        f = Formula.parse('CH3(CH2)nCOOH')

        with self.assertRaises(ValueError):
            f.substitute(lambda v: {'n': -5}.get(v.symbol, v))

        with self.assertRaises(ValueError):
            f.substitute(lambda v: {'n': 0}.get(v.symbol, v))

    def test_formula_is_not_variable(self):
        f = Formula.parse('C6H12O6')
        self.assertFalse(f.is_variable())

    def test_formula_is_variable(self):
        f = Formula.parse('C2H4NO2R(C2H2NOR)n')
        self.assertTrue(f.is_variable())

    def test_formula_balance_missing_on_one_side(self):
        f1, f2 = Formula.balance(Formula.parse('H2O'), Formula.parse('OH'))
        self.assertEqual(f1, Formula())
        self.assertEqual(f2, Formula({Atom.H: 1}))

    def test_formula_balance_missing_on_both_sides(self):
        f1, f2 = Formula.balance(
            Formula.parse('C3H6OH'), Formula.parse('CH6O2'))
        self.assertEqual(f1, Formula({Atom.O: 1}))
        self.assertEqual(f2, Formula({Atom.C: 2, Atom.H: 1}))

    def test_formula_balance_subgroups_cancel_out(self):
        f1, f2 = Formula.balance(
            Formula.parse('H2(CH2)n'), Formula.parse('CH3O(CH2)n'))
        self.assertEqual(f1, Formula({Atom.C: 1, Atom.H: 1, Atom.O: 1}))
        self.assertEqual(f2, Formula())


class TestFormulaParser(unittest.TestCase):
    def test_formula_parse_with_final_digit(self):
        f = Formula.parse('H2O2')
        self.assertEqual(f, Formula({Atom.H: 2, Atom.O: 2}))

    def test_formula_parse_with_implicit_final_digit(self):
        f = Formula.parse('H2O')
        self.assertEqual(f, Formula({Atom.H: 2, Atom.O: 1}))

    def test_formula_parse_with_implicit_digit(self):
        f = Formula.parse('C2H5NO2')
        self.assertEqual(f, Formula(
            {Atom.C: 2, Atom.H: 5, Atom.N: 1, Atom.O: 2}))

    def test_formula_parse_with_wide_element(self):
        f = Formula.parse('ZnO')
        self.assertEqual(f, Formula({Atom.Zn: 1, Atom.O: 1}))

    def test_formula_parse_with_wide_count(self):
        f = Formula.parse('C6H10O2')
        self.assertEqual(f, Formula({Atom.C: 6, Atom.H: 10, Atom.O: 2}))

    def test_formula_parse_with_implicitly_counted_subgroup(self):
        f = Formula.parse('C2H6O2(CH)')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 6, Atom.O: 2,
            Formula({Atom.C: 1, Atom.H: 1}): 1}))

    def test_formula_parse_with_counted_subgroup(self):
        f = Formula.parse('C2H6O2(CH)2')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 6, Atom.O: 2,
            Formula({Atom.C: 1, Atom.H: 1}): 2}))

    def test_formula_parse_with_two_identical_counted_subgroups(self):
        f = Formula.parse('C2H6O2(CH)2(CH)2')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 6, Atom.O: 2,
            Formula({Atom.C: 1, Atom.H: 1}): 4}))

    def test_formula_parse_with_two_distinct_counted_subgroups(self):
        f = Formula.parse('C2H6O2(CH)2(CH2)2')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 6, Atom.O: 2,
            Formula({Atom.C: 1, Atom.H: 1}): 2,
            Formula({Atom.C: 1, Atom.H: 2}): 2}))

    def test_formula_parse_with_wide_counted_subgroup(self):
        f = Formula.parse('C2(CH)10NO2')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.N: 1, Atom.O: 2,
            Formula({Atom.C: 1, Atom.H: 1}): 10}))

    def test_formula_parse_with_radical(self):
        f = Formula.parse('C2H4NO2R')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 4, Atom.N: 1, Atom.O: 2, Radical('R'): 1}))

    def test_formula_parse_with_numbered_radical(self):
        f = Formula.parse('C2H4NO2(R1)')
        self.assertEqual(f, Formula({
            Atom.C: 2, Atom.H: 4, Atom.N: 1, Atom.O: 2, Radical('R1'): 1}))


if __name__ == '__main__':
    unittest.main()
