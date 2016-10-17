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

from psamm.expression.boolean import (Expression, And, Or, Variable,
                                      ParseError, SubstitutionError)


class TestVariable(unittest.TestCase):
    def test_variable_init(self):
        v = Variable('x')

    def test_variable_init_long_symbol(self):
        v = Variable('xyz')

    def test_variable_init_number_symbol(self):
        v = Variable('x2')

    def test_variable_init_underscore_symbol(self):
        v = Variable('x_y')

    def test_variable_init_from_unicode(self):
        symbol = u'\u00c6\u00d8\u00c5'
        v = Variable(symbol)
        self.assertEqual(v.symbol, symbol)

    def test_variable_init_with_dot_symbol(self):
        v = Variable('x12345.6')

    def test_variable_init_with_number(self):
        v = Variable('123')

    def test_variable_symbol(self):
        self.assertEqual(Variable('x').symbol, 'x')

    def test_variables_equals_true(self):
        self.assertEqual(Variable('x'), Variable('x'))

    def test_variables_not_equals(self):
        self.assertNotEqual(Variable('x'), Variable('y'))

    def test_variable_not_equals_bool(self):
        self.assertNotEqual(Variable('x'), True)

    def test_variable_hash(self):
        self.assertEqual(hash(Variable('xyz')), hash(Variable('xyz')))


class TestExpression(unittest.TestCase):
    def test_expression_substitute_existing(self):
        e = Expression('b1')
        e1 = e.substitute(lambda v: {'b1': True}.get(v.symbol, v))
        self.assertTrue(e1.has_value())
        self.assertTrue(e1.value)

    def test_expression_substitute_existing_false(self):
        e = Expression('b1')
        e1 = e.substitute(lambda v: {'b1': False}.get(v.symbol, v))
        self.assertTrue(e1.has_value())
        self.assertFalse(e1.value)

    def test_expression_substitute_unknown_to_variable(self):
        e = Expression('b1')
        e1 = e.substitute(lambda v: v)
        self.assertEqual(e1, Expression(Variable('b1')))

    def test_expression_substitute_unknown_to_expression(self):
        e = Expression('b1 and b2')
        e1 = e.substitute(lambda v: v)
        self.assertEqual(e1, Expression('b1 and b2'))

    def test_expression_substitute_invalid(self):
        e = Expression('b1 and b2')
        with self.assertRaises(SubstitutionError):
            e.substitute(lambda v: {'b1': 17}.get(v.symbol, v))

    def test_expression_substitute_with_short_circuit_and(self):
        e = Expression('a and (b or c) and d')
        e1 = e.substitute(lambda v: {'a': False}.get(v.symbol, v))
        self.assertTrue(e1.has_value())
        self.assertFalse(e1.value)

    def test_expression_substitute_with_short_circuit_or(self):
        e = Expression('(a and b) or c or d')
        e1 = e.substitute(lambda v: {'c': True}.get(v.symbol, v))
        self.assertTrue(e1.has_value())
        self.assertTrue(e1.value)

    def test_expression_substitute_remove_terms_from_and(self):
        e = Expression('a and b')
        e1 = e.substitute(lambda v: {'a': True}.get(v.symbol, v))
        self.assertEqual(e1, Expression(Variable('b')))

    def test_expression_substitute_remove_terms_from_or(self):
        e = Expression('a or b')
        e1 = e.substitute(lambda v: {'a': False}.get(v.symbol, v))
        self.assertEqual(e1, Expression(Variable('b')))

    def test_expression_substitute_evaluate_all_and_terms_to_true(self):
        e = Expression('a and b')
        e1 = e.substitute(lambda v: True)
        self.assertEqual(e1, Expression(True))

    def test_expression_substitute_evaluate_all_or_terms_to_false(self):
        e = Expression('a or b')
        e1 = e.substitute(lambda v: False)
        self.assertEqual(e1, Expression(False))

    def test_expression_iter_variables(self):
        e = Expression(
            Or(
                And(Variable('b1'), Variable('b2')),
                And(Variable('b3'), Variable('b4'))))
        self.assertEqual(list(e.variables), [
            Variable('b1'), Variable('b2'), Variable('b3'), Variable('b4')
        ])

    def test_expression_root(self):
        e = Expression(Or(Variable('b1'), Variable('b2')))
        root = e.root
        self.assertIsInstance(e.root, Or)

    def test_expression_root_boolean(self):
        e = Expression(False)
        self.assertEqual(e.root, False)

    def test_expression_to_string(self):
        e = Expression('(a and b) or (c and d)')
        self.assertEqual(str(e), '(a and b) or (c and d)')

    def test_expression_with_variable_to_string(self):
        e = Expression(Variable('a'))
        self.assertEqual(str(e), 'a')

    def test_expression_with_bool_to_string(self):
        e = Expression(False)
        self.assertEqual(str(e), 'False')

    def test_expression_parse_and(self):
        e = Expression('b1 and b2')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('b2'))))

    def test_expression_parse_or(self):
        e = Expression('b1 or b2')
        self.assertEqual(e, Expression(
            Or(Variable('b1'), Variable('b2'))))

    def test_expression_parse_multiple_with_duplicates(self):
        e = Expression('b1 and b2 and b1 and b4')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('b2'), Variable('b4'))))

    def test_expression_parse_multiple_and(self):
        e = Expression('b1 and b2 and b3 and b4')
        self.assertEqual(e, Expression(
            And(
                Variable('b1'), Variable('b2'), Variable('b3'),
                Variable('b4'))))

    def test_expression_parse_multiple_or(self):
        e = Expression('b1 or b2 or b3 or b4')
        self.assertEqual(e, Expression(
            Or(
                Variable('b1'), Variable('b2'), Variable('b3'),
                Variable('b4'))))

    def test_expression_parse_multiple_parenthesis_and(self):
        e = Expression('b1 and (b2 and b3) and b4')
        self.assertEqual(e, Expression(
            And(
                Variable('b1'), Variable('b2'), Variable('b3'),
                Variable('b4'))))

    def test_expression_parse_multiple_parenthesis_or(self):
        e = Expression('b1 or (b2 or b3) or b4')
        self.assertEqual(e, Expression(
            Or(
                Variable('b1'), Variable('b2'), Variable('b3'),
                Variable('b4'))))

    def test_expression_parse_parentheses_mixed_1(self):
        e = Expression('(b1 and b2) or (b3 and b4)')
        self.assertEqual(e, Expression(
            Or(
                And(Variable('b1'), Variable('b2')),
                And(Variable('b3'), Variable('b4')))))

    def test_expression_parse_parentheses_mixed_2(self):
        e = Expression('(b1 or b2) and (b3 or b4)')
        self.assertEqual(e, Expression(
            And(
                Or(Variable('b1'), Variable('b2')),
                Or(Variable('b3'), Variable('b4')))))

    def test_expression_parse_uppercase_operators(self):
        e = Expression('(b1 AND b2) OR b3')
        self.assertEqual(e, Expression(
            Or(
                Variable('b3'),
                And(Variable('b1'), Variable('b2')))))

    def test_expression_parse_parenthesis_with_space_right(self):
        e = Expression('(b1 and b2 )')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('b2'))))

    def test_expression_parse_parenthesis_with_space_left(self):
        e = Expression('( b1 and b2)')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('b2'))))

    def test_expression_parse_parentheses_right_nested(self):
        e = Expression('(b1 or (b2 and (b3 or (b4 and (b5)))))')
        self.assertEqual(e, Expression(
            Or(
                Variable('b1'),
                And(
                    Variable('b2'),
                    Or(
                        Variable('b3'),
                        And(Variable('b4'), Variable('b5')))))))

    def test_expression_parse_parentheses_left_nested(self):
        e = Expression('(((b1 or b2) and b3) or b4) and (b5)')
        self.assertEqual(e, Expression(
            And(
                Variable('b5'),
                Or(
                    Variable('b4'),
                    And(
                        Variable('b3'),
                        Or(Variable('b1'), Variable('b2')))))))

    def test_expression_parse_implicit_mixed_1(self):
        e = Expression('b1 and b2 or b3 and b4')
        self.assertEqual(e, Expression(
            Or(
                And(Variable('b1'), Variable('b2')),
                And(Variable('b3'), Variable('b4')))))

    def test_expression_parse_implicit_mixed_2(self):
        e = Expression('b1 or b2 and b3 or b4')
        self.assertEqual(e, Expression(
            Or(
                Variable('b1'),
                And(Variable('b2'), Variable('b3')),
                Variable('b4'))))

    def test_expression_parse_longer_names(self):
        e = Expression('b12345 and bxyz and testtesttest')
        self.assertEqual(e, Expression(
            And(
                Variable('b12345'), Variable('bxyz'),
                Variable('testtesttest'))))

    def test_expression_parse_with_missing_variable(self):
        with self.assertRaises(ParseError):
            e = Expression('b1 and and b3')

    def test_expression_parse_with_missing_end_parenthesis(self):
        with self.assertRaises(ParseError):
            e = Expression('b1 and (b2 or b3')

    def test_expression_parse_name_starting_with_or(self):
        e = Expression('b1 and order')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('order'))))

    def test_expression_parse_name_starting_with_and(self):
        e = Expression('b1 or anderson')
        self.assertEqual(e, Expression(
            Or(Variable('b1'), Variable('anderson'))))

    def test_expression_parse_with_extra_space(self):
        e = Expression('b1    and   b2       and b3')
        self.assertEqual(e, Expression(
            And(Variable('b1'), Variable('b2'), Variable('b3'))))

    def test_expression_parse_with_square_mixed_groups(self):
        e = Expression('[(a or b) or (c or d)] or [e or (f and g and h)]')
        self.assertEqual(e, Expression(
            Or(
                Variable('a'), Variable('b'), Variable('c'), Variable('d'),
                Variable('e'),
                And(Variable('f'), Variable('g'), Variable('h')))))

    def test_expression_parse_with_square_groups_unmatched(self):
        with self.assertRaises(ParseError):
            e = Expression('[(a or b) or (c or d])')


class TestParseError(unittest.TestCase):
    def test_parse_error_without_indicator(self):
        e = ParseError('Random error')
        self.assertIsNone(e.indicator)

    def test_parse_error_with_indicator(self):
        e = ParseError('Random error', span=(3, 5))
        self.assertEqual(e.indicator, '   ^^')

    def test_parse_error_with_zero_width_indicator(self):
        e = ParseError('Random error', span=(4, 4))
        self.assertEqual(e.indicator, '    ^')


if __name__ == '__main__':
    unittest.main()
