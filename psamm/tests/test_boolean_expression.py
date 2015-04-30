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

from psamm.expression.boolean import Expression, Variable, And, Or


class TestVariable(unittest.TestCase):
    def test_variable_init(self):
        v = Variable('x')

    def test_variable_init_long_symbol(self):
        v = Variable('xyz')

    def test_variable_init_number_symbol(self):
        v = Variable('x2')

    def test_variable_init_underscore_symbol(self):
        v = Variable('x_y')

    def test_variable_symbol(self):
        self.assertEqual(Variable('x').symbol, 'x')

    def test_variable_substitute(self):
        subst = {'x': True}
        self.assertEqual(Variable('x').substitute(lambda v: subst.get(v.symbol, v)), True)

    def test_variable_substitute_self(self):
        self.assertEqual(Variable('x').substitute(lambda v: v), Variable('x'))

    def test_variable_substitute_multiple(self):
        subst = {'y': True, 'x': False}
        self.assertEqual(Variable('x').substitute(lambda v: subst.get(v.symbol, v)), False)

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
        self.assertTrue(e.substitute(lambda v: {'b1': True}.get(v.symbol, v)))

    def test_expression_substitute_existing_false(self):
        e = Expression('b1')
        self.assertFalse(e.substitute(lambda v: {'b1': False}.get(v.symbol, v)))

    def test_expression_substitute_unknown_to_variable(self):
        e = Expression('b1')
        self.assertEqual(e.substitute(lambda v: v), Variable('b1'))

    def test_expression_substitute_unknown_to_expression(self):
        e = Expression('b1 and b2')
        self.assertEqual(e.substitute(lambda v: v), Expression('b1 and b2'))

    def test_expression_parse_and(self):
        e = Expression('b1 and b2')
        self.assertEqual(e, And(Variable('b1'), Variable('b2')))

    def test_expression_parse_or(self):
        e = Expression('b1 or b2')
        self.assertEqual(e, Or(Variable('b1'), Variable('b2')))

    def test_expression_parse_multiple_with_duplicates(self):
        e = Expression('b1 and b2 and b1 and b4')
        self.assertEqual(e, And(Variable('b1'), Variable('b2'), Variable('b4')))

    def test_expression_parse_multiple_and(self):
        e = Expression('b1 and b2 and b3 and b4')
        self.assertEqual(e, And(Variable('b1'), Variable('b2'), Variable('b3'), Variable('b4')))

    def test_expression_parse_multiple_or(self):
        e = Expression('b1 or b2 or b3 or b4')
        self.assertEqual(e, Or(Variable('b1'), Variable('b2'), Variable('b3'), Variable('b4')))

    def test_expression_parse_multiple_parenthesis_and(self):
        e = Expression('b1 and (b2 and b3) and b4')
        self.assertEqual(e, And(Variable('b1'), Variable('b2'), Variable('b3'), Variable('b4')))

    def test_expression_parse_multiple_parenthesis_or(self):
        e = Expression('b1 or (b2 or b3) or b4')
        self.assertEqual(e, Or(Variable('b1'), Variable('b2'), Variable('b3'), Variable('b4')))

    def test_expression_parse_parentheses_mixed_1(self):
        e = Expression('(b1 and b2) or (b3 and b4)')
        self.assertEqual(e, Or(And(Variable('b1'), Variable('b2')), And(Variable('b3'), Variable('b4'))))

    def test_expression_parse_parentheses_mixed_2(self):
        e = Expression('(b1 or b2) and (b3 or b4)')
        self.assertEqual(e, And(Or(Variable('b1'), Variable('b2')), Or(Variable('b3'), Variable('b4'))))

    def test_expression_parse_parenthesis_with_space_right(self):
        e = Expression('(b1 and b2 )')
        self.assertEqual(e, And(Variable('b1'), Variable('b2')))

    def test_expression_parse_parenthesis_with_space_left(self):
        e = Expression('( b1 and b2)')
        self.assertEqual(e, And(Variable('b1'), Variable('b2')))

    def test_expression_parse_parentheses_right_nested(self):
        e = Expression('(b1 or (b2 and (b3 or (b4 and (b5)))))')
        self.assertEqual(e, Or(Variable('b1'), And(Variable('b2'), Or(Variable('b3'), And(Variable('b4'), Variable('b5'))))))

    def test_expression_parse_parentheses_left_nested(self):
        e = Expression('(((b1 or b2) and b3) or b4) and (b5)')
        self.assertEqual(e, And(Variable('b5'), Or(Variable('b4'), And(Variable('b3'), Or(Variable('b1'), Variable('b2'))))))

    def test_expression_parse_implicit_mixed_1(self):
        e = Expression('b1 and b2 or b3 and b4')
        self.assertEqual(e, Or(And(Variable('b1'), Variable('b2')), And(Variable('b3'), Variable('b4'))))

    def test_expression_parse_implicit_mixed_2(self):
        e = Expression('b1 or b2 and b3 or b4')
        self.assertEqual(e, Or(Variable('b1'), And(Variable('b2'), Variable('b3')), Variable('b4')))

    def test_expression_parse_longer_names(self):
        e = Expression('b12345 and bxyz and testtesttest')
        self.assertEqual(e, And(Variable('b12345'), Variable('bxyz'), Variable('testtesttest')))

    def test_expression_parse_name_starting_with_or(self):
        e = Expression('b1 and order')
        self.assertEqual(e, And(Variable('b1'), Variable('order')))

    def test_expression_parse_name_starting_with_and(self):
        e = Expression('b1 or anderson')
        self.assertEqual(e, Or(Variable('b1'), Variable('anderson')))

    def test_expression_parse_with_extra_space(self):
        e = Expression('b1    and   b2       and b3')
        self.assertEqual(e, And(Variable('b1'), Variable('b2'), Variable('b3')))


if __name__ == '__main__':
    unittest.main()
