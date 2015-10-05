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

from psamm.expression.affine import Expression, Variable


class TestVariable(unittest.TestCase):
    def test_variable_init(self):
        v = Variable('x')

    def test_variable_init_long_symbol(self):
        v = Variable('xyz')

    def test_variable_init_number_symbol(self):
        v = Variable('x2')

    def test_variable_init_underscore_symbol(self):
        v = Variable('x_y')

    def test_variable_init_invalid_number(self):
        with self.assertRaises(ValueError):
            v = Variable('45x')

    def test_variable_symbol(self):
        self.assertEqual(Variable('x').symbol, 'x')

    def test_variable_simplify_returns_self(self):
        self.assertEqual(Variable('y').simplify(), Variable('y'))

    def test_variable_substitute(self):
        self.assertEqual(Variable('x').substitute(lambda v: {'x': 567}.get(v.symbol, v)), 567)

    def test_variable_substitute_unknown(self):
        self.assertEqual(Variable('x').substitute(lambda v: v), Variable('x'))

    def test_variable_substitute_multiple(self):
        self.assertEqual(Variable('x').substitute(lambda v: {'y': 123, 'x': 56}.get(v.symbol, v)), 56)

    def test_variable_add_number(self):
        self.assertEqual(Variable('x') + 1, Expression('x + 1'))

    def test_variable_radd_number(self):
        self.assertEqual(1 + Variable('x'), Expression('x + 1'))

    def test_variable_sub_number(self):
        self.assertEqual(Variable('x') - 4, Expression('x - 4'))

    def test_variable_rsub_number(self):
        self.assertEqual(4 - Variable('x'), Expression('4 - x'))

    def test_variable_mul_zero(self):
        self.assertEqual(Variable('x') * 0, 0)

    def test_variable_mul_one(self):
        self.assertEqual(Variable('x') * 1, Variable('x'))

    def test_variable_mul_number(self):
        self.assertEqual(Variable('x') * 2, Expression('2x'))

    def test_variable_rmul_number(self):
        self.assertEqual(3 * Variable('x'), Expression('3x'))

    def test_variable_div_number(self):
        self.assertEqual(Variable('x')/1, Variable('x'))

    def test_variable_neg(self):
        self.assertEqual(-Variable('x'), Expression('-x'))

    def test_variables_equals_true(self):
        self.assertEqual(Variable('x'), Variable('x'))

    def test_variables_not_equals(self):
        self.assertNotEqual(Variable('x'), Variable('y'))

    def test_variable_not_equals_number(self):
        self.assertNotEqual(Variable('x'), 5)

    def test_variables_ordered_by_symbol(self):
        var_list = sorted([Variable('b'), Variable('cd'), Variable('a'), Variable('cc')])
        self.assertEqual(var_list, [Variable('a'), Variable('b'), Variable('cc'), Variable('cd')])

    def test_variable_hash(self):
        self.assertEqual(hash(Variable('xyz')), hash(Variable('xyz')))

class TestExpression(unittest.TestCase):
    def test_expression_init_empty(self):
        e = Expression()

    def test_expression_init_with_variables(self):
        e = Expression({ Variable('a'): 1, Variable('b'): 2 })

    def test_expression_init_with_offset(self):
        e = Expression({}, 5)

    def test_expression_init_with_variables_and_offset(self):
        e = Expression({ Variable('a'): 4}, -5)

    def test_expression_init_with_zero_variables(self):
        e = Expression({ Variable('a'): 3, Variable('b'): 0})
        self.assertEqual(e, Expression({ Variable('a'): 3 }))

    def test_expression_init_with_non_variables(self):
        with self.assertRaises(ValueError):
            e = Expression({ 'a': 3 })

    def test_expression_simplify_to_number(self):
        result = Expression({}, 5).simplify()
        self.assertEqual(result, 5)
        self.assertIsInstance(result, int)

    def test_expression_simplify_to_variable(self):
        result = Expression({ Variable('x'): 1 }).simplify()
        self.assertEqual(result, Variable('x'))
        self.assertIsInstance(result, Variable)

    def test_expression_simplify_to_expression(self):
        result = Expression({ Variable('x'): 2 }).simplify()
        self.assertEqual(result, Expression('2x'))
        self.assertIsInstance(result, Expression)

    def test_expression_simplify_to_expression_with_offset(self):
        result = Expression({ Variable('x'): 1 }, 2).simplify()
        self.assertEqual(result, Expression('x + 2'))
        self.assertIsInstance(result, Expression)

    def test_expression_substitute_existing(self):
        e = Expression({ Variable('x'): 2 }, 1)
        self.assertEqual(e.substitute(lambda v: {'x': 2}.get(v.symbol, v)), 5)

    def test_expression_substitute_unknown_to_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e.substitute(lambda v: v), Variable('x'))

    def test_expression_substitute_unknown_to_expression(self):
        e = Expression({ Variable('x'): 2 }, 1)
        self.assertEqual(e.substitute(lambda v: v), Expression('2x + 1'))

    def test_expression_substitute_with_variable(self):
        e = Expression({ Variable('x'): 1, Variable('y'): 2 })
        self.assertEqual(e.substitute(lambda v: {'y': Variable('x')}.get(v.symbol, v)), Expression('3x'))

    def test_expression_substitute_with_expression(self):
        e = Expression({ Variable('x'): 3, Variable('y'): -2 })
        es = e.substitute(lambda v: {'x': Expression('y + 2z')}.get(v.symbol, v))
        self.assertEqual(es, Expression('y + 6z'))

    def test_expression_substitute_with_expression_is_atomic(self):
        e = Expression({ Variable('x'): 1, Variable('y'): 1 })
        es = e.substitute(lambda v: {'x': Expression('x + y'),
                                     'y': Expression('x + y') }.get(v.symbol, v))
        self.assertEqual(es, Expression('2x + 2y'))

    def test_expression_variables(self):
        e = Expression({ Variable('x'): 1, Variable('y'): 2 })
        self.assertEqual(sorted(e.variables()), [Variable('x'), Variable('y')])

    def test_expression_add_number(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e + 1, Expression('x + 1'))

    def test_expression_add_other_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e + Variable('y'), Expression('x + y'))

    def test_expression_add_same_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e + Variable('x'), Expression('2x'))

    def test_expression_add_expression(self):
        e1 = Expression({ Variable('x'): 1 })
        e2 = Expression({ Variable('y'): 2 })
        self.assertEqual(e1 + e2, Expression('x + 2y'))

    def test_expression_add_sums_to_zero(self):
        e1 = Expression({ Variable('x'): 2, Variable('y'): 1 })
        e2 = Expression({ Variable('x'): -2, Variable('z'): 3 })
        self.assertEqual(e1 + e2, Expression('y + 3z'))

    def test_expression_sub_number(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e - 1, Expression('x - 1'))

    def text_expression_rsub_number(self):
        e = Expression({ Variable('x'): 2 })
        self.assertEqual(4 - e, Expression('4 - 2x'))

    def test_expression_sub_other_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e - Variable('y'), Expression('x - y'))

    def test_expression_sub_same_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e - Variable('x'), 0)

    def test_expression_sub_expression(self):
        e1 = Expression({ Variable('x'): 1 })
        e2 = Expression({ Variable('y'): 2 })
        self.assertEqual(e1 - e2, Expression('x - 2y'))

    def test_expression_mul_zero(self):
        e = Expression({ Variable('x'): 3, Variable('y'): -2 }, 4)
        self.assertEqual(e * 0, 0)

    def test_expression_mul_one(self):
        e = Expression({ Variable('x'): 3, Variable('y'): -2 }, 4)
        self.assertEqual(e * 1, e)

    def test_expression_mul_number(self):
        e = Expression({ Variable('x'): 3, Variable('y'): -2 }, 4)
        self.assertEqual(e * 4, Expression('12x - 8y + 16'))

    def test_expression_div_number(self):
        e = Expression({ Variable('x'): 8, Variable('y'): -2}, 4)
        self.assertEqual(e/2, Expression('4x - y + 2'))

    def test_expression_neg(self):
        e = Expression({ Variable('x'): 3, Variable('y'): -2 }, 4)
        self.assertEqual(-e, Expression('-3x + 2y - 4'))

    def test_expression_equals_expression(self):
        e1 = Expression({ Variable('x'): 2, Variable('y'): 5 })
        e2 = Expression({ Variable('y'): 5, Variable('x'): 2 })
        self.assertEqual(e1, e2)

    def test_expression_not_equals_same_variables(self):
        e1 = Expression({ Variable('x'): 2, Variable('y'): 5 })
        e2 = Expression({ Variable('y'): 2, Variable('x'): -5 })
        self.assertNotEqual(e1, e2)

    def test_expression_not_equals_diffrent_variables(self):
        e1 = Expression({ Variable('x'): 2, Variable('y'): 5 })
        e2 = Expression({ Variable('y'): 2, Variable('x'): -5, Variable('z'): 1 })
        self.assertNotEqual(e1, e2)

    def test_expression_not_equals_with_offset(self):
        e1 = Expression({ Variable('x'): 2, Variable('y'): 5 }, 4)
        e2 = Expression({ Variable('y'): 5, Variable('x'): 2 })
        self.assertNotEqual(e1, e2)

    def test_expression_equals_variable(self):
        e = Expression({ Variable('x'): 1 })
        self.assertEqual(e, Variable('x'))

    def test_expression_with_coefficient_not_equals_variable(self):
        e = Expression({ Variable('x'): 2 })
        self.assertNotEqual(e, Variable('x'))

    def test_expression_with_offset_not_equals_variable(self):
        e = Expression({ Variable('x'): 1 }, 4)
        self.assertNotEqual(e, Variable('x'))

    def test_expression_with_multiple_not_equals_variable(self):
        e = Expression({ Variable('x'): 1, Variable('y'): 1 })
        self.assertNotEqual(e, Variable('x'))

    def test_expression_with_no_variables_equals_number(self):
        self.assertEqual(Expression({}, 3), 3)

    def test_expression_with_no_variables_not_equals_number(self):
        self.assertNotEqual(Expression({}, -3), 3)

    def test_expression_with_variables_not_equals_number(self):
        self.assertNotEqual(Expression({ Variable('x'): 1 }, 3), 3)

    def test_expression_parse(self):
        e = Expression('2x + 3')
        self.assertEqual(e, Expression({ Variable('x'): 2 }, 3))

    def test_expression_parse_number(self):
        e = Expression('1')
        self.assertEqual(e, Expression({}, 1))
        self.assertEqual(e, 1)

    def test_expression_parse_zero_offset(self):
        e = Expression('x + 4y + 0')
        self.assertEqual(e, Expression({ Variable('x'): 1, Variable('y'): 4 }))

    def test_expression_parse_multiple_occurences(self):
        e = Expression('x + 4y + 2x')
        self.assertEqual(e, Expression({ Variable('x'): 3, Variable('y'): 4 }))

    def test_expression_parse_negative_coefficient(self):
        e = Expression('-x + 4y - 3z')
        self.assertEqual(e, Expression({ Variable('x'): -1, Variable('y'): 4, Variable('z'): -3 }))

    def test_expression_parse_zero_coefficient(self):
        e = Expression('2x + 0y + 4')
        self.assertEqual(e, Expression({ Variable('x'): 2 }, 4))

    def test_expression_parse_long_variable_symbols(self):
        e = Expression('-2x1 + 5pi - 3x2')
        self.assertEqual(e, Expression({ Variable('x1'): -2, Variable('pi'): 5, Variable('x2'): -3 }))

    def test_expression_parse_no_space(self):
        e = Expression('x+3f+6h+10')
        self.assertEqual(e, Expression({ Variable('x'): 1, Variable('f'): 3, Variable('h'): 6 }, 10))

    def test_expression_parse_large_coefficient(self):
        e = Expression('x + 4200y')
        self.assertEqual(e, Expression({ Variable('x'): 1, Variable('y'): 4200 }))

    def test_expression_to_str(self):
        e = Expression({
            Variable('x1'): -2,
            Variable('pi'): 5,
            Variable('x2'): -3
        })
        self.assertEqual(str(e), '5pi - 2x1 - 3x2')


if __name__ == '__main__':
    unittest.main()
