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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import math
import unittest

from psamm.lpsolver import lp


class DummyRangedPropertyContainer(object):
    def __init__(self):
        self._property_rw = 5

    @lp.ranged_property(min=-5, max=10)
    def property_ro(self):
        return 2.0

    @lp.ranged_property(min=-50, max=17)
    def property_rw(self):
        return self._property_rw

    @property_rw.setter
    def property_rw(self, value):
        self._property_rw = value

    @property_rw.deleter
    def property_rw(self):
        self._property_rw = 0


class TestRangedProperty(unittest.TestCase):
    def setUp(self):
        self.c = DummyRangedPropertyContainer()

    def test_read_property_ro(self):
        self.assertEqual(self.c.property_ro.value, 2.0)

    def test_write_property_ro(self):
        with self.assertRaises(AttributeError):
            self.c.property_ro.value = 3.0

    def test_write_property_base(self):
        with self.assertRaises(AttributeError):
            self.c.property_ro = 3.0

    def test_delete_property_base(self):
        with self.assertRaises(AttributeError):
            del self.c.property_ro

    def test_delete_property_ro(self):
        with self.assertRaises(AttributeError):
            del self.c.property_ro.value

    def test_read_property_ro_bounds(self):
        self.assertEqual(self.c.property_ro.min, -5)
        self.assertEqual(self.c.property_ro.max, 10)

    def test_write_property_ro_bounds(self):
        with self.assertRaises(AttributeError):
            self.c.property_ro.min = -50
        with self.assertRaises(AttributeError):
            self.c.property_ro.max = 100

    def test_read_property_rw(self):
        self.assertEqual(self.c.property_rw.value, 5)

    def test_write_property_rw(self):
        self.c.property_rw.value = 10
        self.assertEqual(self.c.property_rw.value, 10)

    def test_delete_property_rw(self):
        del self.c.property_rw.value
        self.assertEqual(self.c.property_rw.value, 0)


class TestExpression(unittest.TestCase):
    def test_create_expression(self):
        e = lp.Expression()
        self.assertEqual(e.offset, 0)
        self.assertEqual(dict(e.values()), {})

    def test_create_expression_with_offset(self):
        e = lp.Expression(offset=42)
        self.assertEqual(e.offset, 42)
        self.assertEqual(dict(e.values()), {})

    def test_create_expression_with_variables(self):
        e = lp.Expression(variables={'x1': -5, 'x2': 100})
        self.assertEqual(e.offset, 0)
        self.assertEqual(dict(e.values()), {'x1': -5, 'x2': 100})

    def test_expression_variables(self):
        e = lp.Expression(variables={'x1': -1, 'x2': 1}, offset=100)
        self.assertEqual(set(e.variables()), {'x1', 'x2'})

    def test_add_expression_and_number(self):
        e = lp.Expression({'x1': -5}, 100)
        e1 = e + 45
        self.assertEqual(e1.offset, 145)
        self.assertEqual(dict(e1.values()), {'x1': -5})

    def test_add_expression_and_expression(self):
        e1 = lp.Expression({'x1': -5, 'x2': 42}, -1000)
        e2 = lp.Expression({'x2': -42, 'x3': -50}, 2)
        e3 = e1 + e2
        self.assertEqual(e3.offset, -998)
        self.assertEqual(dict(e3.values()), {'x1': -5, 'x2': 0, 'x3': -50})

    def test_add_expression_and_infinity(self):
        e = lp.Expression({'x1': -5}, 100)
        e1 = e + float('inf')
        self.assertEqual(e1.offset, float('inf'))

        e2 = e - float('inf')
        self.assertEqual(e2.offset, -float('inf'))

    def test_add_expression_and_infinity_expression(self):
        e = lp.Expression({'x1': -5}, float('inf'))
        e1 = e + lp.Expression({'x2': 4}, 100)
        self.assertEqual(e1.offset, float('inf'))

    def test_add_expression_and_number_in_place(self):
        e = lp.Expression({'x1': -5}, 100)
        e += 45
        self.assertEqual(e.offset, 145)
        self.assertEqual(dict(e.values()), {'x1': -5})

    def test_add_expression_and_expression_in_place(self):
        e1 = lp.Expression({'x1': -5, 'x2': 42}, -1000)
        e2 = lp.Expression({'x2': -42, 'x3': -50}, 2)
        e1 += e2

        self.assertEqual(e1.offset, -998)
        self.assertEqual(dict(e1.values()), {'x1': -5, 'x2': 0, 'x3': -50})

        self.assertEqual(e2.offset, 2)
        self.assertEqual(dict(e2.values()), {'x2': -42, 'x3': -50})

    def test_add_expression_and_infinity_in_place(self):
        e = lp.Expression({'x1': -5}, 100)
        e += float('inf')
        self.assertEqual(e.offset, float('inf'))

    def test_subtract_expressions_in_place(self):
        e = lp.Expression({'x1': 40}, 22)
        e -= lp.Expression({'x1': 5, 'x2': -20}, -2)
        self.assertEqual(e.offset, 24)
        self.assertEqual(dict(e.values()), {'x1': 35, 'x2': 20})

    def test_multiply_expression_and_number(self):
        e = lp.Expression({'x1': -5}, 100)
        e1 = 2 * e
        self.assertEqual(e1.offset, 200)
        self.assertEqual(dict(e1.values()), {'x1': -10})

    def test_multiply_expression_and_infinity(self):
        e = lp.Expression({'x1': -5})
        e1 = float('-inf') * e
        self.assertTrue(math.isnan(e1.offset))

    def test_multiply_expression_and_number_in_place(self):
        e = lp.Expression({'x1': -5}, 100)
        e *= 2
        self.assertEqual(e.offset, 200)
        self.assertEqual(dict(e.values()), {'x1': -10})

    def test_multiply_expression_and_infinity_in_place(self):
        e = lp.Expression({'x1': -5})
        e *= float('inf')
        self.assertTrue(math.isnan(e.offset))

    def test_multiply_expression_and_simple_expression(self):
        e1 = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e2 = lp.Expression({'x1': 1})
        e = e1 * e2
        self.assertEqual(e.offset, 0)
        self.assertEqual(dict(e.values()), {
            lp.Product(['x1', 'x1']): -5,
            lp.Product(['x1', 'x2']): 3,
            'x1': 10
        })

    def test_multiply_expression_and_expression(self):
        e1 = lp.Expression({'x1': 1, 'x2': 2, 'x3': 3}, offset=4)
        e2 = lp.Expression({'x2': 4, 'x4': 5}, offset=-5)
        e = e1 * e2
        self.assertEqual(e.offset, -20)
        self.assertEqual(dict(e.values()), {
            lp.Product(['x1', 'x2']): 4,
            lp.Product(['x2', 'x2']): 8,
            lp.Product(['x2', 'x3']): 12,
            lp.Product(['x1', 'x4']): 5,
            lp.Product(['x2', 'x4']): 10,
            lp.Product(['x3', 'x4']): 15,
            'x1': -5,
            'x2': 6,
            'x3': -15,
            'x4': 20
        })

    def test_multiply_expression_and_expression_in_place(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e *= lp.Expression({'x1': 1})
        self.assertEqual(e.offset, 0)
        self.assertEqual(dict(e.values()), {
            lp.Product(['x1', 'x1']): -5,
            lp.Product(['x1', 'x2']): 3,
            'x1': 10
        })

    def test_expression_pow_negative(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        with self.assertRaises(ValueError):
            e1 = e**-2

    def test_expression_pow_neagtive_in_place(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        with self.assertRaises(ValueError):
            e **= -2

    def test_expression_pow_zero(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e1 = e**0
        self.assertEqual(e1.offset, 1)
        self.assertEqual(dict(e1.values()), {})

    def test_expression_pow_zero_in_place(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e **= 0
        self.assertEqual(e.offset, 1)
        self.assertEqual(dict(e.values()), {})

    def test_expression_pow_one(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e1 = e**1
        self.assertEqual(e1.offset, e.offset)
        self.assertEqual(dict(e1.values()), dict(e.values()))

    def test_expression_pow_one_in_place(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e **= 1
        self.assertEqual(e.offset, 10)
        self.assertEqual(dict(e.values()), {'x1': -5, 'x2': 3})

    def test_expression_pow_two(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e1 = e**2
        e2 = e*e
        self.assertEqual(e1.offset, e2.offset)
        self.assertEqual(dict(e1.values()), dict(e2.values()))

    def test_expression_pow_two_in_place(self):
        e = lp.Expression({'x1': -5, 'x2': 3}, offset=10)
        e **= 2
        self.assertEqual(e.offset, 100)
        self.assertEqual(dict(e.values()), {
            lp.Product(['x1', 'x1']): 25,
            lp.Product(['x1', 'x2']): -30,
            lp.Product(['x2', 'x2']): 9,
            'x1': -100,
            'x2': 60
        })

    def test_negate_expression(self):
        e = lp.Expression({'x1': 52}, -32)
        e1 = -e
        self.assertEqual(e1.offset, 32)
        self.assertEqual(dict(e1.values()), {'x1': -52})

    def test_expression_to_string(self):
        e = lp.Expression({'x1': -4, 'x2': 100, 'x3': 1}, 42)
        self.assertEqual(str(e), '-4*x1 + 100*x2 + x3 + 42')

    def test_expression_with_product_to_string(self):
        e = lp.Expression({'x1': -4, lp.Product(['x1', 'x2']): 2}, 12)
        self.assertEqual(str(e), '-4*x1 + 2*x1*x2 + 12')

    def test_expression_with_tuple_vars_to_string(self):
        e = lp.Expression({('v', 1): 1}, -1)
        self.assertEqual(str(e), "('v', 1) - 1")

    def test_expression_contains(self):
        e = lp.Expression({'x1': 10, 'x2': -5})
        self.assertIn('x1', e)
        self.assertNotIn('x3', e)

    def test_expression_contains_product(self):
        e = lp.Expression({lp.Product(['x1', 'x1']): -1})
        self.assertIn(lp.Product(['x1', 'x1']), e)
        self.assertNotIn('x1', e)


class TestRelation(unittest.TestCase):
    def test_create_relation(self):
        e = lp.Expression({'x1': 4})
        r = lp.Relation(lp.RelationSense.Greater, e)
        self.assertEqual(r.expression, e)
        self.assertEqual(r.sense, lp.RelationSense.Greater)

    def test_relation_with_offset_to_string(self):
        e = lp.Expression({'x1': 4}, -20)
        r = lp.Relation(lp.RelationSense.Less, e)
        self.assertEqual(str(r), '4*x1 <= 20')

    def test_relation_without_offset_to_string(self):
        e = lp.Expression({'x1': 1})
        r = lp.Relation(lp.RelationSense.Equals, e)
        self.assertEqual(str(r), 'x1 == 0')


class MockResult(lp.Result):
    def __init__(self, values):
        self._values = dict(values)

    @property
    def success(self):
        return True

    @property
    def status(self):
        return 'Success'

    @property
    def unbounded(self):
        return False

    def _has_variable(self, var):
        return var in self._values

    def _get_value(self, var):
        return self._values[var]


class TestResult(unittest.TestCase):
    def test_result_to_bool(self):
        result = MockResult({'x': 3})
        self.assertTrue(bool(result))

    def test_result_get_simple_value(self):
        result = MockResult({'x': 3, 'y': 4})
        self.assertEqual(result.get_value('x'), 3)
        self.assertEqual(result.get_value('y'), 4)

    def test_result_get_invalid_value(self):
        result = MockResult({'x': 3})
        with self.assertRaises(ValueError):
            y = result.get_value('y')

    def test_result_get_expression_value(self):
        result = MockResult({'x': 3, 'y': 4})
        expr = lp.Expression({'x': 2, 'y': -6}, offset=3)
        self.assertEqual(result.get_value(expr), -15)

    def test_result_get_expression_product_value(self):
        result = MockResult({'x1': 2, 'x2': 3, 'x3': 4})
        expr = lp.Expression({lp.Product(['x1', 'x2']): 3, 'x3': -1}, offset=4)
        self.assertEqual(result.get_value(expr), 18)
