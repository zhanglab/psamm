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

    def test_negate_expression(self):
        e = lp.Expression({'x1': 52}, -32)
        e1 = -e
        self.assertEqual(e1.offset, 32)
        self.assertEqual(dict(e1.values()), {'x1': -52})

    def test_expression_to_string(self):
        e = lp.Expression({'x1': -4, 'x2': 100, 'x3': 1}, 42)
        self.assertEqual(str(e), '-4*x1 + 100*x2 + x3 + 42')

    def test_expression_with_tuple_vars_to_string(self):
        e = lp.Expression({('v', 1): 1}, -1)
        self.assertEqual(str(e), "('v', 1) - 1")

    def test_expression_contains(self):
        e = lp.Expression({'x1': 10, 'x2': -5})
        self.assertIn('x1', e)
        self.assertNotIn('x3', e)


class TestRelation(unittest.TestCase):
    def test_create_relation(self):
        e = lp.Expression({'x1': 4})
        r = lp.Relation(lp.Relation.Greater, e)
        self.assertEqual(r.expression, e)
        self.assertEqual(r.sense, lp.Relation.Greater)

    def test_relation_with_offset_to_string(self):
        e = lp.Expression({'x1': 4}, -20)
        r = lp.Relation(lp.Relation.Less, e)
        self.assertEqual(str(r), '4*x1 <= 20')

    def test_relation_without_offset_to_string(self):
        e = lp.Expression({'x1': 1})
        r = lp.Relation(lp.Relation.Equals, e)
        self.assertEqual(str(r), 'x1 == 0')
