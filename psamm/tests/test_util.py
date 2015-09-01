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

import unittest

from psamm.util import MaybeRelative


class TestMaybeRelative(unittest.TestCase):
    def test_init_from_float(self):
        arg = MaybeRelative(24.5)
        self.assertFalse(arg.relative)
        self.assertIsNone(arg.reference)
        self.assertAlmostEqual(float(arg), 24.5)

    def test_init_from_float_string(self):
        arg = MaybeRelative('-10')
        self.assertFalse(arg.relative)
        self.assertIsNone(arg.reference)
        self.assertAlmostEqual(float(arg), -10.0)

    def test_init_from_percentage(self):
        arg = MaybeRelative('110%')
        self.assertTrue(arg.relative)
        self.assertIsNone(arg.reference)
        with self.assertRaises(ValueError):
            float(arg)

    def test_resolve_relative(self):
        arg = MaybeRelative('40%')
        arg.reference = 200.0
        self.assertAlmostEqual(float(arg), 80.0)

    def test_invalid(self):
        with self.assertRaises(ValueError):
            arg = MaybeRelative('abc')


if __name__ == '__main__':
    unittest.main()
