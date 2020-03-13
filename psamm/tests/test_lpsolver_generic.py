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
# Copyright 2015-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import os
import unittest
from decimal import Decimal
from fractions import Fraction

from psamm.lpsolver import generic, lp


class TestSolverProblem(unittest.TestCase):
    """Test the current solver chosen by generic.

    To test all solvers, run this test multiple times with the PSAMM_SOLVER
    environment variable set to the solver to test.
    """
    def setUp(self):
        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_set_decimal_variable_bounds(self):
        """Test that Decimal can be used for variable bounds."""
        prob = self.solver.create_problem()
        prob.define('x', lower=Decimal('3.5'), upper=Decimal('500.3'))

    def test_set_decimal_objective(self):
        """Test that Decimal can be used in objective."""
        prob = self.solver.create_problem()
        prob.define('x', lower=3, upper=100)
        prob.set_objective(Decimal('3.4') * prob.var('x'))

    def test_set_fraction_objective(self):
        """Test that Fraction can be used in objective."""
        prob = self.solver.create_problem()
        prob.define('x', lower=-5, upper=5)
        prob.set_objective(Fraction(1, 2) * prob.var('x'))

    def test_set_decimal_constraint(self):
        """Test that Decimal can be used in constraints."""
        prob = self.solver.create_problem()
        prob.define('x', lower=0, upper=5)
        prob.add_linear_constraints(
            Decimal('5.6') * prob.var('x') >= Decimal('8.9'))

    def test_result_evaluate_decimal_expression(self):
        """Test that Decimal in expression can be used to evaluate result."""
        prob = self.solver.create_problem()
        prob.define('x', lower=0, upper=Fraction(5, 2))
        prob.add_linear_constraints(prob.var('x') >= 2)
        prob.set_objective(prob.var('x'))
        result = prob.solve(lp.ObjectiveSense.Maximize)
        value = result.get_value(Decimal('3.4') * prob.var('x'))
        self.assertAlmostEqual(value, Fraction(17, 2))

    def test_redefine_variable(self):
        """Test that redefining a variable fails."""
        prob = self.solver.create_problem()
        prob.define('x', lower=3, upper=100)
        with self.assertRaises(ValueError):
            prob.define('x', lower=3, upper=100)

    def test_has_variable(self):
        """Test whether has_variable() works."""
        prob = self.solver.create_problem()
        self.assertFalse(prob.has_variable('x'))
        prob.define('x', lower=-400)
        self.assertTrue(prob.has_variable('x'))


class TestListSolversCommand(unittest.TestCase):
    def test_list_lpsolvers(self):
        if os.environ.get('PSAMM_SOLVER') in ('nosolver', None):
            with self.assertRaises(SystemExit):
                generic.list_solvers([])
        else:
            generic.list_solvers([])


if __name__ == '__main__':
    unittest.main()
