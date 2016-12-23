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

from psamm.lpsolver import lp

try:
    from psamm.lpsolver import cplex
except ImportError:
    cplex = None

requires_solver = unittest.skipIf(cplex is None, 'Cplex solver not available')


@requires_solver
class TestCplexProblem(unittest.TestCase):
    def setUp(self):
        self.solver = cplex.Solver()

    def test_expression_to_string(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        expr = -3 * prob.var('x') + prob.var('y')
        self.assertEqual(str(expr), '-3*x + y')

    def test_expression_of_tuple_to_string(self):
        prob = self.solver.create_problem()
        prob.define(('x', 1), ('x', 2), lower=0, upper=10)
        expr = -prob.var(('x', 1)) + 2 * prob.var(('x', 2))
        self.assertEqual(str(expr), "-('x', 1) + 2*('x', 2)")

    def test_objective_reset_on_set_objective(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)

        # Solve first time, maximize x
        prob.set_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('x'), 10)

        # Solve second time, maximize y
        # If the objective is not properly cleared,
        # the second solve will still maximize x.
        prob.set_objective(prob.var('y'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('y'), 10)

    def test_quadratic_objective(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)

        prob.set_objective(-prob.var('x')**2 + 5*prob.var('x'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('x'), 2.5)

        # Solve second time, maximize y
        prob.set_objective(prob.var('y'))
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('y'), 10)

    def test_result_to_bool_conversion_on_optimal(self):
        """Run a feasible LP problem and check that the result is True"""
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)

        prob.set_linear_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertTrue(result)

    def test_result_to_bool_conversion_on_infeasible(self):
        """Run an infeasible LP problem and check that the result is False"""
        prob = self.solver.create_problem()
        prob.define('x', 'y', 'z', lower=0, upper=10)
        prob.add_linear_constraints(2*prob.var('x') == -prob.var('y'),
                                    prob.var('x') + prob.var('z') >= 6,
                                    prob.var('z') <= 3)
        prob.set_objective_sense(lp.ObjectiveSense.Maximize)
        prob.set_linear_objective(2*prob.var('x'))
        result = prob.solve()
        self.assertFalse(result)

    def test_constraint_delete(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        c1, = prob.add_linear_constraints(prob.var('x') <= 5)

        prob.set_objective_sense(lp.ObjectiveSense.Maximize)
        prob.set_linear_objective(prob.var('x'))

        result = prob.solve()
        self.assertAlmostEqual(result.get_value('x'), 5)

        # Delete constraint
        c1.delete()
        result = prob.solve()
        self.assertAlmostEqual(result.get_value('x'), 10)


if __name__ == '__main__':
    unittest.main()
