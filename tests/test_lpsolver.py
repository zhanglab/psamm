#!/usr/bin/env python

import unittest

from metnet import lpsolver

class TestCplexProblem(unittest.TestCase):
    def setUp(self):
        self.solver = lpsolver.CplexSolver(None)

    def test_objective_reset_on_set_linear_objective(self):
        prob = self.solver.create_problem()
        prob.define('x', 'y', lower=0, upper=10)
        prob.add_linear_constraints(prob.var('x') + prob.var('y') <= 12)
        prob.set_objective_sense(lpsolver.CplexProblem.Maximize)

        # Solve first time, maximize x
        prob.set_linear_objective(2*prob.var('x'))
        prob.solve()
        self.assertAlmostEqual(prob.get_value('x'), 10)

        # Solve second time, maximize y
        # If the objective is not properly cleared,
        # the second solve will still maximize x.
        prob.set_linear_objective(prob.var('y'))
        prob.solve()
        self.assertAlmostEqual(prob.get_value('y'), 10)

if __name__ == '__main__':
    unittest.main()
