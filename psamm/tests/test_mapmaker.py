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
# Copyright 2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import unittest

from psamm import mapmaker
from psamm.lpsolver import generic
from psamm.formula import Formula, Atom, Radical
from psamm.reaction import Compound
from psamm.datasource.reaction import parse_reaction


class TestMapMakerWeigths(unittest.TestCase):
    def test_default_weight(self):
        self.assertEqual(mapmaker.default_weight(Atom.C), 1)
        self.assertEqual(mapmaker.default_weight(Atom.N), 0.4)
        self.assertEqual(mapmaker.default_weight(Atom.O), 0.4)
        self.assertEqual(mapmaker.default_weight(Atom.P), 0.4)
        self.assertEqual(mapmaker.default_weight(Atom.S), 1)
        self.assertEqual(mapmaker.default_weight(Radical('R')), 40.0)

    def test_weighted_formula(self):
        weighted = mapmaker._weighted_formula(
            Formula.parse('C1H2N3S4R'), mapmaker.default_weight)
        self.assertEqual(set(weighted), {
            (Atom.C, 1, 1),
            # H is discarded
            (Atom.N, 3, 0.4),
            (Atom.S, 4, 1),
            (Radical('R'), 1, 40.0)
        })


class TestMapMaker(unittest.TestCase):
    def setUp(self):
        try:
            self.solver = generic.Solver(integer=True)
        except generic.RequirementsError:
            self.skipTest('Unable to find an MILP solver for tests')

    def test_predict_compound_pairs(self):
        """Test prediction of HEX1 reaction."""
        reaction = parse_reaction(
            'atp[c] + glc-D[c] => adp[c] + g6p[c] + h[c]')
        formulas = {
            'atp': Formula.parse('C10H12N5O13P3'),
            'adp': Formula.parse('C10H12N5O10P2'),
            'glc-D': Formula.parse('C6H12O6'),
            'g6p': Formula.parse('C6H11O9P'),
            'h': Formula.parse('H')
        }

        solutions = list(
            mapmaker.predict_compound_pairs(reaction, formulas, self.solver))

        self.assertEqual(len(solutions), 1)
        self.assertEqual(solutions[0], {
            (Compound('atp', 'c'), Compound('adp', 'c')):
                Formula.parse('C10N5O10P2'),
            (Compound('glc-D', 'c'), Compound('g6p', 'c')):
                Formula.parse('C6O6'),
            (Compound('atp', 'c'), Compound('g6p', 'c')):
                Formula.parse('O3P'),
        })

    def test_predict_compound_pairs_unbalanced(self):
        """Test prediction of (non-sense) unbalanced reaction."""
        reaction = parse_reaction(
            'a[c] <=> b[c] + c[c]')
        formulas = {
            'a': Formula.parse('C10H12'),
            'b': Formula.parse('C9H11'),
            'c': Formula.parse('CO2'),
        }

        with self.assertRaises(mapmaker.UnbalancedReactionError):
            list(mapmaker.predict_compound_pairs(
                reaction, formulas, self.solver))

    def test_predict_compound_pairs_multiple(self):
        """Test prediction of reaction with multiple instances."""
        reaction = parse_reaction(
            'a[c] <=> (2) b[c] + c[c]')
        formulas = {
            'a': Formula.parse('C10H13O6'),
            'b': Formula.parse('C5H6O3'),
            'c': Formula.parse('H'),
        }
        solutions = list(
            mapmaker.predict_compound_pairs(reaction, formulas, self.solver))

        self.assertEqual(len(solutions), 1)
        self.assertEqual(solutions[0], {
            (Compound('a', 'c'), Compound('b', 'c')): Formula.parse('C10O6'),
        })

    def test_predict_with_multiple_solutions(self):
        """Test prediction where multiple solutions exist (TALA)."""
        reaction = parse_reaction(
            'r5p[c] + xu5p-D[c] <=> g3p[c] + s7p[c]')
        formulas = {
            'r5p': Formula.parse('C5H9O8P'),
            'xu5p-D': Formula.parse('C5H9O8P'),
            'g3p': Formula.parse('C3H5O6P'),
            's7p': Formula.parse('C7H13O10P'),
        }

        solutions = list(
            mapmaker.predict_compound_pairs(reaction, formulas, self.solver))
        self.assertEqual(len(solutions), 2)

        # Solution A
        self.assertIn({
            (Compound('r5p', 'c'), Compound('s7p', 'c')):
                Formula.parse('C5O8P'),
            (Compound('xu5p-D', 'c'), Compound('g3p', 'c')):
                Formula.parse('C3O6P'),
            (Compound('xu5p-D', 'c'), Compound('s7p', 'c')):
                Formula.parse('C2O2'),
        }, solutions)

        # Solution B
        self.assertIn({
            (Compound('xu5p-D', 'c'), Compound('s7p', 'c')):
                Formula.parse('C5O8P'),
            (Compound('r5p', 'c'), Compound('g3p', 'c')):
                Formula.parse('C3O6P'),
            (Compound('r5p', 'c'), Compound('s7p', 'c')):
                Formula.parse('C2O2'),
        }, solutions)


if __name__ == '__main__':
    unittest.main()
