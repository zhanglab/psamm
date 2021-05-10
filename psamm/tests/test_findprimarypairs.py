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
# Copyright 2020 Elysha Sameth <esameth1@my.uri.edu>


import unittest

from psamm import findprimarypairs
from psamm.formula import Formula, Atom, Radical
from psamm.reaction import Compound
from psamm.datasource.reaction import parse_reaction
from itertools import permutations


class TestFindPrimaryPairs(unittest.TestCase):
    def test_default_weight(self):
        self.assertEqual(findprimarypairs.element_weight(Atom.C), 1)
        self.assertEqual(findprimarypairs.element_weight(Atom.H), 0)
        self.assertEqual(findprimarypairs.element_weight(Atom.N), 0.82)
        self.assertEqual(findprimarypairs.element_weight(Radical('R')), 0.82)

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
        transfer, balance, ambiguous_pairs = \
            findprimarypairs.predict_compound_pairs(
                reaction, formulas)

        self.assertEqual(balance, {})
        self.assertEqual(transfer, {
            ((Compound('atp', 'c'), 1), (Compound('adp', 'c'), 1)):
                Formula.parse('C10H12N5O10P2'),
            ((Compound('glc-D', 'c'), 1), (Compound('g6p', 'c'), 1)):
                Formula.parse('C6H11O6'),
            ((Compound('atp', 'c'), 1), (Compound('g6p', 'c'), 1)):
                Formula.parse('O3P'),
            ((Compound('glc-D', 'c'), 1), (Compound('h', 'c'), 1)):
                Formula.parse('H')
        })
        self.assertEqual(ambiguous_pairs, set())

    def test_predict_compound_pairs_ambig(self):
        """Test prediction of (non-sense) ambiguous reaction."""
        reaction = parse_reaction(
            'a[c] + b[c] => c[c] + d[c]')
        formulas = {
            'a': Formula.parse('C10H10'),
            'b': Formula.parse('C10H10'),
            'c': Formula.parse('C10H10'),
            'd': Formula.parse('C10H10')
        }
        transfer, balance, ambiguous_pairs = \
            findprimarypairs.predict_compound_pairs(
                reaction, formulas)

        all_perm = set(set(permutations(list(ambiguous_pairs)[0])))
        ambig = list(ambiguous_pairs)[0]

        self.assertEqual(balance, {})
        self.assertEqual(transfer, {
            ((Compound('b', 'c'), 1), (Compound('d', 'c'), 1)):
                Formula.parse('C10H10'),
            ((Compound('a', 'c'), 1), (Compound('c', 'c'), 1)):
                Formula.parse('C10H10')
        })
        self.assertEqual(ambig, ambig if ambig in all_perm else set())

    def test_predict_compound_pairs_unbalanced(self):
        """Test prediction of (non-sense) unbalanced reaction."""
        reaction = parse_reaction(
            'a[c] <=> b[c] + c[c]')
        formulas = {
            'a': Formula.parse('C10H12'),
            'b': Formula.parse('C9H11'),
            'c': Formula.parse('CO2')
        }
        transfer, balance, ambiguous_pairs = \
            findprimarypairs.predict_compound_pairs(
                reaction, formulas)

        self.assertEqual(balance, {
            (Compound('a', 'c'), 1): Formula.parse('H'),
            (Compound('c', 'c'), 1): Formula.parse('O2'),
        })
        self.assertEqual(transfer, {
            ((Compound('a', 'c'), 1), (Compound('b', 'c'), 1)):
                Formula.parse('C9H11'),
            ((Compound('a', 'c'), 1), (Compound('c', 'c'), 1)):
                Formula.parse('C'),
        })
        self.assertEqual(ambiguous_pairs, set())

    def test_predict_compound_pairs_multiple(self):
        """Test prediction of reaction with multiple instances."""
        reaction = parse_reaction(
            'a[c] <=> (2) b[c] + c[c]')
        formulas = {
            'a': Formula.parse('C10H13O6'),
            'b': Formula.parse('C5H6O3'),
            'c': Formula.parse('H'),
        }
        transfer, balance, ambiguous_pairs = \
            findprimarypairs.predict_compound_pairs(
                reaction, formulas)

        self.assertEqual(balance, {})
        self.assertEqual(transfer, {
            ((Compound('a', 'c'), 1), (Compound('b', 'c'), 1)):
                Formula.parse('C5H6O3'),
            ((Compound('a', 'c'), 1), (Compound('b', 'c'), 2)):
                Formula.parse('C5H6O3'),
            ((Compound('a', 'c'), 1), (Compound('c', 'c'), 1)):
                Formula.parse('H')
        })
        self.assertEqual(ambiguous_pairs, set())
