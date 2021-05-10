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
# Copyright 2018-2020  Jing Wang <wjingsjtu@gmail.com>

"""Utility functions."""
from __future__ import division
from future.utils import iteritems
from builtins import range

import re
from itertools import product

from psamm.formula import Atom
from psamm.expression.boolean import Expression


def id_equals(id1, id2):
    """Return True if the two IDs are considered equal."""
    return id1.lower() == id2.lower()


def name_equals(name1, name2):
    """Return True if the two names are considered equal."""
    if name1 is None or name2 is None:
        return False
    # remove special chars, translate to lower case
    pattern = r'[^a-zA-Z0-9]'
    name1 = re.sub(pattern, '', name1.lower())
    name2 = re.sub(pattern, '', name2.lower())
    # unify Coenzyme A and CoA
    pattern = r'coenzymea'
    name1 = re.sub(pattern, 'coa', name1.lower())
    name2 = re.sub(pattern, 'coa', name2.lower())
    return (name1 == name2)


def name_similar(name1, name2):
    """Return the possibility that two names are considered equal."""
    pattern = r'[()\[\]_, -]'
    if name1 is None or name2 is None:
        return 0.01
    name1 = re.sub(pattern, '', name1)
    name2 = re.sub(pattern, '', name2)
    return max(0.01,
               (1 - float(levenshtein(name1, name2)) /
                max(len(name1), len(name2))))


def formula_equals(f1, f2, charge1, charge2):
    """Return True if the two formulas are considered equal."""
    if f1 is None or f2 is None:
        return False

    def calc_formula(f, charge, neutral=False):
        formula = {}
        for e, value in iteritems(f):
            if e != Atom('H') or charge is None or not neutral:
                formula[e] = value
            else:
                formula[e] = value - charge  # No. of H in neutral stat
        return formula

    # in case some model doesn't provide charge for non-neutral formula
    # compare original formula
    formula1 = calc_formula(f1, charge1)
    formula2 = calc_formula(f2, charge2)
    # compare neutral formula
    formula1_neutral = calc_formula(f1, charge1, neutral=True)
    formula2_neutral = calc_formula(f2, charge2, neutral=True)

    return (formula1 == formula2) or (formula1_neutral == formula2_neutral)


def genes_equals(g1, g2, gene_map={}):
    """Return True if the two gene association strings are considered equal.

    Args:
        g1, g2: gene association strings
        gene_map: a dict that maps gene ids in g1 to gene ids in g2 if
                  they have different naming system
    """
    if g1 is None or g2 is None:
        return False

    e1 = Expression(g1)
    e2 = Expression(g2)

    g_list = set([gene_map.get(v.symbol, v.symbol) for v in e2.variables])
    check1 = e1.substitute(lambda v: v.symbol in g_list)

    g_list = set([gene_map.get(v.symbol, v.symbol) for v in e1.variables])
    check2 = e2.substitute(lambda v: v.symbol in g_list)

    return check1.value and check2.value


def formula_exact(f1, f2):
    """Return True if the two formulas are considered equal."""
    if f1 is None or f2 is None:
        return False

    formula1 = {e: value for e, value in iteritems(f1)}
    formula2 = {e: value for e, value in iteritems(f2)}
    return formula1 == formula2


def pairwise_distance(i1, i2, distance, threshold=None):
    """Pairwise distances."""
    for c1, c2 in product(i1, i2):
        score = distance(c1, c2)
        if threshold is None or score < threshold:
            yield (c1, c2), score


def levenshtein(s1, s2):
    """Edit distance function."""
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1.lower()):
        current_row = [i + 1]
        for j, c2 in enumerate(s2.lower()):
            # We use j+1 instead of j since previous_row and current_row
            # are one character longer than s2
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def jaccard(s1, s2):
    """Jaccard similarity function."""
    return len(set(s1) & set(s2)) / len(set(s1) | set(s2))
