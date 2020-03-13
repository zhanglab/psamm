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
# Copyright 2016-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Functions for predicting compound pairs using the MapMaker algorithm.

The MapMaker algorithm is described in [Tervo16]_.
"""

from itertools import product

from six import iteritems

from .lpsolver import lp
from .formula import Atom, Radical, Formula


def default_weight(element):
    """Return weight of formula element.

    This implements the default weight proposed for MapMaker.
    """
    if element in (Atom.N, Atom.O, Atom.P):
        return 0.4
    elif isinstance(element, Radical):
        return 40.0
    return 1.0


def _weighted_formula(form, weight_func):
    """Yield weight of each formula element."""
    for e, mf in form.items():
        if e == Atom.H:
            continue

        yield e, mf, weight_func(e)


def _transfers(reaction, delta, elements, result, epsilon):
    """Yield transfers obtained from result."""
    left = set(c for c, _ in reaction.left)
    right = set(c for c, _ in reaction.right)
    for c1, c2 in product(left, right):
        items = {}
        for e in elements:
            v = result.get_value(delta[c1, c2, e])
            nearest_int = round(v)
            if abs(v - nearest_int) < epsilon:
                v = int(nearest_int)
            if v >= epsilon:
                items[e] = v

        if len(items) > 0:
            yield (c1, c2), Formula(items)


def _reaction_to_dicts(reaction):
    """Convert a reaction to reduced left, right dictionaries."""

    def dict_from_iter_sum(it):
        d = {}
        for k, v in it:
            if k not in d:
                d[k] = 0
            d[k] += v
        return d

    left = dict_from_iter_sum(reaction.left)
    right = dict_from_iter_sum(reaction.right)

    return left, right


class UnbalancedReactionError(ValueError):
    """Raised when an unbalanced reaction is provided."""


def predict_compound_pairs(reaction, compound_formula, solver, epsilon=1e-5,
                           alt_elements=None, weight_func=default_weight):
    """Predict compound pairs of reaction using MapMaker.

    Yields all solutions as dictionaries with compound pairs as keys and
    formula objects as values.

    Args:
        reaction: :class:`psamm.reaction.Reaction` object.
        compound_formula: Dictionary mapping compound IDs to formulas. Formulas
            must be flattened.
        solver: LP solver (MILP).
        epsilon: Threshold for rounding floating-point values to integers in
            the predicted transfers.
        alt_elements: Iterable of elements to consider for alternative
              solutions. Only alternate solutions that have different transfers
              of these elements will be returned (default=all elements).
        weight_func: Optional function that returns a weight for a formula
            element (should handle specific Atom and Radical objects). By
            default, the standard MapMaker weights will be used
            (H=0, R=40, N=0.4, O=0.4, P=0.4, *=1).
    """
    elements = set()
    for compound, value in reaction.compounds:
        if compound.name not in compound_formula:
            return

        f = compound_formula[compound.name]
        elements.update(e for e, _, _ in _weighted_formula(f, weight_func))

    if len(elements) == 0:
        return

    p = solver.create_problem()

    x = p.namespace(name='x')
    y = p.namespace(name='y')
    m = p.namespace(name='m')
    q = p.namespace(name='q')
    omega = p.namespace(name='omega')
    gamma = p.namespace(name='gamma')
    delta = p.namespace(name='delta')

    left, right = _reaction_to_dicts(reaction)

    objective = 0
    for c1, c2 in product(left, right):
        m.define([(c1, c2)], lower=0)
        q.define([(c1, c2)], lower=0)
        omega.define([(c1, c2)], types=lp.VariableType.Binary)

        objective += 100 * m[c1, c2] + 90 * q[c1, c2] - 0.1 * omega[c1, c2]

        gamma.define(((c1, c2, e) for e in elements),
                     types=lp.VariableType.Binary)
        delta.define(((c1, c2, e) for e in elements), lower=0)

        # Eq 12
        x.define([(c1, c2)], types=lp.VariableType.Binary)
        y.define([(c1, c2)], types=lp.VariableType.Binary)
        p.add_linear_constraints(y[c1, c2] <= 1 - x[c1, c2])

        # Eq 6, 9
        delta_wsum = delta.expr(
            ((c1, c2, e), weight_func(e)) for e in elements)
        p.add_linear_constraints(
            m[c1, c2] <= delta_wsum, q[c1, c2] <= delta_wsum)

        objective -= gamma.sum((c1, c2, e) for e in elements)

    p.set_objective(objective)

    # Eq 3
    zs = {}
    for c1, v in iteritems(left):
        f = compound_formula[c1.name]
        for e in elements:
            mf = f.get(e, 0)
            delta_sum = delta.sum((c1, c2, e) for c2 in right)
            try:
                p.add_linear_constraints(delta_sum == v * mf)
            except ValueError:
                raise UnbalancedReactionError('Unable to add constraint')

        # Eq 8, 11
        x_sum = x.sum((c1, c2) for c2 in right)
        y_sum = y.sum((c1, c2) for c2 in right)
        p.add_linear_constraints(x_sum <= 1, y_sum <= 1)

        # Eq 13
        zs[c1] = 0
        for e, mf, w in _weighted_formula(f, weight_func):
            zs[c1] += w * mf

    # Eq 2
    for c2, v in iteritems(right):
        f = compound_formula[c2.name]
        for e in elements:
            mf = f.get(e, 0)
            delta_sum = delta.sum((c1, c2, e) for c1 in left)
            try:
                p.add_linear_constraints(delta_sum == v * mf)
            except ValueError:
                raise UnbalancedReactionError('Unable to add constraint')

    for c1, v1 in iteritems(left):
        for c2 in right:
            f1 = compound_formula[c1.name]
            f2 = compound_formula[c2.name]
            for e in elements:
                mf = f1.get(e, 0)

                # Eq 4
                p.add_linear_constraints(
                    delta[c1, c2, e] <= float(v1) * mf * omega[c1, c2])

                # Eq 5
                if e in f2:
                    p.add_linear_constraints(
                        delta[c1, c2, e] <= float(v1) * mf * gamma[c1, c2, e])

            # Eq 7, 10
            p.add_linear_constraints(
                m[c1, c2] <= float(v1) * zs[c1] * x[c1, c2],
                q[c1, c2] <= float(v1) * zs[c1] * y[c1, c2])

    try:
        result = p.solve(lp.ObjectiveSense.Maximize)
    except lp.SolverError:
        raise UnbalancedReactionError('Unable to solve')

    first_objective = result.get_value(objective)

    # Yield the first result
    yield dict(_transfers(reaction, delta, elements, result, epsilon))

    # Find alternative solutions
    if alt_elements is None:
        alt_elements = elements
    else:
        alt_elements = elements.intersection(alt_elements)

    if len(alt_elements) == 0:
        return

    # Add a constraint that disallows the current transfer solution until no
    # alternative solutions are left.
    add_objective_constraint = True
    while True:
        elem_vars = 0
        elem_count = 0
        for c1, c2 in product(left, right):
            for e in alt_elements:
                if result.get_value(gamma[c1, c2, e]) > 0.5:
                    elem_count += 1
                    elem_vars += gamma[c1, c2, e]

        if elem_count == 0:
            break

        if add_objective_constraint:
            p.add_linear_constraints(objective >= first_objective)
            add_objective_constraint = False

        p.add_linear_constraints(elem_vars <= elem_count - 1)
        try:
            result = p.solve(lp.ObjectiveSense.Maximize)
        except lp.SolverError:
            break

        yield dict(_transfers(reaction, delta, elements, result, epsilon))
