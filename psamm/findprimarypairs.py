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
# Copyright 2018-2020  Ke Zhang <kzhang@my.uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import division

import logging
from itertools import product
from collections import Counter

from six import iteritems
from six.moves import reduce

from .formula import Formula, Atom

from fractions import gcd

logger = logging.getLogger(__name__)


def element_weight(element):
    """Return default element weight.

    This is the default element weight function.
    """
    if element == Atom.H:
        return 0.0
    elif element == Atom.C:
        return 1.0
    return 0.82


def _jaccard_similarity(f1, f2, weight_func):
    """Calculate generalized Jaccard similarity of formulas.

    Returns the weighted similarity value or None if there is no overlap
    at all. If the union of the formulas has a weight of zero (i.e. the
    denominator in the Jaccard similarity is zero), a value of zero is
    returned.
    """
    elements = set(f1)
    elements.update(f2)

    count, w_count, w_total = 0, 0, 0
    for element in elements:
        mi = min(f1.get(element, 0), f2.get(element, 0))
        mx = max(f1.get(element, 0), f2.get(element, 0))
        count += mi
        w = weight_func(element)
        w_count += w * mi
        w_total += w * mx

    if count == 0:
        return None

    return 0.0 if w_total == 0.0 else w_count / w_total


def _reaction_to_dicts(reaction):
    """Convert a reaction to reduced left, right dictionaries.

    Returns a pair of (left, right) dictionaries mapping compounds to
    normalized integer stoichiometric values. If a compound occurs multiple
    times on one side, the occurences are combined into a single entry in the
    dictionary.
    """
    def dict_from_iter_sum(it, div):
        d = {}
        for k, v in it:
            if k not in d:
                d[k] = 0
            d[k] += int(v / div)
        return d

    div = reduce(gcd, (abs(v) for _, v in reaction.compounds), 0)
    if div == 0:
        raise ValueError('Empty reaction')

    left = dict_from_iter_sum(reaction.left, div)
    right = dict_from_iter_sum(reaction.right, div)

    return left, right


class _CompoundInstance(object):
    """Instance of a compound in a reaction during matching.

    A compound in a reaction has a number of instances corresponding to the
    stoichiometric value. The instance is defined by the compound and index
    number. The compound instance also has an attached formula which can be
    updated through the ``formula`` property. The compound instance can be
    created as a copy of an existing instance, or from a
    :psamm:`psamm.reaction.Compound`, number, and
    :psamm:`psamm.formula.Formula`.
    """
    def __init__(self, *args):
        if len(args) == 1 and isinstance(args[0], _CompoundInstance):
            self._compound = args[0]._compound
            self._index = args[0]._index
            self._formula = args[0]._formula
        elif len(args) == 3:
            self._compound = args[0]
            self._index = args[1]
            self._formula = args[2]
        else:
            raise ValueError('Invalid arguments')

    @property
    def compound(self):
        """Return the :class:`psamm.reaction.Compound`."""
        return self._compound

    @property
    def index(self):
        """Return the index."""
        return self._index

    @property
    def formula(self):
        """Return the formula."""
        return self._formula

    @formula.setter
    def formula(self, value):
        self._formula = value

    def __repr__(self):
        return str('<{} {} of {!r}>').format(
            self.__class__.__name__, self._index, self._compound)

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and
                self._compound == other._index and
                self._index == other._index)

    def __hash__(self):
        return hash((self._compound, self._index))


def predict_compound_pairs_iterated(
        reactions, formulas, ambiguous=False, prior=(1, 43),
        max_iterations=None, element_weight=element_weight):
    """Predict reaction pairs using iterated method.

    Returns a tuple containing a dictionary of predictions keyed by the
    reaction IDs, and the final number of iterations. Each reaction prediction
    entry contains a tuple with a dictionary of transfers and a dictionary of
    unbalanced compounds. The dictionary of unbalanced compounds is empty only
    if the reaction is balanced.

    Args:
        reactions: Dictionary or pair-iterable of (id, equation) pairs.
            IDs must be any hashable reaction identifier (e.g. string) and
            equation must be :class:`psamm.reaction.Reaction` objects.
        ambiguous: True or False value to indicate if the ambiguous
            reactions should be printed in the debugging output.
        formulas: Dictionary mapping compound IDs to
            :class:`psamm.formula.Formula`. Formulas must be flattened.
        prior: Tuple of (alpha, beta) parameters for the MAP inference.
            If not provided, the default parameters will be used: (1, 43).
        max_iterations: Maximum iterations to run before stopping. If the
            stopping condition is reached before this number of iterations,
            the procedure also stops. If None, the procedure only stops when
            the stopping condition is reached.
        element_weight: A function providing returning weight value for the
            given :class:`psamm.formula.Atom` or
            :class:`psamm.formula.Radical`. If not provided, the default weight
            will be used (H=0, C=1, *=0.82)
    """
    prior_alpha, prior_beta = prior
    reactions = dict(reactions)

    pair_reactions = {}
    possible_pairs = Counter()
    tie_breakers = {}
    for reaction_id, equation in iteritems(reactions):
        tie_breakers[reaction_id] = {}
        for (c1, _), (c2, _) in product(equation.left, equation.right):
            spair = tuple(sorted([c1.name, c2.name]))
            possible_pairs[spair] += 1
            pair_reactions.setdefault(spair, set()).add(reaction_id)

    final_ambiguous_reactions = {}

    def print_ambiguous_summary():
        '''Function for summarizing and printing ambiguous compound pairs

        '''
        final_ambiguous_count = 0
        final_ambiguous_hydrogen = 0
        for reaction, pairs in sorted(iteritems(final_ambiguous_reactions)):
            all_hydrogen = False
            if len(pairs) > 0:
                final_ambiguous_count += 1
                all_hydrogen = True

                for pair in pairs:
                    for _, _, formula in pair:
                        df = dict(formula.items())
                        if len(df) > 1 or Atom.H not in df:
                            all_hydrogen = False

                if all_hydrogen:
                    final_ambiguous_hydrogen += 1
                    logger.info('Ambiguous Hydrogen '
                                'Transfers in Reaction: {}'.format(reaction))
                else:
                    logger.info('Ambiguous Non-Hydrogen '
                                'Transfers in Reaction: {}'.format(reaction))

        logger.info('{} reactions were decided with ambiguity'.format(
            final_ambiguous_count))
        logger.info('Only hydrogen in decision: {} vs non hydrogen: {}'.format(
            final_ambiguous_hydrogen,
            final_ambiguous_count - final_ambiguous_hydrogen))

    next_reactions = set(reactions)
    pairs_predicted = None
    prediction = {}
    weights = {}
    iteration = 0
    while len(next_reactions) > 0:
        iteration += 1
        if max_iterations is not None and iteration > max_iterations:
            break

        logger.info('Iteration {}: {} reactions...'.format(
            iteration, len(next_reactions)))

        for reaction_id in next_reactions:
            result = predict_compound_pairs(
                reactions[reaction_id], formulas, weights, element_weight)
            if result is None:
                continue

            transfer, balance, ambiguous_pairs = result
            final_ambiguous_reactions[reaction_id] = ambiguous_pairs

            rpairs = {}
            for ((c1, _), (c2, _)), form in iteritems(transfer):
                rpairs.setdefault((c1, c2), []).append(form)

            prediction[reaction_id] = rpairs, balance
        if ambiguous is True:
            logger.info('Ambiguous Reactions:')
            print_ambiguous_summary()

        pairs_predicted = Counter()
        for reaction_id, (rpairs, _) in iteritems(prediction):
            for c1, c2 in rpairs:
                spair = tuple(sorted([c1.name, c2.name]))
                pairs_predicted[spair] += 1

        next_reactions = set()
        for spair, total in sorted(iteritems(possible_pairs)):
            pred = pairs_predicted[spair]

            # The weight is set to the maximum a posteriori (MAP) estimate
            # of the primary pair probability distribution.
            posterior_alpha = prior_alpha + pred
            posterior_beta = prior_beta + total - pred
            pair_weight = ((posterior_alpha - 1) /
                           (posterior_alpha + posterior_beta - 2))

            if (spair not in weights or
                    abs(pair_weight - weights[spair]) > 1e-5):
                next_reactions.update(pair_reactions[spair])

            c1, c2 = spair
            weights[c1, c2] = pair_weight
            weights[c2, c1] = pair_weight

    return prediction, iteration


def _match_greedily(reaction, compound_formula, score_func):
    """Match compounds greedily based on score function.

    Args:
        reaction: Reaction equation :class:`psamm.reaction.Reaction`.
        compound_formula: Dictionary mapping compound IDs to
            :class:`psamm.formula.Formula`. Formulas must be flattened.
        score_func: Function that takes two :class:`_CompoundInstance` and
            returns the score.
    """
    uninstantiated_left, uninstantiated_right = _reaction_to_dicts(reaction)
    ambiguous_pairs = set()

    def compound_instances(uninstantiated):
        instances = []
        for compound, value in iteritems(uninstantiated):
            if value > 0:
                f = compound_formula[compound.name]
                instances.append(_CompoundInstance(compound, value, f))

        for inst in instances:
            uninstantiated[inst.compound] -= 1

        return instances

    def instantiate(uninstantiated, compound):
        n = uninstantiated[compound]
        if n > 0:
            f = compound_formula[compound.name]
            inst = _CompoundInstance(compound, n, f)
            uninstantiated[compound] -= 1
            return inst

        return None

    left = compound_instances(uninstantiated_left)
    right = compound_instances(uninstantiated_right)
    instances = left + right

    pairs = {}
    for inst1, inst2 in product(left, right):
        result = score_func(inst1, inst2)
        if result is not None:
            pairs[inst1, inst2] = result

    def inst_pair_sort_key(entry):
        """Sort key for finding best match among instance pairs.

        Rank by score in general but always match identical compounds first
        (these will always have score equal to one but are handled specially
        to put them ahead of other compounds with score equal to one). Use
        compound names to break ties to produce a deterministic result.
        """
        (inst1, inst2), score = entry
        c1, c2 = inst1.compound, inst2.compound
        same_compound = c1.name == c2.name and c1.compartment != c2.compartment
        return same_compound, score, c1.name, c2.name

    transfer = {}
    while len(pairs) > 0:
        (inst1, inst2), _ = max(iteritems(pairs), key=inst_pair_sort_key)
        common = inst1.formula & inst2.formula
        sorted_pairs = sorted(iteritems(pairs),
                              key=inst_pair_sort_key, reverse=True)
        max_sort_key = inst_pair_sort_key(sorted_pairs[0])

        ambiguous_entry = {(inst1.compound, inst2.compound, common)}
        for (other_inst1, other_inst2), other_score in sorted_pairs[1:]:
            other_sort_key = inst_pair_sort_key(
                ((other_inst1, other_inst2), other_score))
            if other_sort_key[:2] < max_sort_key[:2]:
                break

            other_common = other_inst1.formula & other_inst2.formula
            ambiguous_entry.add(
                (other_inst1.compound, other_inst2.compound, other_common))

        if len(ambiguous_entry) > 1:
            ambiguous_pairs.add(tuple(ambiguous_entry))

        key = (inst1.compound, inst1.index), (inst2.compound, inst2.index)
        if key not in transfer:
            transfer[key] = Formula()
        transfer[key] |= common

        for inst in (inst1, inst2):
            inst.formula -= common

        to_insert = set()

        inst = instantiate(uninstantiated_left, inst1.compound)
        if inst is not None:
            left.append(inst)
            instances.append(inst)
            to_insert.add(inst)

        inst = instantiate(uninstantiated_right, inst2.compound)
        if inst is not None:
            right.append(inst)
            instances.append(inst)
            to_insert.add(inst)

        to_update = {inst1, inst2}

        to_delete = set()
        for inst1, inst2 in pairs:
            if inst1 in to_update or inst2 in to_update:
                if len(inst1.formula) > 0 and len(inst2.formula) > 0:
                    result = score_func(inst1, inst2)
                    if result is None:
                        to_delete.add((inst1, inst2))
                    else:
                        pairs[inst1, inst2] = result
                else:
                    to_delete.add((inst1, inst2))

        for pair in to_delete:
            del pairs[pair]

        for inst1, inst2 in product(left, right):
            if inst1 in to_insert or inst2 in to_insert:
                result = score_func(inst1, inst2)
                if result is not None:
                    pairs[inst1, inst2] = result

    balance = {}
    for inst in instances:
        if len(inst.formula) > 0:
            key = inst.compound, inst.index
            balance[key] = inst.formula

    return transfer, balance, ambiguous_pairs


def predict_compound_pairs(reaction, compound_formula, pair_weights={},
                           weight_func=element_weight):
    """Predict compound pairs for a single reaction.

    Performs greedy matching on reaction compounds using a scoring function
    that uses generalized Jaccard similarity corrected by the weights in the
    given dictionary. Returns a tuple of a transfer dictionary and a dictionary
    of unbalanced compounds. The dictionary of unbalanced compounds is empty
    only if the reaction is balanced.

    Args:
        reaction: :class:`psamm.reaction.Reaction`.
        compound_formula: Dictionary mapping compound IDs to
            :class:`psamm.formula.Formula`. Formulas must be flattened.
        pair_weights: Dictionary mapping pairs of compound IDs to correction
            values. This value is multiplied by the calculated Jaccard
            similarity. If a pair is not in the dictionary, the value 1 is
            used. Pairs are looked up in the weights dictionary as a tuple of
            compound names (``c1``, ``c2``) where ``c1`` is the left-hand side
            and ``c2`` is the right-hand side.
        weight_func: Weight function for caclulating the generalized Jaccard
            similarity. This function will be given an
            :class:`psamm.formula.Atom` or :class:`psamm.formula.Radical` and
            should return a corresponding weight.
    """
    def score_func(inst1, inst2):
        score = _jaccard_similarity(
            inst1.formula, inst2.formula, weight_func)
        if score is None:
            return None
        pair = inst1.compound.name, inst2.compound.name
        pair_weight = pair_weights.get(pair, 1.0)
        return pair_weight * score

    return _match_greedily(reaction, compound_formula, score_func)
