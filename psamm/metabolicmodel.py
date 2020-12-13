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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Representation of metabolic network models."""

from __future__ import unicode_literals

from collections import Mapping

from .database import MetabolicDatabase, StoichiometricMatrixView
from .reaction import Reaction, Direction
from .util import create_unique_id


def create_exchange_id(existing_ids, compound):
    """Create unique ID for exchange of compound."""
    return create_unique_id(u'EX_{}'.format(compound), existing_ids)


def create_transport_id(existing_ids, compound_1, compound_2):
    """Create unique ID for transport reaction of compounds."""
    return create_unique_id(
        u'TP_{}_{}'.format(compound_1, compound_2), existing_ids)


class FluxBounds(object):
    """Represents lower and upper bounds of flux as a mutable object

    This object is used internally in the model representation. Changing
    the state of the object will change the underlying model parameters.
    Deleting a value will reset that value to the defaults.
    """

    def __init__(self, model, reaction):
        self._model = model
        self._reaction = reaction

    def __iter__(self):
        '''Iterator over lower and upper value'''
        yield self.lower
        yield self.upper

    def _assign_lower(self, value):
        self._model._limits_lower[self._reaction] = value

    def _assign_upper(self, value):
        self._model._limits_upper[self._reaction] = value

    def _assign_both(self, lower, upper):
        self._assign_lower(lower)
        self._assign_upper(upper)

    def _check_bounds(self, lower, upper):
        if lower > upper:
            raise ValueError('Lower bound larger than upper bound')

    @property
    def lower(self):
        '''Lower bound'''
        try:
            return self._model._limits_lower[self._reaction]
        except KeyError:
            if self._model.is_reversible(self._reaction):
                return -self._model._v_max
            else:
                return 0

    @lower.setter
    def lower(self, value):
        self._check_bounds(value, self.upper)
        self._assign_lower(value)

    @lower.deleter
    def lower(self):
        self._model._limits_lower.pop(self._reaction, None)

    @property
    def upper(self):
        '''Upper bound'''
        try:
            return self._model._limits_upper[self._reaction]
        except KeyError:
            return self._model._v_max

    @upper.setter
    def upper(self, value):
        self._check_bounds(self.lower, value)
        self._assign_upper(value)

    @upper.deleter
    def upper(self):
        self._model._limits_upper.pop(self._reaction, None)

    @property
    def bounds(self):
        '''Bounds as a tuple'''
        return self.lower, self.upper

    @bounds.setter
    def bounds(self, value):
        lower, upper = value
        self._check_bounds(lower, upper)
        self._assign_both(lower, upper)

    @bounds.deleter
    def bounds(self):
        del self.lower
        del self.upper

    def __eq__(self, other):
        """Equality test"""
        return (isinstance(other, FluxBounds) and
                self.lower == other.lower and
                self.upper == other.upper)

    def __ne__(self, other):
        """Inequality test"""
        return not self == other

    def __repr__(self):
        return u'{}({!r}, {!r})'.format(
            self.__class__.__name__, self.lower, self.upper)


class LimitsView(Mapping):
    """Provides a view of the flux bounds defined in the model

    This object is used internally in MetabolicModel to
    expose a dictonary view of the FluxBounds associated
    with the model reactions.
    """

    def __init__(self, model):
        super(LimitsView, self).__init__()
        self._model = model

    def _create_bounds(self, reaction):
        return FluxBounds(self._model, reaction)

    def __getitem__(self, key):
        if not self._model.has_reaction(key):
            raise KeyError(key)
        return self._create_bounds(key)

    def __iter__(self):
        return self._model.reactions

    def __len__(self):
        return sum(1 for _ in self._model.reactions)


class MetabolicModel(MetabolicDatabase):
    """Represents a metabolic model containing a set of reactions

    The model contains a list of reactions referencing the reactions
    in the associated database.
    """

    def __init__(self, database, v_max=1000):
        self._database = database
        self._limits_lower = {}
        self._limits_upper = {}

        self._reaction_set = set()
        self._compound_set = set()

        self._v_max = v_max

    @property
    def database(self):
        return self._database

    @property
    def reactions(self):
        return iter(self._reaction_set)

    @property
    def compounds(self):
        return iter(self._compound_set)

    @property
    def compartments(self):
        compartment_set = set()
        for compound in self.compounds:
            if compound.compartment not in compartment_set:
                compartment_set.add(compound.compartment)
                yield compound.compartment

    def has_reaction(self, reaction_id):
        return reaction_id in self._reaction_set

    def has_compound(self, compound_id):
        return compound_id in [str(i) for i in self._compound_set]

    def get_reaction(self, reaction_id):
        if reaction_id not in self._reaction_set:
            raise ValueError(u'Reaction not in model: {}'.format(reaction_id))
        return self._database.get_reaction(reaction_id)

    def get_reaction_values(self, reaction_id):
        """Return stoichiometric values of reaction as a dictionary"""
        if reaction_id not in self._reaction_set:
            raise ValueError(u'Unknown reaction: {}'.format(repr(reaction_id)))
        return self._database.get_reaction_values(reaction_id)

    def get_compound_reactions(self, compound_id):
        """Iterate over all reaction ids the includes the given compound"""
        if compound_id not in self._compound_set:
            raise ValueError(u'Compound not in model: {}'.format(compound_id))

        for reaction_id in self._database.get_compound_reactions(compound_id):
            if reaction_id in self._reaction_set:
                yield reaction_id

    def is_reversible(self, reaction_id):
        """Whether the given reaction is reversible"""
        if reaction_id not in self._reaction_set:
            raise ValueError(u'Reaction not in model: {}'.format(reaction_id))
        return self._database.is_reversible(reaction_id)

    def is_exchange(self, reaction_id):
        """Whether the given reaction is an exchange reaction."""
        reaction = self.get_reaction(reaction_id)
        return (len(reaction.left) == 0) != (len(reaction.right) == 0)

    @property
    def limits(self):
        return LimitsView(self)

    def add_reaction(self, reaction_id):
        """Add reaction to model"""

        if reaction_id in self._reaction_set:
            return

        reaction = self._database.get_reaction(reaction_id)
        self._reaction_set.add(reaction_id)
        for compound, _ in reaction.compounds:
            self._compound_set.add(compound)

    def remove_reaction(self, reaction):
        """Remove reaction from model"""

        if reaction not in self._reaction_set:
            return

        self._reaction_set.remove(reaction)
        self._limits_lower.pop(reaction, None)
        self._limits_upper.pop(reaction, None)

        # Remove compound from compound_set if it is not referenced
        # by any other reactions in the model.
        for compound, value in self._database.get_reaction_values(reaction):
            reactions = frozenset(
                self._database.get_compound_reactions(compound))
            if all(other_reaction not in self._reaction_set
                   for other_reaction in reactions):
                self._compound_set.remove(compound)

    def copy(self):
        """Return copy of model"""

        model = self.__class__(self._database)
        model._limits_lower = dict(self._limits_lower)
        model._limits_upper = dict(self._limits_upper)
        model._reaction_set = set(self._reaction_set)
        model._compound_set = set(self._compound_set)
        return model

    @classmethod
    def load_model(cls, database, reaction_iter=None, exchange=None,
                   limits=None, v_max=None):
        """Get model from reaction name iterator.

        The model will contain all reactions of the iterator.
        """

        model_args = {}
        if v_max is not None:
            model_args['v_max'] = v_max

        model = cls(database, **model_args)
        if reaction_iter is None:
            reaction_iter = iter(database.reactions)
        for reaction_id in reaction_iter:
            model.add_reaction(reaction_id)

        # Apply reaction limits
        if limits is not None:
            for reaction_id, lower, upper in limits:
                if lower is not None:
                    model.limits[reaction_id].lower = lower
                if upper is not None:
                    model.limits[reaction_id].upper = upper

        # TODO: Currently we just create a new exchange reaction in the
        # database and add it to the model. Ideally, we should not modify
        # the database. The exchange reaction could be created on the
        # fly when required.
        if exchange is not None:
            for compound, reaction_id, lower, upper in exchange:
                # Create exchange reaction
                if reaction_id is None:
                    reaction_id = create_exchange_id(
                        model.database.reactions, compound)
                model.database.set_reaction(
                    reaction_id, Reaction(Direction.Both, {compound: -1}))
                model.add_reaction(reaction_id)
                if lower is not None:
                    model.limits[reaction_id].lower = lower
                if upper is not None:
                    model.limits[reaction_id].upper = upper

        return model

    def make_irreversible(self, gene_dict={}, exclude_list=[], lumped_rxns={},
                          all_reversible=False):
        """Creates a new metabolic models with only irreversible reactions.

        This function will find every reversible reaction in the
        model and split it into two reactions with the
        {rxnid}_forward or {rxnid}_reverse as the IDs.

        Args:
            self: A metabolicmodel object
            gene_dict: A dictionary stores reaction ids as keys and
                       corresponding gene association as values, it's
                       required only in GIMME function
            exclude_list: list of reactions to exclude in TMFA simulation
            all_reversible: if True make all reactions in model reversible.

        Values:
            mm_irrev: A new metabolic model with only irreversible reactions
            gene_dict_reversible: A dictionary mapping irreversible reaction
                                  ids to gene associations, used only
                                  in GIMME function
            split_rxns: A list of splitted reactions, each element is a tuple
                        of ({rxnid}_forward, {rxnid}_reverse)
        """
        split_reversible = set()
        mm_irrev = self.copy()
        reversible_gene_dict = {}
        new_lump_rxn_dict = {}
        for rxn in self.reactions:
            upper = self.limits[rxn].upper
            lower = self.limits[rxn].lower
            mm_irrev.limits[rxn].upper = upper
            mm_irrev.limits[rxn].lower = lower

            reaction = mm_irrev.get_reaction(rxn)
            if rxn not in exclude_list:
                r = Reaction(Direction.Forward, reaction.left, reaction.right)
                r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
                r_id = u'{}_forward'.format(rxn)
                r2_id = u'{}_reverse'.format(rxn)
                if reaction.direction == Direction.Forward:
                    if all_reversible is False:
                        reversible_gene_dict[rxn] = gene_dict.get(rxn)
                        continue
                    else:
                        mm_irrev.remove_reaction(rxn)
                        mm_irrev.database.set_reaction(r_id, r)
                        mm_irrev.database.set_reaction(r2_id, r2)
                        mm_irrev.add_reaction(r_id)
                        mm_irrev.add_reaction(r2_id)
                        split_reversible.add((r_id, r2_id))
                        reversible_gene_dict[r_id] = gene_dict.get(rxn)
                        reversible_gene_dict[r2_id] = gene_dict.get(rxn)
                elif reaction.direction == Direction.Both:
                    mm_irrev.remove_reaction(rxn)
                    mm_irrev.database.set_reaction(r_id, r)
                    mm_irrev.database.set_reaction(r2_id, r2)
                    mm_irrev.add_reaction(r_id)
                    mm_irrev.add_reaction(r2_id)
                    split_reversible.add((r_id, r2_id))
                    reversible_gene_dict[r_id] = gene_dict.get(rxn)
                    reversible_gene_dict[r2_id] = gene_dict.get(rxn)
                if lower >= 0:
                    mm_irrev.limits[r_id].upper = upper
                    mm_irrev.limits[r_id].lower = lower
                    mm_irrev.limits[r2_id].upper = 0
                    mm_irrev.limits[r2_id].lower = 0
                elif upper <= 0:
                    mm_irrev.limits[r_id].upper = 0
                    mm_irrev.limits[r_id].lower = 0
                    mm_irrev.limits[r2_id].upper = - lower
                    mm_irrev.limits[r2_id].lower = - upper
                else:
                    mm_irrev.limits[r_id].upper = upper
                    mm_irrev.limits[r_id].lower = 0
                    mm_irrev.limits[r2_id].upper = - lower
                    mm_irrev.limits[r2_id].lower = 0

            elif rxn in lumped_rxns.keys():
                final_sub_rxn_list = []
                sub = lumped_rxns[rxn]
                check = 0
                for (x, y) in sub:
                    rn = mm_irrev.get_reaction(x)
                    if rn.direction != Direction.Both:
                        check += 1
                if reaction.direction == Direction.Forward or check != 0:
                    mm_irrev.limits[rxn].upper = 0
                    mm_irrev.limits[rxn].lower = 0
                    sub_rxn_list = lumped_rxns[rxn]
                    for entry in sub_rxn_list:
                        (subrxn, dir) = entry
                        final_sub_rxn_list.append(subrxn)
                    new_lump_rxn_dict[rxn] = final_sub_rxn_list
                elif reaction.direction == Direction.Both:
                    # split the lump reaction itself
                    r = Reaction(Direction.Forward, reaction.left,
                                 reaction.right)
                    r2 = Reaction(Direction.Forward, reaction.right,
                                  reaction.left)
                    r_id = u'{}_forward'.format(rxn)
                    r2_id = u'{}_reverse'.format(rxn)
                    mm_irrev.remove_reaction(rxn)
                    mm_irrev.database.set_reaction(r_id, r)
                    mm_irrev.database.set_reaction(r2_id, r2)
                    mm_irrev.add_reaction(r_id)
                    mm_irrev.add_reaction(r2_id)
                    split_reversible.add((r_id, r2_id))

                    sub_rxn_list = lumped_rxns[rxn]
                    for_sub_rxn_list = []
                    rev_sub_rxn_list = []
                    mm_irrev.limits[r_id].upper = 0
                    mm_irrev.limits[r_id].lower = 0
                    mm_irrev.limits[r2_id].upper = 0
                    mm_irrev.limits[r2_id].lower = 0
                    for entry in sub_rxn_list:
                        (subrxn, dir) = entry
                        dir = int(dir)
                        # subreaction = mm.get_reaction(subrxn)
                        subreaction = self.get_reaction(subrxn)
                        sub_r1 = Reaction(Direction.Forward,
                                          subreaction.left, subreaction.right)
                        sub_r1_id = u'{}_forward'.format(subrxn)
                        sub_r2 = Reaction(Direction.Forward,
                                          subreaction.right, subreaction.left)
                        sub_r2_id = u'{}_reverse'.format(subrxn)

                        split_reversible.add((sub_r1_id, sub_r2_id))
                        if dir == 1:
                            for_sub_rxn_list.append(sub_r1_id)
                            rev_sub_rxn_list.append(sub_r2_id)
                            mm_irrev.database.set_reaction(sub_r1_id, sub_r1)
                            mm_irrev.add_reaction(sub_r1_id)
                            mm_irrev.database.set_reaction(sub_r2_id, sub_r2)
                            mm_irrev.add_reaction(sub_r2_id)
                        elif dir == -1:
                            for_sub_rxn_list.append(sub_r2_id)
                            rev_sub_rxn_list.append(sub_r1_id)
                            mm_irrev.database.set_reaction(sub_r1_id, sub_r1)
                            mm_irrev.add_reaction(sub_r1_id)
                            mm_irrev.database.set_reaction(sub_r2_id, sub_r2)
                            mm_irrev.add_reaction(sub_r2_id)
                        mm_irrev.limits[sub_r1_id].lower = 0
                        mm_irrev.limits[sub_r1_id].upper = 100
                        mm_irrev.limits[sub_r2_id].lower = 0
                        mm_irrev.limits[sub_r2_id].upper = 100
                    new_lump_rxn_dict[r_id] = for_sub_rxn_list
                    new_lump_rxn_dict[r2_id] = rev_sub_rxn_list

        return \
            mm_irrev, reversible_gene_dict, split_reversible, new_lump_rxn_dict


class FlipableFluxBounds(FluxBounds):
    """FluxBounds object for a FlipableModelView.

    This object is used internally in the FlipableModelView to represent
    the bounds of flux on a reaction that can be flipped.
    """

    def __init__(self, view, reaction):
        super(FlipableFluxBounds, self).__init__(view._model, reaction)
        self._view = view

    def _assign_lower(self, value):
        if self._reaction in self._view._flipped:
            super(FlipableFluxBounds, self)._assign_upper(-value)
        else:
            super(FlipableFluxBounds, self)._assign_lower(value)

    def _assign_upper(self, value):
        if self._reaction in self._view._flipped:
            super(FlipableFluxBounds, self)._assign_lower(-value)
        else:
            super(FlipableFluxBounds, self)._assign_upper(value)

    @property
    def lower(self):
        """Lower bound"""
        if self._reaction in self._view._flipped:
            return -super(FlipableFluxBounds, self).upper
        return super(FlipableFluxBounds, self).lower

    @property
    def upper(self):
        """Upper bound"""
        if self._reaction in self._view._flipped:
            return -super(FlipableFluxBounds, self).lower
        return super(FlipableFluxBounds, self).upper


class FlipableStoichiometricMatrixView(StoichiometricMatrixView):
    """Provides a matrix view that flips with the flipable model view.

    This object is used internally in FlipableModelView to
    expose a matrix view that negates the stoichiometric
    values of flipped reactions.
    """

    def __init__(self, view):
        super(FlipableStoichiometricMatrixView, self).__init__(view._model)
        self._view = view

    def _value_mul(self, reaction):
        return -1 if reaction in self._view._flipped else 1

    def __getitem__(self, key):
        if len(key) != 2:
            raise KeyError(key)
        compound, reaction = key
        orig_value = (
            super(FlipableStoichiometricMatrixView, self).__getitem__(key))
        return self._value_mul(reaction) * orig_value


class FlipableLimitsView(LimitsView):
    """Provides a limits view that flips with the flipable model view.

    This object is used internally in FlipableModelView to
    expose a limits view that flips the bounds of all flipped
    reactions.
    """

    def __init__(self, view):
        super(FlipableLimitsView, self).__init__(view._model)
        self._view = view

    def _create_bounds(self, reaction):
        return FlipableFluxBounds(self._view, reaction)


class FlipableModelView(object):
    """Proxy wrapper of model objects allowing a flipped set of reactions.

    The proxy will forward all properties normally except
    that flipped reactions will appear to have stoichiometric
    values negated in the matrix property, and have bounds in
    the limits property flipped. This view is needed for
    some algorithms.
    """

    def __init__(self, model, flipped=set()):
        self._model = model
        self._flipped = set(flipped)

    @property
    def matrix(self):
        return FlipableStoichiometricMatrixView(self)

    @property
    def limits(self):
        return FlipableLimitsView(self)

    def flip(self, subset):
        self._flipped ^= subset

    def __getattr__(self, name):
        return getattr(self._model, name)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
