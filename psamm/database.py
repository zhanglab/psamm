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

"""Representation of metabolic network databases."""

import abc
from collections import defaultdict, Mapping

from six import iteritems, add_metaclass

from .reaction import Reaction, Direction


class StoichiometricMatrixView(Mapping):
    """Provides a sparse matrix view on the stoichiometry of a database.

    This object is used internally in the database to expose a sparse matrix
    view of the model stoichiometry. This class should not be instantied,
    instead use the :attr:`MetabolicDatabase.matrix` property. Any compound,
    reaction-pair can be looked up to obtain the corresponding stoichiometric
    value. If the value is not defined (implicitly zero) a
    :exc:`KeyError <exceptions.KeyError>` will be raised.

    In addition, instances also support the NumPy `__array__` protocol which
    means that a :class:`numpy.array` can by created directly from the matrix.

    >>> model = MetabolicModel()
    >>> matrix = numpy.array(model.matrix)
    """

    def __init__(self, database):
        super(StoichiometricMatrixView, self).__init__()
        self._database = database

    def __getitem__(self, key):
        if len(key) != 2:
            raise KeyError(key)
        compound, reaction = key
        if not self._database.has_reaction(reaction):
            raise KeyError(key)
        reaction_values = dict(self._database.get_reaction_values(reaction))
        if compound not in reaction_values:
            raise KeyError(key)
        return reaction_values[compound]

    def __iter__(self):
        for reaction in self._database.reactions:
            for compound, _ in self._database.get_reaction_values(reaction):
                yield compound, reaction

    def __len__(self):
        return sum(sum(1 for _ in self._database.get_reaction_values(reaction))
                   for reaction in self._database.reactions)

    def __array__(self):
        """Return Numpy ndarray instance of matrix.

        The matrix is indexed by sorted compound, reaction-keys
        """
        import numpy  # NumPy is only required for this method

        compound_list = sorted(self._database.compounds)
        reaction_list = sorted(self._database.reactions)
        matrix = numpy.zeros((len(compound_list), len(reaction_list)))

        compound_map = dict(
            (compound, i) for i, compound in enumerate(compound_list))
        reaction_map = dict(
            (reaction_id, i) for i, reaction_id in enumerate(reaction_list))

        for reaction_id in self._database.reactions:
            for compound, value in self._database.get_reaction_values(
                    reaction_id):
                c_index = compound_map[compound]
                r_index = reaction_map[reaction_id]
                matrix[c_index, r_index] = value

        return matrix


@add_metaclass(abc.ABCMeta)
class MetabolicDatabase(object):
    """Database of metabolic reactions."""

    @abc.abstractproperty
    def reactions(self):
        """Iterator of reactions IDs in the database."""

    @abc.abstractproperty
    def compounds(self):
        """Itertor of :class:`Compounds <psamm.reaction.Compound>` in the
        database."""

    @abc.abstractproperty
    def compartments(self):
        """Iterator of compartment IDs in the database."""

    @abc.abstractmethod
    def has_reaction(self, reaction_id):
        """Whether the given reaction exists in the database."""

    @abc.abstractmethod
    def is_reversible(self, reaction_id):
        """Whether the given reaction is reversible."""

    @abc.abstractmethod
    def get_reaction_values(self, reaction_id):
        """Return an iterator of reaction compounds and stoichiometric values.

        The returned iterator contains
        (:class:`Compound <psamm.reaction.Compound>`, value)-tuples. The value
        is negative for left-hand side compounds and positive for right-hand
        side.
        """

    @abc.abstractmethod
    def get_compound_reactions(self, compound_id):
        """Return an iterator of reactions containing the compound.

        Reactions are returned as IDs.
        """

    @property
    def reversible(self):
        """The set of reversible reactions."""
        return set(reaction_id for reaction_id in self.reactions
                   if self.is_reversible(reaction_id))

    @property
    def matrix(self):
        """Mapping from compound, reaction to stoichiometric value.

        This is an instance of :class:`StoichiometricMatrixView`."""
        return StoichiometricMatrixView(self)

    def get_reaction(self, reaction_id):
        """Return reaction as a :class:`Reaction <psamm.reaction.Reaction>`."""

        direction = Direction.Forward
        if self.is_reversible(reaction_id):
            direction = Direction.Both

        return Reaction(direction, self.get_reaction_values(reaction_id))


class DictDatabase(MetabolicDatabase):
    """Metabolic database backed by in-memory dictionaries

    This is a subclass of :class:`MetabolicDatabase`."""

    def __init__(self):
        super(DictDatabase, self).__init__()
        self._model_name = None
        self._biomass_reaction = None
        self._default_flux_limit = None
        self._extracellular_compartment = None
        # _flux_limits is a map from reaction_id to (upper, lower)
        self._flux_limits = {}
        # _medium is a map from reaction_id to (Compound, upper, lower)
        self._medium = {}
        # reaction IDs in the model
        self._model = set()
        # compound map from compound ID to CompoundEntry
        self._compound_entries = defaultdict(dict)
        # reaction map from reaction ID to ReactionEntry
        self._reaction_entries = defaultdict(dict)
        # original attributes in DictDatabase, consider to be moved
        self._reactions = defaultdict(dict)
        self._compound_reactions = defaultdict(set)
        self._reversible = set()

    @property
    def model_name(self):
        return self._model_name

    @property
    def extracelluar_compartment(self):
        return self._extracellular_compartment

    @property
    def defalut_flux_limit(self):
        return self._default_flux_limit

    @property
    def reactions(self):
        return iter(self._reactions)

    @property
    def compound_entries(self):
        return iter(self._compound_entries)

    @property
    def compounds(self):
        return iter(self._compound_reactions)

    @property
    def model(self):
        return iter(self._model)

    # consider to be moved
    @property
    def compartments(self):
        compartment_set = set()
        for compound in self.compounds:
            if compound.compartment not in compartment_set:
                compartment_set.add(compound.compartment)
                yield compound.compartment

    def has_compound(self, compound_id):
        return compound_id in self._compound_entries

    def has_reaction(self, reaction_id):
        return reaction_id in self._reaction_entries

    def is_reversible(self, reaction_id):
        return reaction_id in self._reversible

    def get_reaction_values(self, reaction_id):
        if reaction_id not in self._reactions:
            raise ValueError('Unknown reaction: {}'.format(repr(reaction_id)))
        return iteritems(self._reactions[reaction_id])

    def get_compound_reactions(self, compound):
        return iter(self._compound_reactions[compound])

    @model_name.setter
    def model_name(self, model_name):
        self._model_name = model_name

    @extracelluar_compartment.setter
    def extracellular_compartment(self, extracelluar_compartment):
        self._extracellular_compartment = extracelluar_compartment

    @defalut_flux_limit.setter
    def default_flux_limit(self, flux_limit):
        self._default_flux_limit = flux_limit

    def add_medium(self, reaction_id, compound, lower, upper):
        self._medium[reaction_id] = (compound, lower, upper)

    def add_compound_entry(self, compound_id, compound_entry):
        self._compound_entries[compound_id] = compound_entry

    def add_reaction_entry(self, reaction_id, reaction_entry):
        self._reaction_entries[reaction_id] = reaction_entry

    def add_model_reaction(self, reaction_id):
        self._model.add(reaction_id)

    def add_flux_limit(self, reaction_id, lower, upper):
        self._flux_limits[reaction_id] = (lower, upper)

    def set_reaction(self, reaction_id, reaction):
        """Set the reaction ID to a reaction given by a
        :class:`Reaction <psamm.reaction.Reaction>`

        If an existing reaction exists with the given reaction ID it will be
        overwritten.
        """

        # Overwrite previous reaction if the same id is used
        if reaction_id in self._reactions:
            # Clean up compound to reaction mapping
            for compound in self._reactions[reaction_id]:
                self._compound_reactions[compound].remove(reaction_id)

            self._reversible.discard(reaction_id)
            del self._reactions[reaction_id]

        # Add values to global (sparse) stoichiometric matrix
        # Compounds that occur on both sides will get a stoichiometric
        # value based on the sum of the signed values on each side.
        for compound, _ in reaction.compounds:
            if compound not in self._reactions[reaction_id]:
                self._reactions[reaction_id][compound] = 0
                self._compound_reactions[compound].add(reaction_id)
        for compound, value in reaction.left:
            self._reactions[reaction_id][compound] -= value
        for compound, value in reaction.right:
            self._reactions[reaction_id][compound] += value

        # Remove reaction from compound reactions if the resulting
        # stoichiometric value turned out to be zero.
        zero_compounds = set()
        for compound, value in iteritems(self._reactions[reaction_id]):
            if value == 0:
                zero_compounds.add(compound)

        for compound in zero_compounds:
            del self._reactions[reaction_id][compound]
            self._compound_reactions[compound].remove(reaction_id)

        if reaction.direction != Direction.Forward:
            self._reversible.add(reaction_id)


class ChainedDatabase(MetabolicDatabase):
    """Links a number of databases so they can be treated a single database

    This is a subclass of :class:`MetabolicDatabase`."""

    def __init__(self, *databases):
        self._databases = list(databases)
        if len(self._databases) == 0:
            self._databases.append(DictDatabase())

    def _is_shadowed(self, reaction_id, database):
        """Whether reaction in database is shadowed by another database"""
        for other_database in self._databases:
            if other_database == database:
                break
            if other_database.has_reaction(reaction_id):
                return True
        return False

    @property
    def reactions(self):
        # Make sure that we only yield each reaction once
        reaction_set = set()
        for database in self._databases:
            for reaction in database.reactions:
                if reaction not in reaction_set:
                    reaction_set.add(reaction)
                    yield reaction

    @property
    def compounds(self):
        # Make sure that we only yield each compound once
        compound_set = set()
        for database in self._databases:
            for compound in database.compounds:
                if compound not in compound_set:
                    compound_set.add(compound)
                    yield compound

    @property
    def compartments(self):
        # Make sure that we yield each compartment once
        compartment_set = set()
        for database in self._databases:
            for compartment in database.compartments:
                if compartment not in compartment_set:
                    compartment_set.add(compartment)
                    yield compartment

    def has_reaction(self, reaction_id):
        return any(database.has_reaction(reaction_id)
                   for database in self._databases)

    def is_reversible(self, reaction_id):
        for database in self._databases:
            if database.has_reaction(reaction_id):
                return database.is_reversible(reaction_id)
        raise ValueError('Unknown reaction: {}'.format(reaction_id))

    def get_reaction_values(self, reaction_id):
        for database in self._databases:
            if database.has_reaction(reaction_id):
                return database.get_reaction_values(reaction_id)
        raise ValueError('Unknown reaction: {}'.format(reaction_id))

    def get_compound_reactions(self, compound):
        # Make sure that we only yield each reaction once
        reaction_set = set()
        for database in self._databases:
            for reaction in database.get_compound_reactions(compound):
                if (reaction not in reaction_set and
                        not self._is_shadowed(reaction, database)):
                    reaction_set.add(reaction)
                    yield reaction

    def set_reaction(self, reaction_id, reaction):
        if hasattr(self._databases[0], 'set_reaction'):
            self._databases[0].set_reaction(reaction_id, reaction)
        else:
            raise ValueError('First database is immutable')
