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
# Copyright 2020-2020  Elysha Sameth <esameth1@my.uri.edu>

"""Representation of compound/reaction entries in models."""

from __future__ import unicode_literals

import abc
from collections.abc import Mapping

from six import add_metaclass


@add_metaclass(abc.ABCMeta)
class ModelEntry(object):
    """Abstract model entry.

    Provdides a base class for model entries which are representations of
    any entity (such as compound, reaction or compartment) in a model. An
    entity has an ID, and may have a name and filemark. The ID is a unique
    string identified within a model. The name is a string identifier for
    human consumption. The filemark indicates where the entry originates from
    (e.g. file name and line number). Any additional properties for an
    entity exist in ``properties`` which is any dict-like object mapping
    from string keys to any value type. The ``name`` entry in the dictionary
    corresponds to the name. Entries can be mutable, where the
    properties can be modified, or immutable, where the properties cannot be
    modified or where modifications are ignored. The ID is always immutable.
    """
    @abc.abstractproperty
    def id(self):
        """Identifier of entry."""

    @property
    def name(self):
        """Name of entry (or None)."""
        return self.properties.get('name')

    @abc.abstractproperty
    def properties(self):
        """Properties of entry as a :class:`Mapping` subclass (e.g. dict).

        Note that the properties are not generally mutable but may be mutable
        for specific subclasses. If the ``id`` exists in this dictionary, it
        must never change the actual entry ID as obtained from the ``id``
        property, even if other properties are mutable.
        """

    @abc.abstractproperty
    def filemark(self):
        """Position of entry in the source file (or None)."""

    def __repr__(self):
        return str('<{} id={!r}>').format(self.__class__.__name__, self.id)


class CompoundEntry(ModelEntry):
    """Abstract compound entry.

    Entry subclass for representing compounds. This standardizes the properties
    ``formula`` and ``charge``.
    """
    @property
    def formula(self):
        """Chemical formula of compound."""
        return self.properties.get('formula')

    @property
    def charge(self):
        """Compound charge value."""
        return self.properties.get('charge')


class ReactionEntry(ModelEntry):
    """Abstract reaction entry.

    Entry subclass for representing compounds. This standardizes the properties
    ``equation`` and ``genes``.
    """
    @property
    def equation(self):
        """Reaction equation."""
        return self.properties.get('equation')

    @property
    def genes(self):
        """Gene association expression."""
        return self.properties.get('genes')


class CompartmentEntry(ModelEntry):
    """Abstract compartment entry.

    Entry subclass for representing compartments.
    """


class _BaseDictEntry(ModelEntry):
    """Base class for concrete entries based on dictionary.

    The properties are mutable for this subclass. If ``id`` is None, the
    value corresponding to the ``id`` key in the dictionary is used. If this is
    not defined, a :class:`ValueError` is raised.
    """
    def __init__(self, abstract_type, properties={}, filemark=None, id=None):
        if isinstance(properties, abstract_type):
            self._id = id
            if self._id is None:
                self._id = properties.id
            self._properties = dict(properties.properties)
            if filemark is None:
                filemark = properties.filemark
        elif isinstance(properties, Mapping):
            self._id = id
            if self._id is None:
                if 'id' not in properties:
                    raise ValueError('id not defined in properties')
                self._id = properties['id']
            self._properties = dict(properties)
        else:
            raise ValueError('Invalid type of properties object')

        self._properties['id'] = self._id
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @ModelEntry.name.setter
    def name(self, value):
        self._properties['name'] = value

    @property
    def properties(self):
        return self._properties

    @property
    def filemark(self):
        return self._filemark


class DictCompoundEntry(CompoundEntry, _BaseDictEntry):
    """Compound entry backed by dictionary.

    The given properties dictionary must contain a key ``id`` with the
    identifier.

    Args:
        properties: dict or :class:`CompoundEntry` to construct from.
        filemark: Where the entry was parsed from (optional)
    """
    def __init__(self, *args, **kwargs):
        super(DictCompoundEntry, self).__init__(CompoundEntry, *args, **kwargs)

    @CompoundEntry.formula.setter
    def formula(self, value):
        self._properties['formula'] = value

    @CompoundEntry.charge.setter
    def charge(self, value):
        self._properties['charge'] = value


class DictReactionEntry(ReactionEntry, _BaseDictEntry):
    """Reaction entry backed by dictionary.

    The given properties dictionary must contain a key ``id`` with the
    identifier.

    Args:
        properties: dict or :class:`ReactionEntry` to construct from.
        filemark: Where the entry was parsed from (optional)
    """
    def __init__(self, *args, **kwargs):
        super(DictReactionEntry, self).__init__(ReactionEntry, *args, **kwargs)

    @ReactionEntry.equation.setter
    def equation(self, value):
        self._properties['equation'] = value

    @ReactionEntry.genes.setter
    def genes(self, value):
        self._properties['genes'] = value


class DictCompartmentEntry(CompartmentEntry, _BaseDictEntry):
    """Compartment entry backed by dictionary.

    The given properties dictionary must contain a key ``id`` with the
    identifier.

    Args:
        properties: dict or :class:`CompartmentEntry` to construct from.
        filemark: Where the entry was parsed from (optional)
    """
    def __init__(self, *args, **kwargs):
        super(DictCompartmentEntry, self).__init__(
            CompartmentEntry, *args, **kwargs)
