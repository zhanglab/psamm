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
# Copyright 2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Representation of compound/reaction entries in models."""

import abc

from six import add_metaclass


@add_metaclass(abc.ABCMeta)
class ModelEntry(object):
    """Abstract model entry."""
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

        Note that the properties are not generally mutable.
        """

    @abc.abstractproperty
    def filemark(self):
        """Position of entry in the source file (or None)."""


class CompoundEntry(ModelEntry):
    """Abstract compound entry."""
    @property
    def formula(self):
        """Chemical formula of compound."""
        return self.properties.get('formula')

    @property
    def charge(self):
        """Compound charge value."""
        return self.properties.get('charge')


class ReactionEntry(ModelEntry):
    """Abstract reaction entry."""
    @property
    def equation(self):
        """Reaction equation."""
        return self.properties.get('equation')

    @property
    def genes(self):
        """Gene association expression."""
        return self.properties.get('genes')


class _BaseDictEntry(ModelEntry):
    """Base class for concrete entries based on dictionary."""
    def __init__(self, abstract_type, properties={}, filemark=None):
        if isinstance(properties, abstract_type):
            self._id = properties.id
            self._properties = dict(properties.properties)
            if filemark is None:
                filemark = properties.filemark
        else:
            if 'id' not in properties:
                raise ValueError('id not defined in properties')
            self._id = properties['id']
            self._properties = dict(properties)

        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def properties(self):
        return self._properties

    @property
    def filemark(self):
        return self._filemark


class DictCompoundEntry(CompoundEntry, _BaseDictEntry):
    """Compound entry backed by dictionary.

    The given properties dictionary must contain a key 'id' with the
    identified.

    Args:
        properties: dict or :class:`CompoundEntry` to construct from.
        filemark: Where the entry was parsed from (optional)
    """
    def __init__(self, *args, **kwargs):
        super(DictCompoundEntry, self).__init__(CompoundEntry, *args, **kwargs)


class DictReactionEntry(ReactionEntry, _BaseDictEntry):
    """Reaction entry backed by dictionary.

    The given properties dictionary must contain a key 'id' with the
    identified.

    Args:
        properties: dict or :class:`ReactionEntry` to construct from.
        filemark: Where the entry was parsed from (optional)
    """
    def __init__(self, *args, **kwargs):
        super(DictReactionEntry, self).__init__(ReactionEntry, *args, **kwargs)
