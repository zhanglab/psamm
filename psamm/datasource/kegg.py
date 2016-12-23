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
# Copyright 2014-2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Module related to loading KEGG database files."""

import re
from types import FunctionType
from collections import Mapping

from .context import FileMark
from .entry import (CompoundEntry as BaseCompoundEntry,
                    ReactionEntry as BaseReactionEntry)
from ..reaction import Reaction, Compound, Direction
from ..expression.affine import Expression
from ..util import DictView

from six import iteritems, add_metaclass


class ParseError(Exception):
    """Exception used to signal errors while parsing"""


class KEGGEntry(object):
    """Base class for KEGG entry with raw values from KEGG."""

    def __init__(self, properties, filemark=None):
        self._properties = DictView(dict(properties))
        self._filemark = filemark

    @property
    def properties(self):
        return self._properties

    @property
    def filemark(self):
        return self._filemark


class _MappedEntry(object):
    """Base class for KEGG entry with values mapped to standard keys."""

    def __init__(self, mapper, entry):
        self._entry = entry
        self._properties = mapper(self._entry)

    @property
    def id(self):
        return self._properties['id']

    @property
    def properties(self):
        return self._properties

    @property
    def raw_entry(self):
        return self._entry

    @property
    def filemark(self):
        return self._entry.filemark


class CompoundEntry(_MappedEntry, BaseCompoundEntry):
    """KEGG entry with mapped compound properties."""
    def __init__(self, entry):
        super(CompoundEntry, self).__init__(CompoundMapper, entry)


class ReactionEntry(_MappedEntry, BaseReactionEntry):
    """KEGG entry with mapped reaction properties."""
    def __init__(self, entry):
        super(ReactionEntry, self).__init__(ReactionMapper, entry)


class _MapperMeta(Mapping.__class__):
    """Metaclass for mapper that turns methods into cached keys."""
    def __new__(cls, name, bases, namespace):
        property_accessors = {}
        new_namespace = {}
        for attr_name, attr in iteritems(namespace):
            if (isinstance(attr, FunctionType) and
                    not attr_name.startswith('_')):
                property_accessors[attr_name] = attr
            else:
                new_namespace[attr_name] = attr

        def getitem_props(self, key):
            if key in property_accessors:
                if key not in self._cache:
                    self._cache[key] = property_accessors[key](
                        self, self._entry.properties)
                return self._cache[key]
            raise KeyError('key not found: {!r}'.format(key))

        def iter_props(self):
            return iter(property_accessors)

        def len_props(self):
            return len(property_accessors)

        new_namespace['__getitem__'] = getitem_props
        new_namespace['__iter__'] = iter_props
        new_namespace['__len__'] = len_props

        return super(_MapperMeta, cls).__new__(
            cls, name, bases, new_namespace)

    def __call__(self, entry, *args, **kwargs):
        inst = super(_MapperMeta, self).__call__(*args, **kwargs)
        inst._entry = entry
        inst._cache = {}
        return inst


@add_metaclass(_MapperMeta)
class CompoundMapper(Mapping):
    """Mapper for raw KEGG compound properties to standard properties.

    Public methods are automatically translated into cached properties by the
    metaclass.
    """
    def id(self, raw):
        return raw['entry'][0].split(None, 1)[0]

    def name(self, raw):
        try:
            return next(self._iter_names(raw))
        except StopIteration:
            return None

    def names(self, raw):
        return list(self._iter_names(raw))

    def reactions(self, raw):
        return list(self._iter_reactions(raw))

    def enzymes(self, raw):
        return list(self._iter_enzymes(raw))

    def formula(self, raw):
        if 'formula' in raw:
            return raw['formula'][0]
        return None

    def exact_mass(self, raw):
        if 'exact_mass' in raw:
            return float(raw['exact_mass'][0])
        return None

    def mol_weight(self, raw):
        if 'mol_weight' in raw:
            return float(raw['mol_weight'][0])
        return None

    def pathways(self, raw):
        return list(self._iter_pathways(raw))

    def dblinks(self, raw):
        return list(self._iter_dblinks(raw))

    def comment(self, raw):
        if 'comment' in raw:
            return '\n'.join(raw['comment'])
        return None

    def _iter_names(self, raw):
        if 'name' in raw:
            for line in raw['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    def _iter_reactions(self, raw):
        if 'reaction' in raw:
            for line in raw['reaction']:
                for rxnid in line.split():
                    yield rxnid

    def _iter_enzymes(self, raw):
        if 'enzyme' in raw:
            for line in raw['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    def _iter_pathways(self, raw):
        if 'pathway' in raw:
            for line in raw['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    def _iter_dblinks(self, raw):
        if 'dblinks' in raw:
            for line in raw['dblinks']:
                database, entry = line.split(':', 1)
                yield database.strip(), entry.strip()


@add_metaclass(_MapperMeta)
class ReactionMapper(Mapping):
    """Mapper for raw KEGG reaction properties to standard properties.

    Methods are automatically translated into cached properties by the
    metaclass.
    """
    def id(self, raw):
        return raw['entry'][0].split(None, 1)[0]

    def name(self, raw):
        try:
            return next(self._iter_names(raw))
        except StopIteration:
            return None

    def names(self, raw):
        return list(self._iter_names(raw))

    def definition(self, raw):
        if 'definition' in raw:
            return raw['definition'][0]
        return None

    def equation(self, raw):
        if 'equation' in raw:
            return parse_reaction(raw['equation'][0])
        return None

    def enzymes(self, raw):
        return list(self._iter_enzymes(raw))

    def pathways(self, raw):
        return list(self._iter_pathways(raw))

    def comment(self, raw):
        if 'comment' in raw:
            return '\n'.join(raw['comment'])
        return None

    def rpairs(self, raw):
        return list(self._iter_rpairs(raw))

    def _iter_names(self, raw):
        if 'name' in raw:
            for line in raw['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    def _iter_enzymes(self, raw):
        if 'enzyme' in raw:
            for line in raw['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    def _iter_pathways(self, raw):
        if 'pathway' in raw:
            for line in raw['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    def _iter_rpairs(self, raw):
        if 'rpair' in raw:
            for line in raw['rpair']:
                pair, compounds, rp_type = line.split(None, 2)
                compounds = tuple(compounds.split('_', 1))
                yield pair, compounds, rp_type


def parse_kegg_entries(f, context=None):
    """Iterate over entries in KEGG file."""

    section_id = None
    entry_line = None
    properties = {}
    for lineno, line in enumerate(f):
        if line.strip() == '///':
            # End of entry
            mark = FileMark(context, entry_line, 0)
            yield KEGGEntry(properties, filemark=mark)
            properties = {}
            section_id = None
            entry_line = None
        else:
            if entry_line is None:
                entry_line = lineno

            # Look for beginning of section
            m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
            if m is not None:
                section_id = m.group(1).lower()
                properties[section_id] = [m.group(2)]
            elif section_id is not None:
                properties[section_id].append(line.strip())
            else:
                raise ParseError(
                    'Missing section identifier at line {}'.format(lineno))


def parse_compound_file(f, context=None):
    """Iterate over the compound entries in the given file."""
    for entry in parse_kegg_entries(f, context):
        yield CompoundEntry(entry)


def parse_reaction_file(f, context=None):
    """Iterate over the reaction entries in the given file."""
    for entry in parse_kegg_entries(f, context):
        yield ReactionEntry(entry)


def parse_reaction(s):
    """Parse a KEGG reaction string"""

    def parse_count(s):
        m = re.match(r'^\((.+)\)$', s)
        if m is not None:
            s = m.group(1)

        m = re.match(r'^\d+$', s)
        if m is not None:
            return int(m.group(0))

        return Expression(s)

    def parse_compound(s):
        m = re.match(r'(.+)\((.+)\)', s)
        if m is not None:
            return Compound(m.group(1), arguments=[Expression(m.group(2))])
        return Compound(s)

    def parse_compound_list(s):
        for cpd in s.split(' + '):
            if cpd == '':
                continue

            fields = cpd.strip().split(' ')
            if len(fields) > 2:
                raise ParseError(
                    'Malformed compound specification: {}'.format(cpd))
            if len(fields) == 1:
                count = 1
                compound = parse_compound(fields[0])
            else:
                count = parse_count(fields[0])
                compound = parse_compound(fields[1])

            yield compound, count

    cpd_left, cpd_right = s.split('<=>')
    left = parse_compound_list(cpd_left.strip())
    right = parse_compound_list(cpd_right.strip())

    return Reaction(Direction.Both, left, right)
