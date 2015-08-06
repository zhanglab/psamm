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

"""Module related to loading KEGG database files."""

import re

from .context import FileMark
from ..reaction import Reaction, Compound
from ..expression.affine import Expression


class ParseError(Exception):
    """Exception used to signal errors while parsing"""


class CompoundEntry(object):
    """Representation of entry in KEGG compound file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing compound identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def reactions(self):
        if 'reaction' in self.values:
            for line in self.values['reaction']:
                for rxnid in line.split():
                    yield rxnid

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def formula(self):
        if 'formula' not in self.values:
            return None
        return self.values['formula'][0]

    @property
    def exact_mass(self):
        if 'exact_mass' not in self.values:
            return None
        return float(self.values['exact_mass'][0])

    @property
    def mol_weight(self):
        if 'mol_weight' not in self.values:
            return None
        return float(self.values['mol_weight'][0])

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def dblinks(self):
        if 'dblinks' in self.values:
            for line in self.values['dblinks']:
                database, entry = line.split(':', 1)
                yield database.strip(), entry.strip()

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<CompoundEntry "{}">'.format(self.id)


class ReactionEntry(object):
    """Representation of entry in KEGG reaction file"""

    def __init__(self, values, filemark=None):
        self.values = dict(values)
        if 'entry' not in values:
            raise ParseError('Missing reaction identifier')
        self._id, _ = values['entry'][0].split(None, 1)
        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        try:
            return next(self.names)
        except StopIteration:
            return None

    @property
    def names(self):
        if 'name' in self.values:
            for line in self.values['name']:
                for name in line.rstrip(';').split(';'):
                    yield name.strip()

    @property
    def definition(self):
        if 'definition' not in self.values:
            return None
        return self.values['definition'][0]

    @property
    def equation(self):
        if 'equation' not in self.values:
            return None
        return self.values['equation'][0]

    @property
    def enzymes(self):
        if 'enzyme' in self.values:
            for line in self.values['enzyme']:
                for enzyme in line.split():
                    yield enzyme

    @property
    def pathways(self):
        if 'pathway' in self.values:
            for line in self.values['pathway']:
                pathway, name = line.split(None, 1)
                yield pathway, name

    @property
    def comment(self):
        if 'comment' not in self.values:
            return None
        return '\n'.join(self.values['comment'])

    @property
    def rpairs(self):
        if 'rpair' in self.values:
            for line in self.values['rpair']:
                pair, compounds, rp_type = line.split(None, 2)
                compounds = tuple(compounds.split('_', 1))
                yield pair, compounds, rp_type

    def __getitem__(self, name):
        if name not in self.values:
            raise AttributeError('Attribute does not exist: {}'.format(name))
        return self._values[name]

    @property
    def filemark(self):
        return self._filemark

    def __repr__(self):
        return '<ReactionEntry "{}">'.format(self.id)


def _parse_kegg_entries(f, entry_class, context=None):
    """Iterate over entries in KEGG file."""

    section_id = None
    entry_line = None
    compound = {}
    for lineno, line in enumerate(f):
        if line.strip() == '///':
            # End of compound
            mark = FileMark(context, entry_line, 0)
            yield entry_class(compound, filemark=mark)
            compound = {}
            section_id = None
            entry_line = None
        else:
            if entry_line is None:
                entry_line = lineno

            # Look for beginning of section
            m = re.match(r'([A-Z_]+)\s+(.*)', line.rstrip())
            if m is not None:
                section_id = m.group(1).lower()
                compound[section_id] = [m.group(2)]
            elif section_id is not None:
                compound[section_id].append(line.strip())
            else:
                raise ParseError(
                    'Missing section identifier at line {}'.format(lineno))


def parse_compound_file(f, context=None):
    """Iterate over the compound entries in the given file."""
    return _parse_kegg_entries(f, CompoundEntry, context)


def parse_reaction_file(f, context=None):
    """Iterate over the reaction entries in the given file."""
    return _parse_kegg_entries(f, ReactionEntry, context)


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

    return Reaction('<=>', left, right)
