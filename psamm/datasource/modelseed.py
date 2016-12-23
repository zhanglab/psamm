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

"""Module related to loading ModelSEED database files."""

import csv
import re

from .context import FileMark
from .entry import CompoundEntry as BaseCompoundEntry


class ParseError(Exception):
    """Exception used to signal errors while parsing"""


def decode_name(s):
    """Decode names in ModelSEED files"""
    # Some names contain XML-like entity codes
    return re.sub(r'&#(\d+);', lambda x: chr(int(x.group(1))), s)


class CompoundEntry(BaseCompoundEntry):
    """Representation of entry in a ModelSEED compound table"""

    def __init__(self, id, names, formula, filemark=None):
        self._id = id
        self._properties = {
            'id': self._id,
            'names': list(names),
            'formula': formula
        }

        # Find shortest name
        # Usually the best name is the shortest one,
        # but in the case of ties, the best one is
        # usually the last one to appear in the list.
        # Since the sort is stable we can obtain this
        # by first reversing the list, then sorting
        # by length.
        names = self._properties['names']
        self._properties['name'] = None
        if len(self._properties['names']) > 0:
            name = sorted(reversed(names), key=lambda x: len(x))[0]
            self._properties['name'] = name

        self._filemark = filemark

    @property
    def id(self):
        return self._id

    @property
    def name(self):
        return self._properties.get('name')

    @property
    def names(self):
        return iter(self._properties.get('names'))

    @property
    def formula(self):
        return self._properties.get('formula')

    @property
    def properties(self):
        return dict(self._properties)

    @property
    def filemark(self):
        return self._filemark


def parse_compound_file(f, context=None):
    """Iterate over the compound entries in the given file"""

    f.readline()  # Skip header
    for lineno, row in enumerate(csv.reader(f, delimiter='\t')):
        compound_id, names, formula = row[:3]
        names = (decode_name(name) for name in names.split(',<br>'))

        # ModelSEED sometimes uses an asterisk and number at
        # the end of formulas. This seems to have a similar
        # meaning as '(...)n'.
        m = re.match(r'^(.*)\*(\d*)$', formula)
        if m is not None:
            if m.group(2) != '':
                formula = '({}){}'.format(m.group(1), m.group(2))
            else:
                formula = '({})n'.format(m.group(1))

        formula = formula.strip()
        if formula == '' or formula == 'noformula':
            formula = None

        mark = FileMark(context, lineno, 0)
        yield CompoundEntry(compound_id, names, formula, filemark=mark)
