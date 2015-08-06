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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Utilities for keeping track of parsing context."""


import os

from six import string_types


class FilePathContext(object):
    """File context that keeps track of contextual information.

    When a file is loaded, all files specified in that file must be loaded
    relative to the first file. This is made possible by keeping a context
    that remembers where a file was loaded so that other files can be loaded
    relatively.
    """

    def __init__(self, arg):
        """Create new context from a path or existing context"""

        if isinstance(arg, string_types):
            self._filepath = arg
        else:
            self._filepath = arg.filepath
        self._basepath = os.path.dirname(self._filepath)

    @property
    def filepath(self):
        return self._filepath

    @property
    def basepath(self):
        return self._basepath

    def resolve(self, relpath):
        return FilePathContext(os.path.join(self._basepath, relpath))

    def open(self, mode='r'):
        return open(self.filepath, mode)

    def __str__(self):
        return self._filepath


class FileMark(object):
    """Marks a position in a file.

    This is used when parsing input files, to keep track of the position that
    generates an entry.
    """

    def __init__(self, filecontext, line, column):
        self._filecontext = filecontext
        self._line = line
        self._column = column

    @property
    def filecontext(self):
        return self._filecontext

    @property
    def line(self):
        return self._line

    @property
    def column(self):
        return self._column

    def __str__(self):
        return '{}:{}:{}'.format(
            str(self._filecontext), '?' if self._line is None else self._line,
            '?' if self._column is None else self._column)

    def __repr__(self):
        return '{}({}, {}, {})'.format(
            self.__class__.__name__, repr(self._filecontext), repr(self._line),
            repr(self._column))
