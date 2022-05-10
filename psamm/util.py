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

"""Various utilities."""

from __future__ import unicode_literals

import os
import errno
import re
import math
import subprocess
import collections

import six
from six import iteritems, text_type


def mkdir_p(path):
    """Make directory path if it does not already exist."""
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST or not os.path.isdir(path):
            raise


class LoggerFile(object):
    """File-like object that forwards to a logger.

    The Cplex API takes a file-like object for writing log output.
    This class allows us to forward the Cplex messages to the
    Python logging system.
    """

    def __init__(self, logger, level):
        self._logger = logger
        self._level = level

    def write(self, s):
        """Write message to logger."""
        for line in re.split(r'\n+', s):
            if line != '':
                self._logger.log(self._level, line)

    def flush(self):
        """Flush stream.

        This is a noop."""


@six.python_2_unicode_compatible
class MaybeRelative(object):
    """Helper type for parsing possibly relative parameters.

    >>> arg = MaybeRelative('40%')
    >>> arg.reference = 200.0
    >>> float(arg)
    80.0

    >>> arg = MaybeRelative('24.5')
    >>> arg.reference = 150.0
    >>> float(arg)
    24.5
    """

    def __init__(self, s):
        try:
            self._value = float(s)
            self._relative = False
        except ValueError:
            self._value = self._parse_percentage(s)
            self._relative = self._value != 0.0

        self._reference = None

    @classmethod
    def _parse_percentage(cls, s):
        """Parse string as a percentage (e.g. '42.0%') and return as float."""
        m = re.match('^(.+)%$', s)
        if not m:
            raise ValueError('Unable to parse as percentage: {}'.format(
                repr(s)))

        return float(m.group(1)) / 100.0

    @property
    def relative(self):
        """Whether the parsed number was relative."""
        return self._relative

    @property
    def reference(self):
        """The reference used for converting to absolute value."""
        return self._reference

    @reference.setter
    def reference(self, value):
        self._reference = value

    def __float__(self):
        if self._relative:
            if self._reference is None:
                raise ValueError('Reference not set!')
            return self._reference * self._value
        else:
            return self._value

    def __str__(self):
        if self._relative:
            f = None if self._reference is None else float(self)
            return '{:.1%} of {} = {}'.format(
                self._value, self._reference, f)
        else:
            return text_type(self._value)

    def __repr__(self):
        return str('<{}, {}>').format(self.__class__.__name__, str(self))


class FrozenOrderedSet(collections.abc.Set, collections.abc.Hashable):
    """An immutable set that retains insertion order."""

    def __init__(self, seq=[]):
        self.__map = collections.OrderedDict()
        for e in seq:
            self.__map[e] = None

    def __contains__(self, element):
        return element in self.__map

    def __iter__(self):
        return iter(self.__map)

    def __len__(self):
        return len(self.__map)

    def __hash__(self):
        h = 0
        for e in self:
            h ^= 31 * hash(e)
        return h

    def __repr__(self):
        return str('{}({})').format(self.__class__.__name__, list(self))


class DictView(collections.abc.Mapping):
    """An immutable wrapper around another dict-like object."""

    def __init__(self, d):
        self.__d = d

    def __getitem__(self, key):
        return self.__d[key]

    def __iter__(self):
        return iter(self.__d)

    def __len__(self):
        return len(self.__d)


def create_unique_id(prefix, existing_ids):
    """Return a unique string ID from the prefix.

    First check if the prefix is itself a unique ID in the set-like parameter
    existing_ids. If not, try integers in ascending order appended to the
    prefix until a unique ID is found.
    """
    if prefix in existing_ids:
        suffix = 1
        while True:
            new_id = '{}_{}'.format(prefix, suffix)
            if new_id not in existing_ids:
                return new_id
            suffix += 1

    return prefix


def git_try_describe(repo_path):
    """Try to describe the current commit of a Git repository.

    Return a string containing a string with the commit ID and/or a base tag,
    if successful. Otherwise, return None.
    """
    try:
        p = subprocess.Popen(['git', 'describe', '--always', '--dirty'],
                             cwd=repo_path, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output, _ = p.communicate()
    except:
        return None
    else:
        if p.returncode == 0:
            return output.strip()

    return None


def convex_cardinality_relaxed(f, epsilon=1e-5):
    """Transform L1-norm optimization function into cardinality optimization.

    The given function must optimize a convex problem with
    a weighted L1-norm as the objective. The transformed function
    will apply the iterated weighted L1 heuristic to approximately
    optimize the cardinality of the solution. This method is
    described by S. Boyd, "L1-norm norm methods for convex cardinality
    problems." Lecture Notes for EE364b, Stanford University, 2007.
    Available online at www.stanford.edu/class/ee364b/.

    The given function must take an optional keyword parameter weights
    (dictionary), and the weights must be set to one if not specified.
    The function must return the non-weighted solution as an iterator
    over (identifier, value)-tuples, either directly or as the first
    element of a tuple.
    """

    def convex_cardinality_wrapper(*args, **kwargs):
        def dict_result(r):
            if isinstance(r, tuple):
                return dict(r[0])
            return dict(r)

        # Initial run with default weights
        full_result = f(*args, **kwargs)
        result = dict_result(full_result)

        def update_weight(value):
            return 1/(epsilon + abs(value))

        # Iterate until the difference from one iteration to
        # the next is less than epsilon.
        while True:
            weights = {identifier: update_weight(value)
                       for identifier, value in iteritems(result)}
            kwargs['weights'] = weights

            last_result = result
            full_result = f(*args, **kwargs)
            result = dict_result(full_result)

            delta = math.sqrt(sum(pow(value - last_result[identifier], 2)
                                  for identifier, value in iteritems(result)))
            if delta < epsilon:
                break

        if isinstance(full_result, tuple):
            return (iteritems(result),) + full_result[1:]
        return iteritems(result)

    return convex_cardinality_wrapper
