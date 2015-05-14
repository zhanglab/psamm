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

"""Various utilities"""

import re


class LoggerFile(object):
    """File-like object that forwards to a logger

    The Cplex API takes a file-like object for writing log output.
    This class allows us to forward the Cplex messages to the
    Python logging system.
    """

    def __init__(self, logger, level):
        self._logger = logger
        self._level = level

    def write(self, s):
        """Write message to logger"""
        for line in re.split(r'\n+', s):
            if line != '':
                self._logger.log(self._level, line)

    def flush(self):
        """Flush stream

        This is a noop."""


def convex_cardinality_relaxed(f, epsilon=1e-5):
    """Transform L1-norm optimization function into approximate cardinality optimization

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
            weights = { identifier: update_weight(value) for identifier, value in result.iteritems() }
            kwargs['weights'] = weights

            last_result = result
            full_result = f(*args, **kwargs)
            result = dict_result(full_result)

            delta = math.sqrt(sum(pow(value - last_result[identifier], 2) for identifier, value in result.iteritems()))
            if delta < epsilon:
                break

        if isinstance(full_result, tuple):
            return (result.iteritems(),) + full_result[1:]
        return result.iteritems()

    return convex_cardinality_wrapper
