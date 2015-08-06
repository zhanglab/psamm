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

"""Generic interface to LP solver instantiation."""

from __future__ import absolute_import

import os
import logging

from .lp import Solver as BaseSolver

from six import iteritems

logger = logging.getLogger(__name__)
_solvers = []


# Try to load Cplex solver
try:
    from . import cplex
    _solvers.append({
        'class': cplex.Solver,
        'name': 'cplex',
        'integer': True,
        'rational': False,
        'priority': 10
    })
except ImportError:
    pass

# Try to load QSopt_ex solver
try:
    from . import qsoptex
    _solvers.append({
        'class': qsoptex.Solver,
        'name': 'qsoptex',
        'integer': False,
        'rational': True,
        'priority': 5
    })
except ImportError:
    pass


class RequirementsError(Exception):
    """Error resolving solver requirements"""


class Solver(BaseSolver):
    """Generic solver interface based on requirements

    Use the any of the following keyword arguments to restrict which
    underlying solver is used:

    - `integer`: Use a solver that support integer variables (MILP)
    - `rational`: Use a solver that returns rational results
    - `name`: Select a specific solver based on the name
    """

    def __init__(self, **kwargs):
        solvers = _solvers
        if len(solvers) == 0:
            raise RequirementsError('No solvers available')

        self._requirements = {key: value for key, value in iteritems(kwargs)
                              if value is not None}
        for req, value in iteritems(self._requirements):
            if req in ('integer', 'rational', 'name'):
                solvers = [s for s in solvers if req in s and s[req] == value]

        # Obtain solver priority from environment variable, if specified.
        priority = {}
        if 'PSAMM_SOLVER' in os.environ:
            solver_names = os.environ['PSAMM_SOLVER'].split(',')
            for i, solver_name in enumerate(solver_names):
                priority[solver_name] = len(solver_names) - i
            solvers = [s for s in solvers if s['name'] in priority]
        else:
            # Use built-in priorities
            for solver in solvers:
                priority[solver['name']] = solver['priority']

        if len(solvers) == 0:
            raise RequirementsError(
                'Unable to find a solver matching the specified requirements:'
                ' {}'.format(self._requirements))

        solver = max(solvers, key=lambda s: priority.get(s['name'], 0))
        logger.debug('Using solver {}'.format(solver['name']))

        self._solver = solver['class']()

    def create_problem(self):
        """Create a :class:`Problem <psamm.lpsolver.lp.Problem>` instance"""
        return self._solver.create_problem(**self._requirements)


def parse_solver_setting(s):
    """Parse a string containing a solver setting"""

    try:
        key, value = s.split('=', 1)
    except ValueError:
        key, value = s, 'yes'

    if key in ('rational', 'integer'):
        value = value.lower() in ('1', 'yes', 'true', 'on')
    elif key in ('threads',):
        value = int(value)
    elif key in ('feasibility_tolerance',):
        value = float(value)

    return key, value
