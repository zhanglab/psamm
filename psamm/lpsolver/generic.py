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

from __future__ import absolute_import, print_function

import argparse
import os
import logging

from .lp import Solver as BaseSolver

from six import iteritems

logger = logging.getLogger(__name__)
_solvers = []
_solver_import_errors = {}


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
except ImportError as e:
    _solver_import_errors['cplex'] = str(e)

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
except ImportError as e:
    _solver_import_errors['qsoptex'] = str(e)

# Try to load Gurobi solver
try:
    from . import gurobi
    _solvers.append({
        'class': gurobi.Solver,
        'name': 'gurobi',
        'integer': True,
        'rational': False,
        'priority': 9
    })
except ImportError as e:
    _solver_import_errors['gurobi'] = str(e)


class RequirementsError(Exception):
    """Error resolving solver requirements"""


def filter_solvers(solvers, requirements):
    """Yield solvers that fullfil the requirements."""
    for solver in solvers:
        for req, value in iteritems(requirements):
            if (req in ('integer', 'rational', 'name') and
                    (req not in solver or solver[req] != value)):
                break
        else:
            yield solver


class Solver(BaseSolver):
    """Generic solver interface based on requirements

    Use the any of the following keyword arguments to restrict which
    underlying solver is used:

    - `integer`: Use a solver that support integer variables (MILP)
    - `rational`: Use a solver that returns rational results
    - `name`: Select a specific solver based on the name
    """

    def __init__(self, **kwargs):
        if len(_solvers) == 0:
            raise RequirementsError('No solvers available')

        self._requirements = {key: value for key, value in iteritems(kwargs)
                              if value is not None}
        solvers = list(filter_solvers(_solvers, self._requirements))

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


def list_solvers():
    """Entry point for listing available solvers."""
    parser = argparse.ArgumentParser(
        description='''List LP solver available in PSAMM. This will produce a
                       list of all of the available LP solvers in prioritized
                       order. Addtional requirements can be imposed with the
                       arguments (e.g. integer=yes to select only solvers that
                       support MILP problems). The list will also be influenced
                       by the PSAMM_SOLVER environment variable which can be
                       used to only allow specific solvers (e.g.
                       PSAMM_SOLVER=cplex).''')
    parser.add_argument(
        'requirement', nargs='*', type=str,
        help='Additional requirements on the selected solvers')
    args = parser.parse_args()

    requirements = {}
    for arg in args.requirement:
        try:
            key, value = parse_solver_setting(arg)
        except ValueError as e:
            parser.error(str(e))
        else:
            requirements[key] = value

    solvers = list(filter_solvers(_solvers, requirements))
    solver_names = set(solver['name'] for solver in solvers)

    # Obtain solver priority from environment variable, if specified.
    priority = {}
    if 'PSAMM_SOLVER' in os.environ:
        names = os.environ['PSAMM_SOLVER'].split(',')
        for i, solver_name in enumerate(names):
            priority[solver_name] = len(names) - i
        solvers = [s for s in solvers if s['name'] in priority]
        solver_names = set(priority)
    else:
        # Use built-in priorities
        for solver in solvers:
            priority[solver['name']] = solver['priority']

    solvers = sorted(solvers, key=lambda s: priority.get(s['name'], 0),
                     reverse=True)

    status = 0

    if len(solvers) > 0:
        print('Prioritized solvers:')
        for solver in solvers:
            print('Name: {}'.format(solver['name']))
            print('Priority: {}'.format(solver['priority']))
            print('MILP (integer) problem support: {}'.format(
                solver['integer']))
            print('Rational solution: {}'.format(solver['rational']))
            print('Class: {}'.format(solver['class']))
            print()
    else:
        status = 1
        print('No solvers fullfil the requirements!')
        print()

    filtered_solvers_count = len(_solvers) - len(solvers)
    if filtered_solvers_count > 0 or len(_solver_import_errors) > 0:
        print('Unavailable solvers:')
        for solver in _solvers:
            if solver['name'] not in solver_names:
                print('{}: Does not fullfil the specified requirements'.format(
                    solver['name']))

        for solver, error in iteritems(_solver_import_errors):
            print('{}: Error loading solver: {}'.format(solver, error))

    if status != 0:
        parser.exit(status)
