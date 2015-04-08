
"""Generic interface to LP solver instantiation"""

from __future__ import absolute_import

import operator
import logging

from .lp import Solver as BaseSolver

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

        requirements = {key: value for key, value in kwargs.iteritems()
                        if value is not None}
        for req, value in requirements.iteritems():
            solvers = [s for s in solvers if req in s and s[req] == value]

        if len(solvers) == 0:
            raise RequirementsError('Unable to find a solver matching the'
                                    ' specified requirements')

        solver = max(solvers, key=operator.itemgetter('priority'))
        logger.debug('Using solver {}'.format(solver['name']))

        self._solver = solver['class']()

    def create_problem(self):
        """Create a :class:`Problem <metnet.lpsolver.lp.Problem>` instance"""
        return self._solver.create_problem()
