
'''Wrapper around Cplex solver with better control of log output'''

import sys
import cplex

class CplexSolver(object):
    def __init__(self, stream=sys.stderr):
        self._stream = stream

    def create_problem(self):
        prob = cplex.Cplex()
        prob.set_results_stream(self._stream)
        prob.set_warning_stream(self._stream)
        prob.set_error_stream(self._stream)
        prob.set_log_stream(self._stream)
        return prob
