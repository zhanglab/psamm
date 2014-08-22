
'''Wrapper around Cplex solver with better control of log output'''

import sys
import math
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


def convex_cardinality_relaxed(f, epsilon=1e-5):
    '''Transform L1-norm optimization function into approximate cardinality optimization

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
    element of a tuple.'''

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