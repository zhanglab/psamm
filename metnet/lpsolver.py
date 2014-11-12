
'''Wrapper around Cplex solver'''

import sys
import math
import numbers
from itertools import repeat, count, izip
from collections import Counter

import cplex as cp

class CplexSolver(object):
    '''Represents an LP-solver using Cplex'''
    def __init__(self, stream=sys.stderr):
        self._stream = stream

    def create_problem(self):
        '''Create a new LP-problem using the solver'''
        return CplexProblem(stream=self._stream)

class CplexProblem(object):
    '''Represents an LP-problem of a CplexSolver'''

    Maximize = 0
    Minimize = 1

    Continuous = 'C'
    Binary = 'B'
    Integer = 'I'

    def __init__(self, stream=sys.stderr):
        self._cp = cp.Cplex()
        self._cp.set_results_stream(stream)
        self._cp.set_warning_stream(stream)
        self._cp.set_error_stream(stream)
        self._cp.set_log_stream(stream)

        self._variables = {}
        self._var_names = ('x'+str(i) for i in count(1))

        self._result = None

    @property
    def cplex(self):
        '''The underlying Cplex object'''
        return self._cp

    def define(self, *names, **kwargs):
        '''Define variable in the problem

        Variables must be defined before they can be accessed by var() or set().
        This function takes keyword arguments lower and upper to define the
        bounds of the variable (default: -inf to inf, where inf is a very large number
        defined by Cplex). The keyword argument types can be used to select the type of
        the variable (Continuous (default), Binary or Integer). Setting any variables
        different than Continuous will turn the problem into an MILP problem.'''

        names = tuple(names)
        lower = kwargs.get('lower', None)
        upper = kwargs.get('upper', None)
        vartype = kwargs.get('types', None)

        # Repeat values if a scalar is given
        if lower is None or isinstance(lower, numbers.Number):
            lower = repeat(lower, len(names))
        if upper is None or isinstance(upper, numbers.Number):
            upper = repeat(upper, len(names))
        if vartype is None or vartype in (CplexProblem.Continuous, CplexProblem.Binary, CplexProblem.Integer):
            vartype = repeat(vartype, len(names))

        lp_names = tuple(next(self._var_names) for name in names)

        # Assign default values
        lower = (-cp.infinity if value is None else value for value in lower)
        upper = (cp.infinity if value is None else value for value in upper)
        vartype = tuple(CplexProblem.Continuous if value is None else value for value in vartype)

        args = { 'names': lp_names, 'lb': tuple(lower), 'ub': tuple(upper) }
        if any(value != CplexProblem.Continuous for value in vartype):
            # Set types only if some are integer (otherwise Cplex will change
            # the solver to MILP).
            args['types'] = vartype

        self._variables.update(izip(names, lp_names))
        self._cp.variables.add(**args)

    def var(self, name):
        '''Return the variable as an expression'''
        if name not in self._variables:
            raise ValueError('Undefined variable: {}'.format(name))
        return Expression({ name: 1 })

    def set(self, names):
        '''Return the set of variables as an expression'''
        names = tuple(names)
        if any(name not in self._variables for name in names):
            raise ValueError('Undefined variables: {}'.format(set(names) - set(self._variables)))
        return Expression({ VariableSet(names): 1 })

    def add_linear_constraints(self, *relations):
        '''Add constraints to the problem

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.'''
        for relation in relations:
            if isinstance(relation, bool):
                # A bool in place of a relation is accepted to mean
                # a relation that does not involve any variables and
                # has therefore been evaluated to a truth-value (e.g
                # '0 == 0' or '2 >= 3').
                if not relation:
                    raise ValueError('Unsatisfiable relation added')
            else:
                if relation.sense in (Relation.StrictlyGreater, Relation.StrictlyLess):
                    raise ValueError('Strict relations are invalid in LP-problems: {}'.format(relation))

                expression = relation.expression
                pairs = []
                for value_set in expression.value_sets():
                    ind, val = zip(*((self._variables[variable], float(value)) for variable, value in value_set))
                    pairs.append(cp.SparsePair(ind=ind, val=val))
                self._cp.linear_constraints.add(lin_expr=pairs, senses=tuple(repeat(relation.sense, len(pairs))),
                                                rhs=tuple(repeat(float(-expression.offset), len(pairs))))

    def set_linear_objective(self, expression):
        '''Set linear objective of problem'''

        if isinstance(expression, numbers.Number):
            # Allow expressions with no variables as objective,
            # represented as a number
            expression = Expression()

        self._cp.objective.set_linear((lp_name, expression.value(var)) for var, lp_name in self._variables.iteritems())

    def set_objective_sense(self, sense):
        '''Set type of problem (maximize or minimize)'''
        if sense == CplexProblem.Minimize:
            self._cp.objective.set_sense(self._cp.objective.sense.minimize)
        elif sense == CplexProblem.Maximize:
            self._cp.objective.set_sense(self._cp.objective.sense.maximize)
        else:
            raise ValueError('Invalid objective sense')

    def solve(self, sense=None):
        '''Solve problem'''
        if sense is not None:
            self.set_objective_sense(sense)
        self._cp.solve()

        self._result = CplexResult(self)
        return self._result

    @property
    def result(self):
        return self._result

class CplexResult(object):
    '''Represents the solution to a CplexProblem

    This object will be returned from the CplexProblem.solve() method or by
    accessing the CplexProblem.result property after solving a problem. This
    class should not be instantiated manually.

    CplexResult will evaluate to a boolean according to the success of the
    solution, so checking the truth value of the result will immediately
    indicate whether solving was successful.'''

    def __init__(self, prob):
        self._problem = prob

    def _check_valid(self):
        if self._problem.result != self:
            raise Exception('Previous result is no longer valid as problem has been solved again')

    @property
    def success(self):
        '''Return boolean indicating whether a solution was found'''
        self._check_valid()
        return self._problem._cp.solution.get_status() in (1, 101)

    def __bool__(self):
        return self.success

    @property
    def status(self):
        '''Return string indicating the error encountered on failure'''
        self._check_valid()
        return self._problem._cp.solution.get_status_string()

    def get_value(self, expression):
        '''Return value of expression'''

        self._check_valid()
        if isinstance(expression, Expression):
            return sum(self._problem._cp.solution.get_values(self._problem._variables[var])*value for var, value in expression.values())
        elif expression not in self._problem._variables:
            raise ValueError('Unknown expression: {}'.format(expression))
        return self._problem._cp.solution.get_values(self._problem._variables[expression])

class VariableSet(tuple):
    '''A tuple used to represent sets of variables'''

class Expression(object):
    '''Represents a linear expression

    The variables can be any hashable objects. If one or more variables
    are instead VariableSets, this will be taken to represent a set
    of expressions separately using a different element of the
    VariableSet.'''

    def __init__(self, variables={}, offset=0):
        self._variables = Counter(variables)
        self._offset = offset

    @property
    def offset(self):
        '''Value of the offset'''
        return self._offset

    def variables(self):
        '''Iterator of variables in expression'''
        return self._variables.iterkeys()

    def values(self):
        '''Iterator of variable, value-pairs in expression'''
        return self._variables.iteritems()

    def value(self, variable):
        return self._variables.get(variable, 0)

    def value_sets(self):
        '''Iterator of expression sets

        This will yield an iterator of variable, value-pairs for
        each expression in the expression set (each equivalent to
        values()). If none of the variables is a set variable then
        a single iterator will be yielded.'''
        count = max(1 if not isinstance(var, VariableSet) else len(var) for var in self._variables)
        def value_set(n):
            for variable, value in self._variables.iteritems():
                if isinstance(variable, VariableSet):
                    yield variable[n], value
                else:
                    yield variable, value
        for i in xrange(count):
            yield value_set(i)

    def __add__(self, other):
        '''Add expression with a number or another expression'''
        if isinstance(other, numbers.Number):
            return self.__class__(self._variables, self._offset + other)
        elif isinstance(other, self.__class__):
            variables = Counter(self._variables)
            variables.update(other._variables)
            return self.__class__(variables, self._offset + other._offset)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + -other

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        return self.__class__({ var: value*other for var, value in self._variables.iteritems()}, self._offset*other)

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return self * -1

    def __eq__(self, other):
        '''Return equality relation (equation): self == other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Equals, self - other)

    def __ge__(self, other):
        '''Return greater-than relation (inequality): self >= other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Greater, self - other)

    def __le__(self, other):
        '''Return less-than relation (inequality): self <= other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.Less, self - other)

    def __gt__(self, other):
        '''Return strictly greater-than relation (inequality): self > other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.StrictlyGreater, self - other)

    def __lt__(self, other):
        '''Return strictly less-than relation (inequality): self < other

        This method is overloaded so that relations can be
        formed using a natural syntax.'''
        return Relation(Relation.StrictlyLess, self - other)

    def __str__(self):
        '''Return string representation of expression'''

        def all_terms():
            count_vars = 0
            for name, value in sorted(self._variables.iteritems()):
                if value != 0:
                    count_vars += 1
                    if isinstance(name, VariableSet):
                        yield '<set>', value
                    else:
                        yield name, value
            if self._offset != 0 or count_vars == 0:
                yield None, self._offset

        terms = []
        for i, spec in enumerate(all_terms()):
            name, value = spec
            if i == 0:
                # First term is special
                if name is None:
                    terms.append('{}'.format(value))
                elif abs(value) == 1:
                    terms.append(name if value > 0 else '-'+name)
                else:
                    terms.append('{}*{}'.format(value, name))
            else:
                prefix = '+' if value >= 0 else '-'
                if name is None:
                    terms.append('{} {}'.format(prefix, abs(value)))
                elif abs(value) == 1:
                    terms.append('{} {}'.format(prefix, name))
                else:
                    terms.append('{} {}*{}'.format(prefix, abs(value), name))
        return ' '.join(terms)

    def __repr__(self):
        return '<Expression \'{}\'>'.format(str(self))

class Relation(object):
    '''Represents a binary relation (equation or inequality)

    Relations can be equalities or inequalities. All relations
    of this type can be represented as a left-hand side expression
    and the type of relation. In this representation, the right-hand
    side is always zero.'''

    Equals = 'E'
    Greater = 'G'
    Less = 'L'
    StrictlyGreater = 'SG'
    StrictlyLess = 'SL'

    SYMBOL = {
        Equals: '==',
        Greater: '>=',
        Less: '<=',
        StrictlyGreater: '>',
        StrictlyLess: '<'
    }

    def __init__(self, sense, expression):
        self._sense = sense
        self._expression = expression

    @property
    def sense(self):
        '''Type of relation (equality or inequality)

        Can be one of Equal, Greater or Less, or one of the
        strict relations, StrictlyGreater or StrictlyLess.'''
        return self._sense

    @property
    def expression(self):
        '''Left-hand side expression'''
        return self._expression

    def __str__(self):
        '''Convert relation to string representation'''
        return '{} {} 0'.format(str(self._expression), Relation.SYMBOL[self._sense])

    def __repr__(self):
        return '<Relation \'{}\'>'.format(str(self))
