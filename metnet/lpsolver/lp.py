
'''Base objects for representation of LP problems'''

import numbers
from collections import Counter
import abc

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


class ObjectiveSense(object):
    '''Enumeration of objective sense values'''
    Minimize = object()
    Maximize = object()

class VariableType(object):
    '''Enumeration of variable types'''
    Continuous = object()
    Integer = object()
    Binary = object()


class Solver(object):
    '''Factory for LP Problem instances'''

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def create_problem(self):
        pass

class Problem(object):
    '''Representation of LP Problem instance

    Names in the problem can be any hashable object. It is the responsibility
    of the implementation to translate the object into a unique string, if
    required by the underlying LP solver.'''

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def define(self, *names, **kwargs):
        '''Define a variable in the problem'''
        pass

    @abc.abstractmethod
    def var(self, name):
        '''Return variable as an expression'''
        pass

    @abc.abstractmethod
    def set(self, names):
        '''Return variables as a set expression'''
        pass

    @abc.abstractmethod
    def add_linear_constraints(self, *relations):
        '''Add constraints to the problem

        Each constraint is represented by a Relation, and the
        expression in that relation can be a set expression.'''
        pass

    @abc.abstractmethod
    def set_linear_objective(self, expression):
        '''Set linear objective of the problem'''
        pass

    @abc.abstractmethod
    def set_objective_sense(self, sense):
        '''Set type of problem (minimize or maximize)'''
        pass

    @abc.abstractmethod
    def solve(self, sense=None):
        '''Solve problem and return result'''
        pass

    @abc.abstractproperty
    def result(self):
        '''Result of solved problem'''
        pass

class Result(object):
    '''Result of solving an LP problem'''

    @abc.abstractproperty
    def success(self):
        '''Whether solution was optimal'''
        pass

    def __bool__(self):
        '''Whether solution was optimal'''
        return self.success

    @abc.abstractproperty
    def status(self):
        '''String indicating the status of the problem result'''
        pass

    @abc.abstractmethod
    def get_value(self, expression):
        '''Get value of variable or expression in result

        Expression can be an object defined as a name in the problem, in which case
        the corresponding value is simply returned. If expression is an actual
        Expression object, it will be evaluated using the values from the result.'''
        pass
