
'''Representation of metabolic networks'''

from collections import defaultdict
from itertools import chain
import fastcore
import fluxanalysis

class FluxBounds(object):
    '''Represents lower and upper bounds of flux as a mutable object

    >>> FluxBounds(-5, 1000)
    FluxBounds(-5, 1000)

    >>> FluxBounds(100, 0)
    Traceback (most recent call last):
        ...
    ValueError: Lower bound larger than upper bound
    '''

    def __init__(self, lower=0, upper=0):
        self._check_bounds(lower, upper)
        self._lower = lower
        self._upper = upper

    def __iter__(self):
        '''Iterator over lower and upper value

        >>> it = iter(FluxBounds(-5, 1000))
        >>> next(it)
        -5
        >>> next(it)
        1000
        >>> next(it)
        Traceback (most recent call last):
            ...
        StopIteration
        '''
        yield self._lower
        yield self._upper

    def _check_bounds(self, lower, upper):
        if lower > upper:
            raise ValueError('Lower bound larger than upper bound')

    @property
    def lower(self):
        '''Lower bound

        >>> f = FluxBounds(-5, 1000)
        >>> f.lower
        -5
        >>> f.lower = -20
        >>> f
        FluxBounds(-20, 1000)
        >>> f.lower = 1200
        Traceback (most recent call last):
            ...
        ValueError: Lower bound larger than upper bound'''
        return self._lower

    @lower.setter
    def lower(self, value):
        self._check_bounds(value, self._upper)
        self._lower = value

    @property
    def upper(self):
        '''Upper bound

        >>> f = FluxBounds(-5, 1000)
        >>> f.upper
        1000
        >>> f.upper = 0
        >>> f
        FluxBounds(-5, 0)
        >>> f.upper = -6
        Traceback (most recent call last):
            ...
        ValueError: Lower bound larger than upper bound'''
        return self._upper

    @upper.setter
    def upper(self, value):
        self._check_bounds(self._lower, value)
        self._upper = value

    @property
    def bounds(self):
        '''Bounds as a tuple

        >>> f = FluxBounds(-5, 1000)
        >>> f.bounds
        (-5, 1000)
        >>> f.bounds = 10, 20
        >>> f
        FluxBounds(10, 20)
        >>> f.bounds = 20, -10
        Traceback (most recent call last):
            ...
        ValueError: Lower bound larger than upper bound'''
        return self._lower, self._upper

    @bounds.setter
    def bounds(self, value):
        lower, upper = value
        self._check_bounds(lower, upper)
        self._lower, self._upper = value

    def flipped(self):
        '''New FluxBounds with bounds mirrored around zero

        >>> FluxBounds(-5, 1000).flipped()
        FluxBounds(-1000, 5)'''
        return self.__class__(-self._upper, -self._lower)

    def __eq__(self, other):
        '''Equality test

        >>> f = FluxBounds(-5, 10)
        >>> f == FluxBounds(-5, 10)
        True
        >>> f == FluxBounds(0, 10)
        False'''
        return isinstance(other, FluxBounds) and self.bounds == other.bounds

    def __ne__(self, other):
        '''Inequality test

        >>> f = FluxBounds(-5, 10)
        >>> f != FluxBounds(-5, 10)
        False'''
        return not self == other

    def __repr__(self):
        return 'FluxBounds({}, {})'.format(repr(self._lower), repr(self._upper))

class MetabolicReaction(object):
    def __init__(self, reversible, metabolites):
        self.reversible = bool(reversible)
        self.metabolites = dict(metabolites)

    @classmethod
    def from_reaction(cls, reaction):
        if reaction.direction not in ('=>', '<=>'):
            raise ValueError('Invalid direction in reaction: {}'.format(reaction.direction))
        reversible = reaction.direction == '<=>'
        return cls(reversible, chain((((cpdid, comp), -value) for cpdid, value, comp in reaction.left),
                                        (((cpdid, comp), value) for cpdid, value, comp in reaction.right)))

    def __repr__(self):
        return 'MetabolicReaction({}, {})'.format(repr(self.reversible), repr(self.metabolites))

class MetabolicDatabase(object):
    '''Database of metabolic reactions'''

    def __init__(self):
        self.reactions = defaultdict(dict)
        self.compound_reactions = defaultdict(set)
        self.reversible = set()

    def load_model_from_file(self, v_max=1000):
        '''Load model from definition in current directory'''

        model = MetabolicModel(self)

        with open('modelrxn.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                model.add_reaction(rxnid)

        for rxnid in model.reaction_set:
            model.reset_flux_bounds(rxnid, v_max)

        return model

    @classmethod
    def load_from_file(cls):
        '''Load database from definition in current directory'''

        database = cls()

        with open('mat.txt', 'r') as f:
            for line in f:
                fields = line.strip().split()
                cpdid, rxnid = fields[0].split('.', 2)
                value = float(fields[1])

                database.reactions[rxnid][cpdid] = value
                database.compound_reactions[cpdid].add(rxnid)

        with open('rev.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                if rxnid in database.reactions:
                    database.reversible.add(rxnid)

        return database

class MetabolicModel(object):
    '''Represents a metabolic model containing a set of reactions

    The model contains a list of reactions referencing the reactions
    in the associated database. To support certain optimization
    algorithms that operate on the stoichiometric matrix it is possible
    to flip a subset of the reactions such that both flux limits and
    stoichiometric values are reversed.'''

    def __init__(self, database):
        self.database = database
        self._limits = defaultdict(FluxBounds)
        self.reaction_set = set()
        self.compound_set = set()
        self._flipped = set()

    def add_reaction(self, reaction, v_max=1000):
        '''Add reaction to model'''

        if reaction in self.reaction_set:
            return

        if reaction not in self.database.reactions:
            raise Exception('Model reaction does not reference a database reaction: {}'.format(reaction))

        self.reaction_set.add(reaction)
        self._limits[reaction] = FluxBounds(-v_max, v_max) if reaction in self.database.reversible else FluxBounds(0, v_max)

        for compound, value in self.database.reactions[reaction].iteritems():
            self.compound_set.add(compound)

    def reset_flux_bounds(self, reaction, v_max=1000):
        '''Reset flux bounds of model reaction

        A reversible reaction will be reset to (-v_max, v_max) and
        an irreversible reaction will be set to (0, v_max).'''

        self._limits[reaction].bounds = (-v_max, v_max) if reaction in self.database.reversible else (0, v_max)

    def load_exchange_limits(self, v_max=1000):
        '''Load exchange limits from external file'''

        with open('exchangerxn.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                if rxnid in self.reaction_set:
                    self._limits[rxnid].bounds = 0, v_max

        with open('exchangelimit.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if line == '' or line[0] == '*':
                    continue
                rxnid, value = line.split()

                if rxnid in self.reaction_set:
                    self._limits[rxnid].bounds = float(value), v_max

    @property
    def reversible(self):
        '''The set of reversible reactions'''
        return self.database.reversible & self.reaction_set

    @property
    def matrix(self):
        '''Mapping from compound, reaction to stoichiometric value'''

        class MatrixView(object):
            def __init__(self, model):
                self._model = model

            def _value_mul(self, rxnid):
                return -1 if rxnid in self._model._flipped else 1

            def __getitem__(self, key):
                if len(key) != 2:
                    raise KeyError(repr(key))
                cpdid, rxnid = key
                if rxnid not in self._model.reaction_set:
                    raise KeyError(repr(key))
                value = self._model.database.reactions[rxnid][cpdid]
                return value * self._value_mul(rxnid)

            def __contains__(self, key):
                try:
                    self.__getitem__(key)
                except KeyError:
                    return False
                return True

            def __iter__(self):
                return self.iterkeys()

            def iterkeys(self):
                for rxnid in self._model.reaction_set:
                    for cpdid in self._model.database.reactions[rxnid]:
                        yield cpdid, rxnid

            def iteritems(self):
                for rxnid in self._model.reaction_set:
                    for cpdid, value in self._model.database.reactions[rxnid].iteritems():
                        value = value * self._value_mul(rxnid)
                        yield (cpdid, rxnid), value

        return MatrixView(self)

    @property
    def limits(self):
        '''Bounds for reactions fluxes'''

        class LimitsView(object):
            def __init__(self, model):
                self._model = model

            def __getitem__(self, key):
                if key in self._model._flipped:
                    return self._model._limits[key].flipped()
                return self._model._limits[key]

            def __contains__(self, key):
                return key in self._model._limits

            def __iter__(self):
                return self.iterkeys()

            def iterkeys(self):
                return self._model._limits.iterkeys()

            def iteritems(self):
                for key in iter(self):
                    yield key, self[key]

        return LimitsView(self)

    def copy(self):
        '''Return copy of model'''

        model = self.__class__(self.database)
        model._limits = dict(self._limits)
        model.reaction_set = set(self.reaction_set)
        model.compound_set = set(self.compound_set)
        model._flipped = set(self._flipped)
        return model

    def flip(self, subset):
        '''Flip stoichiometry and limits of reactions in subset'''
        self._flipped ^= subset

    def flipped(self, subset):
        '''Return model with flipped stoichiometry of the subset reactions

        Both stoichiometric values and flux bounds are reversed in the
        returned model.'''

        model = self.copy()
        model.flip(subset)
        return model


if __name__ == '__main__':
    import doctest
    doctest.testmod()