
'''Representation of metabolic networks'''

from collections import defaultdict
from itertools import chain

from .reaction import ModelSEED, Reaction, Compound

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

    @property
    def compounds(self):
        return set(self.compound_reactions)

    def get_reaction(self, rxnid):
        if rxnid not in self.reactions:
            raise ValueError('Unknown reaction {}'.format(repr(rxnid)))

        direction = '=>'
        if rxnid in self.reversible:
            direction = '<=>'

        left = ((Compound(compound[0]), -value, compound[1]) for compound, value in self.reactions[rxnid].iteritems() if value < 0)
        right = ((Compound(compound[0]), value, compound[1]) for compound, value in self.reactions[rxnid].iteritems() if value > 0)
        return Reaction(direction, left, right)

    def set_reaction(self, rxnid, reaction):
        # Overwrite previous reaction if the same id is used
        if rxnid in self.reactions:
            # Clean up compound to reaction mapping
            for compound in self.reactions[rxnid]:
                self.compound_reactions[compound].remove(rxnid)

            self.reversible.discard(rxnid)
            del self.reactions[rxnid]

        # Add values to global (sparse) stoichiometric matrix
        for compound, value, comp in reaction.left:
            self.reactions[rxnid][(compound.name, comp)] = -value
            self.compound_reactions[(compound.name, comp)].add(rxnid)
        for compound, value, comp in reaction.right:
            self.reactions[rxnid][(compound.name, comp)] = value
            self.compound_reactions[(compound.name, comp)].add(rxnid)

        if reaction.direction != '=>':
            self.reversible.add(rxnid)

    def load_model_from_file(self, reaction_list, v_max=1000):
        '''Load model defined by given reaction list file'''

        model = MetabolicModel(self)
        for line in reaction_list:
            rxnid = line.strip()
            model.add_reaction(rxnid)

        for rxnid in model.reaction_set:
            model.reset_flux_bounds(rxnid, v_max)

        return model

    @classmethod
    def load_from_files(cls, *files):
        '''Load database from given reactions definition lists'''

        database = cls()

        for file in files:
            for line in file:
                rxnid, equation = line.strip().split(None, 1)
                rx = ModelSEED.parse(equation).normalized()
                database.set_reaction(rxnid, rx)

        return database

class MetabolicModel(object):
    '''Represents a metabolic model containing a set of reactions

    The model contains a list of reactions referencing the reactions
    in the associated database.'''

    def __init__(self, database):
        self._database = database
        self._limits = defaultdict(FluxBounds)
        self._reaction_set = set()
        self._compound_set = set()

    @property
    def database(self):
        return self._database

    @property
    def reaction_set(self):
        return set(self._reaction_set)

    @property
    def compound_set(self):
        return set(self._compound_set)

    def add_reaction(self, reaction, v_max=1000):
        '''Add reaction to model'''

        if reaction in self._reaction_set:
            return

        if reaction not in self._database.reactions:
            raise Exception('Model reaction does not reference a database reaction: {}'.format(reaction))

        self._reaction_set.add(reaction)
        self._limits[reaction] = FluxBounds(-v_max, v_max) if reaction in self._database.reversible else FluxBounds(0, v_max)

        for compound, value in self._database.reactions[reaction].iteritems():
            self._compound_set.add(compound)

    def remove_reaction(self, reaction):
        '''Remove reaction from model'''

        if reaction not in self._reaction_set:
            return

        self._reaction_set.remove(reaction)
        del self._limits[reaction]

        # Remove compound from compound_set if it is not referenced
        # by any other reactions in the model.
        for compound, value in self._database.reactions[reaction].iteritems():
            if all(other_reaction not in self._reaction_set for other_reaction in self._database.compound_reactions[compound]):
                self._compound_set.remove(compound)

    def add_all_database_reactions(self, compartments={None, 'e'}):
        '''Add all reactions from database that occur in given compartments'''

        for rxnid in self._database.reactions:
            reaction = self._database.get_reaction(rxnid)
            if all(comp in compartments for compound, value, comp in reaction.compounds):
                self.add_reaction(rxnid)

    def add_all_exchange_reactions(self, allow_duplicates=False):
        '''Add all exchange reactions to database and to model'''

        all_reactions = {}
        if not allow_duplicates:
            # TODO: Avoid adding reactions that already exist in the database.
            # This should be integrated in the database.
            for rxnid in self._database.reactions:
                rx = self._database.get_reaction(rxnid)
                all_reactions[rx] = rxnid

        for cpdid, comp in sorted(self.compound_set):
            rxnid_ex = 'rxnex_'+cpdid
            if rxnid_ex not in self._database.reactions:
                reaction_ex = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [])
                if reaction_ex not in all_reactions:
                    self._database.set_reaction(rxnid_ex, reaction_ex)
                    self.add_reaction(rxnid_ex)
                else:
                    self.add_reaction(all_reactions[reaction_ex])
            else:
                self.add_reaction(rxnid_ex)

    def add_all_transport_reactions(self, allow_duplicates=False):
        '''Add all transport reactions to database and to model'''

        all_reactions = {}
        if not allow_duplicates:
            # TODO: Avoid adding reactions that already exist in the database.
            # This should be integrated in the database.
            for rxnid in self._database.reactions:
                rx = self._database.get_reaction(rxnid)
                all_reactions[rx] = rxnid

        for cpdid, comp in sorted(self.compound_set):
            rxnid_tp = 'rxntp_'+cpdid
            if rxnid_tp not in self._database.reactions:
                reaction_tp = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [(Compound(cpdid), 1, None)])
                if reaction_tp not in all_reactions:
                    self._database.set_reaction(rxnid_tp, reaction_tp)
                    self.add_reaction(rxnid_tp)
                else:
                    self.add_reaction(all_reactions[reaction_tp])
            else:
                self.add_reaction(rxnid_tp)

    def reset_flux_bounds(self, reaction, v_max=1000):
        '''Reset flux bounds of model reaction

        A reversible reaction will be reset to (-v_max, v_max) and
        an irreversible reaction will be set to (0, v_max).'''

        self._limits[reaction].bounds = (-v_max, v_max) if reaction in self._database.reversible else (0, v_max)

    def load_exchange_limits(self, v_max=1000):
        '''Load exchange limits from external file'''

        with open('exchangerxn.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                if rxnid in self._reaction_set:
                    self._limits[rxnid].bounds = 0, v_max

        with open('exchangelimit.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if line == '' or line[0] == '*':
                    continue
                rxnid, value = line.split()

                if rxnid in self._reaction_set:
                    self._limits[rxnid].bounds = float(value), v_max

    @property
    def reversible(self):
        '''The set of reversible reactions'''
        return self._database.reversible & self._reaction_set

    class MatrixView(object):
        def __init__(self, model):
            self._model = model

        def __getitem__(self, key):
            if len(key) != 2:
                raise TypeError(repr(key))
            cpdid, rxnid = key
            if rxnid not in self._model._reaction_set:
                raise KeyError(repr(key))
            value = self._model._database.reactions[rxnid][cpdid]
            return value

        def __contains__(self, key):
            try:
                self.__getitem__(key)
            except KeyError:
                return False
            return True

        def __iter__(self):
            return self.iterkeys()

        def iterkeys(self):
            for rxnid in self._model._reaction_set:
                for cpdid in self._model._database.reactions[rxnid]:
                    yield cpdid, rxnid

        def iteritems(self):
            for rxnid in self._model._reaction_set:
                for cpdid, value in self._model._database.reactions[rxnid].iteritems():
                    yield (cpdid, rxnid), value

    @property
    def matrix(self):
        '''Mapping from compound, reaction to stoichiometric value'''
        return MetabolicModel.MatrixView(self)

    @property
    def limits(self):
        return self._limits

    def copy(self):
        '''Return copy of model'''

        model = self.__class__(self._database)
        model._limits = dict(self._limits)
        model._reaction_set = set(self._reaction_set)
        model._compound_set = set(self._compound_set)
        return model

class FlipableModelView(object):
    '''Proxy wrapper of model objects allowing a flipped set of reactions

    The proxy will forward all properties normally except
    that flipped reactions will appear to have stoichiometric
    values negated in the matrix property, and have bounds in
    the limits property flipped. This view is needed for the
    some algorithms.'''

    def __init__(self, model, flipped=set()):
        self._model = model
        self._flipped = set(flipped)

    class MatrixView(object):
        def __init__(self, view):
            self._view = view

        def _value_mul(self, rxnid):
            return -1 if rxnid in self._view._flipped else 1

        def __getitem__(self, key):
            if len(key) != 2:
                raise TypeError(repr(key))
            cpdid, rxnid = key
            value = self._view._model.matrix[key]
            return value * self._value_mul(rxnid)

        def __contains__(self, key):
            return key in self._view._model.matrix

        def __iter__(self):
            return self.iterkeys()

        def iterkeys(self):
            return self._view._model.matrix.iterkeys()

        def iteritems(self):
            for key, value in self._view._model.matrix.iteritems():
                cpdid, rxnid = key
                yield key, value * self._value_mul(rxnid)

    class LimitsView(object):
        def __init__(self, view):
            self._view = view

        def __getitem__(self, key):
            if key in self._view._flipped:
                return self._view._model.limits[key].flipped()
            return self._view._model.limits[key]

        def __contains__(self, key):
            return key in self._view._model.limits

        def __iter__(self):
            return self.iterkeys()

        def iterkeys(self):
            return self._view._model.limits.iterkeys()

        def iteritems(self):
            for key in iter(self):
                yield key, self[key]

    @property
    def matrix(self):
        return FlipableModelView.MatrixView(self)

    @property
    def limits(self):
        return FlipableModelView.LimitsView(self)

    def flip(self, subset):
        self._flipped ^= subset

    def __getattr__(self, name):
        return getattr(self._model, name)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
