
'''Representation of metabolic networks'''

from collections import defaultdict, Mapping
from itertools import chain

from .reaction import ModelSEED, Reaction, Compound

class FluxBounds(object):
    '''Represents lower and upper bounds of flux as a mutable object

    This object is used internally in the model representation. Changing
    the state of the object will change the underlying model parameters.
    Deleting a value will reset that value to the defaults.'''

    def __init__(self, model, reaction):
        self._model = model
        self._reaction = reaction

    def __iter__(self):
        '''Iterator over lower and upper value'''
        yield self.lower
        yield self.upper

    def _assign_lower(self, value):
        self._model._limits_lower[self._reaction] = value

    def _assign_upper(self, value):
        self._model._limits_upper[self._reaction] = value

    def _assign_both(self, lower, upper):
        self._assign_lower(lower)
        self._assign_upper(upper)

    def _check_bounds(self, lower, upper):
        if lower > upper:
            raise ValueError('Lower bound larger than upper bound')

    @property
    def lower(self):
        '''Lower bound'''
        try:
            return self._model._limits_lower[self._reaction]
        except KeyError:
            return -self._model._v_max if self._reaction in self._model._database.reversible else 0

    @lower.setter
    def lower(self, value):
        self._check_bounds(value, self.upper)
        self._assign_lower(value)

    @lower.deleter
    def lower(self):
        self._model._limits_lower.pop(self._reaction, None)

    @property
    def upper(self):
        '''Upper bound'''
        try:
            return self._model._limits_upper[self._reaction]
        except KeyError:
            return self._model._v_max

    @upper.setter
    def upper(self, value):
        self._check_bounds(self.lower, value)
        self._assign_upper(value)

    @upper.deleter
    def upper(self):
        self._model._limits_upper.pop(self._reaction, None)

    @property
    def bounds(self):
        '''Bounds as a tuple'''
        return self.lower, self.upper

    @bounds.setter
    def bounds(self, value):
        lower, upper = value
        self._check_bounds(lower, upper)
        self._assign_both(lower, upper)

    @bounds.deleter
    def bounds(self):
        del self.lower
        del self.upper

    def __eq__(self, other):
        '''Equality test'''
        return isinstance(other, FluxBounds) and self.lower == other.lower and self.upper == other.upper

    def __ne__(self, other):
        '''Inequality test'''
        return not self == other

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, repr(self.lower), repr(self.upper))

class StoichiometricMatrixView(Mapping):
    '''Provides a sparse matrix view on the stoichiometry of a model

    This object is used internally in the MetabolicModel to expose
    a sparse matrix view of the model stoichiometry.

    Any compound, reaction-pair can be looked up to obtain the
    corresponding stoichiometric value. If the value is not
    defined (implicitly zero) a KeyError will be raised.'''

    def __init__(self, model):
        super(StoichiometricMatrixView, self).__init__()
        self._model = model

    def __getitem__(self, key):
        if len(key) != 2:
            raise KeyError(key)
        compound, reaction = key
        if reaction not in self._model._reaction_set:
            raise KeyError(key)
        value = self._model._database.reactions[reaction][compound]
        return value

    def __iter__(self):
        for reaction in self._model._reaction_set:
            for compound in self._model._database.reactions[reaction]:
                yield compound, reaction

    def __len__(self):
        return sum(len(self._model._database.reactions[reaction]) for reaction in self._model._reaction_set)

class LimitsView(Mapping):
    '''Provides a view of the flux bounds defined in the model

    This object is used internally in MetabolicModel to
    expose a dictonary view of the FluxBounds associated
    with the model reactions.'''

    def __init__(self, model):
        super(LimitsView, self).__init__()
        self._model = model

    def _create_bounds(self, reaction):
        return FluxBounds(self._model, reaction)

    def __getitem__(self, key):
        if key not in self._model._reaction_set:
            raise KeyError(key)
        return self._create_bounds(key)

    def __iter__(self):
        return iter(self._model._reaction_set)

    def __len__(self):
        return len(self._model._reaction_set)

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
        return iter(self.compound_reactions)

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
        # Compounds that occur on both sides will get a stoichiometric
        # value based on the sum of the signed values on each side.
        for compound, value, comp in reaction.compounds:
            if (compound.name, comp) not in self.reactions[rxnid] and value != 0:
                self.reactions[rxnid][compound.name, comp] = 0
                self.compound_reactions[compound.name, comp].add(rxnid)
        for compound, value, comp in reaction.left:
            self.reactions[rxnid][(compound.name, comp)] -= value
        for compound, value, comp in reaction.right:
            self.reactions[rxnid][(compound.name, comp)] += value

        if reaction.direction != '=>':
            self.reversible.add(rxnid)

    def load_model_from_file(self, file):
        '''Load model defined by given reaction list file

        Comments, indicated by pound sign (#), are skipped.'''

        def reaction_file_iter(f):
            for line in f:
                line, _, comment = line.partition('#')
                line = line.strip()
                if line == '':
                    continue
                yield line

        return self.get_model(reaction_file_iter(file))

    def get_model(self, reaction_iter):
        '''Get model from reaction name iterator

        The model will contain all reactions of the iterator.'''
        model = MetabolicModel(self)
        for rxnid in reaction_iter:
            model.add_reaction(rxnid)

        for rxnid in model.reaction_set:
            del model.limits[rxnid].bounds

        return model

    @classmethod
    def load_from_files(cls, *files):
        '''Load database from given reactions definition lists'''

        database = cls()

        for file in files:
            for line in file:
                line, _, comment = line.partition('#')
                line = line.strip()
                if line == '':
                    continue
                rxnid, equation = line.split(None, 1)
                rx = ModelSEED.parse(equation).normalized()
                database.set_reaction(rxnid, rx)

        return database

class MetabolicModel(object):
    '''Represents a metabolic model containing a set of reactions

    The model contains a list of reactions referencing the reactions
    in the associated database.'''

    def __init__(self, database, v_max=1000):
        self._database = database
        self._limits_lower = {}
        self._limits_upper = {}
        self._reaction_set = set()
        self._compound_set = set()
        self._v_max = v_max

    @property
    def database(self):
        return self._database

    @property
    def reaction_set(self):
        return iter(self._reaction_set)

    @property
    def compound_set(self):
        return iter(self._compound_set)

    def add_reaction(self, reaction):
        '''Add reaction to model'''

        if reaction in self._reaction_set:
            return

        if reaction not in self._database.reactions:
            raise Exception('Model reaction does not reference a database reaction: {}'.format(reaction))

        self._reaction_set.add(reaction)
        for compound, value in self._database.reactions[reaction].iteritems():
            self._compound_set.add(compound)

    def remove_reaction(self, reaction):
        '''Remove reaction from model'''

        if reaction not in self._reaction_set:
            return

        self._reaction_set.remove(reaction)
        self._limits_lower.pop(reaction, None)
        self._limits_upper.pop(reaction, None)

        # Remove compound from compound_set if it is not referenced
        # by any other reactions in the model.
        for compound, value in self._database.reactions[reaction].iteritems():
            if all(other_reaction not in self._reaction_set for other_reaction in self._database.compound_reactions[compound]):
                self._compound_set.remove(compound)

    def add_all_database_reactions(self, compartments={None, 'e'}):
        '''Add all reactions from database that occur in given compartments'''

        added = set()
        for rxnid in self._database.reactions:
            reaction = self._database.get_reaction(rxnid)
            if all(comp in compartments for compound, value, comp in reaction.compounds):
                if rxnid not in self.reaction_set:
                    added.add(rxnid)
                self.add_reaction(rxnid)

        return added

    def add_all_exchange_reactions(self, allow_duplicates=False):
        '''Add all exchange reactions to database and to model'''

        all_reactions = {}
        if not allow_duplicates:
            # TODO: Avoid adding reactions that already exist in the database.
            # This should be integrated in the database.
            for rxnid in self._database.reactions:
                rx = self._database.get_reaction(rxnid)
                all_reactions[rx] = rxnid

        added = set()
        for cpdid, comp in sorted(self.compound_set):
            rxnid_ex = 'rxnex_'+cpdid
            if rxnid_ex not in self._database.reactions:
                reaction_ex = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [])
                if reaction_ex not in all_reactions:
                    self._database.set_reaction(rxnid_ex, reaction_ex)
                else:
                    rxnid_ex = all_reactions[reaction_ex]

            if rxnid_ex not in self.reaction_set:
                added.add(rxnid_ex)
            self.add_reaction(rxnid_ex)

        return added

    def add_all_transport_reactions(self, allow_duplicates=False):
        '''Add all transport reactions to database and to model'''

        all_reactions = {}
        if not allow_duplicates:
            # TODO: Avoid adding reactions that already exist in the database.
            # This should be integrated in the database.
            for rxnid in self._database.reactions:
                rx = self._database.get_reaction(rxnid)
                all_reactions[rx] = rxnid

        def tp_id(cpdid, comp):
            if comp is None:
                return 'rxntp_'+cpdid
            return 'rxntp_'+cpdid+'_'+comp

        added = set()
        for cpdid, comp in sorted(self.compound_set):
            if comp == 'e':
                # A transport reaction with exchange would not be valid
                continue

            rxnid_tp = tp_id(cpdid, comp)
            if rxnid_tp not in self._database.reactions:
                reaction_tp = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [(Compound(cpdid), 1, comp)])
                if reaction_tp not in all_reactions:
                    self._database.set_reaction(rxnid_tp, reaction_tp)
                else:
                    rxnid_tp = all_reactions[reaction_tp]

            if rxnid_tp not in self.reaction_set:
                added.add(rxnid_tp)
            self.add_reaction(rxnid_tp)

        return added

    def load_reaction_limits(self, limits_file):
        '''Load reaction limits from external file'''

        for line in limits_file:
            line, _, comment = line.partition('#')
            line = line.strip()
            # TODO Comments can start with an asterisk to remain
            # compatible with GAMS files. Can be removed when
            # compatibility is no longer needed.
            if line == '' or line[0] == '*':
                continue

            # A line can specify lower limit only (useful for
            # exchange reactions), or both lower and upper limit.
            fields = line.split(None)
            if len(fields) == 2:
                reaction_id, lower = fields
                if reaction_id in self._reaction_set:
                    self.limits[reaction_id].lower = float(lower)
            elif len(fields) == 3:
                reaction_id, lower, upper = fields
                if reaction_id in self._reaction_set:
                    self.limits[reaction_id].bounds = float(lower), float(upper)
            else:
                raise ValueError('Malformed reaction limit: {}'.format(fields))

    @property
    def reversible(self):
        '''The set of reversible reactions'''
        return self._database.reversible & self._reaction_set

    @property
    def matrix(self):
        '''Mapping from compound, reaction to stoichiometric value'''
        return StoichiometricMatrixView(self)

    @property
    def limits(self):
        return LimitsView(self)

    def copy(self):
        '''Return copy of model'''

        model = self.__class__(self._database)
        model._limits_lower = dict(self._limits_lower)
        model._limits_upper = dict(self._limits_upper)
        model._reaction_set = set(self._reaction_set)
        model._compound_set = set(self._compound_set)
        return model


class FlipableFluxBounds(FluxBounds):
    '''FluxBounds object for a FlipableModelView

    This object is used internally in the FlipableModelView to represent
    the bounds of flux on a reaction that can be flipped.'''

    def __init__(self, view, reaction):
        super(FlipableFluxBounds, self).__init__(view._model, reaction)
        self._view = view

    def _assign_lower(self, value):
        if self._reaction in self._view._flipped:
            super(FlipableFluxBounds, self)._assign_upper(-value)
        else:
            super(FlipableFluxBounds, self)._assign_lower(value)

    def _assign_upper(self, value):
        if self._reaction in self._view._flipped:
            super(FlipableFluxBounds, self)._assign_lower(-value)
        else:
            super(FlipableFluxBounds, self)._assign_upper(value)

    @property
    def lower(self):
        '''Lower bound'''
        if self._reaction in self._view._flipped:
            return -super(FlipableFluxBounds, self).upper
        return super(FlipableFluxBounds, self).lower

    @property
    def upper(self):
        '''Upper bound'''
        if self._reaction in self._view._flipped:
            return -super(FlipableFluxBounds, self).lower
        return super(FlipableFluxBounds, self).upper

class FlipableStoichiometricMatrixView(StoichiometricMatrixView):
    '''Provides a matrix view that flips with the underlying flipable model view

    This object is used internally in FlipableModelView to
    expose a matrix view that negates the stoichiometric
    values of flipped reactions.'''

    def __init__(self, view):
        super(FlipableStoichiometricMatrixView, self).__init__(view._model)
        self._view = view

    def _value_mul(self, reaction):
        return -1 if reaction in self._view._flipped else 1

    def __getitem__(self, key):
        if len(key) != 2:
            raise KeyError(key)
        compound, reaction = key
        return self._value_mul(reaction) * super(FlipableStoichiometricMatrixView, self).__getitem__(key)

class FlipableLimitsView(LimitsView):
    '''Provides a limits view that flips with the underlying flipable model view

    This object is used internally in FlipableModelView to
    expose a limits view that flips the bounds of all flipped
    reactions.'''

    def __init__(self, view):
        super(FlipableLimitsView, self).__init__(view._model)
        self._view = view

    def _create_bounds(self, reaction):
        return FlipableFluxBounds(self._view, reaction)

class FlipableModelView(object):
    '''Proxy wrapper of model objects allowing a flipped set of reactions

    The proxy will forward all properties normally except
    that flipped reactions will appear to have stoichiometric
    values negated in the matrix property, and have bounds in
    the limits property flipped. This view is needed for
    some algorithms.'''

    def __init__(self, model, flipped=set()):
        self._model = model
        self._flipped = set(flipped)

    @property
    def matrix(self):
        return FlipableStoichiometricMatrixView(self)

    @property
    def limits(self):
        return FlipableLimitsView(self)

    def flip(self, subset):
        self._flipped ^= subset

    def __getattr__(self, name):
        return getattr(self._model, name)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
