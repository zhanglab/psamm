
'''Representation of metabolic network databases'''

import abc
from collections import defaultdict, Mapping

from .reaction import ModelSEED, Reaction

class StoichiometricMatrixView(Mapping):
    '''Provides a sparse matrix view on the stoichiometry of a database

    This object is used internally in the database to expose
    a sparse matrix view of the model stoichiometry.

    Any compound, reaction-pair can be looked up to obtain the
    corresponding stoichiometric value. If the value is not
    defined (implicitly zero) a KeyError will be raised.'''

    def __init__(self, database):
        super(StoichiometricMatrixView, self).__init__()
        self._database = database

    def __getitem__(self, key):
        if len(key) != 2:
            raise KeyError(key)
        compound, reaction = key
        if not self._database.has_reaction(reaction):
            raise KeyError(key)
        reaction_values = dict(self._database.get_reaction_values(reaction))
        if compound not in reaction_values:
            raise KeyError(key)
        return reaction_values[compound]

    def __iter__(self):
        for reaction in self._database.reactions:
            for compound, _ in self._database.get_reaction_values(reaction):
                yield compound, reaction

    def __len__(self):
        return sum(sum(1 for _ in self._database.get_reaction_values(reaction)) for reaction in self._database.reactions)

class MetabolicReaction(object):
    def __init__(self, reversible, metabolites):
        self.reversible = bool(reversible)
        self.metabolites = dict(metabolites)

    @classmethod
    def from_reaction(cls, reaction):
        if reaction.direction not in (Reaction.Right, Reaction.Bidir):
            raise ValueError('Invalid direction in reaction: {}'.format(reaction.direction))
        reversible = reaction.direction == Reaction.Bidir

        metabolites = {}
        for compound, value in reaction.left:
            metabolites[compound] -= value
        for compound, value in reaction.right:
            metabolites[compound] += value

        return cls(reversible, metabolites)

    def __repr__(self):
        return 'MetabolicReaction({}, {})'.format(repr(self.reversible), repr(self.metabolites))

class MetabolicDatabase(object):
    '''Database of metabolic reactions'''

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def reactions(self):
        pass

    @abc.abstractproperty
    def compounds(self):
        pass

    @abc.abstractmethod
    def has_reaction(self, reaction_id):
        pass

    @abc.abstractmethod
    def is_reversible(self, reaction_id):
        pass

    @abc.abstractmethod
    def get_reaction_values(self, reaction_id):
        pass

    @abc.abstractmethod
    def get_compound_reactions(self, compound_id):
        pass

    @property
    def reversible(self):
        '''The set of reversible reactions'''
        return set(reaction_id for reaction_id in self.reactions if self.is_reversible(reaction_id))

    @property
    def matrix(self):
        '''Mapping from compound, reaction to stoichiometric value'''
        return StoichiometricMatrixView(self)

    def get_reaction(self, reaction_id):
        direction = Reaction.Bidir if self.is_reversible(reaction_id) else Reaction.Right
        left = ((compound, -value) for compound, value in self.get_reaction_values(reaction_id) if value < 0)
        right = ((compound, value) for compound, value in self.get_reaction_values(reaction_id) if value > 0)
        return Reaction(direction, left, right)

class DictDatabase(MetabolicDatabase):
    '''Metabolic database backed by in-memory dictionaries'''

    def __init__(self):
        super(DictDatabase, self).__init__()
        self._reactions = defaultdict(dict)
        self._compound_reactions = defaultdict(set)
        self._reversible = set()

    @property
    def reactions(self):
        return iter(self._reactions)

    @property
    def compounds(self):
        return iter(self._compound_reactions)

    def has_reaction(self, reaction_id):
        return reaction_id in self._reactions

    def is_reversible(self, reaction_id):
        '''Whether reaction is reversible'''
        return reaction_id in self._reversible

    def get_reaction_values(self, reaction_id):
        '''Yield compound and stoichiometric values of reaction as tuples'''
        if reaction_id not in self._reactions:
            raise ValueError('Unknown reaction: {}'.format(repr(reaction_id)))
        return self._reactions[reaction_id].iteritems()

    def get_compound_reactions(self, compound_id):
        '''Iterator of reactions that include the given compound'''
        return iter(self._compound_reactions[compound_id])

    def set_reaction(self, reaction_id, reaction):
        # Overwrite previous reaction if the same id is used
        if reaction_id in self._reactions:
            # Clean up compound to reaction mapping
            for compound in self._reactions[reaction_id]:
                self._compound_reactions[compound].remove(reaction_id)

            self._reversible.discard(reaction_id)
            del self._reactions[reaction_id]

        # Add values to global (sparse) stoichiometric matrix
        # Compounds that occur on both sides will get a stoichiometric
        # value based on the sum of the signed values on each side.
        for compound, value in reaction.compounds:
            if compound not in self._reactions[reaction_id] and value != 0:
                self._reactions[reaction_id][compound] = 0
                self._compound_reactions[compound].add(reaction_id)
        for compound, value in reaction.left:
            self._reactions[reaction_id][compound] -= value
        for compound, value in reaction.right:
            self._reactions[reaction_id][compound] += value

        if reaction.direction != '=>':
            self._reversible.add(reaction_id)

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
                reaction_id, equation = line.split(None, 1)
                reaction = ModelSEED.parse(equation).normalized()
                database.set_reaction(reaction_id, reaction)
        return database

class ChainedDatabase(MetabolicDatabase):
    '''Links a number of databases so they can be treated a single database'''

    def __init__(self, *databases):
        self._databases = list(databases)
        if len(self._databases) == 0:
            self._databases.append(DictDatabase())

    def _is_shadowed(self, reaction_id, database):
        '''Whether reaction in database is shadowed by another database'''
        for other_database in self._databases:
            if other_database == database:
                break
            if other_database.has_reaction(reaction_id):
                return True
        return False

    @property
    def reactions(self):
        # Make sure that we only yield each reaction once
        reaction_set = set()
        for database in self._databases:
            for reaction in database.reactions:
                if reaction not in reaction_set:
                    reaction_set.add(reaction)
                    yield reaction

    @property
    def compounds(self):
        # Make sure that we only yield each reaction once
        compound_set = set()
        for database in self._databases:
            for compound in database.compounds:
                if compound not in compound_set:
                    compound_set.add(compound)
                    yield compound

    def has_reaction(self, reaction_id):
        return any(database.has_reaction(reaction_id) for database in self._databases)

    def is_reversible(self, reaction_id):
        for database in self._databases:
            if database.has_reaction(reaction_id):
                return database.is_reversible(reaction_id)
        raise ValueError('Unknown reaction: {}'.format(reaction_id))

    def get_reaction_values(self, reaction_id):
        for database in self._databases:
            if database.has_reaction(reaction_id):
                return database.get_reaction_values(reaction_id)
        raise ValueError('Unknown reaction: {}'.format(reaction_id))

    def get_compound_reactions(self, compound):
        # Make sure that we only yield each reaction once
        reaction_set = set()
        for database in self._databases:
            for reaction in database.get_compound_reactions(compound):
                if reaction not in reaction_set and not self._is_shadowed(reaction, database):
                    reaction_set.add(reaction)
                    yield reaction

    def set_reaction(self, reaction_id, reaction):
        if hasattr(self._databases[0], 'set_reaction'):
            self._databases[0].set_reaction(reaction_id, reaction)
        else:
            raise ValueError('First database is immutable')
