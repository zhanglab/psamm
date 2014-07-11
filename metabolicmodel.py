#!/usr/bin/env python

import pulp
from itertools import izip
from collections import namedtuple, defaultdict
import fastcore

Bounds = namedtuple('Bounds', ('lower', 'upper'))

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
            if rxnid in self.reversible:
                model._limits[rxnid] = Bounds(lower=-v_max, upper=v_max)
            else:
                model._limits[rxnid] = Bounds(lower=0, upper=v_max)

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
        self._limits = {}
        self.reaction_set = set()
        self.compound_set = set()
        self._flipped = set()

    def add_reaction(self, reaction):
        '''Add reaction to model'''

        if reaction in self.reaction_set:
            return

        if reaction not in self.database.reactions:
            raise Exception('Model reaction does not reference a database reaction: {}'.format(rxnid))

        self.reaction_set.add(reaction)
        for compound, value in self.database.reactions[reaction].iteritems():
            self.compound_set.add(compound)

    def load_exchange_limits(self, v_max=1000):
        '''Load exchange limits from external file'''

        with open('exchangelimit.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if line == '' or line[0] == '*':
                    continue
                rxnid, value = line.split()

                if rxnid in model.reaction_set:
                    model._limits[rxnid] = Bounds(lower=float(value), upper=v_max)

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
                    lower, upper = self._model._limits[key]
                    return Bounds(lower=-upper, upper=-lower)
                return self._model._limits[key]

            def __setitem__(self, key, value):
                self._model._limits[key] = value

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

def flux_balance(model, reaction='Biomass'):
    # Create Flux balance problem
    prob = pulp.LpProblem('FluxBalance', pulp.LpMaximize)

    reaction_list = sorted(model.reaction_set)
    compound_list = sorted(model.compound_set)
    reaction_map = dict((rxnid, i) for i, rxnid in enumerate(reaction_list))
    compound_map = dict((cpdid, i) for i, cpdid in enumerate(compound_list))

    # Add variables to problem
    fluxes = []
    for rxnid in reaction_list:
        lower, upper = model.limits[rxnid]
        var = pulp.LpVariable('v_'+rxnid, lower, upper)
        fluxes.append(var)

    # Add constraints to problem
    massbalance_lhs = [0 for cpdid in compound_list]
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        cpdindex = compound_map[cpdid]
        rxnindex = reaction_map[rxnid]
        massbalance_lhs[cpdindex] += value*fluxes[rxnindex]

    for lhs in massbalance_lhs:
        prob += lhs == 0

    # Add objective
    objective = fluxes[reaction_map[reaction]]
    prob += objective

    # Solve
    status = prob.solve()
    if pulp.LpStatus[status] != 'Optimal':
        raise Exception('Non-optimal solution: {}'.format(pulp.LpStatus[status]))

    for rxnid, flux in izip(reaction_list, fluxes):
        yield rxnid, flux.value()

if __name__ == '__main__':
    database = MetabolicDatabase.load_from_file()
    model = database.load_model_from_file()

    # Run Fastcc and print the inconsistent set
    print 'Fastcc inconsistent set'
    print model.reaction_set - fastcore.fastcc(model, 0.001)

    # Run Fastcore and print the induced reaction set
    model_complete = model.copy()
    core = set(model.reaction_set)
    for rxnid in database.reactions:
        model_complete.add_reaction(rxnid)
        model_complete.limits[rxnid] = Bounds(lower=-1000, upper=1000) if rxnid in database.reversible else Bounds(lower=0, upper=1000)
    print 'Fastcore induced set with core = {}'.format(core)
    induced = fastcore.fastcore(model_complete, core, 0.001)
    print '|A| = {}, A = {}'.format(len(induced), induced)

    # Load bounds on exchange reactions
    model.load_exchange_limits()

    print 'Flux balance maximizing Biomass'
    for rxnid, flux in flux_balance(model, 'Biomass'):
        print '{}\t{}'.format(rxnid, flux)
