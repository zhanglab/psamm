#!/usr/bin/env python

import pulp
from itertools import izip
from collections import namedtuple
import fastcore

Bounds = namedtuple('Bounds', ('lower', 'upper'))

class MetabolicModel(object):
    '''Represents a metabolic model containing a set of reactions

    The model contains a list of reactions as a stoichiometric matrix,
    as well as the reaction class (reversible, exchange) and flux
    bounds.'''

    def __init__(self):
        self.matrix = {}
        self.reversible = set()
        self.exchange = set()
        self.limits = {}
        self.compound_set = set()
        self.reaction_set = set()

    def flipped(self, subset):
        '''Return model with flipped stoichiometry of reversible reactions

        Both stoichiometric values and flux bounds are reversed in the
        returned model.'''

        model = self.__class__()
        model.matrix = dict((spec, -value if spec[1] in subset else value) for spec, value in self.matrix.iteritems())
        model.reversible = set(self.reversible)
        model.limits = dict((rxnid, Bounds(-b.upper, -b.lower) if rxnid in subset else Bounds(b.lower, b.upper))
                                for rxnid, b in self.limits.iteritems())
        model.compound_set = set(self.compound_set)
        model.reaction_set = set(self.reaction_set)
        return model

    @classmethod
    def load_model(cls, v_max=1000):
        '''Load model from definition in current directory'''

        model = cls()
        model_reactions = set()

        with open('modelrxn.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                model_reactions.add(rxnid)

        with open('exchangerxn.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                model.exchange.add(rxnid)

        with open('rev.txt', 'r') as f:
            for line in f:
                rxnid = line.strip()
                model.reversible.add(rxnid)

        with open('exchangelimit.txt', 'r') as f:
            for line in f:
                line = line.strip()
                if len(line) == 0 or line[0] == '*':
                    continue
                rxnid, value = line.split()

                if rxnid not in model.exchange:
                    raise Exception('Invalid exchange reaction in exchangelimit.txt')
                model.limits[rxnid] = Bounds(lower=float(value), upper=v_max)

        for rxnid in model.exchange:
            if rxnid not in model.limits:
                model.limits[rxnid] = Bounds(lower=0, upper=v_max)

        with open('mat.txt', 'r') as f:
            for line in f:
                fields = line.strip().split()
                cpdid, rxnid = fields[0].split('.', 2)
                value = float(fields[1])

                if rxnid in model_reactions:
                    model.compound_set.add(cpdid)
                    model.reaction_set.add(rxnid)
                    model.matrix[(cpdid, rxnid)] = value

                    if rxnid not in model.limits:
                        if rxnid in model.reversible:
                            model.limits[rxnid] = Bounds(lower=-v_max, upper=v_max)
                        else:
                            model.limits[rxnid] = Bounds(lower=0, upper=v_max)

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
    model = MetabolicModel.load_model()
    for rxnid, flux in flux_balance(model, 'Biomass'):
        print '{}\t{}'.format(rxnid, flux)

    # Remove bounds on exchange reactions
    for rxnid in model.exchange:
        model.limits[rxnid] = Bounds(lower=-1000, upper=1000)

    #print model.reaction_set - fastcore.fastcc(model, 0.001)

    core = { 'Biomass' }
    print fastcore.find_sparse_mode(model, core, model.reaction_set - core, False, 0.001)