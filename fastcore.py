
'''Fastcore module implementing the fastcore algorithm

This is an implementation of the algorithms described in
Vlassis, Nikos, Maria Pires Pacheco, and Thomas Sauter.
"Fast reconstruction of compact context-specific metabolic
network models." PLoS computational biology 10.1 (2014):
e1003424.'''

import pulp
from itertools import izip

def fastcore_lp3(model, reaction_subset):
    # Create LP-3 problem of Fastcore
    prob = pulp.LpProblem('Fastcore-LP-3', pulp.LpMaximize)

    reaction_list = sorted(model.reaction_set)
    compound_list = sorted(model.compound_set)
    reaction_map = dict((rxnid, i) for i, rxnid in enumerate(reaction_list))
    compound_map = dict((cpdid, i) for i, cpdid in enumerate(compound_list))

    # Define flux variables
    fluxes = []
    for rxnid in reaction_list:
        lower, upper = model.limits[rxnid]
        var = pulp.LpVariable('v_'+rxnid, lower, upper)
        fluxes.append(var)

    # Define constraints
    massbalance_lhs = [0 for cpdid in compound_list]
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        cpdindex = compound_map[cpdid]
        rxnindex = reaction_map[rxnid]
        massbalance_lhs[cpdindex] += value*fluxes[rxnindex]

    for lhs in massbalance_lhs:
        prob += lhs == 0

    # Define objective
    objective = sum(fluxes[reaction_map[rxnid]] for rxnid in reaction_subset)
    prob += objective

    # Solve
    status = prob.solve()
    if pulp.LpStatus[status] != 'Optimal':
        raise Exception('Non-optimal solution: {}'.format(pulp.LpStatus[status]))

    for rxnid, flux in izip(reaction_list, fluxes):
        #print '{}\t{}'.format(rxnid, flux.value())
        yield rxnid, flux.value()

def fastcore_lp7(model, reaction_subset, epsilon):
    # Create LP-7 problem of Fastcore
    prob = pulp.LpProblem('Fastcore-LP-7', pulp.LpMaximize)

    reaction_list = sorted(model.reaction_set)
    compound_list = sorted(model.compound_set)
    reaction_map = dict((rxnid, i) for i, rxnid in enumerate(reaction_list))
    compound_map = dict((cpdid, i) for i, cpdid in enumerate(compound_list))

    # Define flux variables
    fluxes = []
    for rxnid in reaction_list:
        lower, upper = model.limits[rxnid]
        var = pulp.LpVariable('v_'+rxnid, lower, upper)
        fluxes.append(var)

    # Define z variables and constraints
    zs = []
    for rxnid in reaction_subset:
        var = pulp.LpVariable('z_'+rxnid, 0, epsilon)
        zs.append(var)

        rxnindex = reaction_map[rxnid]
        prob += fluxes[rxnindex] >= var

    # Define constraints
    massbalance_lhs = [0 for cpdid in compound_list]
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        cpdindex = compound_map[cpdid]
        rxnindex = reaction_map[rxnid]
        massbalance_lhs[cpdindex] += value*fluxes[rxnindex]

    for lhs in massbalance_lhs:
        prob += lhs == 0

    # Define objective
    objective = sum(zs)
    prob += objective

    # Solve
    status = prob.solve()
    if pulp.LpStatus[status] != 'Optimal':
        raise Exception('Non-optimal solution: {}'.format(pulp.LpStatus[status]))

    for rxnid, flux in izip(reaction_list, fluxes):
        #print '{}\t{}'.format(rxnid, flux.value())
        yield rxnid, flux.value()

def fastcore_lp10(model, subset_k, subset_p, epsilon):
    if len(subset_p) == 0 or len(subset_k) == 0:
        return

    scaling = 1e5

    # Create LP-10 problem of Fastcore
    prob = pulp.LpProblem('Fastcore-LP-10', pulp.LpMinimize)

    reaction_list = sorted(model.reaction_set)
    compound_list = sorted(model.compound_set)
    reaction_map = dict((rxnid, i) for i, rxnid in enumerate(reaction_list))
    compound_map = dict((cpdid, i) for i, cpdid in enumerate(compound_list))

    # Define flux variables
    max_bound = 0
    fluxes = []
    for rxnid in reaction_list:
        lower, upper = model.limits[rxnid]
        for b in (lower, upper):
            if abs(b) >= max_bound:
                max_bound = abs(b)
        var = pulp.LpVariable('v_'+rxnid, lower*scaling, upper*scaling)
        fluxes.append(var)

    # Define z variables and constraints
    zs = []
    for rxnid in subset_p:
        var = pulp.LpVariable('z_'+rxnid, 0, max_bound*scaling)
        zs.append(var)

        rxnindex = reaction_map[rxnid]
        prob += fluxes[rxnindex] >= -var
        prob += fluxes[rxnindex] <= var

    # Define constraints
    for rxnid in subset_k:
        rxnindex = reaction_map[rxnid]
        prob += fluxes[rxnindex] >= epsilon*scaling

    massbalance_lhs = [0 for cpdid in compound_list]
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        cpdindex = compound_map[cpdid]
        rxnindex = reaction_map[rxnid]
        massbalance_lhs[cpdindex] += value*fluxes[rxnindex]

    for lhs in massbalance_lhs:
        prob += lhs == 0

    # Define objective
    objective = sum(zs)
    prob += objective
    print prob

    # Solve
    status = prob.solve()
    if pulp.LpStatus[status] != 'Optimal':
        raise Exception('Non-optimal solution: {}'.format(pulp.LpStatus[status]))

    for rxnid, flux in izip(reaction_list, fluxes):
        print '{}\t{}'.format(rxnid, flux.value())
        yield rxnid, flux.value()

def support_set(fluxiter, threshold):
    '''Return the support set of the fluxes

    The fluxiter argument must be an iterable returning
    reaction, flux-pairs. The support set is the set of
    reactions that have a flux value above a certain
    threshold.'''
    return set(rxnid for rxnid, v in fluxiter if abs(v) >= threshold)

def fastcc(model, epsilon):
    '''Find the largest consistent subset of the model

    Returns a set containing the reactions in this subset
    or an empty set if no such subset can be found.'''

    consistent_subset = set()

    subset = model.reaction_set - model.reversible
    print '|J| = {}'.format(len(subset))
    print 'J = {}'.format(subset)

    support = support_set(fastcore_lp7(model, subset, epsilon), 0.99*epsilon)
    consistent_subset |= support
    print '|A| = {}'.format(len(consistent_subset))
    print 'A = {}'.format(consistent_subset)

    inconsistent_subset = subset - support

    if len(inconsistent_subset) > 0:
        print 'Subset inconsistent: {}'.format(inconsistent_subset)

    subset = (model.reaction_set - support) - inconsistent_subset
    print '|J| = {}'.format(len(subset))
    print 'J = {}'.format(subset)

    flipped = False
    singleton = False
    while len(subset) > 0:
        if singleton:
            subset_i = set((next(iter(subset)),))
            print 'LP3 on {}'.format(subset_i)
            support = support_set(fastcore_lp3(model, subset_i), 0.99*epsilon)
        else:
            subset_i = subset
            print 'LP7 on {}'.format(subset_i)
            support = support_set(fastcore_lp7(model, subset_i, epsilon), 0.99*epsilon)
        consistent_subset |= support
        print '|A| = {}'.format(len(consistent_subset))
        print 'A = {}'.format(consistent_subset)

        if not subset.isdisjoint(consistent_subset):
            subset -= consistent_subset
            print '|J| = {}'.format(len(subset))
            print 'J = {}'.format(subset)
            flipped = False
        else:
            subset_rev_i = subset_i & model.reversible
            if flipped or len(subset_rev_i) == 0:
                flipped = False
                if singleton:
                    subset -= subset_rev_i
                    print 'Reversible reaction inconsistent: {}'.format(subset_rev_i)
                else:
                    singleton = True
            else:
                model = model.flipped(subset_rev_i)
                flipped = True
                print 'Flip'

    return consistent_subset

def find_sparse_mode(model, subset, additional, singleton, epsilon):
    if len(subset) == 0:
        return set()

    if singleton:
        subset_i = set((next(iter(subset)),))
        print 'LP7 on {}'.format(subset_i)
        support = support_set(fastcore_lp7(model, subset_i, epsilon), 0.99*epsilon)
    else:
        print 'LP7 on {}'.format(subset)
        support = support_set(fastcore_lp7(model, subset, epsilon), 0.99*epsilon)

    print 'Support = {}'.format(support)
    k = subset & support

    if len(k) == 0:
        return set()

    print 'K = {}'.format(k)
    return support_set(fastcore_lp10(model, k, additional, epsilon), 0.99*epsilon)
