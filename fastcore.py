
'''Fastcore module implementing the fastcore algorithm

This is an implementation of the algorithms described in
Vlassis, Nikos, Maria Pires Pacheco, and Thomas Sauter.
"Fast reconstruction of compact context-specific metabolic
network models." PLoS computational biology 10.1 (2014):
e1003424.'''

import pulp
import cplex
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
        yield rxnid, flux.value()

def fastcore_lp3_cplex(model, reaction_subset):
    # Create LP-3 problem of Fastcore
    prob = cplex.Cplex()

    # Define flux variables
    flux_names = []
    flux_lower = []
    flux_upper = []
    flux_obj = []
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        flux_names.append('v_'+rxnid)
        flux_lower.append(lower)
        flux_upper.append(upper)
        flux_obj.append(1 if rxnid in reaction_subset else 0)
    prob.variables.add(names=flux_names, lb=flux_lower, ub=flux_upper, obj=flux_obj)

    # Define constraints
    massbalance_lhs = { cpdid: [] for cpdid in model.compound_set }
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        massbalance_lhs[cpdid].append(('v_'+rxnid, value))
    for cpdid, lhs in massbalance_lhs.iteritems():
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'], rhs=[0], names=['massbalance_'+cpdid])

    # Solve
    prob.objective.set_sense(prob.objective.sense.maximize)
    prob.solve()
    status = prob.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

    for rxnid in model.reaction_set:
        yield rxnid, prob.solution.get_values('v_'+rxnid)

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
        yield rxnid, flux.value()

def fastcore_lp7_cplex(model, reaction_subset, epsilon):
    # Create LP-7 problem of Fastcore
    prob = cplex.Cplex()

    # Define flux variables
    flux_names = []
    flux_lower = []
    flux_upper = []
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        flux_names.append('v_'+rxnid)
        flux_lower.append(lower)
        flux_upper.append(upper)
    prob.variables.add(names=flux_names, lb=flux_lower, ub=flux_upper)

    # Define z variables
    zs_names = []
    for rxnid in reaction_subset:
        zs_names.append('z_'+rxnid)
    prob.variables.add(names=zs_names, lb=[0]*len(zs_names), ub=[epsilon]*len(zs_names), obj=[1]*len(zs_names))

    # Define constraints
    for rxnid in reaction_subset:
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(1, -1))], senses=['G'], rhs=[0])

    massbalance_lhs = { cpdid: [] for cpdid in model.compound_set }
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        massbalance_lhs[cpdid].append(('v_'+rxnid, value))
    for cpdid, lhs in massbalance_lhs.iteritems():
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'], rhs=[0], names=['massbalance_'+cpdid])

    # Solve
    prob.objective.set_sense(prob.objective.sense.maximize)
    prob.solve()
    status = prob.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

    for rxnid in model.reaction_set:
        yield rxnid, prob.solution.get_values('v_'+rxnid)

def fastcore_lp10(model, subset_k, subset_p, epsilon):
    if len(subset_k) == 0:
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

    # Solve
    status = prob.solve()
    if pulp.LpStatus[status] != 'Optimal':
        raise Exception('Non-optimal solution: {}'.format(pulp.LpStatus[status]))

    for rxnid, flux in izip(reaction_list, fluxes):
        yield rxnid, flux.value()

def fastcore_lp10_cplex(model, subset_k, subset_p, epsilon):
    if len(subset_k) == 0:
        return

    scaling = 1e5

    # Create LP-10 problem of Fastcore
    prob = cplex.Cplex()

    # Define flux variables
    flux_names = []
    flux_lower = []
    flux_upper = []
    max_bound = 0
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        for b in (lower, upper):
            if abs(b) >= max_bound:
                max_bound = abs(b)
        flux_names.append('v_'+rxnid)
        flux_lower.append(lower*scaling)
        flux_upper.append(upper*scaling)
    prob.variables.add(names=flux_names, lb=flux_lower, ub=flux_upper)

    # Define z variables
    zs_names = []
    for rxnid in subset_p:
        zs_names.append('z_'+rxnid)
    prob.variables.add(names=zs_names, lb=[0]*len(zs_names), ub=[cplex.infinity]*len(zs_names), obj=[1]*len(zs_names))

    # Define constraints
    for rxnid in subset_p:
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(1, 1)), # v >= -z
                                                cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(-1, 1))], # z >= v
                                    senses=['G', 'G'], rhs=[0, 0])

    prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=['v_'+rxnid], val=[1]) for rxnid in subset_k],
                                senses=['G'] * len(subset_k), rhs=[epsilon*scaling] * len(subset_k))

    massbalance_lhs = { cpdid: [] for cpdid in model.compound_set }
    for spec, value in model.matrix.iteritems():
        cpdid, rxnid = spec
        massbalance_lhs[cpdid].append(('v_'+rxnid, value))
    for cpdid, lhs in massbalance_lhs.iteritems():
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'], rhs=[0], names=['massbalance_'+cpdid])

    # Solve
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.solve()
    status = prob.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

    for rxnid in model.reaction_set:
        yield rxnid, prob.solution.get_values('v_'+rxnid)

def support_set(fluxiter, threshold):
    '''Return the support set of the fluxes

    The fluxiter argument must be an iterable returning
    reaction, flux-pairs. The support set is the set of
    reactions that have a flux value above a certain
    threshold.'''
    return set(rxnid for rxnid, v in fluxiter if abs(v) >= threshold)

def fastcc(model, epsilon):
    '''Check consistency of model reactions

    Returns a set containing the reactions in the largest
    consistent subset of the model or an empty set if no such
    subset can be found.'''

    consistent_subset = set()

    subset = model.reaction_set - model.reversible
    #print '|J| = {}'.format(len(subset))
    #print 'J = {}'.format(subset)

    support = support_set(fastcore_lp7_cplex(model, subset, epsilon), 0.99*epsilon)
    consistent_subset |= support
    #print '|A| = {}'.format(len(consistent_subset))
    #print 'A = {}'.format(consistent_subset)

    inconsistent_subset = subset - support

    if len(inconsistent_subset) > 0:
        print 'Subset inconsistent: {}'.format(inconsistent_subset)

    subset = (model.reaction_set - support) - inconsistent_subset
    #print '|J| = {}'.format(len(subset))
    #print 'J = {}'.format(subset)

    flipped = False
    singleton = False
    while len(subset) > 0:
        if singleton:
            subset_i = set((next(iter(subset)),))
            #print 'LP3 on {}'.format(subset_i)
            support = support_set(fastcore_lp3_cplex(model, subset_i), 0.99*epsilon)
        else:
            subset_i = subset
            #print 'LP7 on {}'.format(subset_i)
            support = support_set(fastcore_lp7_cplex(model, subset_i, epsilon), 0.99*epsilon)
        consistent_subset |= support
        #print '|A| = {}'.format(len(consistent_subset))
        #print 'A = {}'.format(consistent_subset)

        if not subset.isdisjoint(consistent_subset):
            subset -= consistent_subset
            #print '|J| = {}'.format(len(subset))
            #print 'J = {}'.format(subset)
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
                #print 'Flip'

    return consistent_subset

def find_sparse_mode(model, core, additional, singleton, epsilon):
    '''Find the support of a sparse mode containing the core subset.'''

    if len(core) == 0:
        return set()

    if singleton:
        subset_i = set((next(iter(core)),))
        #print 'LP7 on {}'.format(subset_i)
        support = support_set(fastcore_lp7_cplex(model, subset_i, epsilon), 0.99*epsilon)
    else:
        #print 'LP7 on {}'.format(core)
        support = support_set(fastcore_lp7_cplex(model, core, epsilon), 0.99*epsilon)

    #print 'Support = {}'.format(support)
    k = core & support

    if len(k) == 0:
        return set()

    #print 'K = {}'.format(k)
    #print 'LP10 on K={} and P={}'.format(k, additional)
    return support_set(fastcore_lp10_cplex(model, k, additional, epsilon), 0.99*epsilon)

def fastcore(model, core, epsilon):
    '''Find a flux consistent subnetwork containing the core subset

    The result will contain the core subset and as few of the additional
    reactions as possible.'''

    consistent_subset = set()

    subset = core - model.reversible
    #print '|J| = {}, J = {}'.format(len(subset), subset)

    penalty_set = model.reaction_set - core
    #print '|P| = {}, P = {}'.format(len(penalty_set), penalty_set)

    mode = find_sparse_mode(model, subset, penalty_set, False, epsilon)
    if len(subset - mode) > 0:
        raise Exception('Inconsistent irreversible core reactions')

    consistent_subset |= mode
    #print '|A| = {}, A = {}'.format(len(consistent_subset), consistent_subset)

    subset = core - mode
    #print '|J| = {}, J = {}'.format(len(subset), subset)

    flipped = False
    singleton = False
    while len(subset) > 0:
        penalty_set -= consistent_subset
        mode = find_sparse_mode(model, subset, penalty_set, singleton, epsilon)
        consistent_subset |= mode
        #print '|A| = {}, A = {}'.format(len(consistent_subset), consistent_subset)

        if not subset.isdisjoint(consistent_subset):
            subset -= consistent_subset
            #print '|J| = {}, J = {}'.format(len(subset), subset)
            flipped = False
        else:
            if singleton:
                subset_i = set((next(iter(core)),))
                subset_rev_i = subset_i & model.reversible
            else:
                subset_rev_i = subset & model.reversible

            if flipped or len(subset_rev_i) == 0:
                flipped = False
                if singleton:
                    raise Exception('Global network inconsistent: {}'.format(subset_rev_i))
                else:
                    singleton = True
            else:
                model = model.flipped(subset_rev_i)
                flipped = True
                #print 'Flip'

    return consistent_subset