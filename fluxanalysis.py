
'''Implementation of Flux Balance Analysis'''

import cplex
import lpsolver

def cpdid_str(compound):
    cpdid, comp = compound
    if comp is None:
        return cpdid
    return cpdid+'_'+comp

def flux_balance(model, reaction='Biomass', solver=lpsolver.CplexSolver()):
    # Create Flux balance problem
    prob = solver.create_problem()

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
        flux_obj.append(1 if rxnid == reaction else 0)
    prob.variables.add(names=flux_names, lb=flux_lower, ub=flux_upper, obj=flux_obj)

    # Define constraints
    massbalance_lhs = { compound: [] for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound].append(('v_'+rxnid, float(value)))
    for compound, lhs in massbalance_lhs.iteritems():
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))],
                                    senses=['E'], rhs=[0], names=['massbalance_'+cpdid_str(compound)])

    # Solve
    prob.objective.set_sense(prob.objective.sense.maximize)
    prob.solve()
    status = prob.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

    for rxnid in model.reaction_set:
        yield rxnid, prob.solution.get_values('v_'+rxnid)

def naive_consistency_check(model, subset, epsilon, solver=lpsolver.CplexSolver()):
    subset = set(subset)
    while len(subset) > 0:
        reaction = next(iter(subset))
        print '{} left, checking {}...'.format(len(subset), reaction)
        fluxiter = flux_balance(model, reaction, solver)
        support = set(rxnid for rxnid, v in fluxiter if abs(v) >= epsilon)
        subset -= support
        if reaction in support:
            continue
        elif reaction in model.reversible:
            model2 = model.flipped({ reaction })
            fluxiter = flux_balance(model2, reaction, solver)
            support = set(rxnid for rxnid, v in fluxiter if abs(v) >= epsilon)
            subset -= support
            if reaction in support:
                continue
        print '{} not consistent!'.format(reaction)
        yield reaction