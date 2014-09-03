
'''Implementation of Flux Balance Analysis'''

from .metabolicmodel import FlipableModelView
from . import lpsolver

def cpdid_str(compound):
    cpdid, comp = compound
    if comp is None:
        return cpdid
    return cpdid+'_'+comp

def flux_balance(model, reaction='Biomass', solver=lpsolver.CplexSolver()):
    '''Maximize the flux of a specific reaction

    Yields tuples of reaction id and flux for all model reactions.'''

    # Create Flux balance problem
    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        prob.define('v_'+rxnid, lower=lower, upper=upper)

    objective = prob.var('v_'+reaction)
    prob.set_linear_objective(objective)

    # Define constraints
    massbalance_lhs = { compound: 0 for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += prob.var('v_'+rxnid) * value
    for compound, lhs in massbalance_lhs.iteritems():
        prob.add_linear_constraints(lhs == 0)

    # Solve
    prob.solve(lpsolver.CplexProblem.Maximize)
    status = prob.cplex.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    for rxnid in model.reaction_set:
        yield rxnid, prob.get_value('v_'+rxnid)

def flux_minimization(model, fixed, weights={}, solver=lpsolver.CplexSolver()):
    '''Minimize flux of all reactions while keeping certain fluxes fixed

    The fixed reactions are given in a dictionary as reaction id
    to value mapping. The weighted L1-norm of the fluxes is minimized.'''

    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        if rxnid in fixed:
            lower = max(lower, fixed[rxnid])
        prob.define('v_'+rxnid, lower=lower, upper=upper)

    # Define z variables
    prob.define(*('z_'+rxnid for rxnid in model.reaction_set), lower=0)
    objective = sum(prob.var('z_'+rxnid) * weights.get(rxnid, 1) for rxnid in model.reaction_set)
    prob.set_linear_objective(objective)

    # Define constraints
    v = prob.set('v_'+rxnid for rxnid in model.reaction_set)
    z = prob.set('z_'+rxnid for rxnid in model.reaction_set)
    prob.add_linear_constraints(z >= v, v >= -z)

    massbalance_lhs = { compound: 0 for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += value * prob.var('v_'+rxnid)
    for compound, lhs in massbalance_lhs.iteritems():
        prob.add_linear_constraints(lhs == 0)

    # Solve
    prob.solve(lpsolver.CplexProblem.Minimize)
    status = prob.cplex.solution.get_status()
    if status != 1:
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    return ((rxnid, prob.get_value('v_'+rxnid)) for rxnid in model.reaction_set)

def naive_consistency_check(model, subset, epsilon, solver=lpsolver.CplexSolver()):
    '''Check that reaction subset of model is consistent using FBA

    A reaction is consistent if there is at least one flux solution
    to the model that both respects the compound balance of the
    stoichiometric matrix and also allows the reaction to have non-zero
    flux.

    The naive way to check this is to run FBA on each reaction in turn
    and checking whether the flux in the solution is non-zero. Since FBA
    only tries to maximize the flux (and the flux can be negative for
    reversible reactions), we have to try to both maximize and minimize
    the flux. An optimization to this method is implemented such that if
    checking one reaction results in flux in another unchecked reaction,
    that reaction will immediately be marked flux consistent.'''

    # Wrap model in flipable proxy so reactions can be flipped
    model = FlipableModelView(model)

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
            model.flip({ reaction })
            fluxiter = flux_balance(model, reaction, solver)
            support = set(rxnid for rxnid, v in fluxiter if abs(v) >= epsilon)
            subset -= support
            if reaction in support:
                continue
        print '{} not consistent!'.format(reaction)
        yield reaction
        subset.remove(reaction)
