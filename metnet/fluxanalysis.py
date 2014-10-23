
'''Implementation of Flux Balance Analysis'''

from .metabolicmodel import FlipableModelView
from . import lpsolver

def flux_balance(model, reaction='Biomass', solver=lpsolver.CplexSolver()):
    '''Maximize the flux of a specific reaction

    Yields tuples of reaction id and flux for all model reactions.'''

    # Create Flux balance problem
    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reactions:
        lower, upper = model.limits[rxnid]
        prob.define('v_'+rxnid, lower=lower, upper=upper)

    objective = prob.var('v_'+reaction)
    prob.set_linear_objective(objective)

    # Define constraints
    massbalance_lhs = { compound: 0 for compound in model.compounds }
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

    for rxnid in model.reactions:
        yield rxnid, prob.get_value('v_'+rxnid)

def flux_balance_td(model, reaction, solver=lpsolver.CplexSolver()):
    '''Maximize the flux of a specific reaction with thermodynamic constraints

    Like FBA, but with the additional constraint that the flux must be
    thermodynamically feasible. This is solved as a MILP problem and the
    problem has been shown to be NP-hard in general.

    Described by Muller, Arne C., and Alexander Bockmayr. "Fast thermodynamically
    constrained flux variability analysis." Bioinformatics (2013): btt059.'''

    em = 1e5
    epsilon = 1e-5

    # Create Flux balance problem
    prob = solver.create_problem()

    # Define flux variables
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        prob.define('v_'+reaction_id, lower=lower, upper=upper)

        # Constrain internal reactions to a direction determined
        # by alpha.
        if not model.is_exchange(reaction_id):
            prob.define('alpha_'+reaction_id, types=lpsolver.CplexProblem.Binary)
            prob.define('dmu_'+reaction_id) # Delta mu

            flux = prob.var('v_'+reaction_id)
            alpha = prob.var('alpha_'+reaction_id)
            dmu = prob.var('dmu_'+reaction_id)

            prob.add_linear_constraints(flux >= lower*(1 - alpha),
                                        flux <= upper*alpha,
                                        dmu >= -em*alpha + epsilon,
                                        dmu <= em*(1 - alpha) - epsilon)

    objective = prob.var('v_'+reaction)
    prob.set_linear_objective(objective)

    # Define mu variables
    prob.define(*('mu_'+compound.id for compound in model.compounds))

    # Define constraints
    massbalance_lhs = { compound: 0 for compound in model.compounds }
    tdbalance_lhs = { reaction_id: 0 for reaction_id in model.reactions }
    for spec, value in model.matrix.iteritems():
        compound, reaction_id = spec
        massbalance_lhs[compound] += prob.var('v_'+reaction_id) * value
        if not model.is_exchange(reaction_id):
            tdbalance_lhs[reaction_id] += prob.var('mu_'+compound.id) * value
    for compound, lhs in massbalance_lhs.iteritems():
        prob.add_linear_constraints(lhs == 0)
    for reaction_id, lhs in tdbalance_lhs.iteritems():
        if not model.is_exchange(reaction_id):
            prob.add_linear_constraints(lhs == prob.var('dmu_'+reaction_id))

    # Solve
    prob.solve(lpsolver.CplexProblem.Maximize)
    status = prob.cplex.solution.get_status()
    if status != 101:
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    for reaction_id in model.reactions:
        yield reaction_id, prob.get_value('v_'+reaction_id)

def flux_variability(model, reactions, fixed, solver=lpsolver.CplexSolver()):
    '''Find the variability of each reaction while fixing certain fluxes

    Yields the reaction id, and a tuple of minimum and maximum value for each
    of the given reactions. The fixed reactions are given in a dictionary as
    a reaction id to value mapping.'''

    prob = solver.create_problem()

    # Define flux variables
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        if reaction_id in fixed:
            lower = max(lower, fixed[reaction_id])
        prob.define('v_'+reaction_id, lower=lower, upper=upper)

    massbalance_lhs = { compound: 0 for compound in model.compounds }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += value * prob.var('v_'+rxnid)
    for compound, lhs in massbalance_lhs.iteritems():
        prob.add_linear_constraints(lhs == 0)

    def min_max_solve(reaction_id):
        for direction in (lpsolver.CplexProblem.Minimize, lpsolver.CplexProblem.Maximize):
            prob.solve(direction)
            status = prob.cplex.solution.get_status()
            if status != 1:
                raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))
            yield prob.get_value('v_'+reaction_id)

    # Solve for each reaction
    for reaction_id in reactions:
        prob.set_linear_objective(prob.var('v_'+reaction_id))
        yield reaction_id, tuple(min_max_solve(reaction_id))

def flux_minimization(model, fixed, weights={}, solver=lpsolver.CplexSolver()):
    '''Minimize flux of all reactions while keeping certain fluxes fixed

    The fixed reactions are given in a dictionary as reaction id
    to value mapping. The weighted L1-norm of the fluxes is minimized.'''

    prob = solver.create_problem()

    # Define flux variables
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        if reaction_id in fixed:
            lower = max(lower, fixed[reaction_id])
        prob.define('v_'+reaction_id, lower=lower, upper=upper)

    # Define z variables
    prob.define(*('z_'+reaction_id for reaction_id in model.reactions), lower=0)
    objective = sum(prob.var('z_'+reaction_id) * weights.get(reaction_id, 1) for reaction_id in model.reactions)
    prob.set_linear_objective(objective)

    # Define constraints
    v = prob.set('v_'+rxnid for rxnid in model.reactions)
    z = prob.set('z_'+rxnid for rxnid in model.reactions)
    prob.add_linear_constraints(z >= v, v >= -z)

    massbalance_lhs = { compound: 0 for compound in model.compounds }
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

    return ((reaction_id, prob.get_value('v_'+reaction_id)) for reaction_id in model.reactions)

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
