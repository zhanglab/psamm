
'''Identify blocked metabolites and possible reconstructions

This implements a variant of the algorithms described by Kumar,
Vinay Satish, Madhukar S. Dasika, and Costas D. Maranas.
"Optimization based automated curation of metabolic reconstructions."'''

from . import lpsolver

def gapfind(model, epsilon=1e-5, v_max=1000, solver=lpsolver.CplexSolver()):
    '''Identify compounds in the model that cannot be produced

    Yields all compounds that cannot be produced. This method
    assumes implicit sinks for all compounds in the model so
    the only factor that influences whether a compound can be
    produced is the presence of the compounds needed to produce it.

    Epsilon indicates the threshold amount of a compound produced
    for it to not be considered blocked. V_max indicates the
    maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.'''

    prob = solver.create_problem()

    # Define flux variables
    for rxnid in model.reaction_set:
        lower, upper = model.limits[rxnid]
        prob.define('v_'+rxnid, lower=lower, upper=upper)

    # Define constraints on production of metabolites in reaction
    binary_cons_lhs = { compound: 0 for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        if value != 0 and (rxnid in model.reversible or value > 0):
            prob.define('w_'+rxnid+'_'+compound.id, types=lpsolver.CplexProblem.Binary)

            w = prob.var('w_'+rxnid+'_'+compound.id)
            sv = float(value) * prob.var('v_'+rxnid)

            prob.add_linear_constraints(sv <= v_max*w)
            if rxnid in model.reversible:
                prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            else:
                prob.add_linear_constraints(sv >= epsilon*w)

            binary_cons_lhs[compound] += w

    prob.define(*('xp_'+compound.id for compound in model.compound_set), types=lpsolver.CplexProblem.Binary)

    objective = sum(prob.var('xp_'+compound.id) for compound in model.compound_set)
    prob.set_linear_objective(objective)

    for compound, lhs in binary_cons_lhs.iteritems():
        prob.add_linear_constraints(lhs >= prob.var('xp_'+compound.id))

    # Define mass balance constraints
    massbalance_lhs = { compound: 0 for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += prob.var('v_'+rxnid) * value
    for compound, lhs in massbalance_lhs.iteritems():
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    prob.solve(lpsolver.CplexProblem.Maximize)
    status = prob.cplex.solution.get_status()
    if status != 101: # MILP solution
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    for compound in model.compound_set:
        if prob.get_value('xp_'+compound.id) == 0:
            yield compound

def gapfill(model, core, blocked, epsilon=1e-5, v_max=1000, solver=lpsolver.CplexSolver()):
    '''Find a set of reactions to add such that no compounds are blocked

    Returns two iterators: first an iterator of reactions not in
    core, that were added to resolve the model. Second, an
    iterator of reaction in core, that were reversed. Similarly to
    GapFind, this method assumes implicit sinks for all compounds in
    the model so the only factor that influences whether a compound
    can be produced is the presence of the compounds needed to produce
    it. This means that the resulting model will not necessarily be
    flux consistent.

    Core indicates the core set of reactions in the model. GapFill will
    minimize the number of added reactions that are not in core. Blocked
    indicates the set of compounds to be resolved. Epsilon indicates the
    threshold amount of a compound produced for it to not be considered
    blocked. V_max indicates the maximum flux.

    This method is implemented as a MILP-program. Therefore it may
    not be efficient for larger models.'''

    prob = solver.create_problem()

    # Define flux variables
    prob.define(*('v_'+rxnid for rxnid in model.reaction_set), lower=-v_max, upper=v_max)

    # Add binary indicator variables
    database_reactions = set(model.reaction_set).difference(core)
    prob.define(*('ym_'+rxnid for rxnid in core), types=lpsolver.CplexProblem.Binary)
    prob.define(*('yd_'+rxnid for rxnid in database_reactions), types=lpsolver.CplexProblem.Binary)

    objective = sum(prob.var('ym_'+rxnid) for rxnid in core) + sum(prob.var('yd_'+rxnid) for rxnid in database_reactions)
    prob.set_linear_objective(objective)

    # Add constraints on core reactions
    for rxnid in core:
        v = prob.var('v_'+rxnid)
        if rxnid in model.reversible:
            prob.add_linear_constraints(v >= model.limits[rxnid].lower)
        else:
            prob.add_linear_constraints(v >= -v_max*prob.var('ym_'+rxnid))
        prob.add_linear_constraints(v <= model.limits[rxnid].upper)

    # Add constraints on database reactions
    for rxnid in database_reactions:
        v = prob.var('v_'+rxnid)
        yd = prob.var('yd_'+rxnid)
        prob.add_linear_constraints(v >= yd * model.limits[rxnid].lower)
        prob.add_linear_constraints(v <= yd * model.limits[rxnid].upper)

    # Define constraints on production of blocked metabolites in reaction
    binary_cons_lhs = { compound: 0 for compound in blocked }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        if compound in blocked and value != 0:
            prob.define('w_'+rxnid+'_'+compound.id, types=lpsolver.CplexProblem.Binary)

            w = prob.var('w_'+rxnid+'_'+compound.id)
            sv = float(value) * prob.var('v_'+rxnid)

            prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            prob.add_linear_constraints(sv <= v_max*w)

            if rxnid in model.reversible or value > 0:
                binary_cons_lhs[compound] += w

    for compound, lhs in binary_cons_lhs.iteritems():
        if compound in blocked:
            prob.add_linear_constraints(lhs >= 1)

    # Define mass balance constraints
    massbalance_lhs = { compound: 0 for compound in model.compound_set }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += prob.var('v_'+rxnid) * value
    for compound, lhs in massbalance_lhs.iteritems():
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    prob.solve(lpsolver.CplexProblem.Minimize)
    status = prob.cplex.solution.get_status()
    if status != 101: # MILP solution
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    def added_iter():
        for rxnid in database_reactions:
            if prob.get_value('yd_'+rxnid) > 0:
                yield rxnid

    def reversed_iter():
        for rxnid in core:
            if prob.get_value('ym_'+rxnid) > 0:
                yield rxnid

    return added_iter(), reversed_iter()
