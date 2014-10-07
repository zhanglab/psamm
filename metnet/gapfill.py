
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
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        prob.define('v_'+reaction_id, lower=lower, upper=upper)

    # Define constraints on production of metabolites in reaction
    binary_cons_lhs = { compound: 0 for compound in model.compounds }
    for spec, value in model.matrix.iteritems():
        compound, reaction_id = spec
        if value != 0 and (reaction_id in model.reversible or value > 0):
            prob.define('w_'+reaction_id+'_'+compound.id, types=lpsolver.CplexProblem.Binary)

            w = prob.var('w_'+reaction_id+'_'+compound.id)
            sv = float(value) * prob.var('v_'+reaction_id)

            prob.add_linear_constraints(sv <= v_max*w)
            if model.is_reversible(reaction_id):
                prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            else:
                prob.add_linear_constraints(sv >= epsilon*w)

            binary_cons_lhs[compound] += w

    prob.define(*('xp_'+compound.id for compound in model.compounds), types=lpsolver.CplexProblem.Binary)

    objective = sum(prob.var('xp_'+compound.id) for compound in model.compounds)
    prob.set_linear_objective(objective)

    for compound, lhs in binary_cons_lhs.iteritems():
        prob.add_linear_constraints(lhs >= prob.var('xp_'+compound.id))

    # Define mass balance constraints
    massbalance_lhs = { compound: 0 for compound in model.compounds }
    for spec, value in model.matrix.iteritems():
        compound, reaction_id = spec
        massbalance_lhs[compound] += prob.var('v_'+reaction_id) * value
    for compound, lhs in massbalance_lhs.iteritems():
        # The constraint is merely >0 meaning that we have implicit sinks
        # for all compounds.
        prob.add_linear_constraints(lhs >= 0)

    # Solve
    prob.solve(lpsolver.CplexProblem.Maximize)
    status = prob.cplex.solution.get_status()
    if status != 101: # MILP solution
        raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

    for compound in model.compounds:
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
    prob.define(*('v_'+reaction_id for reaction_id in model.reactions), lower=-v_max, upper=v_max)

    # Add binary indicator variables
    database_reactions = set(model.reactions).difference(core)
    prob.define(*('ym_'+reaction_id for reaction_id in core), types=lpsolver.CplexProblem.Binary)
    prob.define(*('yd_'+reaction_id for reaction_id in database_reactions), types=lpsolver.CplexProblem.Binary)

    objective = (sum(prob.var('ym_'+reaction_id) for reaction_id in core) +
                    sum(prob.var('yd_'+reaction_id) for reaction_id in database_reactions))
    prob.set_linear_objective(objective)

    # Add constraints on core reactions
    for reaction_id in core:
        v = prob.var('v_'+reaction_id)
        if model.is_reversible(reaction_id):
            prob.add_linear_constraints(v >= model.limits[reaction_id].lower)
        else:
            prob.add_linear_constraints(v >= -v_max*prob.var('ym_'+reaction_id))
        prob.add_linear_constraints(v <= model.limits[reaction_id].upper)

    # Add constraints on database reactions
    for reaction_id in database_reactions:
        v = prob.var('v_'+reaction_id)
        yd = prob.var('yd_'+reaction_id)
        prob.add_linear_constraints(v >= yd * model.limits[reaction_id].lower)
        prob.add_linear_constraints(v <= yd * model.limits[reaction_id].upper)

    # Define constraints on production of blocked metabolites in reaction
    binary_cons_lhs = { compound: 0 for compound in blocked }
    for spec, value in model.matrix.iteritems():
        compound, reaction_id = spec
        if compound in blocked and value != 0:
            prob.define('w_'+reaction_id+'_'+compound.id, types=lpsolver.CplexProblem.Binary)

            w = prob.var('w_'+reaction_id+'_'+compound.id)
            sv = float(value) * prob.var('v_'+reaction_id)

            prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            prob.add_linear_constraints(sv <= v_max*w)

            if reaction_id in model.reversible or value > 0:
                binary_cons_lhs[compound] += w

    for compound, lhs in binary_cons_lhs.iteritems():
        if compound in blocked:
            prob.add_linear_constraints(lhs >= 1)

    # Define mass balance constraints
    massbalance_lhs = { compound: 0 for compound in model.compounds }
    for spec, value in model.matrix.iteritems():
        compound, reaction_id = spec
        massbalance_lhs[compound] += prob.var('v_'+reaction_id) * value
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
        for reaction_id in database_reactions:
            if prob.get_value('yd_'+reaction_id) > 0:
                yield reaction_id

    def reversed_iter():
        for reaction_id in core:
            if prob.get_value('ym_'+reaction_id) > 0:
                yield reaction_id

    return added_iter(), reversed_iter()
