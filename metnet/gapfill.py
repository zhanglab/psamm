
'''Identify blocked metabolites and possible reconstructions

This implements a variant of the algorithms described by Kumar,
Vinay Satish, Madhukar S. Dasika, and Costas D. Maranas.
"Optimization based automated curation of metabolic reconstructions."'''

from . import lpsolver

def cpdid_str(compound):
    cpdid, comp = compound
    if comp is None:
        return cpdid
    return cpdid+'_'+comp

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
            prob.define('w_'+rxnid+'_'+cpdid_str(compound), types=lpsolver.CplexProblem.Binary)

            w = prob.var('w_'+rxnid+'_'+cpdid_str(compound))
            sv = float(value) * prob.var('v_'+rxnid)

            prob.add_linear_constraints(sv <= v_max*w)
            if rxnid in model.reversible:
                prob.add_linear_constraints(sv >= epsilon-v_max*(1 - w))
            else:
                prob.add_linear_constraints(sv >= epsilon*w)

            binary_cons_lhs[compound] += w

    prob.define(*('xp_'+cpdid_str(compound) for compound in model.compound_set), types=lpsolver.CplexProblem.Binary)

    objective = sum(prob.var('xp_'+cpdid_str(compound)) for compound in model.compound_set)
    prob.set_linear_objective(objective)

    for compound, lhs in binary_cons_lhs.iteritems():
        prob.add_linear_constraints(lhs >= prob.var('xp_'+cpdid_str(compound)))

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
        if prob.get_value('xp_'+cpdid_str(compound)) == 0:
            yield compound

def gapfill(model, core, blocked, epsilon=1e-5, v_max=1000, solver=lpsolver.CplexSolver()):
    pass
