
'''Mass consistency analysis of metabolic databases

A stoichiometric matrix, S, is said to be mass-consistent if
S^Tm = 0 has a positive solution (m_i > 0). This corresponds to
assigning a positive mass to each compound in the stoichiometric
matrix and having each reaction preserve mass. Exchange reactions
will have to be excluded from this check, as they are not able to
preserve mass (by definition). In addition some databases may contain
pseudo-compounds (e.g. "photon") that also has to be excluded.'''

import cplex

from . import lpsolver

class MassConsistencyCheck(object):
    def __init__(self, solver=lpsolver.CplexSolver()):
        self._solver = solver

    def _non_localized_compounds(self, database):
        '''Return set of non-localized compounds in database (i.e. in compartment None)'''
        return set(compound.in_compartment(None) for compound in database.compounds)

    def is_consistent(self, database, exchange=set(), zeromass=set()):
        '''Try to assign a positive mass to each compound and return True on success

        The masses are simply constrained by m_i > 1 and finding a solution
        under these conditions proves that the database is mass consistent.'''

        prob = self._solver.create_problem()
        compound_set = self._non_localized_compounds(database)

        # Define mass variables
        for compound in compound_set:
            prob.define(('m', compound), lower=(0 if compound in zeromass else 1))
        prob.set_linear_objective(sum(prob.var(('m', compound)) for compound in compound_set))

        # Define constraints
        massbalance_lhs = { reaction: 0 for reaction in database.reactions }
        for spec, value in database.matrix.iteritems():
            compound, reaction = spec
            massbalance_lhs[reaction] += prob.var(('m', compound.in_compartment(None))) * value
        for reaction, lhs in massbalance_lhs.iteritems():
            if reaction not in exchange:
                prob.add_linear_constraints(lhs == 0)

        prob.solve(lpsolver.CplexProblem.Minimize)
        status = prob.cplex.solution.get_status()

        return status == 1

    def check_reaction_consistency(self, database, exchange=set(), zeromass=set(), weights={}):
        '''Check inconsistent reactions by minimizing mass residuals for each reaction

        Returns a reaction iterable, and compound iterable. The
        reaction iterable yields reaction ids and mass residuals.
        The compound iterable yields compound ids and mass assignments.

        Each compound is assigned a mass of at least one, and the
        masses are balanced using the stoichiometric matrix. In addition,
        each reaction has a residual mass that is included in the mass
        balance equations. The L1-norm of the residuals is minimized.'''

        # Create Flux balance problem
        prob = self._solver.create_problem()
        compound_set = self._non_localized_compounds(database)

        # Define mass variables
        for compound in compound_set:
            prob.define(('m', compound), lower=(0 if compound in zeromass else 1))

        # Define residual mass variables and objective constriants
        prob.define(*(('z', reaction_id) for reaction_id in database.reactions), lower=0)
        prob.define(*(('r', reaction_id) for reaction_id in database.reactions))

        objective = sum(prob.var(('z', reaction_id)) * weights.get(reaction_id, 1) for reaction_id in database.reactions)
        prob.set_linear_objective(objective)

        r = prob.set(('r', reaction_id) for reaction_id in database.reactions)
        z = prob.set(('z', reaction_id) for reaction_id in database.reactions)
        prob.add_linear_constraints(z >= r, r >= -z)

        massbalance_lhs = { reaction_id: 0 for reaction_id in database.reactions }
        for spec, value in database.matrix.iteritems():
            compound, reaction_id = spec
            massbalance_lhs[reaction_id] += prob.var(('m', compound.in_compartment(None))) * value
        for reaction_id, lhs in massbalance_lhs.iteritems():
            if reaction_id not in exchange:
                r = prob.var(('r', reaction_id)) # residual
                prob.add_linear_constraints(lhs + r == 0)

        # Solve
        prob.solve(lpsolver.CplexProblem.Minimize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        def iterate_reactions():
            for reaction_id in database.reactions:
                residual = prob.get_value(('r', reaction_id))
                yield reaction_id, residual

        def iterate_compounds():
            for compound in compounds:
                yield compound, prob.get_value(('m', compound))

        return iterate_reactions(), iterate_compounds()

    def check_compound_consistency(self, database, exchange=set(), zeromass=set()):
        '''Yield each compound in the database with assigned mass

        Each compound will be assigned a mass and the number of compounds
        having a positive mass will be approximately maximized.

        This is an implementation of the solution originally proposed by
        Albert Gevorgyan, Mark G. Poolman and David A. Fell, "Detection of
        stoichiometric inconsistencies in biomolecular models", but using the
        new method proposed by Nikos Vlassis and Ronan Fleming,
        "fastGapFill: efficient gap filling in metabolic networks" to avoid
        MILP constraints. This is similar to the way Fastcore avoids MILP
        contraints.'''

        # Create mass balance problem
        prob = self._solver.create_problem()
        compound_set = self._non_localized_compounds(database)

        # Define mass variables
        prob.define(*(('m', compound) for compound in compound_set), lower=0)

        # Define z variables
        prob.define(*(('z', compound) for compound in compound_set), lower=0, upper=1)
        prob.set_linear_objective(sum(prob.var(('z', compound)) for compound in compound_set))

        z = prob.set(('z', compound) for compound in compound_set)
        m = prob.set(('m', compound) for compound in compound_set)
        prob.add_linear_constraints(m >= z)

        massbalance_lhs = { reaction_id: 0 for reaction_id in database.reactions }
        for spec, value in database.matrix.iteritems():
            compound, reaction_id = spec
            massbalance_lhs[reaction_id] += prob.var(('m', compound.in_compartment(None))) * value
        for reaction_id, lhs in massbalance_lhs.iteritems():
            if reaction_id not in exchange:
                prob.add_linear_constraints(lhs == 0)

        # Solve
        prob.solve(lpsolver.CplexProblem.Maximize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        for compound in compound_set:
            yield compound, prob.get_value(('m', compound))
