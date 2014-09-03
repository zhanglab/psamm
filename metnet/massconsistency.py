
'''Mass consistency analysis of metabolic models

A stoichiometric matrix, S, is said to be mass-consistent if
S^Tm = 0 has a positive solution (m_i > 0). This corresponds to
assigning a positive mass to each compound in the stoichiometric
matrix and having each reaction preserve mass. Exchange reactions
will have to be excluded from this check, as they are not able to
preserve mass (by definition). In addition some models may contain
pseudo-compounds (e.g. "photon") that also has to be excluded.'''

import cplex

from . import lpsolver

class MassConsistencyCheck(object):
    def __init__(self, solver=lpsolver.CplexSolver()):
        self._solver = solver

    def _cpdid_str(self, compound):
        cpdid, comp = compound
        if comp is None:
            return cpdid
        return cpdid+'_'+comp

    def _cplex_add_compound_mass(self, prob, model, zeromass=set(), lower=1):
        '''Add variables for compound mass'''
        for compound in model.compound_set:
            prob.define('m_'+self._cpdid_str(compound), lower=(0 if compound[0] in zeromass else lower))

    def _cplex_constrain_identical(self, prob, model):
        '''Constrain identical compounds in different compartments to the same mass'''
        compound_id_constr = []
        for compound in model.compound_set:
            if compound[1] is not None and (compound[0], None) in model.compound_set:
                mass_c = prob.var('m_'+compound[0])
                mass_other = prob.var('m_'+self._cpdid_str(compound))
                prob.add_linear_constraints(mass_c == mass_other)

    def is_consistent(self, model, exchange=set(), zeromass=set()):
        '''Try to assign a positive mass to each compound and return True on success

        The masses are simply constrained by m_i > 1 and finding a solution
        under these conditions proves that the model is mass consistent.'''

        prob = self._solver.create_problem()

        # Define mass variables
        self._cplex_add_compound_mass(prob, model, zeromass)
        self._cplex_constrain_identical(prob, model)
        prob.set_linear_objective(sum(prob.var('m_'+self._cpdid_str(compound)) for compound in model.compound_set))

        # Define constraints
        massbalance_lhs = { rxnid: 0 for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            cpdid, rxnid = spec
            massbalance_lhs[rxnid] += prob.var('m_'+self._cpdid_str(cpdid)) * value
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                prob.add_linear_constraints(lhs == 0)

        prob.solve(lpsolver.CplexProblem.Minimize)
        status = prob.cplex.solution.get_status()

        return status == 1

    def check_reaction_consistency(self, model, exchange=set(), zeromass=set(), weights={}):
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

        # Define mass variables
        self._cplex_add_compound_mass(prob, model, zeromass)
        self._cplex_constrain_identical(prob, model)

        # Define residual mass variables and objective constriants
        prob.define(*('z_'+rxnid for rxnid in model.reaction_set), lower=0)
        prob.define(*('r_'+rxnid for rxnid in model.reaction_set))

        objective = sum(prob.var('z_'+rxnid) * weights.get(rxnid, 1) for rxnid in model.reaction_set)
        prob.set_linear_objective(objective)

        r = prob.set('r_'+rxnid for rxnid in model.reaction_set)
        z = prob.set('z_'+rxnid for rxnid in model.reaction_set)
        prob.add_linear_constraints(z >= r, r >= -z)

        massbalance_lhs = { rxnid: 0 for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[rxnid] += prob.var('m_'+self._cpdid_str(compound)) * value
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                r = prob.var('r_'+rxnid)
                prob.add_linear_constraints(lhs + r == 0)

        # Solve
        prob.solve(lpsolver.CplexProblem.Minimize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        def iterate_reactions():
            for rxnid in model.reaction_set:
                residual = prob.get_value('r_'+rxnid)
                yield rxnid, residual

        def iterate_compounds():
            for compound in model.compound_set:
                yield self._cpdid_str(compound), prob.get_value('m_'+self._cpdid_str(compound))

        return iterate_reactions(), iterate_compounds()

    def check_compound_consistency(self, model, exchange=set(), zeromass=set()):
        '''Yield each compound in the model with assigned mass

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

        # Define mass variables
        self._cplex_add_compound_mass(prob, model, zeromass, 0)
        self._cplex_constrain_identical(prob, model)

        # Define z variables
        prob.define(*('z_'+self._cpdid_str(compound) for compound in model.compound_set), lower=0, upper=1)
        prob.set_linear_objective(sum(prob.var('z_'+self._cpdid_str(compound)) for compound in model.compound_set))

        z = prob.set('z_'+self._cpdid_str(compound) for compound in model.compound_set)
        m = prob.set('m_'+self._cpdid_str(compound) for compound in model.compound_set)
        prob.add_linear_constraints(m >= z)

        massbalance_lhs = { rxnid: 0 for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            cpdid, rxnid = spec
            massbalance_lhs[rxnid] += prob.var('m_'+self._cpdid_str(cpdid)) * value
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                prob.add_linear_constraints(lhs == 0)

        # Solve
        prob.solve(lpsolver.CplexProblem.Maximize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        for compound in model.compound_set:
            yield self._cpdid_str(compound), prob.get_value('m_'+self._cpdid_str(compound))
