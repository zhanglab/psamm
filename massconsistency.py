
'''Mass consistency analysis of metabolic models

A stoichiometric matrix, S, is said to be mass-consistent if
S^Tm = 0 has a positive solution (m_i > 0). This corresponds to
assigning a positive mass to each compound in the stoichiometric
matrix and having each reaction preserve mass. Exchange reactions
will have to be excluded from this check, as they are not able to
preserve mass (by definition). In addition some models may contain
pseudo-compounds (e.g. "photon") that also has to be excluded.'''

import lpsolver
import cplex

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
        mass_names = []
        mass_lower = []
        for compound in model.compound_set:
            mass_names.append('m_'+self._cpdid_str(compound))
            mass_lower.append(0 if compound[0] in zeromass else lower)
        prob.variables.add(names=mass_names, lb=mass_lower, ub=[cplex.infinity]*len(mass_names))

    def _cplex_constrain_identical(self, prob, model):
        '''Constrain identical compounds in different compartments to the same mass'''
        compound_id_constr = []
        for compound in model.compound_set:
            if compound[1] is not None and (compound[0], None) in model.compound_set:
                compound_id_constr.append(('m_'+compound[0], 'm_'+self._cpdid_str(compound)))
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=constr, val=(1, -1)) for constr in compound_id_constr],
                                    senses=['E']*len(compound_id_constr), rhs=[0]*len(compound_id_constr))

    def is_consistent(self, model, exchange=set(), zeromass=set()):
        '''Try to assign a positive mass to each compound and return True on success

        The masses are simply constrained by m_i > 1 and finding a solution
        under these conditions proves that the model is mass consistent.'''

        prob = self._solver.create_problem()

        # Define mass variables
        self._cplex_add_compound_mass(prob, model, zeromass)
        self._cplex_constrain_identical(prob, model)
        prob.objective.set_linear(('m_'+self._cpdid_str(compound), 1) for compound in model.compound_set)

        # Define constraints
        massbalance_lhs = { rxnid: [] for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            cpdid, rxnid = spec
            massbalance_lhs[rxnid].append(('m_'+self._cpdid_str(cpdid), float(value)))
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                ind, val = zip(*lhs)
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=ind, val=val)],
                                            senses=['E'], rhs=[0], names=['massbalance_'+rxnid])

        # Solve
        prob.objective.set_sense(prob.objective.sense.minimize)
        try:
            prob.solve()
            status = prob.solution.get_status()
        except cplex.exceptions.CplexSolveError as e:
            status = e.args[2]

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

        # Initialize default weights of the weighted L1-norm
        weights = { rxnid: weights.get(rxnid, 1) for rxnid in model.reaction_set }

        # Define residual mass variables and objective constriants
        rs_names = []
        zs_names = []
        for rxnid in model.reaction_set:
            rs_names.append('r_'+rxnid) # reaction mass residual
            zs_names.append('z_'+rxnid) # objective variable
        prob.variables.add(names=rs_names, lb=[-cplex.infinity]*len(rs_names), ub=[cplex.infinity]*len(rs_names))
        prob.variables.add(names=zs_names, lb=[0]*len(zs_names), ub=[cplex.infinity]*len(zs_names))
        prob.objective.set_linear(('z_'+rxnid, weights[rxnid]) for rxnid in model.reaction_set)

        # Define constraints
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('r_'+rxnid, 'z_'+rxnid),
                                                                val=(1, 1)) for rxnid in model.reaction_set],
                                    senses=['G']*len(model.reaction_set), rhs=[0]*len(model.reaction_set))
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('r_'+rxnid, 'z_'+rxnid),
                                                                val=(-1, 1)) for rxnid in model.reaction_set],
                                    senses=['G']*len(model.reaction_set), rhs=[0]*len(model.reaction_set))

        massbalance_lhs = { rxnid: [] for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[rxnid].append(('m_'+self._cpdid_str(compound), float(value)))
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                ind, val = zip(*(lhs + [('r_'+rxnid, 1)]))
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=ind, val=val)],
                                            senses=['E'], rhs=[0], names=['massbalance_'+rxnid])

        # Solve
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.solve()
        status = prob.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

        def iterate_reactions():
            for rxnid in model.reaction_set:
                residual = prob.solution.get_values('r_'+rxnid)
                yield rxnid, residual

        def iterate_compounds():
            for compound in model.compound_set:
                yield self._cpdid_str(compound), prob.solution.get_values('m_'+self._cpdid_str(compound))

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
        zs_names = []
        for cpdid in model.compound_set:
            zs_names.append('z_'+self._cpdid_str(cpdid))
        prob.variables.add(names=zs_names, lb=[0]*len(zs_names), ub=[1]*len(zs_names),
                            obj=[1]*len(zs_names))

        # Define constraints
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('m_'+self._cpdid_str(compound),
                                                                    'z_'+self._cpdid_str(compound)),
                                                                val=(1, -1)) for compound in model.compound_set],
                                    senses=['G']*len(model.compound_set), rhs=[0]*len(model.compound_set))

        massbalance_lhs = { rxnid: [] for rxnid in model.reaction_set }
        for spec, value in model.matrix.iteritems():
            cpdid, rxnid = spec
            massbalance_lhs[rxnid].append(('m_'+self._cpdid_str(cpdid), float(value)))
        for rxnid, lhs in massbalance_lhs.iteritems():
            if rxnid not in exchange:
                ind, val = zip(*lhs)
                prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=ind, val=val)],
                                            senses=['E'], rhs=[0], names=['massbalance_'+rxnid])

        # Solve
        prob.objective.set_sense(prob.objective.sense.maximize)
        prob.solve()
        status = prob.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

        for compound in model.compound_set:
            yield self._cpdid_str(compound), prob.solution.get_values('m_'+self._cpdid_str(compound))
