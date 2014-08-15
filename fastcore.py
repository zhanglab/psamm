
'''Fastcore module implementing the fastcore algorithm

This is an implementation of the algorithms described in
Vlassis, Nikos, Maria Pires Pacheco, and Thomas Sauter.
"Fast reconstruction of compact context-specific metabolic
network models." PLoS computational biology 10.1 (2014):
e1003424.'''

import sys
import lpsolver
import cplex

def support_set(fluxiter, threshold):
    '''Return the support set of the fluxes

    The fluxiter argument must be an iterable returning
    reaction, flux-pairs. The support set is the set of
    reactions that have an absolute flux value above a
    certain threshold.'''
    return set(rxnid for rxnid, v in fluxiter if abs(v) >= threshold)

def support_set_positive(fluxiter, threshold):
    '''Return the support set of the fluxes (only positive fluxes)

    The fluxiter argument must be an iterable returning
    reaction, flux-pairs. The support set is the set of
    reactions that have a flux value above a certain
    threshold.'''
    return set(rxnid for rxnid, v in fluxiter if v >= threshold)


class Fastcore(object):
    '''Fastcore computation object containing reference to the solver'''

    def __init__(self, solver=lpsolver.CplexSolver()):
        self._solver = solver

    def _cpdid_str(self, compound):
        '''Create string identifier for compound with compartment'''
        cpdid, comp = compound
        if comp is None:
            return cpdid
        return cpdid+'_'+comp

    def lp3(self, model, reaction_subset):
        # This is in fact just a standard FBA, maximizing the flux of
        # the reactions in the subset.

        # Create LP-3 problem of Fastcore
        prob = self._solver.create_problem()

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
        massbalance_lhs = { compound: [] for compound in model.compound_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[compound].append(('v_'+rxnid, float(value)))
        for compound, lhs in massbalance_lhs.iteritems():
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'],
                                        rhs=[0], names=['massbalance_'+self._cpdid_str(compound)])

        # Solve
        prob.objective.set_sense(prob.objective.sense.maximize)
        prob.solve()
        status = prob.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

        for rxnid in model.reaction_set:
            yield rxnid, prob.solution.get_values('v_'+rxnid)

    def lp7(self, model, reaction_subset, epsilon):
        # This is similar to FBA but approximately maximizing the number
        # of reactions in subset with flux > epsilon, instead of just
        # maximizing the flux of one particular reaction. LP7 prefers
        # "flux splitting" over "flux concentrating".

        # Create LP-7 problem of Fastcore
        prob = self._solver.create_problem()

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
        prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(1, -1)) for rxnid in reaction_subset],
                                    senses=['G']*len(reaction_subset), rhs=[0]*len(reaction_subset))

        massbalance_lhs = { compound: [] for compound in model.compound_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[compound].append(('v_'+rxnid, float(value)))
        for compound, lhs in massbalance_lhs.iteritems():
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'],
                                        rhs=[0], names=['massbalance_'+self._cpdid_str(compound)])

        # Solve
        prob.objective.set_sense(prob.objective.sense.maximize)
        prob.solve()
        status = prob.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

        for rxnid in sorted(model.reaction_set):
            yield rxnid, prob.solution.get_values('v_'+rxnid)

    def lp10(self, model, subset_k, subset_p, epsilon):
        # This program forces reactions in subset K to attain flux > epsilon
        # while minimizing the sum of absolute flux values for reactions
        # in subset P (L1-regularization).

        if len(subset_k) == 0:
            return

        scaling = 1e5

        # Create LP-10 problem of Fastcore
        prob = self._solver.create_problem()

        # Define flux variables
        flux_names = []
        flux_lower = []
        flux_upper = []
        max_penalty_bound = 0
        for rxnid in model.reaction_set:
            lower, upper = model.limits[rxnid]
            if rxnid in subset_k:
                lower = max(lower, epsilon)
            if rxnid in subset_p:
                for b in (lower, upper):
                    if abs(b) > max_penalty_bound:
                        max_penalty_bound = abs(b)
            flux_names.append('v_'+rxnid)
            flux_lower.append(lower*scaling)
            flux_upper.append(upper*scaling)
        prob.variables.add(names=flux_names, lb=flux_lower, ub=flux_upper)

        # Define z variables
        zs_names = []
        for rxnid in subset_p:
            zs_names.append('z_'+rxnid)
        prob.variables.add(names=zs_names, lb=[0]*len(zs_names), ub=[max_penalty_bound*scaling]*len(zs_names), obj=[1]*len(zs_names))

        # Define constraints
        for rxnid in subset_p:
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(1, 1)), # v >= -z
                                                  cplex.SparsePair(ind=('v_'+rxnid, 'z_'+rxnid), val=(-1, 1))], # z >= v
                                        senses=['G', 'G'], rhs=[0, 0])

        massbalance_lhs = { compound: [] for compound in model.compound_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[compound].append(('v_'+rxnid, float(value)))
        for compound, lhs in massbalance_lhs.iteritems():
            prob.linear_constraints.add(lin_expr=[cplex.SparsePair(*zip(*lhs))], senses=['E'],
                                        rhs=[0], names=['massbalance_'+self._cpdid_str(compound)])

        # Solve
        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.solve()
        status = prob.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.solution.get_status_string()))

        for rxnid in sorted(model.reaction_set):
            yield rxnid, prob.solution.get_values('v_'+rxnid)

    def fastcc(self, model, epsilon):
        '''Check consistency of model reactions

        Yields all reactions in the model that are not part
        of the consistent subset.'''

        consistent_subset = set()

        subset = model.reaction_set - model.reversible
        #print '|J| = {}'.format(len(subset))
        #print 'J = {}'.format(subset)

        support = support_set(self.lp7(model, subset, epsilon), epsilon)
        consistent_subset |= support
        #print '|A| = {}'.format(len(consistent_subset))
        #print 'A = {}'.format(consistent_subset)

        inconsistent_subset = subset - support
        for reaction in inconsistent_subset:
            yield reaction

        # Check remaining reactions
        subset = (model.reaction_set - support) - inconsistent_subset
        #print '|J| = {}'.format(len(subset))
        #print 'J = {}'.format(subset)

        flipped = False
        singleton = False
        while len(subset) > 0:
            if singleton:
                subset_i = set((next(iter(subset)),))
                #print 'LP3 on {}'.format(subset_i)
                support = support_set(self.lp3(model, subset_i), epsilon)
            else:
                subset_i = subset
                #print 'LP7 on {}'.format(subset_i)
                support = support_set(self.lp7(model, subset_i, epsilon), epsilon)
            consistent_subset |= support
            #print '|A| = {}'.format(len(consistent_subset))
            #print 'A = {}'.format(consistent_subset)

            if not subset.isdisjoint(consistent_subset):
                subset -= consistent_subset
                #print '|J| = {}'.format(len(subset))
                #print 'J = {}'.format(subset)
                flipped = False
            else:
                # TODO: irreversible reactions are taken care of before the
                # loop so at this point all reactions in subset_i are reversble(?).
                subset_rev_i = subset_i & model.reversible
                if flipped or len(subset_rev_i) == 0:
                    flipped = False
                    if singleton:
                        subset -= subset_rev_i
                        for reaction in subset_rev_i:
                            yield reaction
                    else:
                        singleton = True
                else:
                    model = model.flipped(subset_rev_i)
                    flipped = True
                    #print 'Flip'

    def fastcc_is_consistent(self, model, epsilon):
        '''Quickly check whether model is consistent

        Returns true if the model is consistent. If it is only necessary
        to know whether a model is consistent, this function is faster
        as it will return the result as soon as it finds a single
        inconsistent reaction.'''
        for reaction in self.fastcc(model, epsilon):
            return False
        return True

    def fastcc_consistent_subset(self, model, epsilon):
        '''Return consistent subset of model

        The largest consistent subset is returned as
        a set of reaction names.'''
        return model.reaction_set - set(self.fastcc(model, epsilon))

    def find_sparse_mode(self, model, core, additional, epsilon):
        '''Find the support of a sparse mode containing the core subset.'''

        if len(core) == 0:
            return set()

        #print 'LP7 on {}'.format(core)
        support = support_set_positive(self.lp7(model, core, epsilon), epsilon)

        #print '|Support| = {}'.format(len(support))
        k = core & support

        if len(k) == 0:
            return set()

        #print '|K| = {}'.format(len(k))
        #print 'LP10 on K ({}) and P ({})'.format(len(k), len(additional))
        return support_set(self.lp10(model, k, additional, epsilon), epsilon)

    def fastcore(self, model, core, epsilon):
        '''Find a flux consistent subnetwork containing the core subset

        The result will contain the core subset and as few of the additional
        reactions as possible.'''

        consistent_subset = set()

        subset = core - model.reversible
        #print '|J| = {}, J = {}'.format(len(subset), subset)
        #print '|J| = {}'.format(len(subset))

        penalty_set = model.reaction_set - core
        #print '|P| = {}, P = {}'.format(len(penalty_set), penalty_set)
        #print '|P| = {}'.format(len(penalty_set))

        mode = self.find_sparse_mode(model, subset, penalty_set, epsilon)
        if len(subset - mode) > 0:
            raise Exception('Inconsistent irreversible core reactions: {}'.format(subset - mode))

        consistent_subset |= mode
        #print '|A| = {}, A = {}'.format(len(consistent_subset), consistent_subset)
        #print '|A| = {}'.format(len(consistent_subset))

        subset = core - mode
        #print '|J| = {}, J = {}'.format(len(subset), subset)
        #print '|J| = {}'.format(len(subset))

        flipped = False
        singleton = False
        while len(subset) > 0:
            penalty_set -= consistent_subset
            if singleton:
                subset_i = set((next(iter(subset)),))
            else:
                subset_i = subset

            mode = self.find_sparse_mode(model, subset_i, penalty_set, epsilon)
            consistent_subset |= mode
            #print '|A| = {}, A = {}'.format(len(consistent_subset), consistent_subset)
            #print '|A| = {}'.format(len(consistent_subset))

            if not subset.isdisjoint(consistent_subset):
                #print 'Subset improved {} -> {}'.format(len(subset), len(subset - consistent_subset))
                subset -= consistent_subset
                #print '|J| = {}, J = {}'.format(len(subset), subset)
                flipped = False
            else:
                #print 'Nothing found, changing state...'
                subset_rev_i = subset_i & model.reversible
                if flipped or len(subset_rev_i) == 0:
                    if singleton:
                        raise Exception('Global network inconsistent: {}'.format(subset_rev_i))

                    #print 'Going to non-flipped, singleton state...'
                    singleton = True
                    flipped = False
                else:
                    model = model.flipped(subset_rev_i)
                    flipped = True
                    #print 'Going to flipped state... {}'.format(subset_rev_i)

        return consistent_subset
