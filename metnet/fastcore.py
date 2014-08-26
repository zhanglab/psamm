
'''Fastcore module implementing the fastcore algorithm

This is an implementation of the algorithms described in
Vlassis, Nikos, Maria Pires Pacheco, and Thomas Sauter.
"Fast reconstruction of compact context-specific metabolic
network models." PLoS computational biology 10.1 (2014):
e1003424.'''

import sys
import cplex

from . import lpsolver
from .fluxanalysis import flux_balance
from .metabolicmodel import FlipableModelView


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

    def lp7(self, model, reaction_subset, epsilon):
        # This is similar to FBA but approximately maximizing the number
        # of reactions in subset with flux > epsilon, instead of just
        # maximizing the flux of one particular reaction. LP7 prefers
        # "flux splitting" over "flux concentrating".

        # Create LP-7 problem of Fastcore
        prob = self._solver.create_problem()

        # Define flux variables
        for rxnid in model.reaction_set:
            lower, upper = model.limits[rxnid]
            prob.define('v_'+rxnid, lower=lower, upper=upper)

        # Define z variables
        prob.define(*('z_'+rxnid for rxnid in reaction_subset), lower=0, upper=epsilon)
        prob.set_linear_objective(sum(prob.var('z_'+rxnid) for rxnid in reaction_subset))
        v = prob.set('v_'+rxnid for rxnid in reaction_subset)
        z = prob.set('z_'+rxnid for rxnid in reaction_subset)
        prob.add_linear_constraints(v >= z)

        massbalance_lhs = { compound: 0 for compound in model.compound_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[compound] += prob.var('v_'+rxnid) * value
        prob.add_linear_constraints(*(lhs == 0 for compound, lhs in massbalance_lhs.iteritems()))

        # Solve
        prob.solve(lpsolver.CplexProblem.Maximize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        for rxnid in sorted(model.reaction_set):
            yield rxnid, prob.get_value('v_'+rxnid)

    def lp10(self, model, subset_k, subset_p, epsilon):
        # This program forces reactions in subset K to attain flux > epsilon
        # while minimizing the sum of absolute flux values for reactions
        # in subset P (L1-regularization).

        if len(subset_k) == 0:
            return

        scaling = 1e10

        # Create LP-10 problem of Fastcore
        prob = self._solver.create_problem()

        # Define flux variables
        max_penalty_bound = 0
        for rxnid in model.reaction_set:
            lower, upper = model.limits[rxnid]
            if rxnid in subset_k:
                lower = max(lower, epsilon)
            if rxnid in subset_p:
                for b in (lower, upper):
                    if abs(b) > max_penalty_bound:
                        max_penalty_bound = abs(b)
            prob.define('v_'+rxnid, lower=lower*scaling, upper=upper*scaling)

        # Define z variables
        prob.define(*('z_'+rxnid for rxnid in subset_p), lower=0, upper=max_penalty_bound*scaling)
        prob.set_linear_objective(sum(prob.var('z_'+rxnid) for rxnid in subset_p))
        z = prob.set('z_'+rxnid for rxnid in subset_p)
        v = prob.set('v_'+rxnid for rxnid in subset_p)
        prob.add_linear_constraints(z >= v, v >= -z)

        massbalance_lhs = { compound: 0 for compound in model.compound_set }
        for spec, value in model.matrix.iteritems():
            compound, rxnid = spec
            massbalance_lhs[compound] += prob.var('v_'+rxnid) * value
        prob.add_linear_constraints(*(lhs == 0 for compound, lhs in massbalance_lhs.iteritems()))

        # Solve
        prob.solve(lpsolver.CplexProblem.Minimize)
        status = prob.cplex.solution.get_status()
        if status != 1:
            raise Exception('Non-optimal solution: {}'.format(prob.cplex.solution.get_status_string()))

        for rxnid in sorted(model.reaction_set):
            yield rxnid, prob.get_value('v_'+rxnid)

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

        # Wrap model in flipable proxy so reactions can be flipped
        model = FlipableModelView(model)

        flipped = False
        singleton = False
        while len(subset) > 0:
            if singleton:
                reaction = next(iter(subset))
                subset_i = { reaction }
                #print 'LP3 on {}'.format(subset_i)
                support = support_set(flux_balance(model, reaction, self._solver), epsilon)
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
                    model.flip(subset_rev_i)
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

        # Wrap model in flipable proxy so reactions can be flipped
        model = FlipableModelView(model)

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
                    model.flip(subset_rev_i)
                    flipped = True
                    #print 'Going to flipped state... {}'.format(subset_rev_i)

        return consistent_subset
