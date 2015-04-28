
"""Implementation of Flux Balance Analysis"""

import logging
import random

from .lpsolver import lp

# Module-level logging
logger = logging.getLogger(__name__)

class FluxBalanceError(Exception):
    '''Error indicating that a flux balance cannot be solved'''

class FluxBalanceProblem(object):
    '''Maximize the flux of a specific reaction

    The solve method solves the problem for a specific parameter. After solving, the
    flux can be obtained using the get_flux method.'''

    def __init__(self, model, solver):
        self._prob = solver.create_problem()

        # Define flux variables
        for reaction_id in model.reactions:
            lower, upper = model.limits[reaction_id]
            self._prob.define(('v', reaction_id), lower=lower, upper=upper)

        # Define constraints
        massbalance_lhs = { compound: 0 for compound in model.compounds }
        for spec, value in model.matrix.iteritems():
            compound, reaction_id = spec
            massbalance_lhs[compound] += self.get_flux_var(reaction_id) * value
        for compound, lhs in massbalance_lhs.iteritems():
            self._prob.add_linear_constraints(lhs == 0)

    @property
    def prob(self):
        """Return the underlying LP problem"""
        return self._prob

    def solve(self, reaction):
        '''Solve problem maximizing the given reaction

        If reaction is a dictionary object, each entry is interpreted as a weight on
        the objective for that reaction (non-existent reaction will have zero weight).'''

        if isinstance(reaction, dict):
            objective = sum(v * self.get_flux_var(r) for r, v in reaction.iteritems())
        else:
            objective = self.get_flux_var(reaction)

        # Set objective and solve
        self._prob.set_linear_objective(objective)
        result = self._prob.solve(lp.ObjectiveSense.Maximize)
        if not result:
            raise FluxBalanceError('Non-optimal solution: {}'.format(result.status))

    def get_flux_var(self, reaction):
        """Get LP variable representing the reaction flux"""
        return self._prob.var(('v', reaction))

    def get_flux(self, reaction):
        '''Get resulting flux value for reaction'''
        return self._prob.result.get_value(('v', reaction))

class FluxBalanceTDProblem(FluxBalanceProblem):
    """Maximize the flux of a specific reaction with thermodynamic constraints

    Like FBA, but with the additional constraint that the flux must be
    thermodynamically feasible. This is solved as a MILP problem and the
    problem has been shown to be NP-hard in general.

    Described in [Muller13]_.
    """

    def __init__(self, model, solver):
        super(FluxBalanceTDProblem, self).__init__(model, solver)
        p = self.prob

        em = 1e5
        epsilon = 1e-5

        for reaction_id in model.reactions:
            # Constrain internal reactions to a direction determined
            # by alpha.
            if not model.is_exchange(reaction_id):
                p.define(('alpha', reaction_id), types=lp.VariableType.Binary)
                p.define(('dmu', reaction_id)) # Delta mu

                flux = self.get_flux_var(reaction_id)
                alpha = p.var(('alpha', reaction_id))
                dmu = p.var(('dmu', reaction_id))

                lower, upper = model.limits[reaction_id]
                p.add_linear_constraints(flux >= lower*(1 - alpha),
                                            flux <= upper*alpha,
                                            dmu >= -em*alpha + epsilon,
                                            dmu <= em*(1 - alpha) - epsilon)

        # Define mu variables
        p.define(*(('mu', compound) for compound in model.compounds))

        tdbalance_lhs = { reaction_id: 0 for reaction_id in model.reactions }
        for spec, value in model.matrix.iteritems():
            compound, reaction_id = spec
            if not model.is_exchange(reaction_id):
                tdbalance_lhs[reaction_id] += p.var(('mu', compound)) * value
        for reaction_id, lhs in tdbalance_lhs.iteritems():
            if not model.is_exchange(reaction_id):
                p.add_linear_constraints(lhs == p.var(('dmu', reaction_id)))

def flux_balance(model, reaction, solver):
    '''Run flux balance analysis on the given model

    Yields the reaction id and flux value for each reaction in the model.

    This is a convenience function for sertting up and running the
    FluxBalanceProblem. If the FBA is solved for more than one parameter
    it is recommended to setup and reuse the FluxBalanceProblem manually
    for a speed up.'''

    fba = FluxBalanceProblem(model, solver=solver)
    fba.solve(reaction)
    for reaction in model.reactions:
        yield reaction, fba.get_flux(reaction)

def flux_balance_td(model, reaction, solver):
    '''Run flux balance analysis with thermodynamic constraints

    Yields the reaction id and flux value for each reaction in the model.

    This is a convenience function for sertting up and running the
    FluxBalanceTDProblem. If the tFBA is solved for more than one parameter
    it is recommended to setup and reuse the FluxBalanceTDProblem manually
    for a speed up.'''

    tfba = FluxBalanceTDProblem(model, solver=solver)
    tfba.solve(reaction)
    for reaction in model.reactions:
        yield reaction, tfba.get_flux(reaction)


def flux_variability(model, reactions, fixed, solver):
    """Find the variability of each reaction while fixing certain fluxes

    Yields the reaction id, and a tuple of minimum and maximum value for each
    of the given reactions. The fixed reactions are given in a dictionary as
    a reaction id to value mapping.
    """

    fba = FluxBalanceProblem(model, solver=solver)
    for reaction_id, value in fixed.iteritems():
        flux = fba.get_flux_var(reaction_id)
        fba.prob.add_linear_constraints(flux >= value)

    def min_max_solve(reaction_id):
        for direction in (-1, 1):
            fba.solve({ reaction_id: direction })
            yield fba.get_flux(reaction_id)

    # Solve for each reaction
    for reaction_id in reactions:
        yield reaction_id, tuple(min_max_solve(reaction_id))


def flux_minimization(model, fixed, solver, weights={}):
    '''Minimize flux of all reactions while keeping certain fluxes fixed

    The fixed reactions are given in a dictionary as reaction id
    to value mapping. The weighted L1-norm of the fluxes is minimized.'''

    prob = solver.create_problem()

    # Define flux variables
    for reaction_id in model.reactions:
        lower, upper = model.limits[reaction_id]
        if reaction_id in fixed:
            lower = max(lower, fixed[reaction_id])
        prob.define(('v', reaction_id), lower=lower, upper=upper)

    # Define z variables
    prob.define(*(('z', reaction_id) for reaction_id in model.reactions), lower=0)
    objective = sum(prob.var(('z', reaction_id)) * weights.get(reaction_id, 1) for reaction_id in model.reactions)
    prob.set_linear_objective(objective)

    # Define constraints
    v = prob.set(('v', rxnid) for rxnid in model.reactions)
    z = prob.set(('z', rxnid) for rxnid in model.reactions)
    prob.add_linear_constraints(z >= v, v >= -z)

    massbalance_lhs = { compound: 0 for compound in model.compounds }
    for spec, value in model.matrix.iteritems():
        compound, rxnid = spec
        massbalance_lhs[compound] += value * prob.var(('v', rxnid))
    for compound, lhs in massbalance_lhs.iteritems():
        prob.add_linear_constraints(lhs == 0)

    # Solve
    result = prob.solve(lp.ObjectiveSense.Minimize)
    if not result:
        raise FluxBalanceError('Non-optimal solution: {}'.format(result.status))

    return ((reaction_id, result.get_value(('v', reaction_id))) for reaction_id in model.reactions)

def flux_randomization(model, fixed, solver):
    """Find a uniformly random flux mode in the solution space

    This can be used to generate a random sample from the solution
    space. If many random samples are needed more efficient methods
    (e.g. random walk) are useful.

    The reactions in the fixed dictionary are constrained with the
    associated lower bound."""

    optimize = {}
    for reaction_id in model.reactions:
        if model.is_reversible(reaction_id):
            optimize[reaction_id] = 2*random.random() - 1.0
        else:
            optimize[reaction_id] = random.random()

    fba = FluxBalanceProblem(model, solver=solver)
    for reaction_id, value in fixed.iteritems():
        fba.prob.add_linear_constraints(fba.get_flux_var(reaction_id) >= value)

    fba.solve(optimize)
    for reaction_id in model.reactions:
        yield reaction_id, fba.get_flux(reaction_id)

def consistency_check(model, subset, epsilon, solver):
    """Check that reaction subset of model is consistent using FBA

    Yields all reactions that are *not* flux consistent. A reaction is
    consistent if there is at least one flux solution to the model that both
    respects the model constraints and also allows the reaction in question to
    have non-zero flux.

    This can be determined by running FBA on each reaction in turn
    and checking whether the flux in the solution is non-zero. Since FBA
    only tries to maximize the flux (and the flux can be negative for
    reversible reactions), we have to try to both maximize and minimize
    the flux. An optimization to this method is implemented such that if
    checking one reaction results in flux in another unchecked reaction,
    that reaction will immediately be marked flux consistent.
    """

    fba = FluxBalanceProblem(model, solver=solver)

    subset = set(subset)
    while len(subset) > 0:
        reaction = next(iter(subset))

        logger.debug('{} left, checking {}...'.format(len(subset), reaction))

        fba.solve(reaction)
        support = set(rxnid for rxnid in model.reactions if abs(fba.get_flux(rxnid)) >= epsilon)
        subset -= support
        if reaction in support:
            continue
        elif model.is_reversible(reaction):
            fba.solve({ reaction: -1 })
            support = set(rxnid for rxnid in model.reactions if abs(fba.get_flux(rxnid)) >= epsilon)
            subset -= support
            if reaction in support:
                continue

        logger.debug('{} not consistent!'.format(reaction))

        yield reaction
        subset.remove(reaction)
