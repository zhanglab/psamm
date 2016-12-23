# This file is part of PSAMM.
#
# PSAMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PSAMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PSAMM.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014-2016  Jon Lund Steffensen <jon_steffensen@uri.edu>

from __future__ import unicode_literals

import logging
from six import iteritems

from ..command import Command, MetabolicMixin, SolverCommandMixin
from ..gapfill import gapfind, GapFillError
from ..lpsolver import lp

logger = logging.getLogger(__name__)


class GapCheckCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Check for compound production gaps in model."""
    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method', choices=['gapfind', 'prodcheck', 'pccheck'],
            help='Adds implicit sinks for all compounds',type=str, required=True)
        parser.add_argument(
            '--epsilon', type=float, default=1e-5,
            help='Threshold for compound production')

        super(GapCheckCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        compound_name = {}
        compound_list = []
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id
        for compound in self._mm.compounds:
            compound_list.append(compound)

        solver = self._get_solver(integer=True)
        extracellular_comp = self._model.extracellular_compartment
        epsilon = self._args.epsilon
        v_max = float(self._model.default_flux_limit)
        model = self._mm

        if self._args.method == 'gapfind':
            # Run GapFind on model
            logger.info('Searching for blocked compounds...')
            result = gapfind(self._mm, solver=solver, epsilon=epsilon, v_max=v_max)

            try:
                blocked = set(compound for compound in result
                              if compound.compartment != extracellular_comp)
            except GapFillError as e:
                self._log_epsilon_and_fail(epsilon, e)

            if len(blocked) > 0:
                logger.info('Blocked compounds: {}'.format(len(blocked)))
                for compound in sorted(blocked):
                    name = compound_name.get(compound.name, compound.name)
                    print('{}\t{}'.format(compound, name))
            else:
                logger.info('No blocked compounds found')
        if self._args.method == 'prodcheck':
            self.run_prodcheck(model, solver, compound_list, self._args.epsilon)
        if self._args.method == 'pccheck':
            self.run_prodconcheck(model, solver, compound_list, self._args.epsilon)

    def run_prodcheck(self, model, solver, compound_name, threshold):
        prob = solver.create_problem()

        v = prob.namespace()
        for reaction_id in model.reactions:
            lower, upper = model.limits[reaction_id]
            v.define([reaction_id], lower=lower, upper=upper)

        massbalance_lhs = {compound: 0 for compound in compound_name}
        for spec, value in iteritems(model.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value

        for compound, lhs in iteritems(massbalance_lhs):
            #This constraint results in implicit sinks being present for each compound
            prob.add_linear_constraints(lhs >= 0)

        blocked_dict = {}
        for cpd_id in compound_name:
            prob_testing = prob
            prob_testing.set_objective(massbalance_lhs.get(cpd_id))
            prob_testing.solve(lp.ObjectiveSense.Maximize)
            if prob_testing.result.get_value(massbalance_lhs.get(cpd_id)) >= 0 + threshold:
                blocked_dict[cpd_id] = 1
            else:
                blocked_dict[cpd_id] = 0
        for i, j in blocked_dict.iteritems():
            if j == 0:
                print(i)

    def run_prodconcheck(self, model, solver, compound_names, threshold):
        prob = solver.create_problem()

        v = prob.namespace()
        for reaction_id in model.reactions:
            lower, upper = model.limits[reaction_id]
            v.define([reaction_id], lower=lower, upper=upper)

        massbalance_lhs = {compound: 0 for compound in compound_names}
        for spec, value in iteritems(model.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value

        mass_balance_constrs = {}
        for compound, lhs in iteritems(massbalance_lhs):
            # The constraint is merely >0 meaning that we have implicit sinks
            # for all compounds.
            c, = prob.add_linear_constraints(lhs == 0)
            mass_balance_constrs[compound] = c

        blocked_dict = {}
        for cpd_id in compound_names:
            prob_testing = prob
            mass_balance_constrs.get(cpd_id).delete()

            prob_testing.set_objective(massbalance_lhs.get(cpd_id))
            prob_testing.solve(lp.ObjectiveSense.Maximize)
            if prob_testing.result.get_value(massbalance_lhs.get(cpd_id)) >= 0+threshold:
                blocked_dict[cpd_id] = 1
            else:
                blocked_dict[cpd_id] = 0
            prob.add_linear_constraints(massbalance_lhs.get(cpd_id) == 0)

        for i, j in blocked_dict.iteritems():
            if j == 0:
                print(i)

    def _log_epsilon_and_fail(self, epsilon, exc):
        msg = ('Finding blocked compounds failed with epsilon set to {}. Try'
               ' lowering the epsilon value to reduce artifical constraints'
               ' on the model.'.format(epsilon))
        self.fail(msg, exc)
