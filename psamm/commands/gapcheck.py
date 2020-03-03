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
# Copyright 2014-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

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
            '--method', choices=['gapfind', 'prodcheck', 'sinkcheck'],
            default='prodcheck',
            help='Method to use for gap checking (default: prodcheck)')
        parser.add_argument(
            '--no-implicit-sinks', action='store_true',
            help='Do not include implicit sinks when gap-filling')
        parser.add_argument(
            '--unrestricted-exchange', action='store_true',
            help='Remove all limits on exchange reactions while gap checking')
        parser.add_argument(
            '--exclude-extracellular', action='store_true',
            help='Exclude extracellular compounds from the result')
        parser.add_argument(
            '--epsilon', type=float, default=1e-5,
            help='Threshold for compound production')

        super(GapCheckCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        def compound_name(id):
            if id not in self._model.compounds:
                return id
            return self._model.compounds[id].properties.get('name', id)

        extracellular_comp = self._model.extracellular_compartment
        epsilon = self._args.epsilon
        v_max = float(self._model.default_flux_limit)

        implicit_sinks = not self._args.no_implicit_sinks

        if self._args.unrestricted_exchange:
            # Allow all exchange reactions with no flux limits
            for reaction in self._mm.reactions:
                if self._mm.is_exchange(reaction):
                    del self._mm.limits[reaction].bounds

        logger.info('Searching for blocked compounds...')

        if self._args.method == 'gapfind':
            # Run GapFind on model
            solver = self._get_solver(integer=True)
            try:
                blocked = sorted(gapfind(
                    self._mm, solver=solver, epsilon=epsilon, v_max=v_max,
                    implicit_sinks=implicit_sinks))
            except GapFillError as e:
                self._log_epsilon_and_fail(epsilon, e)
        elif self._args.method == 'prodcheck':
            # Run production check on model
            solver = self._get_solver()
            blocked = self.run_reaction_production_check(
                self._mm, solver=solver, threshold=self._args.epsilon,
                implicit_sinks=implicit_sinks)
        elif self._args.method == 'sinkcheck':
            # Run sink check on model
            solver = self._get_solver()
            blocked = self.run_sink_check(
                self._mm, solver=solver, threshold=self._args.epsilon,
                implicit_sinks=implicit_sinks)
        else:
            self.argument_error('Invalid method: {}'.format(self._args.method))

        # Show result
        count = 0
        for compound in blocked:
            if (self._args.exclude_extracellular and
                    compound.compartment == extracellular_comp):
                continue

            count += 1
            print('{}\t{}'.format(compound, compound_name(compound.name)))

        logger.info('Blocked compounds: {}'.format(count))

    def run_sink_check(self, model, solver, threshold, implicit_sinks=True):
        """Run sink production check method."""
        prob = solver.create_problem()

        # Create flux variables
        v = prob.namespace()
        for reaction_id in model.reactions:
            lower, upper = model.limits[reaction_id]
            v.define([reaction_id], lower=lower, upper=upper)

        # Build mass balance constraints
        massbalance_lhs = {compound: 0 for compound in model.compounds}
        for spec, value in iteritems(model.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value

        mass_balance_constrs = {}
        for compound, lhs in iteritems(massbalance_lhs):
            if implicit_sinks:
                # The constraint is merely >0 meaning that we have implicit
                # sinks for all compounds.
                prob.add_linear_constraints(lhs >= 0)
            else:
                # Save these constraints so we can temporarily remove them
                # to create a sink.
                c, = prob.add_linear_constraints(lhs == 0)
                mass_balance_constrs[compound] = c

        for compound, lhs in sorted(iteritems(massbalance_lhs)):
            if not implicit_sinks:
                mass_balance_constrs[compound].delete()

            prob.set_objective(lhs)
            try:
                result = prob.solve(lp.ObjectiveSense.Maximize)
            except lp.SolverError as e:
                logger.warning('Failed to solve for compound: {} ({})'.format(
                    compound, e))

            if result.get_value(lhs) < threshold:
                yield compound

            if not implicit_sinks:
                # Restore mass balance constraint.
                c, = prob.add_linear_constraints(lhs == 0)
                mass_balance_constrs[compound] = c

    def run_reaction_production_check(self, model, solver, threshold,
                                      implicit_sinks=True):
        """Run reaction production check method."""
        prob = solver.create_problem()

        # Create flux variables
        v = prob.namespace()
        for reaction_id in model.reactions:
            lower, upper = model.limits[reaction_id]
            v.define([reaction_id], lower=lower, upper=upper)

        # Build mass balance constraints
        massbalance_lhs = {compound: 0 for compound in model.compounds}
        for spec, value in iteritems(model.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value

        # Create production variables and apply constraints
        for compound, lhs in iteritems(massbalance_lhs):
            if implicit_sinks:
                # The constraint is merely >0 meaning that we have implicit
                # sinks for all compounds.
                prob.add_linear_constraints(lhs >= 0)
            else:
                prob.add_linear_constraints(lhs == 0)

        confirmed_production = set()
        for reaction in model.reactions:
            if all(c in confirmed_production for c, _ in
                   model.get_reaction_values(reaction)):
                continue

            prob.set_objective(v(reaction))
            for sense in (lp.ObjectiveSense.Maximize,
                          lp.ObjectiveSense.Minimize):
                try:
                    result = prob.solve(sense)
                except lp.SolverError as e:
                    self.fail(
                        'Failed to solve for compound, reaction: {}, {}:'
                        ' {}'.format(compound, reaction, e))

                flux = result.get_value(v(reaction))
                for compound, value in model.get_reaction_values(reaction):
                    if compound in confirmed_production:
                        continue

                    production = 0
                    if sense == lp.ObjectiveSense.Maximize and flux > 0:
                        production = float(value) * flux
                    elif sense == lp.ObjectiveSense.Minimize and flux < 0:
                        production = float(value) * flux

                    if production >= threshold:
                        confirmed_production.add(compound)

        for compound in sorted(model.compounds):
            if compound not in confirmed_production:
                yield compound

    def _log_epsilon_and_fail(self, epsilon, exc):
        msg = ('Finding blocked compounds failed with epsilon set to {}. Try'
               ' lowering the epsilon value to reduce artifical constraints'
               ' on the model.'.format(epsilon))
        self.fail(msg, exc)
