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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import logging

from ..command import Command, SolverCommandMixin
from .. import massconsistency

logger = logging.getLogger(__name__)


class MassConsistencyCommand(SolverCommandMixin, Command):
    """Check whether the model is mass consistent."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', action='append', type=str,
            default=[], help='Exclude reaction from mass consistency')
        parser.add_argument(
            '--epsilon', type=float, help='Mass threshold',
            default=1e-5)
        parser.add_argument(
            '--checked', metavar='reaction', action='append', type=str,
            default=[], help='Mark reaction as already checked (no residual)')
        super(MassConsistencyCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        compound_name = {}
        zeromass = set()
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']
            elif compound.id not in compound_name:
                compound_name[compound.id] = compound.id

            if compound.properties.get('zeromass', False):
                zeromass.add(compound.id)

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        # Other reactions to exclude from consistency check
        exclude = set(self._args.exclude)

        # Add biomass reaction to be excluded
        biomass_reaction = self._model.get_biomass_reaction()
        if biomass_reaction is not None:
            exclude.add(biomass_reaction)

        # Create set of checked reactions
        checked = set(self._args.checked)

        solver = self._get_solver()

        known_inconsistent = exclude | exchange

        logger.info('Mass consistency on database')
        epsilon = self._args.epsilon
        compound_iter = massconsistency.check_compound_consistency(
            self._mm, solver, known_inconsistent, zeromass)

        good = 0
        total = 0
        for compound, mass in sorted(compound_iter, key=lambda x: (x[1], x[0]),
                                     reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print('{}: {}'.format(compound.translate(
                lambda x: compound_name.get(x, x)), mass))
        logger.info('Consistent compounds: {}/{}'.format(good, total))

        good = 0
        total = 0
        reaction_iter, compound_iter = (
            massconsistency.check_reaction_consistency(
                self._mm, solver=solver, exchange=known_inconsistent,
                checked=checked, zeromass=zeromass))
        for reaction_id, residual in sorted(
                reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            total += 1
            if abs(residual) >= epsilon:
                reaction = self._mm.get_reaction(reaction_id)
                rxt = reaction.translated_compounds(
                    lambda x: compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction_id, residual, rxt))
            else:
                good += 1
        logger.info('Consistent reactions: {}/{}'.format(good, total))
