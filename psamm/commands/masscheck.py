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

from __future__ import unicode_literals

import time
import logging

from six import iteritems

from ..command import (Command, SolverCommandMixin,
                       MetabolicMixin, FilePrefixAppendAction)
from .. import massconsistency
from ..reaction import Compound

logger = logging.getLogger(__name__)


class MassConsistencyCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Check whether the model is mass consistent."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', action=FilePrefixAppendAction,
            type=str, default=[],
            help='Exclude reaction from mass consistency')
        parser.add_argument(
            '--epsilon', type=float, help='Mass threshold',
            default=1e-5)
        parser.add_argument(
            '--checked', metavar='reaction', action=FilePrefixAppendAction,
            type=str, default=[],
            help='Mark reaction as already checked (no residual)')
        parser.add_argument(
            '--type', choices=['compound', 'reaction'],
            default='compound', help='Type of check to perform')
        super(MassConsistencyCommand, cls).init_parser(parser)

    def run(self):
        # Load compound information
        self._compound_name = {}
        zeromass = set()
        for compound in self._model.parse_compounds():
            if 'name' in compound.properties:
                self._compound_name[compound.id] = compound.properties['name']
            elif compound.id not in self._compound_name:
                self._compound_name[compound.id] = compound.id

            if compound.properties.get('zeromass', False):
                zeromass.add(Compound(compound.id))

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in self._mm.reactions:
            if self._mm.is_exchange(reaction_id):
                exchange.add(reaction_id)

        # Create a set of excluded reactions
        exclude = set(self._args.exclude)

        # Add biomass reaction to be excluded
        biomass_reaction = self._model.biomass_reaction
        if biomass_reaction is not None:
            exclude.add(biomass_reaction)

        solver = self._get_solver()

        known_inconsistent = exclude | exchange

        if self._args.type == 'compound':
            self._check_compounds(known_inconsistent, zeromass, solver)
        elif self._args.type == 'reaction':
            self._check_reactions(known_inconsistent, zeromass, solver)
        else:
            self.argument_error(
                'Invalid type of check: {}'.format(self._args.type))

    def _check_compounds(self, known_inconsistent, zeromass, solver):
        logger.info('Checking stoichiometric consistency of compounds...')

        epsilon = self._args.epsilon

        start_time = time.time()

        masses = dict(massconsistency.check_compound_consistency(
            self._mm, solver, known_inconsistent, zeromass))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        good = 0
        total = 0
        for compound, mass in sorted(
                iteritems(masses), key=lambda x: (x[1], x[0]), reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print('{}\t{}\t{}'.format(
                compound, mass, compound.translate(
                    lambda x: self._compound_name.get(x, x))))
        logger.info('Consistent compounds: {}/{}'.format(good, total))

    def _check_reactions(self, known_inconsistent, zeromass, solver):
        logger.info('Checking stoichiometric consistency of reactions...')

        # Create a set of checked reactions
        checked = set(self._args.checked)

        epsilon = self._args.epsilon

        start_time = time.time()

        reaction_iter, compound_iter = (
            massconsistency.check_reaction_consistency(
                self._mm, solver=solver, exchange=known_inconsistent,
                checked=checked, zeromass=zeromass))

        logger.info('Solving took {:.2f} seconds'.format(
            time.time() - start_time))

        good = 0
        total = 0
        for reaction_id, residual in sorted(
                reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            total += 1
            if abs(residual) >= epsilon:
                reaction = self._mm.get_reaction(reaction_id)
                rxt = reaction.translated_compounds(
                    lambda x: self._compound_name.get(x, x))
                print('{}\t{}\t{}'.format(reaction_id, residual, rxt))
            else:
                good += 1
        logger.info('Consistent reactions: {}/{}'.format(good, total))
