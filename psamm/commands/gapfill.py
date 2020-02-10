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
import argparse

from six import text_type

from ..command import (Command, MetabolicMixin, SolverCommandMixin,
                       FilePrefixAppendAction)
from ..gapfill import gapfill, GapFillError
from ..datasource.reaction import parse_compound
from ..gapfilling import create_extended_model

logger = logging.getLogger(__name__)


class GapFillCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Run the GapFill algorithms on the model."""
    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--compound', metavar='compound', action=FilePrefixAppendAction,
            type=parse_compound, default=[],
            help='Select compounds to try to unblock')
        parser.add_argument(
            '--penalty', metavar='file', type=argparse.FileType('r'),
            help='List of penalty scores for database reactions')
        parser.add_argument(
            '--db-penalty', metavar='penalty', type=float,
            help='Default penalty for database reactions')
        parser.add_argument(
            '--tp-penalty', metavar='penalty', type=float,
            help='Default penalty for transport reactions')
        parser.add_argument(
            '--ex-penalty', metavar='penalty', type=float,
            help='Default penalty for exchange reactions')
        parser.add_argument(
            '--epsilon', type=float, default=1e-5,
            help='Threshold for reaction flux')
        parser.add_argument(
            '--no-implicit-sinks', action='store_true',
            help='Do not include implicit sinks when gap-filling')
        parser.add_argument(
            '--allow-bounds-expansion', action='store_true',
            help=('Allow GapFill to propose expansion of flux bounds. This'
                  ' includes turning irreversible reactions reversible.'))
        super(GapFillCommand, cls).init_parser(parser)

    def run(self):
        """Run GapFill command"""

        # Load compound information
        def compound_name(id):
            if id not in self._model.compounds:
                return id
            return self._model.compounds[id].properties.get('name', id)

        # Calculate penalty if penalty file exists
        penalties = {}
        if self._args.penalty is not None:
            for line in self._args.penalty:
                line, _, comment = line.partition('#')
                line = line.strip()
                if line == '':
                    continue
                rxnid, penalty = line.split(None, 1)
                penalties[rxnid] = float(penalty)

        core = set(self._mm.reactions)

        solver = self._get_solver(integer=True)
        default_comp = self._model.default_compartment
        epsilon = self._args.epsilon
        v_max = float(self._model.default_flux_limit)

        blocked = set()
        for compound in self._args.compound:
            if compound.compartment is None:
                compound = compound.in_compartment(default_comp)
            blocked.add(compound)

        if len(blocked) > 0:
            logger.info('Unblocking compounds: {}...'.format(
                ', '.join(text_type(c) for c in sorted(blocked))))
        else:
            logger.info(
                'Unblocking all compounds in model. Use --compound option to'
                ' unblock specific compounds.')
            blocked = set(self._mm.compounds)

        exclude = set()
        if self._model.biomass_reaction is not None:
            exclude.add(self._model.biomass_reaction)

        # Add exchange and transport reactions to database
        model_complete, weights = create_extended_model(
            self._model,
            db_penalty=self._args.db_penalty,
            ex_penalty=self._args.ex_penalty,
            tp_penalty=self._args.tp_penalty,
            penalties=penalties)

        implicit_sinks = not self._args.no_implicit_sinks

        logger.info('Searching for reactions to fill gaps')
        try:
            added_reactions, no_bounds_reactions = gapfill(
                model_complete, core, blocked, exclude, solver=solver,
                epsilon=epsilon, v_max=v_max, weights=weights,
                implicit_sinks=implicit_sinks,
                allow_bounds_expansion=self._args.allow_bounds_expansion)
        except GapFillError as e:
            self._log_epsilon_and_fail(epsilon, e)

        for reaction_id in sorted(self._mm.reactions):
            rx = self._mm.get_reaction(reaction_id)
            rxt = rx.translated_compounds(compound_name)
            print('{}\t{}\t{}\t{}'.format(reaction_id, 'Model', 0, rxt))

        for rxnid in sorted(added_reactions):
            rx = model_complete.get_reaction(rxnid)
            rxt = rx.translated_compounds(compound_name)
            print('{}\t{}\t{}\t{}'.format(
                rxnid, 'Add', weights.get(rxnid, 1), rxt))

        for rxnid in sorted(no_bounds_reactions):
            rx = model_complete.get_reaction(rxnid)
            rxt = rx.translated_compounds(compound_name)
            print('{}\t{}\t{}\t{}'.format(
                rxnid, 'Remove bounds', weights.get(rxnid, 1), rxt))

    def _log_epsilon_and_fail(self, epsilon, exc):
        msg = ('Finding blocked compounds failed with epsilon set to {}. Try'
               ' lowering the epsilon value to reduce artifical constraints on'
               ' the model.'.format(epsilon))
        self.fail(msg, exc)
