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

import argparse
import logging

from ..command import (Command, MetabolicMixin, LoopRemovalMixin,
                       SolverCommandMixin)
from .. import fastcore, fluxanalysis
from ..fastgapfill import create_extended_model, fastgapfill

logger = logging.getLogger(__name__)


class FastGapFillCommand(MetabolicMixin, LoopRemovalMixin, SolverCommandMixin,
                         Command):
    """Run the FastGapFill gap-filling algorithm on model."""

    _supported_loop_removal = ['none', 'tfba']

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--subset', metavar='file', type=argparse.FileType('r'),
            help='specify core reaction subset to use')
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
            '--epsilon', type=float, help='Threshold for Fastcore',
            default=1e-5)
        parser.add_argument(
            'reaction', help='Reaction to maximize', nargs='?')
        super(FastGapFillCommand, cls).init_parser(parser)

    def run(self):
        """Run FastGapFill command"""

        loop_removal = self._get_loop_removal_option()

        # Create solver
        enable_tfba = loop_removal == 'tfba'
        if enable_tfba:
            solver = self._get_solver(integer=True)
        else:
            solver = self._get_solver()

        # Load compound information
        compound_name = {}
        for compound in self._model.parse_compounds():
            compound_name[compound.id] = (
                compound.name if compound.name is not None else compound.id)

        # TODO: The exchange and transport reactions have tuple names. This
        # means that in Python 3 the reactions can no longer be directly
        # compared (e.g. while sorting) so define this helper function as a
        # workaround.
        def reaction_key(r):
            return r if isinstance(r, tuple) else (r,)

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

        model_extended, weights = create_extended_model(
            self._model,
            db_penalty=self._args.db_penalty,
            ex_penalty=self._args.ex_penalty,
            tp_penalty=self._args.tp_penalty,
            penalties=penalties)

        epsilon = self._args.epsilon
        core = set()
        if self._args.subset is None:
            for r in self._mm.reactions:
                if not self._mm.is_exchange(r):
                    core.add(r)
        else:
            for line in self._args.subset:
                line = line.strip()
                if line == '':
                    continue
                core.add(line)
        induced = fastgapfill(model_extended, core, weights=weights,
                              epsilon=epsilon, solver=solver)

        if self._args.reaction is not None:
            maximized_reaction = self._args.reaction
        else:
            maximized_reaction = self._model.biomass_reaction
            if maximized_reaction is None:
                self.argument_error(
                    'The maximized reaction was not specified')

        if not self._mm.has_reaction(maximized_reaction):
            self.fail('The biomass reaction is not a valid model'
                      ' reaction: {}'.format(maximized_reaction))

        logger.info('Flux balance on induced model maximizing {}'.format(
            maximized_reaction))
        model_induced = model_extended.copy()
        for rxnid in model_extended.reactions:
            if not self._mm.has_reaction(rxnid) and rxnid not in induced:
                model_induced.remove_reaction(rxnid)
        for rxnid, flux in sorted(fluxanalysis.flux_balance(
                model_induced, maximized_reaction, tfba=enable_tfba,
                solver=solver), key=lambda x: (reaction_key(x[0]), x[1])):
            reaction_class = 'Dbase'
            penalty = weights.get(rxnid, 1)
            if self._mm.has_reaction(rxnid):
                reaction_class = 'Model'
                penalty = 0
            rx = model_induced.get_reaction(rxnid)
            rxt = rx.translated_compounds(lambda x: compound_name.get(x, x))
            print('{}\t{}\t{}\t{}\t{}'.format(
                rxnid, reaction_class, penalty, flux, rxt))

        logger.info('Calculating Fastcc consistent subset of induced model')
        consistent_core = fastcore.fastcc_consistent_subset(
            model_induced, epsilon, solver=solver)
        logger.info('Result: |A| = {}, A = {}'.format(
            len(consistent_core), consistent_core))
        removed_reactions = set(model_induced.reactions) - consistent_core
        logger.info('Removed: |R| = {}, R = {}'.format(
            len(removed_reactions), removed_reactions))
