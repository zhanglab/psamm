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

import argparse
import logging

from ..command import Command, MetabolicMixin, SolverCommandMixin
from ..fastgapfill import fastgapfill
from ..gapfilling import create_extended_model

logger = logging.getLogger(__name__)


class FastGapFillCommand(MetabolicMixin, SolverCommandMixin, Command):
    """Run the FastGapFill gap-filling algorithm on model."""

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
        super(FastGapFillCommand, cls).init_parser(parser)

    def run(self):
        """Run FastGapFill command"""

        # Create solver
        solver = self._get_solver()

        # Load compound information
        def compound_name(id):
            if id not in self._model.compounds:
                return id
            return self._model.compounds[id].properties.get('name', id)

        # TODO: The exchange and transport reactions have tuple names. This
        # means that in Python 3 the reactions can no longer be directly
        # compared (e.g. while sorting) so define this helper function as a
        # workaround.
        def reaction_key(r):
            return r if isinstance(r, tuple) else (r,)

        # Calculate penalty if penalty file exists
        penalties = {}
        for id in sorted(self._mm.reactions):
            penalties[id] = 0
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

        for reaction_id in sorted(self._mm.reactions):
            rx = self._mm.get_reaction(reaction_id)
            rxt = rx.translated_compounds(compound_name)
            print('{}\t{}\t{}\t{}'.format(
                reaction_id, 'Model', weights[reaction_id], rxt))

        for rxnid in sorted(induced, key=reaction_key):
            if self._mm.has_reaction(rxnid):
                continue
            rx = model_extended.get_reaction(rxnid)
            rxt = rx.translated_compounds(compound_name)
            print('{}\t{}\t{}\t{}'.format(
                rxnid, 'Add', weights.get(rxnid), rxt))
