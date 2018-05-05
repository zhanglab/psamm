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
# Copyright 2018-2018  Ke Zhang <kzhang@my.uri.edu>

from __future__ import unicode_literals

import time
import logging
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command)
import csv
from ..reaction import Direction
form six import text_type, iteritems
from .. import findprimatypairs
from ..formula import Formula, Atom, ParseError
from .. import graph
from collections import Counter


logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#b3cde3'
ACTIVE_COLOR = '#fbb4ae'
ALT_COLOR = '#ccb460'

class VisualizationCommand(MetabolicMixin, ObjectiveMixin,
                         SolverCommandMixin, Command):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help=('Reaction to exclude (e.g. biomass reactions or'
                  ' macromolecule synthesis)'))
        parser.add_argument(
            '--edge-values', type=text_type, default=None,
            help='Values for edges')
        parser.add_argument(
            '--element', type=text_type, default=None,
            help='primary element flow')
        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""

        # parse compound id and formula
        compound_formula = {}
        for compound in self._model.compounds:
            if compound.formula is not None:
                try:
                    f = Formula.parse(compound.formula).flattened()
                    if not f.is_variable():
                        compound_formula[compound.id] = f
                    else:
                        logger.warning(
                            'Skipping variable formula {}: {}'.format(
                                compound.id, compound.formula))    #skip generic compounds
                except ParseError as e:
                    msg = (
                        'Error parsing formula'
                        ' for compound {}:\n{}\n{}'.format(
                            compound.id, e, compound.formula))
                    if e.indicator is not None:
                        msg += '\n{}'.format(e.indicator)
                    logger.warning(msg)


        #set edge_values
        edge_values = None
        if self._args.edge_values is not None:
            raw_values = {}
            with open(self._args.edge_values, 'r') as f:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    raw_values[row[0]] = float(row[1])

            edge_values = {}
            for reaction in self._mm.reactions:
                rx = self._mm.get_reaction(reaction)
                if reaction in raw_values:
                    flux = raw_values[reaction]
                    if abs(flux) < 1e-9:
                        continue

                    if flux > 0:
                        for compound, value in rx.right:  #value = stoichiometry
                            edge_values[reaction, compound] = (
                                    flux * float(value))
                        for compound, value in rx.left:
                            edge_values[compound, reaction] = (
                                    flux * float(value))
                    else:
                        for compound, value in rx.left:
                            edge_values[reaction, compound] = (
                                    -flux * float(value))
                        for compound, value in rx.right:
                            edge_values[compound, reaction] = (
                                    -flux * float(value))