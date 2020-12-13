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
# Copyright 2015-2017  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import logging
import random

from ..command import Command, SolverCommandMixin, \
    FilePrefixAppendAction, convert_to_unicode
from ..formula import Formula, Atom, Radical, ParseError
from .. import findprimarypairs
from .. import mapmaker

from six import text_type, iteritems

logger = logging.getLogger(__name__)


class PrimaryPairsCommand(SolverCommandMixin, Command):
    """Predict primary pairs of reactions.

    This command is used to predict element-transferring reactant/product pairs
    in the reactions of the model. This can be used to determine the flow of
    elements through reactions.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method',
            choices=['fpp', 'mapmaker'],
            default='fpp', help='Primary pair prediction method')
        parser.add_argument(
            '--exclude', metavar='reaction', type=convert_to_unicode,
            default=[], action=FilePrefixAppendAction,
            help=('Reaction to exclude (e.g. biomass reactions or'
                  ' macromolecule synthesis)'))
        parser.add_argument(
            '--report-element', metavar='element', type=text_type, default=[],
            action='append', help=('Only output pairs predicted to transfer'
                                   ' this element (e.g. C, N, P)'))
        parser.add_argument(
            '--report-all-transfers', action='store_true',
            help=('Report all transfers instead of combining elements for'
                  ' each pair (only available for fpp).'))
        parser.add_argument(
            '--ambiguous', action='store_true',
            help=('Report additional information about reactions where'
                  'primary pair identification was ambiguous'))
        parser.add_argument(
            '--weights', action='append', default=[], type=text_type,
            help=('Set weights for elements for inferring compound'
                  ' similarities (e.g. --weights N=0.4,C=1,H=0,R=20,*=0.6).'
                  ' The default value depends on the prediction method.'))
        super(PrimaryPairsCommand, cls).init_parser(parser)

    def run(self):
        # Check that elements are valid
        elements = set()
        for element in self._args.report_element:
            if not hasattr(Atom, element):
                self.argument_error('Invalid element: {}'.format(element))
            elements.add(Atom(element))

        # Check that method can report all transfers if requested
        if self._args.report_all_transfers and self._args.method != 'fpp':
            self.argument_error(
                'Reporting all transfers is not available for this method:'
                ' {}'.format(self._args.method))

        if len(self._args.weights) > 0:
            weights_dict, r_group_weight, default_weight = _parse_weights(
                self._args.weights, default_weight=0.6)

            def element_weight(element):
                if isinstance(element, Radical):
                    return r_group_weight
                return weights_dict.get(element, default_weight)
        else:
            if self._args.method == 'fpp':
                logger.info(
                    'Using default element weights for {}:'
                    ' C=1, H=0, *=0.82'.format(self._args.method))
                element_weight = findprimarypairs.element_weight
            elif self._args.method == 'mapmaker':
                logger.info(
                    'Using default element weights for {}:'
                    ' H=0, N=0.4, O=0.4, P=0.4, R=40, *=1'.format(
                        self._args.method))
                element_weight = mapmaker.default_weight

        # Mapping from compound id to formula
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
                                compound.id, compound.formula))
                except ParseError as e:
                    msg = (
                        'Error parsing formula'
                        ' for compound {}:\n{}\n{}'.format(
                            compound.id, e, compound.formula))
                    if e.indicator is not None:
                        msg += '\n{}'.format(e.indicator)
                    logger.warning(msg)

        # Set of excluded reactions
        exclude = set(self._args.exclude)

        def iter_reactions():
            for reaction in self._model.reactions:
                if (reaction.id not in self._model.model or
                        reaction.id in exclude):
                    continue

                if reaction.equation is None:
                    logger.warning(
                        'Reaction {} has no reaction equation'.format(
                            reaction.id))
                    continue

                if any(c.name not in compound_formula
                       for c, _ in reaction.equation.compounds):
                    logger.warning(
                        'Reaction {} has compounds with undefined'
                        ' formula'.format(reaction.id))
                    continue

                yield reaction

        if self._args.method == 'fpp':
            result = self._run_find_primary_pairs(
                compound_formula, iter_reactions(), element_weight)
        elif self._args.method == 'mapmaker':
            result = self._run_mapmaker(
                compound_formula, iter_reactions(), element_weight)
        else:
            self.argument_error('Unknown method: {}'.format(self._args.method))

        if not self._args.report_all_transfers:
            result = self._combine_transfers(result)

        for reaction_id, c1, c2, form in result:
            if len(elements) == 0 or any(e in form for e in elements):
                print('{}\t{}\t{}\t{}'.format(reaction_id, c1, c2, form))

    def _combine_transfers(self, result):
        """Combine multiple pair transfers into one."""
        transfers = {}
        for reaction_id, c1, c2, form in result:
            key = reaction_id, c1, c2
            combined_form = transfers.setdefault(key, Formula())
            transfers[key] = combined_form | form

        for (reaction_id, c1, c2), form in iteritems(transfers):
            yield reaction_id, c1, c2, form

    def _run_find_primary_pairs(
            self, compound_formula, reactions, element_weight):
        reaction_pairs = [(r.id, r.equation) for r in reactions]
        ambig = self._args.ambiguous
        prediction, _ = findprimarypairs.predict_compound_pairs_iterated(
            reaction_pairs, compound_formula,
            ambig, element_weight=element_weight)

        for reaction_id, _ in reaction_pairs:
            if reaction_id not in prediction:
                logger.warning('Failed to predict pairs for {}'.format(
                    reaction_id))
                continue

            pairs, balance = prediction[reaction_id]
            if len(balance) > 0:
                logger.warning('Reaction {} is not balanced!'.format(
                    reaction_id))

            for (c1, c2), forms in sorted(iteritems(pairs)):
                for form in forms:
                    yield reaction_id, c1, c2, form

    def _run_mapmaker(self, compound_formula, reactions, element_weight):
        solver = self._get_solver(integer=True)
        ambiguous = set()
        for reaction in reactions:
            logger.info('Predicting reaction {}...'.format(reaction.id))
            try:
                transfer = list(mapmaker.predict_compound_pairs(
                                reaction.equation, compound_formula, solver,
                                weight_func=element_weight))
                if self._args.ambiguous:
                    if len(transfer) > 1:
                        ambiguous.add(reaction.id)
            except mapmaker.UnbalancedReactionError:
                logger.info('Reaction {} is not balanced! Skipping...'.format(
                    reaction.id))
                continue
            if len(transfer) == 0:
                logger.warning('Reaction {} has no predicted transfer!'.format(
                    reaction.id))
                continue

            for (c1, c2), form in sorted(
                    iteritems(random.sample(transfer, 1)[0])):
                yield reaction.id, c1, c2, form
        if self._args.ambiguous:
            for rxn_ambig in ambiguous:
                logger.info('Ambiguous Transfers '
                            'in Reaction: {}'.format(rxn_ambig))


def _parse_weights(weight_args, default_weight=0.6):
    """Parse list of weight assignments."""
    weights_dict = {}
    r_group_weight = default_weight
    for weight_arg in weight_args:
        for weight_assignment in weight_arg.split(','):
            if '=' not in weight_assignment:
                raise ValueError(
                    'Invalid weight assignment: {}'.format(weight_assignment))

            key, value = weight_assignment.split('=', 1)
            value = float(value)
            if key == 'R':
                r_group_weight = value
            elif key == '*':
                default_weight = value
            elif hasattr(Atom, key):
                weights_dict[Atom(key)] = value
            else:
                raise ValueError('Invalid element: {}'.format(key))

    return weights_dict, r_group_weight, default_weight
