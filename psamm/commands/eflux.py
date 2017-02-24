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
# Copyright 2016  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2016  Julie Cuddigan <julie_cuddigan@my.uri.edu>
# Copyright 2016  Matthew Gentry <mgentry@umass.edu>

from __future__ import unicode_literals

import math
import csv
import logging
from itertools import chain

from six import iteritems, itervalues, text_type

from ..command import SolverCommandMixin, MetabolicMixin, Command
from ..expression import boolean
from ..util import MaybeRelative
from .. import fluxanalysis

logger = logging.getLogger(__name__)


class EFluxBalance(MetabolicMixin, SolverCommandMixin, Command):
    """Run eflux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--threshold',
            help='Minimum objective flux as a decimal or percent',
            type=MaybeRelative, default=MaybeRelative('100%'))
        parser.add_argument(
            '--condition', help='Condition values to use (default: first)',
            type=text_type)
        parser.add_argument(
            '--data', help='Path to transcriptomics data file', metavar='FILE')
        parser.add_argument(
            '--fva', help='Enable FVA', action='store_true')
        parser.add_argument(
            '--transform',
            help=('Transformation to apply to reaction expression values'
                  ' (default: none)'),
            default='none', choices=['none', 'logistic'])
        super(EFluxBalance, cls).init_parser(parser)

    def run(self):
        """Run E-Flux implementation."""
        condition_names = []
        conditions = {}
        if self._args.data is not None:
            condition_names, conditions = parse_transcription_file(
                self._args.data)

        if len(condition_names) > 0:
            logger.info('Available conditions: {}'.format(
                ', '.join(condition_names)))

        condition = None
        if self._args.condition is not None:
            if self._args.condition not in conditions:
                self.argument_error('Invalid condition: {}'.format(
                    self._args.condition))
            condition = self._args.condition
        elif len(condition_names) > 0:
            condition = condition_names[0]

        if condition is not None:
            logger.info('Using condition: {}'.format(condition))
        else:
            logger.warning(
                'No condition selected! Specify a transcriptomics data file'
                ' using the --data option.')

        rxn_genelogic = self.parse_gene_logic()

        biomass = self._model.biomass_reaction

        solver = self._get_solver()
        prob = fluxanalysis.FluxBalanceProblem(self._mm, solver)

        # Determine combined expression level for each reaction
        reaction_levels = {}
        for cond, values in iteritems(conditions):
            reaction_levels[cond] = dict(reaction_expression_level(
                rxn_genelogic, conditions[cond]))

        # Optionally transform reaction levels
        if self._args.transform == 'logistic':
            for cond, levels in iteritems(reaction_levels):
                reaction_levels[cond] = logistic_transform(
                    reaction_levels[cond])

        condition_levels = {}
        if condition is not None:
            condition_levels = reaction_levels[condition]

        max_level = 0
        levels = list(chain(*(
            itervalues(levels) for levels in itervalues(reaction_levels))))
        if len(levels) > 0:
            max_level = max(levels)

        # Create flux bounds based on expression level and solve
        create_expression_bounds(prob, self._mm, condition_levels, max_level)
        prob.maximize(biomass)

        if self._args.fva:
            thresh = self._args.threshold
            fba_biomass_flux = prob.get_flux(biomass)
            thresh.reference = fba_biomass_flux
            prob.prob.add_linear_constraints(
                prob.get_flux_var(biomass) >= float(thresh))

            for rxn in sorted(self._mm.reactions):
                print('{}\t{}\t{}\t{}\t{}'.format(
                    rxn, prob.flux_bound(rxn, -1),
                    prob.flux_bound(rxn, 1), self._mm.get_reaction(rxn),
                    rxn_genelogic.get(rxn, '')))
        else:
            for rxn in sorted(self._mm.reactions):
                print('{}\t{}\t{}\t{}'.format(
                    rxn, prob.get_flux(rxn), self._mm.get_reaction(rxn),
                    rxn_genelogic.get(rxn, '')))

    def parse_gene_logic(self):
        """Return dictionary of gene logic with reaction ID keys."""
        gene_dict = {}
        for reaction in self._model.parse_reactions():
            if self._mm.has_reaction(reaction.id):
                if reaction.genes is None:
                    continue
                gene_dict[reaction.id] = boolean.Expression(reaction.genes)

        return gene_dict


def combined_eflux_expression(term, data):
    """Return combined eflux expression level for boolean term.

    Recursively unpacks gene logic and returns the expression level of the
    boolean term according to E-Flux rules: Genes are translated directly to
    the expression level in data; AND-terms are translated to the minimum
    value of its subterms; and OR-terms are translated to the sum of its
    subterms.

    Args:
        term: boolean term.
        data: Dictionary of expression level for genes.
    """
    if isinstance(term, boolean.Variable):
        return data.get(term.symbol)
    else:
        values = []
        for subterm in term:
            value = combined_eflux_expression(subterm, data)
            if value is not None:
                values.append(value)

        if len(values) == 0:
            return None

        if isinstance(term, boolean.Or):
            return sum(values)
        else:
            return min(values)


def create_expression_bounds(prob, mm, reaction_level, max_level):
    """Create bounds on reaction fluxes based on expression.

    Alters the bounds on metabolic reactions by a factor equal to the
    reaction expression, relative to the maximum expression across all
    conditions.

    Args:
        prob: Linear programming problem to apply bounds to.
        mm: :class:`MetabolicModel`.
        rxn_exp: Dictionary mapping reactions to expression level.
        max_level: Maximum expression level to normalize by.
    """

    for reaction, level in sorted(iteritems(reaction_level)):
        factor = level / float(max_level)

        lower = float(mm.limits[reaction].lower) * factor
        prob.prob.add_linear_constraints(
            prob.get_flux_var(reaction) >= lower)

        upper = float(mm.limits[reaction].upper) * factor
        prob.prob.add_linear_constraints(
            prob.get_flux_var(reaction) <= upper)


def reaction_expression_level(gene_logic, condition):
    """Yield pairs of reaction ID and expression level for condition.

    Args:
        gene_logic: Dictionary of reaction IDs mapped to gene logic.
        condition: Dictionary of gene expression levels.
    """
    for reaction, genes in iteritems(gene_logic):
        level = combined_eflux_expression(genes.root, condition)
        if level is None:
            continue

        yield reaction, level


def logistic(x, mean=0.0, k=0.005):
    """Logistic function centered about the midpoint (default zero).

    Args:
        mean: The midpoint of the sigmoid.
        k: Steepness of the curve.
    """
    return 1.0/(1.0 + math.exp(-k*(x-mean)))


def logistic_transform(condition):
    """Apply logistic transform on transcription level dictionaries.

    Args:
        condition: Dictionary mapping reaction IDs to expression level.
    """
    if len(condition) == 0:
        return {}

    mean = sum(itervalues(condition)) / float(len(condition))
    return {k: logistic(v, mean) for k, v in iteritems(condition)}


def parse_transcription_file(path):
    """Returns the contents of transcription file in a tuple of dictionaries.

    First column is the gene ID. Additional columns are conditions. A header
    must be included.
    """
    condition_names = None
    conditions = {}
    with open(path) as file1:
        for i, row in enumerate(csv.reader(file1, delimiter=str('\t'))):
            if i == 0:
                condition_names = row[1:]
                for name in condition_names:
                    conditions[name] = {}
            else:
                gene = row[0]
                for col_index, value in enumerate(row[1:]):
                    col = condition_names[col_index]
                    conditions[col][gene] = float(value)

    return condition_names, conditions
