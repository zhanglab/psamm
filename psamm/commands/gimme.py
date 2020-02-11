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
# Copyright 2019-2019  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Implementation of Flux Balance Analysis."""

from __future__ import unicode_literals

import logging
import csv
from ..lpsolver import lp
import argparse
from six import iteritems
from ..expression import boolean
from ..lpsolver.lp import Expression, ObjectiveSense
from ..fluxanalysis import FluxBalanceProblem
from ..reaction import Direction, Reaction
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command)

logger = logging.getLogger(__name__)


class GimmeCommand(MetabolicMixin, LoopRemovalMixin, ObjectiveMixin,
                   SolverCommandMixin, Command):
    """Run flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--biomass-threshold', type=float, default=None,
            help='Threshold for biomass reaction flux')
        parser.add_argument(
            '--expression-threshold', type=float,
            help='Threshold for gene expression')
        parser.add_argument(
            '--transcriptome-file', type=argparse.FileType('r'),
            help='Two column file of gene ID to expression')
        super(GimmeCommand, cls).init_parser(parser)

    def run(self):
        solver = self._get_solver()
        model = self._model
        mm = model.create_metabolic_model()

        base_gene_dict = {}
        for rxn in model.reactions:
            base_gene_dict[rxn.id] = rxn.genes

        threshold_value = parse_transcriptome_file(
            self._args.transcriptome_file,
            self._args.expression_threshold)

        gene_dict = {}

        exchange_exclude = [rxn for rxn in mm.reactions
                            if mm.is_exchange(rxn)]

        mm_irreversible, reversible_gene_assoc, split_rxns = make_irreversible(
            mm, base_gene_dict, exclude_list=exchange_exclude)

        p = FluxBalanceProblem(mm_irreversible, solver)

        final_model, below_threshold_ids, incon_score = solve_gimme_problem(p,
                        mm_irreversible, self._get_objective(),
                        reversible_gene_assoc, split_rxns, threshold_value,
                        self._args.biomass_threshold)

        for reaction in sorted(final_model):
            print(reaction)

        logger.info('Used {} below threshold reacctions out of {} total '
                    'below threshold reactions'.format(str(
            len([x for x in below_threshold_ids if x in final_model])),
            str(len(below_threshold_ids))))
        logger.info('Inconsistency Score: {}'.format(incon_score))


def solve_gimme_problem(problem, mm, biomass, reversible_gene_assoc,
                        split_rxns, transcript_values, threshold):
    """Formulates and Solves a GIMME model.

    Implementation of the GIMME algorithm (Becker and Pallson 2008).
    Accepts an irreversible metabolic model, LP problem, and GIMME specific
    data. Will add in relevant GIMME constraints to the LP problem and
    generates a contextualized model based on the media constraints and the
    transcriptome data provided.

    Args:
        problem: :class:`FluxBalanceProblem` to solve.
        mm: An irreversible metabolic model.
        biomass: Biomass reaction ID.
        reversible_gene_assoc: A dictionary of gene IDs from make_irreversible.
        split_rxns: A set of tuples of reaction IDs from make_irreversible.
        transcript_values: A dictionary returned from parse_transcriptome_file.
        threshold: A threshold that the biomass flux needs to stay above.
    """
    ci_dict = {}
    for reaction in mm.reactions:
        problem.prob.define('zi_{}'.format(reaction),
                            types=lp.VariableType.Binary)
        gene_string = reversible_gene_assoc.get(reaction)
        if gene_string is None:
            continue
        else:
            e = boolean.Expression(gene_string)
            ci = get_rxn_value(e._root, transcript_values)
            if ci is not None:
                ci_dict[reaction] = ci

    gimme_objective = Expression()

    for (forward, reverse) in split_rxns:
        problem.prob.add_linear_constraints(
            problem.prob.var('zi_{}'.format(forward)) + problem.prob.var(
                'zi_{}'.format(reverse)) <= 1)

    if threshold is None:
        problem.maximize(biomass)
        problem.prob.add_linear_constraints(
            problem.get_flux_var(biomass) >= problem.get_flux(biomass))
    else:
        problem.maximize(biomass)
        if problem.get_flux(biomass) < threshold:
            logger.warning('Input threshold greater than '
            'maximum biomass: {}'.format(problem.get_flux(biomass)))
        else:
            problem.prob.add_linear_constraints(
                problem.get_flux_var(biomass) >= threshold)

    for key, value in iteritems(ci_dict):
        gimme_objective += problem.get_flux_var(key) * value
    problem.prob.define('gimme_objective', types=lp.VariableType.Continuous)
    obj = problem.prob.var('gimme_objective')
    problem.prob.add_linear_constraints(obj == gimme_objective)
    problem.prob.set_objective(obj)
    problem.prob.set_objective_sense(ObjectiveSense.Minimize)
    problem.prob.solve()
    used_rxns = []
    for reaction in mm.reactions:
        if abs(problem.get_flux(reaction)) > 1e-12:
            original_id = reaction.replace('_forward', '')
            original_id = original_id.replace('_reverse', '')
            used_rxns.append(original_id)
    used_below = 0
    used_above = 0
    off_below = 0
    off_above = 0
    final_model = set()
    original_reaction_set = set()
    for reaction in mm.reactions:
        if mm.is_exchange(reaction):
            continue
        else:
            reaction = reaction.replace('_forward', '')
            reaction = reaction.replace('_reverse', '')
            original_reaction_set.add(reaction)
    below_threshold_ids = set()
    for reaction in ci_dict.keys():
        reaction = reaction.replace('_forward', '')
        reaction = reaction.replace('_reverse', '')
        below_threshold_ids.add(reaction)

    used_below_list = []
    for reaction in original_reaction_set:
        if reaction in used_rxns:
            if reaction in below_threshold_ids:
                used_below += 1
                used_below_list.append(reaction)
                final_model.add(reaction)
            else:
                used_above += 1
                final_model.add(reaction)
        else:
            if reaction in below_threshold_ids:
                off_below += 1
            else:
                off_above += 1
                final_model.add(reaction)
    return final_model, below_threshold_ids, problem.prob.result.get_value(obj)


def parse_transcriptome_file(f, threshold):
    """Parses a file containing a gene to expression mapping.

    Parses a tab separated two column table. The first column contains gene
    IDs while the second column contains expression values. Will compare the
    expression values to the given threshold and return a dict with any
    genes that fall under the expression threshold, and the amount that they
    fall under the threshold.

    Args:
        f: a file containing a two column gene to expression table.
        threshold: Expression threshold value.
    """
    threshold_value = {}
    for row in csv.reader(f, delimiter=str('\t')):
        try:
            gene = row[0]
            gene = gene.replace('"', '')
            if float(row[1]) < threshold:
                threshold_value[gene] = threshold - float(row[1])
            elif float(row[1]) >= threshold:
                continue
        except ValueError:
            logger.warning('Invalid expression value provided: '
            'gene: {} value: {}'.format(row[0], row[1]))
    return threshold_value


def get_rxn_value(root, gene_dict):
    """Gets overall expression value for a reaction gene association.

    Recursive function designed to parse a gene expression and return
    an expression value to use in the GIMME algorithm. This function is
    designed to return the value directly if the expression only has
    one gene. If the expression has multiple genes related by 'OR'
    associations, then it will return the highest of the set of values.
    If the genes are associated with 'AND' logic then the function
    will return the lowest of the set of values.

    Args:
        root: object of boolean.Expression()._root
        gene_dict: dict of gene expression from parse_transcriptome_file
    """
    if type(root) == boolean.Variable:
        return gene_dict.get(root.symbol)
    elif type(root) is boolean.And:
        val_list = [x for x in
                    [get_rxn_value(i, gene_dict) for i in root._terms]
                    if x is not None]
        if len(val_list) == 0:
            return None
        else:
            return min(x for x in val_list if x is not None)
        # get min value
    elif type(root) is boolean.Or:
        val_list = [x for x in
                    [get_rxn_value(i, gene_dict)
                     for i in root._terms]]
        if None in val_list:

            return None
        else:
            return max(x for x in val_list if x is not None)


def make_irreversible(mm, gene_dict, exclude_list=[],
                      all_reversible=False):
    """Creates a new metabolic models with only irreversible reactions.

    This function will find every reversible reaction in the
    model and split it into two reactions with the
    {rxnid}_forward or {rxnid}_reverse as the IDs.

    Args:
        mm: A metabolicmodel object
        exclude_list: list of reactions to exclude in TMFA simulation
        all_reversible: if True make all reactions in model reversible.
    """
    split_reversible = set()
    mm_irrev = mm.copy()
    reversible_gene_dict = {}
    for rxn in mm.reactions:
        upper = mm.limits[rxn].upper
        lower = mm.limits[rxn].lower
        mm_irrev.limits[rxn].upper = upper
        mm_irrev.limits[rxn].lower = lower

        reaction = mm_irrev.get_reaction(rxn)
        if rxn not in exclude_list:
            r = Reaction(Direction.Forward, reaction.left, reaction.right)
            r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
            r_id = str('{}_forward'.format(rxn))
            r2_id = str('{}_reverse'.format(rxn))
            if reaction.direction == Direction.Forward:
                if all_reversible is False:
                    reversible_gene_dict[rxn] = gene_dict.get(rxn)
                    continue
                else:
                    mm_irrev.remove_reaction(rxn)
                    mm_irrev.database.set_reaction(r_id, r)
                    mm_irrev.database.set_reaction(r2_id, r2)
                    mm_irrev.add_reaction(r_id)
                    mm_irrev.add_reaction(r2_id)
                    split_reversible.add((r_id, r2_id))
                    reversible_gene_dict[r_id] = gene_dict.get(rxn)
                    reversible_gene_dict[r2_id] = gene_dict.get(rxn)
            elif reaction.direction == Direction.Both:
                mm_irrev.remove_reaction(rxn)
                mm_irrev.database.set_reaction(r_id, r)
                mm_irrev.database.set_reaction(r2_id, r2)
                mm_irrev.add_reaction(r_id)
                mm_irrev.add_reaction(r2_id)
                split_reversible.add((r_id, r2_id))
                reversible_gene_dict[r_id] = gene_dict.get(rxn)
                reversible_gene_dict[r2_id] = gene_dict.get(rxn)
            if upper == lower:
                mm_irrev.limits[r_id].upper = upper
                mm_irrev.limits[r_id].lower = lower
                mm_irrev.limits[r2_id].upper = 0
                mm_irrev.limits[r2_id].lower = 0
            else:
                mm_irrev.limits[r_id].upper = upper
                mm_irrev.limits[r_id].lower = 0
                if lower == 0:
                    mm_irrev.limits[r2_id].upper = upper
                else:
                    mm_irrev.limits[r2_id].upper = -lower
                mm_irrev.limits[r2_id].lower = 0

    return mm_irrev, reversible_gene_dict, split_reversible
