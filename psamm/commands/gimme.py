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
# Copyright 2019-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

"""Implementation of Flux Balance Analysis."""

from __future__ import unicode_literals

import logging
import csv
from os import mkdir, path
from ..lpsolver import lp
import argparse
from six import iteritems
from ..expression import boolean
from ..lpsolver.lp import Expression, ObjectiveSense
from ..fluxanalysis import FluxBalanceProblem
from ..command import (ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command)
from psamm.importer import write_yaml_model
from ..util import MaybeRelative

logger = logging.getLogger(__name__)


class GimmeCommand(MetabolicMixin, ObjectiveMixin,
                   SolverCommandMixin, Command):
    """Subset a metabolic model based on gene expression data."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--biomass-threshold', help='Threshold of objective reaction '
                                        'flux. Can be an absolute flux value '
                                        'or percentage with percent sign.',
            type=MaybeRelative, default=MaybeRelative('100%'))
        parser.add_argument(
            '--expression-threshold', type=float,
            help='Threshold for gene expression')
        parser.add_argument(
            '--transcriptome-file', type=argparse.FileType('r'),
            help='Two column file of gene ID to expression')
        parser.add_argument('--detail', action='store_true',
                            help='Show model statistics.')
        parser.add_argument(
            '--export-model', type=str, default=None,
            help='Path to directory for full model export.')
        super(GimmeCommand, cls).init_parser(parser)

    def run(self):
        solver = self._get_solver()
        model = self._model
        mm = model.create_metabolic_model()

        if self._args.export_model is not None:
            if path.exists('{}'.format(self._args.export_model)):
                logger.warning('Output directory {} already exists.'.format(
                    self._args.export_model))
                quit()

        base_gene_dict = {}
        for rxn in model.reactions:
            base_gene_dict[rxn.id] = rxn.genes

        threshold_value = parse_transcriptome_file(
            self._args.transcriptome_file,
            self._args.expression_threshold)

        exchange_exclude = [rxn for rxn in mm.reactions
                            if mm.is_exchange(rxn)]
        exchange_exclude.append(self._get_objective())
        mm_irreversible, reversible_gene_assoc, split_rxns, _ = \
            mm.make_irreversible(base_gene_dict,
                                 exclude_list=exchange_exclude)

        p = FluxBalanceProblem(mm_irreversible, solver)

        final_model, used_exchange, below_threshold_ids, \
            incon_score = \
            solve_gimme_problem(p, mm_irreversible,
                                self._get_objective(),
                                reversible_gene_assoc,
                                split_rxns, threshold_value,
                                self._args.biomass_threshold)

        print('# internal Reactions')
        for reaction in sorted(final_model):
            print(reaction)
        print('# Exchange Reactions')
        for reaction in used_exchange:
            print(reaction)
        used_below = str(len([x for x in
                              below_threshold_ids if x in final_model]))
        below = str(len(below_threshold_ids))
        if self._args.detail:
            logger.info('Used {} below threshold reactions out of {} total '
                        'below threshold reactions'.format(used_below, below))
            logger.info('Inconsistency Score: {}'.format(incon_score))
        if self._args.export_model:
            mkdir('{}'.format(self._args.export_model))
            reactions_to_discard = []
            for reaction in self._model.reactions:
                if reaction.id not in final_model:
                    reactions_to_discard.append(reaction.id)
            for rid in reactions_to_discard:
                self._model.reactions.discard(rid)
            compound_set = set()
            for reaction in self._model.reactions:
                for compound in reaction.equation.compounds:
                    compound_set.add(compound[0].name)
            compounds_to_discard = []
            for compound in self._model.compounds:
                if compound.id not in compound_set:
                    compounds_to_discard.append(compound.id)
            for cid in compounds_to_discard:
                self._model.compounds.discard(cid)
            exchanges_to_discard = []
            for key in self._model.exchange:
                if key.name not in compound_set:
                    exchanges_to_discard.append(key)
            for key in exchanges_to_discard:
                del self._model.exchange[key]
            limits_to_discard = []
            for key in self._model.limits:
                if key not in final_model:
                    limits_to_discard.append(key)
            for key in limits_to_discard:
                del self._model.limits[key]
            write_yaml_model(self._model, dest=self._args.export_model,
                             split_subsystem=False)


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
        gene_string = reversible_gene_assoc.get(reaction)
        if gene_string is None:
            continue
        else:
            e = boolean.Expression(gene_string)
            ci = get_rxn_value(e._root, transcript_values)
            if ci is not None:
                ci_dict[reaction] = ci

    gimme_objective = Expression()

    threshold = threshold
    problem.maximize(biomass)
    if threshold.relative:
        threshold.reference = problem.get_flux(biomass)

    logger.info('Setting objective threshold to {}'.format(
        threshold))

    if problem.get_flux(biomass) < float(threshold):
        logger.warning('Input threshold '
                       'greater than maximum '
                       'biomass: {}'.format(problem.get_flux(biomass)))
        quit()
    problem.prob.add_linear_constraints(
        problem.get_flux_var(biomass) >= float(threshold))

    for key, value in iteritems(ci_dict):
        gimme_objective += problem.get_flux_var(key) * value
    problem.prob.define('gimme_objective', types=lp.VariableType.Continuous)
    obj = problem.prob.var('gimme_objective')
    problem.prob.add_linear_constraints(obj == gimme_objective)
    problem.prob.set_objective(obj)
    problem.prob.set_objective_sense(ObjectiveSense.Minimize)
    problem.prob.solve()
    used_rxns = set()
    sum_fluxes = 0
    for reaction in mm.reactions:
        if abs(problem.get_flux(reaction)) > 1e-12:
            original_id = reaction.replace('_forward', '')
            original_id = original_id.replace('_reverse', '')
            used_rxns.add(original_id)
            sum_fluxes += abs(problem.get_flux(reaction))
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
    used_exchange = used_rxns - original_reaction_set

    return final_model, used_exchange, below_threshold_ids, \
        problem.prob.result.get_value(obj)


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
            logger.warning('Invalid expression value '
                           'provided: gene: {} '
                           'value: {}'.format(row[0], row[1]))
    return threshold_value


def get_rxn_value(root, gene_dict):
    """Gets overall expression value for a reaction gene association.

    Recursive function designed to parse a gene expression and return
    a penalty value to use in the GIMME algorithm. This function is
    designed to return the value directly if the expression only has
    one gene. If the expression has multiple genes related by 'OR'
    associations, then it will return the highest lowest penalty value.
    If the genes are associated with 'AND' logic then the function
    will return the highest penalty value of the set of values.

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
            return max(x for x in val_list if x is not None)
    elif type(root) is boolean.Or:
        val_list = [x for x in
                    [get_rxn_value(i, gene_dict)
                     for i in root._terms]]
        if None in val_list:
            return None
        else:
            return min(x for x in val_list if x is not None)
