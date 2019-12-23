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
# Copyright 2019-2019  Jon Lund Steffensen <jon_steffensen@uri.edu>

"""Implementation of Flux Balance Analysis."""

from __future__ import unicode_literals

import logging
import csv
from ..lpsolver import lp, cplex

from ..expression import boolean
from ..lpsolver.lp import Expression, ObjectiveSense
from ..fluxanalysis import FluxBalanceProblem
from ..reaction import Direction, Reaction
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command)

# Module-level logging
logger = logging.getLogger(__name__)


class GimmeCommand(MetabolicMixin, LoopRemovalMixin, ObjectiveMixin,
                         SolverCommandMixin, Command):
    """Run flux balance analysis on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--biomass-threshold', type=float, default=None)
        parser.add_argument(
            '--expression-threshold', type=float)
        parser.add_argument(
            '--transcriptome-file', type=file)
        super(GimmeCommand, cls).init_parser(parser)

    def run(self):
        """Run flux analysis command."""
        threshold_dict, threshold_value = parse_transcriptome_file(
	        self._args.transcriptome_file,
	        self._args.expression_threshold)
        make_gimme_model(self._model,
            self._get_objective(),
            threshold_dict, threshold_value, self._args.biomass_threshold)


def make_gimme_model(model, biomass, transcriptome_dict,
                     trasncript_values, objective_threshold):
    base_gene_dict = {}
    for rxn in model.reactions:
        base_gene_dict[rxn.id] = rxn.genes
    base_mm = model.create_metabolic_model()
    exchange_exclude = [rxn for rxn in base_mm.reactions
                        if base_mm.is_exchange(rxn)]

    mm, reversible_gene_assoc, split_rxns = make_irreversible(
	    model.create_metabolic_model(), base_gene_dict,
	    exclude_list=exchange_exclude)
    ci_dict = {}

    solver = cplex.Solver()
    p = FluxBalanceProblem(mm, solver)

    for reaction in mm.reactions:
        p.prob.define('zi_{}'.format(reaction), types=lp.VariableType.Binary)
        gene_string = reversible_gene_assoc.get(reaction)
        if gene_string is None:
            continue
        else:
            e = boolean.Expression(gene_string)
            ci = get_rxn_value(e._root, trasncript_values)
            if ci is not None:
                ci_dict[reaction] = ci

    gimme_objective = Expression()

    for (forward, reverse) in split_rxns:
        p.prob.add_linear_constraints(
	        p.prob.var('zi_{}'.format(forward)) + p.prob.var(
		        'zi_{}'.format(reverse)) <= 1)

    if objective_threshold is None:
        p.maximize(biomass)
        p.prob.add_linear_constraints(
	        p.get_flux_var(biomass) >= p.get_flux(biomass))
    else:
        p.maximize(biomass)
        logger.info('{}'.format(p.get_flux(biomass)))
        p.prob.add_linear_constraints(
	        p.get_flux_var(biomass) >= objective_threshold)

    for key, value in ci_dict.iteritems():
        gimme_objective += p.get_flux_var(key) * value
    p.prob.define('gimme_objective', types=lp.VariableType.Continuous)
    obj = p.prob.var('gimme_objective')
    p.prob.add_linear_constraints(obj == gimme_objective)
    p.prob.set_objective(obj)
    p.prob.set_objective_sense(ObjectiveSense.Minimize)
    p.prob.solve()
    used_rxns = []
    for reaction in mm.reactions:
        if abs(p.get_flux(reaction)) > 0.00000000001:
            original_id = reaction.replace('_forward', '')
            original_id = original_id.replace('_reverse', '')
            used_rxns.append(original_id)
    used_below = 0
    used_above = 0
    off_below = 0
    off_above = 0
    final_model = []
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

    for reaction in original_reaction_set:
        if reaction in used_rxns:
            if reaction in below_threshold_ids:
                used_below += 1
                final_model.append(reaction)
            else:
                used_above += 1
                final_model.append(reaction)
        else:
            if reaction in below_threshold_ids:
                off_below += 1
            else:
                off_above += 1
                final_model.append(reaction)

    logger.info('Inconsistency Score\t{}'.format(p.prob.result.get_value(obj)))
    logger.info('final model size\t{}'.format(len(final_model)))
    logger.info('reactions below expression threshold used/total {}/{}'.format(
	    str(used_below), str(used_below+off_below)))
    logger.info('Biomass flux: {}'.format(p.get_flux(biomass)))
    for i in sorted(final_model):
        print(i)


def parse_transcriptome_file(f, threshold):
    threshold_value = {}
    threshold_dict = {}
    for row in csv.reader(f, delimiter=str('\t')):
        try:
            gene = row[0]
            gene = gene.replace('"', '')
            if float(row[1]) < threshold:
                threshold_dict[gene] = False
                threshold_value[gene] = threshold - float(row[1])
            elif float(gene) >= threshold:
                threshold_dict[gene] = True
        except ValueError:
            continue
    return threshold_dict, threshold_value


def get_rxn_value(root, gene_dict):
    if type(root) == boolean.Variable:
        # return (gene_dict.get(var.symbol))
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


def make_irreversible(mm, gene_dict, exclude_list = [],
                      all_reversible = False):
    """Creates a new metabolic models with only irreversible reactions.

    This function will find every reversible reaction in the
    model and split it into two reactions with the
    {rxnid}_forward or {rxnid}_reverse as the IDs.

    Args:
        mm: A metabolicmodel object
        exclude_list: list of reactions to exclude in TMFA simulation
        all_reversible: if True make all reactions in model reversible.
    """
    split_reversible = []
    mm_irrev = mm.copy()
    lumped_rxns = []
    new_lump_rxn_dict = {}
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
                    split_reversible.append((r_id, r2_id))
                    reversible_gene_dict[r_id] = gene_dict.get(rxn)
                    reversible_gene_dict[r2_id] = gene_dict.get(rxn)
            elif reaction.direction == Direction.Both:
                mm_irrev.remove_reaction(rxn)
                mm_irrev.database.set_reaction(r_id, r)
                mm_irrev.database.set_reaction(r2_id, r2)
                mm_irrev.add_reaction(r_id)
                mm_irrev.add_reaction(r2_id)
                split_reversible.append((r_id, r2_id))
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
