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
# Copyright 2017-2018  Keith Dufault-Thompson <keitht547@uri.edu>

from __future__ import unicode_literals
import logging
from ..command import SolverCommandMixin, MetabolicMixin, Command, ObjectiveMixin
from ..reaction import Direction, Reaction, Compound
from .. import fluxanalysis
from ..lpsolver import lp
from .. import lpsolver
from ..datasource.reaction import parse_reaction
from decimal import Decimal
import copy
import csv
import math
import random
import yaml
import argparse
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                    DictReactionEntry as ReactionEntry,
                    DictCompartmentEntry as CompartmentEntry)
from collections import Counter, defaultdict
from psamm.randomsparse import  GeneDeletionStrategy, random_sparse_return_all, get_gene_associations
from six import iteritems
logger = logging.getLogger(__name__)


class TMFACommand(MetabolicMixin, SolverCommandMixin, ObjectiveMixin, Command):
    """Run TMFA on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', help='single column list of reactions to exclude from the analysis', type=file)
        parser.add_argument(
            '--lump-file', help='tab separate file, lump ID, lump dG_r, comma sep list of reactions in lump.',
            type=file)
        parser.add_argument(
            '--dgf-file', help='tab separated file containing estimated deltaGf values for all metabolites',
            type=argparse.FileType('rU'))
        parser.add_argument(
            '--set-concentrations', help='Tab seperated file with Reaction ID, [lower], [upper]', type=file)
        parser.add_argument('--transport-parameters', help='file containing parameters for transport rxn dgr', type=file)
        parser.add_argument('--rxn-conc-only', help='file containing reactions where deltaGf of intracellular/extracellular compounds should only be determined by ph difference', type=file)
        parser.add_argument('--scaled-compounds', help='compounds to scale deltaGf', type=file)
        parser.add_argument('--dgr-file', type=file, help='file containing estimated deltarG values for reactions in model.')
        parser.add_argument('--err', action='store_true', help='use error estimates when running TMFA')
        parser.add_argument('--random-addition', action='store_true', help='perform random reaction constraint addition in model')
        parser.add_argument('--hamilton', action='store_true', help='run model using Hamilton TMFA method')
        parser.add_argument('--conc-testing', action='store_true', help='peform random compound constraint addition in model')
        parser.add_argument('--temp', help='Temperature in Celsius')
        parser.add_argument('--tfba', action='store_true', help = 'run TMFA with tFBA like constraints (testing only)')
        parser.add_argument('--threshold', default=None, type=Decimal, help='value to fix biomass flux to during tmfa simulations (default = max biomass)')
        parser.add_argument('--verbose', action='store_true', help='print out all linear constraints and equalities from LP problem. ')
        parser.add_argument('--randomsparse', action='store_true', help='run randomsparse on reactions in model with TMFA constraints applied')
        parser.add_argument('--randomsparse_genes', action='store_true', help='run randomsparse on genes in model with TMFA constraints applied')
        parser.add_argument('--water', type=str, help='water compound ID without compartment (cpd_h2o, C00001)')
        parser.add_argument('--proton-in', type=str, help='id of proton compound with compartment (in the cell)')
        parser.add_argument('--proton-out', type=str, help='id of proton compound with compartment (outside of the cell)')
        parser.add_argument('--single-solution', type=str, choices=['fba', 'l1min','random'], help='get a single TMFA solution', default=None)
        parser.add_argument('--min-max-energy', type=file, help='', default=None)
        parser.add_argument('--config', type=file, help='Config file for TMFA settings')
        super(TMFACommand, cls).init_parser(parser)

    def run(self):
        """Run TMFA command."""
        solver = self._get_solver()
        objective = self._get_objective()
        mm = self._mm


        if self._args.config is None:
            logger.warning('No configuration file was provided with --config option')
            quit()
        else:
            config_dict = yaml.safe_load(self._args.config)


        # Parse exclude file provided through command line.
        base_exclude_list = []
        if config_dict.get('exclude') is not None:
            for line in open(config_dict.get('exclude')).readlines():
                base_exclude_list.append(line.rstrip())

        # Parse set of reactions that will have their deltaG constrainted only by deltaG of the transport component.
        ph_difference_rxn = []

        if config_dict.get('rxn_conc_only') is not None:
            for line in open(config_dict.get('rxn_conc_only')).readlines():
                rxn = line.rstrip()

                ph_difference_rxn.append(rxn)
                ph_difference_rxn.append('{}_forward'.format(rxn))
                ph_difference_rxn.append('{}_reverse'.format(rxn))

        # Parse out the file of lumped reactions or set the lumps dictionaries to be empty
        if config_dict.get('lumped_reactions') is not None:
            lumpid_to_rxn, rxn_lumpid, lump_to_rxnids, lump_to_rxnids_dir = lump_parser(open(config_dict.get('lumped_reactions')))
        else:
            logger.info('No lumped reactions were provided')
            lumpid_to_rxn = {}
            lump_to_rxnids = {}
            lump_to_rxnids_dir = {}

        # Add exchange reactions to the exclude list
        for reaction in mm.reactions:
            if mm.is_exchange(reaction):
                base_exclude_list.append(reaction)

        # Remove following two lines in future. these lines will add ph_difference list to the base exclude list.
        base_exclude_list = base_exclude_list + ph_difference_rxn
        ph_difference_rxn = []

        # exclude lump_list is a list of all excluded reactions and lump reactions.
        # exclude_lump_unkown list is a list of all excluded reactions, lump reactions and lumped reactions.
        exclude_lump_list = set(base_exclude_list)
        exclude_unkown_list = set(base_exclude_list)
        for lump_id, sub_rxn_list in lump_to_rxnids.iteritems():
            exclude_lump_list.add(lump_id)
            for sub_rxn in sub_rxn_list:
                exclude_unkown_list.add(sub_rxn)
        # Make list of all excluded, lumped, and unknown dgr.
        exclude_lump_unkown = exclude_unkown_list.union(exclude_lump_list)

        cpd_conc_dict = {}
        # Parse file containing set concentrations and ranges.
        if config_dict.get('concentrations') is not None:
            for row in csv.reader(open(config_dict.get('concentrations')), delimiter=str('\t')):
                cpd, lower, upper = row
                cpd_conc_dict[cpd] = [lower, upper]

        # Add lump reactions to the mm_irreversible model
        for lump_id, rxn in lumpid_to_rxn.iteritems():
            reaction = parse_reaction(rxn)
            mm.database.set_reaction(lump_id, reaction)
            mm.add_reaction(lump_id)
            mm.limits[lump_id].lower = 0
            mm.limits[lump_id].upper = 0




        # make an irreversible version of the metabolic model:
        mm_irreversible, split_reversible, reversible_lump_to_rxn_dict = make_irreversible(mm, exclude_lump_unkown, lump_to_rxnids_dir, self._args.hamilton)
        split_reactions = []
        for_rev_reactions = []
        for (i,j) in split_reversible:
            for_rev_reactions.append(i)
            for_rev_reactions.append(j)
            split_reactions.append(i[:-8])

        model = self._model

        gene_dictionary = defaultdict(list)
        # for rx in split_reactions:
        # 	print(rx)
        for rx in model.reactions:
            if rx.id in split_reactions:
                f = ReactionEntry(dict(id='{}_forward'.format(rx.id), genes=rx.genes))
                r = ReactionEntry(dict(id='{}_reverse'.format(rx.id), genes=rx.genes))
                model.reactions.add_entry(f)
                model.reactions.add_entry(r)

        # update these lists for the new reversible lumps and constituent reactions.
        # exclude lump_list is a list of all excluded reactions and lump reactions.
        # exclude_lump_unkown list is a list of all excluded reactions, lump reactions and lumped reactions.
        for lump_id, sub_rxn_list in reversible_lump_to_rxn_dict.iteritems():
            exclude_lump_list.add(lump_id)
            for sub_rxn in sub_rxn_list:
                exclude_unkown_list.add(sub_rxn)
        # Make list of all excluded, lumped, and unknown dgr.
        exclude_lump_unkown = exclude_unkown_list.union(exclude_lump_list)



        if config_dict.get('deltaGf') is not None:
            # Parse the deltaGf values for all metaobolites from a supplied file.
            dgf_dict = parse_dgf(mm_irreversible, open(config_dict.get('deltaGf')))
            # Calculate the deltaGr values for all reactions
            dgr_dict = calculate_dgr(mm_irreversible, dgf_dict, exclude_unkown_list, open(config_dict.get('transporters')), ph_difference_rxn, open(config_dict.get('scaled_compounds')))
            print('using dgf file')
        if config_dict.get('deltaG') is not None:
            dgr_dict = parse_dgr_file(open(config_dict.get('deltaG')), mm_irreversible)
            print('using dgr file')

        if config_dict.get('transporters') is not None:
            transport_parameters = parse_tparam_file(open(config_dict.get('transporters')))



        prob = solver.create_problem()
        prob.cplex.parameters.threads.set(1)
        prob.integrality_tolerance.value = 0
        prob.cplex.parameters.emphasis.numerical.value = 1

        self._v = v = prob.namespace(name='flux')
        self._zi = zi = prob.namespace(name='zi')
        self._dgri = dgri = prob.namespace(name='dgri')
        self._xij = xij = prob.namespace(name='xij')
        for reaction in mm_irreversible.reactions:
            lower, upper = mm_irreversible.limits[reaction]
            v.define([reaction], lower=lower, upper=upper, types=lp.VariableType.Continuous)
            zi.define([reaction], lower=int(0), upper=int(1), types=lp.VariableType.Binary)
            dgri.define([reaction], lower=-1000, upper=1000, types=lp.VariableType.Continuous)

        massbalance_lhs = {compound: 0 for compound in mm_irreversible.compounds}
        for spec, value in iteritems(mm_irreversible.matrix):
            compound, reaction_id = spec
            massbalance_lhs[compound] += v(reaction_id) * value
        for compound, lhs in iteritems(massbalance_lhs):
            prob.add_linear_constraints(lhs == 0)

        def get_var_bound(var, objective_sense):
            prob.set_objective(var)
            result = prob.solve_unchecked(objective_sense)
            if not result.success:
                logger.error('Solution not optimal: {}'.format(result.status))
                quit()
            return result.get_value(var)

        testing_list_tmp = list(mm_irreversible.reactions)
        cp_list = []
        for cp in mm_irreversible.compounds:
            cp_list.append(str(cp))
            xij.define([str(cp)], lower=-50, upper=50, types=lp.VariableType.Continuous)

        if self._args.random_addition:
            testing_list = list(mm_irreversible.reactions)
            random.shuffle(testing_list)
            testing_list_tmp = []
            for rx in testing_list:
                testing_list_iter = testing_list_tmp + [rx]
                logger.info('current testing list: {}'.format(testing_list_tmp))
                logger.info('testing reaction: {}'.format(rx))

                prob = solver.create_problem()
                prob.cplex.parameters.threads.set(1)
                prob.integrality_tolerance.value = 0
                prob.cplex.parameters.emphasis.numerical.value = 1

                self._v = v = prob.namespace(name='flux')
                self._zi = zi = prob.namespace(name='zi')
                self._dgri = dgri = prob.namespace(name='dgri')
                self._xij = xij = prob.namespace(name='xij')
                for reaction in mm_irreversible.reactions:
                    lower, upper = mm_irreversible.limits[reaction]
                    v.define([reaction], lower=lower, upper=upper, types=lp.VariableType.Continuous)
                    zi.define([reaction], lower=int(0), upper=int(1), types=lp.VariableType.Binary)
                    dgri.define([reaction], lower=-1000, upper=1000, types=lp.VariableType.Continuous)

                massbalance_lhs = {compound: 0 for compound in mm_irreversible.compounds}
                for spec, value in iteritems(mm_irreversible.matrix):
                    compound, reaction_id = spec
                    massbalance_lhs[compound] += v(reaction_id) * value
                for compound, lhs in iteritems(massbalance_lhs):
                    prob.add_linear_constraints(lhs == 0)

                cp_list = []
                for cp in mm_irreversible.compounds:
                    cp_list.append(str(cp))
                    xij.define([str(cp)], lower=-50, upper=50, types=lp.VariableType.Continuous)

                prob, cpd_xij_dict = add_conc_constraints(self, prob, cpd_conc_dict, cp_list)

                prob = add_reaction_constraints(self, prob, mm_irreversible, exclude_lump_list,
                                                        exclude_unkown_list,
                                                        exclude_lump_unkown, dgr_dict, reversible_lump_to_rxn_dict,
                                                        split_reversible, transport_parameters, testing_list_iter,
                                                        self._args.scaled_compounds, self._args.water,
                                                        self._args.proton_in,
                                                        self._args.proton_out, self._args.temp, self._args.err)
                try:
                    biomass = get_var_bound(v(self._get_objective()), lp.ObjectiveSense.Maximize)
                    logger.info('Current Biomass: {}'.format(biomass))
                    if biomass < self._args.threshold:
                        continue
                    else:
                        testing_list_tmp.append(rx)
                except:
                    continue

            for rx in testing_list:
                if rx not in testing_list_tmp:
                    print('{}\tBad Constraint'.format(rx))
                else:
                    print('{}\tGood Constraint'.format(rx))
            quit()


        prob, cpd_xij_dict = add_conc_constraints(self, prob, cpd_conc_dict, cp_list)

        prob = add_reaction_constraints(self, prob, mm_irreversible, exclude_lump_list,
                                                exclude_unkown_list,
                                                exclude_lump_unkown, dgr_dict, reversible_lump_to_rxn_dict,
                                                split_reversible, transport_parameters, testing_list_tmp,
                                                self._args.scaled_compounds, self._args.water, self._args.proton_in,
                                                self._args.proton_out, self._args.temp, self._args.err)

        if self._args.threshold is not None:
            prob.add_linear_constraints(v(self._get_objective()) == self._args.threshold)
            logger.info('Set biomass based on threshold to {}'.format(self._args.threshold))
        else:
            max_biomass = get_var_bound(v(self._get_objective()), lp.ObjectiveSense.Maximize)
            prob.add_linear_constraints(v(self._get_objective()) == max_biomass)
            logger.info('Set biomass based on max biomass to {}'.format(max_biomass))

        if self._args.single_solution is not None:
            print_fba(self, prob, mm_irreversible, cp_list)
        elif self.min_max_energy is not None:
            energy_min_max(self, prob, mm_irreversible, self.min_max_energy)
        else:
            rand_reactions = [m for m in mm_irreversible.reactions]
            random.shuffle(rand_reactions)
            for step, reaction in enumerate(sorted(rand_reactions)):
                logger.info('Testing Reaction {}/{}'.format(step, len(rand_reactions)))
                min_flux = get_var_bound(v(reaction), lp.ObjectiveSense.Minimize)
                max_flux = get_var_bound(v(reaction), lp.ObjectiveSense.Maximize)
                print('Flux\t{}\t{}\t{}'.format(reaction, min_flux, max_flux))
                if reaction not in exclude_unkown_list:
                    min_dgr = get_var_bound(dgri(reaction), lp.ObjectiveSense.Minimize)
                    max_dgr = get_var_bound(dgri(reaction), lp.ObjectiveSense.Maximize)
                    min_zi = get_var_bound(zi(reaction), lp.ObjectiveSense.Minimize)
                    max_zi = get_var_bound(zi(reaction), lp.ObjectiveSense.Maximize)
                    print('DGR\t{}\t{}\t{}'.format(reaction, min_dgr, max_dgr))
                    print('Zi\t{}\t{}\t{}'.format(reaction, min_zi, max_zi))
                else:
                    print('DGR\t{}\t{}\t{}'.format(reaction, 'NA', 'NA'))
                    print('Zi\t{}\t{}\t{}'.format(reaction, 'NA', 'NA'))
            for step, cpd in enumerate(sorted(cp_list)):
                logger.info('Testing Compound {}\{}'.format(step, len(cp_list)))
                min_cpd = get_var_bound(xij(cpd), lp.ObjectiveSense.Minimize)
                max_cpd = get_var_bound(xij(cpd), lp.ObjectiveSense.Maximize)
                print('CONC\t{}\t{}\t{}'.format(cpd, min_cpd, max_cpd))

            if self._args.verbose:
                index_dict_vars = {}
                for i, j in prob._variables.iteritems():
                    index_dict_vars[j] = str(i)
                for key, value in index_dict_vars.iteritems():
                    print('## LP variable name, lp var lower bound, lp var upper bound, var type')
                    print(value, key, prob.cplex.variables.get_lower_bounds(key),
                          prob.cplex.variables.get_upper_bounds(key),
                          prob.cplex.variables.get_types(key))

                for i in prob.cplex.linear_constraints.get_names():
                    linear_constraint = prob.cplex.linear_constraints.get_rows(i)
                    vars = linear_constraint.ind
                    tmp_vars = []
                    for var in vars:
                        tmp_vars.append(index_dict_vars[var])
                    print('## Raw sparse pair from cplex')
                    print(linear_constraint)
                    print('## lhs variables, coefficients')
                    print(tmp_vars, linear_constraint.val)
                    print('## rhs value')
                    print(TMFA_Problem.prob.cplex.linear_constraints.get_rhs(i))
                    print('## rhs equation sense (L = less than or equal to, G = greater than or equal to, E = equal to')
                    print(TMFA_Problem.prob.cplex.linear_constraints.get_senses(i))
                    print('condensed LP constraint')
                    equation = []
                    for j in range(0, len(vars), 1):
                        equation.append('{}*{}'.format(tmp_vars[j], linear_constraint.val[j]))
                    sense = prob.cplex.linear_constraints.get_senses(i)
                    if sense == 'L':
                        sign = '<='
                    elif sense == 'G':
                        sign = '>='
                    elif sense == 'E':
                        sign = '=='
                    print('{} {} {}'.format(' + '.join(equation), sign,
                                            prob.cplex.linear_constraints.get_rhs(i)))
                    print('-------------------------------------------------------------------')
        quit()


def energy_min_max(self, problem, mm, enfile):
    stoichs = {}
    for row in csv.reader(enfile, delimiter='\t'):
        stoichs['{}_reverse'.formrat(row[0])] = -1 * float(row[1])
        stoichs['{}_forward'.formrat(row[0])] = float(row[1])
        stoichs['{}'.formrat(row[0])] = float(row[1])

    objective = None
    for rxn in stoichs:
        if rxn in mm.reactions:
            if objective is None:
                obejctive = (self._v(rxn) * stoichs)
            else:
                objective += (self._v(rxn) * stoichs)
    problem.set_objective(objective)
    result = problem.solve_unchecked(lp.ObjectiveSense.Maximize)
    if not result.success:
        logger.error('Solution not optimal: {}'.format(result.status))
        quit()
    max = result.get_value(objective)
    result = problem.solve_unchecked(lp.ObjectiveSense.Minimize)
    if not result.success:
        logger.error('Solution not optimal: {}'.format(result.status))
        quit()
    min = result.get_value(objective)
    print('min_sum_flux\t{}'.format(min))
    print('max_sum_flux\t{}'.format(max))


def print_fba(self, prob, mm, cp_list):
    prob.set_objective(self._v(self._get_objective()))
    result = prob.solve_unchecked(lp.ObjectiveSense.Maximize)
    if not result.success:
        logger.error('Solution not optimal: {}'.format(result.status))
        quit()
    if self._args.single_solution == 'fba':
        for rxn in mm.reactions:
            print('Flux\t{}\t{}'.format(rxn, result.get_value(self._v(rxn))))
            print('DGR\t{}\t{}'.format(rxn, result.get_value(self._dgri(rxn))))
            print('Zi\t{}\t{}'.format(rxn, result.get_value(self._zi(rxn))))
        for cp in cp_list:
            print('CPD\t{}\t{}'.format(cp, result.get_value(self._xij(cp))))
    elif self._args.single_solution == 'l1min':
        self._z = prob.namespace(mm.reactions, lower=0)
        z = self._z.set(mm.reactions)
        v = self._v.set(mm.reactions)
        prob.add_linear_constraints(z >= v, v >= -z)

        objective = self._z.expr(
            (reaction_id, -1)
            for reaction_id in mm.reactions)
        prob.set_objective(objective)
        result = prob.solve_unchecked(lp.ObjectiveSense.Maximize)
        if not result.success:
            logger.error('Solution not optimal: {}'.format(result.status))
            quit()
        for rxn in mm.reactions:
            print('Flux\t{}\t{}'.format(rxn, result.get_value(self._v(rxn))))
            print('DGR\t{}\t{}'.format(rxn, result.get_value(self._dgri(rxn))))
            print('Zi\t{}\t{}'.format(rxn, result.get_value(self._zi(rxn))))
        for cp in cp_list:
            print('CPD\t{}\t{}'.format(cp, result.get_value(self._xij(cp))))
        print(result.get_value(objective))
    elif self._args.single_solution == 'random':
        optimize = None
        for rxn_id in mm.reactions:
            if optimize is None:
                optimize = self._v(rxn_id) * random.random()
            else:
                optimize += self._v(rxn_id) * random.random()
        prob.add_linear_constraints(self._v(self._get_objective()) == result.get_value(self._v(self._get_objective())))

        prob.set_objective(optimize)

        result = prob.solve_unchecked()
        for rxn in mm.reactions:
            print('Flux\t{}\t{}'.format(rxn, result.get_value(self._v(rxn))))
            print('DGR\t{}\t{}'.format(rxn, result.get_value(self._dgri(rxn))))
            print('Zi\t{}\t{}'.format(rxn, result.get_value(self._zi(rxn))))
        for cp in cp_list:
            print('CPD\t{}\t{}'.format(cp, result.get_value(self._xij(cp))))


def lump_parser(lump_file):
    """Parses a supplied file and returns dictionaries containing lump information.

    The supplied file should be in a tab separated table in the format of
    lumpID  lump_deltaG rxn1,rxn2,rxn3  lumpRXN
    returns two dictionaries, lump_to_dg where the keys are lump reaction ids
    and the values are deltaG values, and the rxn_to_lump_id dict where the keys
    are reaction ids and the values are what lump reaction they belong to.
    """
    rxn_to_lump_id = {}
    lump_to_rxn = {}
    lump_to_rxnids = {}
    lump_to_rxnids_dir = {}
    for row in csv.reader(lump_file, delimiter=str('\t')):
        lump_id, lump_rxn_list, lump_rxn = row
        rx_dir_list = []
        rx_list = []
        for i in lump_rxn_list.split(','):
            subrxn, dir = i.split(':')
            rx_dir_list.append((subrxn, dir))
            rx_list.append(subrxn)
            rxn_to_lump_id[i] = lump_id
        lump_to_rxn[lump_id] = lump_rxn
        lump_to_rxnids[lump_id] = rx_list
        lump_to_rxnids_dir[lump_id] = rx_dir_list
    return lump_to_rxn, rxn_to_lump_id, lump_to_rxnids, lump_to_rxnids_dir


def parse_tparam_file(file):
    t_param = {}
    if file is not None:
        for row in csv.reader(file, delimiter=str('\t')):
            rxn, c, h = row
            t_param[rxn] = (c, h)
            t_param['{}_forward'.format(rxn)] = (c, h)
            t_param['{}_reverse'.format(rxn)] = (-Decimal(c), -Decimal(h))
    return t_param


def parse_dgf(mm, dgf_file):
    """A function that will parse a supplied deltaG of formation file.

    Compound IDs in this file do not need to contain the compartments.
    compound deltaGf values should be in Kcal/mol

    Args:
        mm: a metabolic model object
        dgf_file: a file that containing 2 columns of compound ids and deltaGf values
    """

    cpd_dgf_dict = {}
    for row in csv.reader(dgf_file, delimiter=str('\t')):
        for cpt in mm.compartments:
            try:
                dg = Decimal(row[1])
                derr = Decimal(row[2])
                cpd_dgf_dict[Compound(row[0], cpt)] = (dg, derr)
            except:
                logger.info('Compound {} has an assigned detltGf value of {}. This is not an number and will be treated\
                            as a missing value.'.format(Compound(row[0], cpt), row[1]))
    return cpd_dgf_dict


def add_conc_constraints(self, problem, cpd_conc_dict, cp_list):
    # Water needs to be excluded from these concentration constraints.
    # excluded_compounds = ['h2o[c]', 'h2o[e]', 'h[c]', 'h[e]']
    excluded_compounds = ['cpd_h2o[c]', 'cpd_h2o[e]', 'cpd_h[c]', 'cpd_h[e]', 'cpd_h2o[p]', 'cpd_h[p]', 'C00080[c]', 'C00080[e]', 'C00001[c]', 'C00001[e]']

    # excluded_compounds = ['h2o[c]', 'h2o[e]', 'pe_ec[c]', '12dgr_ec[c]', 'agpc_EC[c]', 'agpe_EC[c]', 'agpg_EC[c]', 'cdpdag_EC[c]',
    #                       'clpn_EC[c]', 'pa_EC[c]', 'pc_EC[c]', 'pg_EC[c]', 'pe_EC[c]', 'pgp_EC[c]', 'ps_EC[c]']
    cpdid_xij_dict = {}
    # print(cp_list)
    # for cp in problem._model.compounds:
    for cp in cp_list:
        # print(str(cp))
        # define concentration variable for compound.
        # problem.prob.define(str(cp))
        # var = problem.prob.var(str(cp))
        cpdid_xij_dict[str(cp)] = self._xij(str(cp))
        var = self._xij(str(cp))
        # Define default constraints for anything not set in the conc file
        # if str(cp) in cp_list:
        if str(cp) not in cpd_conc_dict.keys():
                if str(cp) not in excluded_compounds:
                    # print('Default Conc Constraint Applied\t{}\t{}\t{}'.format(str(cp), math.log(0.0001), math.log(0.02)))
                    # Add concentration constraints as the ln of the concentration (M).
                    # problem._prob.add_linear_constraints(var >= math.log(0.00001))
                    # problem._prob.add_linear_constraints(var >= math.log(0.0000001))
                    problem.add_linear_constraints(var >= math.log(0.00001))
                    problem.add_linear_constraints(var <= math.log(0.02))
                    # print('default', str(cp))
        elif str(cp) in cpd_conc_dict.keys():
            if str(cp) not in excluded_compounds:
                conc_limits = cpd_conc_dict[str(cp)]
                if conc_limits[0] > conc_limits[1]:
                    logger.error('lower bound for {} concentration higher than upper bound'.format(conc_limits))
                    quit()
                if Decimal(conc_limits[0]) == Decimal(conc_limits[1]):
                    problem.add_linear_constraints(var == math.log(Decimal(conc_limits[0])))
                    # Constraints to allow the concentration constraints on the problem to be more flexible
                    # problem._prob.add_linear_constraints(var >= math.log(Decimal(conc_limits[0]))-1)
                    # problem._prob.add_linear_constraints(var <= math.log(Decimal(conc_limits[1])+1))
                else:
                    problem.add_linear_constraints(var >= math.log(Decimal(conc_limits[0])))
                    problem.add_linear_constraints(var <= math.log(Decimal(conc_limits[1])))
                    # Constraints to allow the concentration constraints to be more flexible
                    # problem._prob.add_linear_constraints(var >= math.log(Decimal(conc_limits[0]))-1)
                    # problem._prob.add_linear_constraints(var <= math.log(Decimal(conc_limits[1]))+1)
                lower = math.log(Decimal(conc_limits[0]))
                upper = math.log(Decimal(conc_limits[1]))
                # logger.info('Non default Conc Constraints Applied\t{}\t{}\t{}'.format(str(cp), lower, upper))
    return problem, cpdid_xij_dict


def parse_dgr_file(dgr_file, mm):
    def is_number(val):
        try:
            float(val)
            return True
        except ValueError:
            return False
    dgr_dict = {}
    for row in csv.reader(dgr_file, delimiter=str('\t')):
        rxn, dgr, err = row
        if is_number(dgr):

            if is_number(err):
                err = Decimal(err)
            else:
                err = Decimal(2)
            if rxn in mm.reactions:
                dgr_dict[rxn] = (Decimal(dgr), err)
            elif '{}_forward'.format(rxn) in mm.reactions:
                dgr_dict['{}_forward'.format(rxn)] = (Decimal(dgr), err)
                dgr_dict['{}_reverse'.format(rxn)] = (-Decimal(dgr), err)
        else:
            logger.info('Reaction {} was provided with dgr value of {}. Treating as an unknown value.'.format(rxn, dgr))
    return dgr_dict


def calculate_dgr(mm, dgf_dict, excluded_reactions, transport_parameters, ph_difference_rxn, scaled_compounds):
    F = Decimal(0.02306)
    dpsi = Decimal(-130)
    R = Decimal(1.9858775 / 1000)
    # T = Decimal(303.15)


    dph = Decimal(0.4)
    t_param = {}
    # for row in csv.reader(transport_parameters, delimiter=str('\t')):
    # 	rxn, c, h = row
    # 	t_param[rxn] = (c, h)
    # 	t_param['{}_forward'.format(rxn)] = (c, h)
    # 	t_param['{}_reverse'.format(rxn)] = (-Decimal(c), -Decimal(h))
    dgf_scaling = {}
    if scaled_compounds is not None:
        scaled_compounds.seek(0)
        for row in csv.reader(scaled_compounds, delimiter=str('\t')):
            dgf_scaling[row[0]] = Decimal(row[1])

    dgr_dict = {}
    for reaction in mm.reactions:
        if reaction not in excluded_reactions:
            dgr = 'NA'
            dgerr = 0
            rxn = mm.get_reaction(reaction)
            if any(dgf_dict.get(j[0]) is None for j in rxn.compounds):
                # print('Reaction DGR Bad\t{}\t{}'.format(reaction, 'NA'))
                for j in rxn.compounds:
                    if j[0] not in dgf_dict.keys():
                        print(j[0])
                if reaction not in ph_difference_rxn:
                    logger.error('Reaction {} contains at least 1 compound with an unknown deltaGf value'.format(reaction))
                    # print('UnknownDGF\t{}'.format(reaction))
                    quit()
            else:
                dgr = 0
                # Make a variable dgf_sum that represents the sum of sij *  (stoichiometry *
                # deltaGf for reaction j.
                for cpd in rxn.compounds:
                    # print('dgr current', dgr)
                    if str(cpd[0]) in dgf_scaling.keys():
                        dgscale = dgf_scaling[str(cpd[0])]
                    else:
                        dgscale = 1
                    (dg, dge) = dgf_dict[cpd[0]]
                    # print('{}'.format(dgscale))
                    dgs = Decimal(dg) * (Decimal(cpd[1])*dgscale)
                    # print(cpd[0], dg, cpd[1], dgs)
                    # print(cpd, dg, Decimal(cpd[1]), dgs)
                    dgr += dgs
                    dgerr += Decimal(cpd[1])*Decimal(dgscale) * dge
            # if reaction in t_param.keys():
            # 	(c, h) = t_param[reaction]
            # 	dgr += (Decimal(c) * F * dpsi) - (2.3 * Decimal(h) * R * T * dph)
            dgr_dict[reaction] = (dgr, dgerr)
            # dgr_dict['{}_forward'.format(reaction)] = (dgr)#, dgerr)
            # dgr_dict['{}_reverse'.format(reaction)] = (-dgr)#, dgerr)
            # print('Reaction DGR\t{}\t{}'.format(reaction, dgr))
    # quit()
    return dgr_dict


def detect_transporter(reaction):
    compound_counter = Counter()
    for compound in reaction.compounds:
        compound_counter[str(compound)] += 1
    over_list = [cp for cp, value in compound_counter.iteritems() if value > 1]
    if len(over_list) > 0:
        return True, over_list
    else:
        return False, []


def make_irreversible(mm, exclude_list, lump_rxn_dir, all_reversible):
    """Returns a new metabolic model with only irreversible reactions and list of split reactions.

    This function will find every reversible reaction in the model and split it into two reactions
    with the {rxnid}_forward or {rxnid}_reverse as the IDs.

    Args:
        mm: A metabolicmodel object
        exclude_list: list of reactions to exclude in TMFA simulation
        reversible: if True make all reactions in model reversible.
    """
    split_reversible = []
    mm_irrev = mm.copy()
    lumped_rxns = []
    new_lump_rxn_dict = {}
    for rxn in mm.reactions:


        upper = mm.limits[rxn].upper
        lower = mm.limits[rxn].lower
        mm_irrev.limits[rxn].upper = upper
        mm_irrev.limits[rxn].lower = lower

        # print(rxn)
        # print(mm.limits[rxn])
        # print(mm_irrev.limits[rxn])

        reaction = mm_irrev.get_reaction(rxn)
        if rxn not in exclude_list:
            # Allow for making all reversible reactions into a _forward and _reverse split reaction pair
            r = Reaction(Direction.Forward, reaction.left, reaction.right)
            r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
            r_id = str('{}_forward'.format(rxn))
            r2_id = str('{}_reverse'.format(rxn))
            if reaction.direction == Direction.Forward:
                if all_reversible is False:
                    continue
                else:
                    mm_irrev.remove_reaction(rxn)
                    mm_irrev.database.set_reaction(r_id, r)
                    mm_irrev.database.set_reaction(r2_id, r2)
                    mm_irrev.add_reaction(r_id)
                    mm_irrev.add_reaction(r2_id)
                    split_reversible.append((r_id, r2_id))
                    # Reapply the limits here from the original reaction limits. (for and rev)
            elif reaction.direction == Direction.Both:
                mm_irrev.remove_reaction(rxn)
                mm_irrev.database.set_reaction(r_id, r)
                mm_irrev.database.set_reaction(r2_id, r2)
                mm_irrev.add_reaction(r_id)
                mm_irrev.add_reaction(r2_id)
                split_reversible.append((r_id, r2_id))
            # print (upper, lower)
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

        elif rxn in lump_rxn_dir.keys():
            final_sub_rxn_list = []
            sub = lump_rxn_dir[rxn]
            check = 0
            for (x,y) in sub:
                rn = mm_irrev.get_reaction(x)
                if rn.direction != Direction.Both:
                    check += 1
            if reaction.direction == Direction.Forward or check != 0:
                mm_irrev.limits[rxn].upper = 0
                mm_irrev.limits[rxn].lower = 0
                sub_rxn_list = lump_rxn_dir[rxn]
                for entry in sub_rxn_list:
                    (subrxn, dir)= entry
                    final_sub_rxn_list.append(subrxn)
                new_lump_rxn_dict[rxn] = final_sub_rxn_list
            elif reaction.direction == Direction.Both:
                # split the lump reaction itself
                lumped_rxns.append(rxn)
                r = Reaction(Direction.Forward, reaction.left, reaction.right)
                r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
                r_id = str('{}_forward'.format(rxn))
                r2_id = str('{}_reverse'.format(rxn))
                mm_irrev.remove_reaction(rxn)
                mm_irrev.database.set_reaction(r_id, r)
                mm_irrev.database.set_reaction(r2_id, r2)
                mm_irrev.add_reaction(r_id)
                mm_irrev.add_reaction(r2_id)
                split_reversible.append((r_id, r2_id))

                sub_rxn_list = lump_rxn_dir[rxn]
                for_sub_rxn_list = []
                rev_sub_rxn_list = []
                mm_irrev.limits[r_id].upper = 0
                mm_irrev.limits[r_id].lower = 0
                mm_irrev.limits[r2_id].upper = 0
                mm_irrev.limits[r2_id].lower = 0
                for entry in sub_rxn_list:
                    (subrxn, dir) = entry
                    dir = int(dir)
                    lumped_rxns.append(subrxn)
                    subreaction = mm.get_reaction(subrxn)
                    sub_r1 = Reaction(Direction.Forward, subreaction.left, subreaction.right)
                    sub_r1_id = '{}_forward'.format(subrxn)
                    sub_r2 = Reaction(Direction.Forward, subreaction.right, subreaction.left)
                    sub_r2_id = '{}_reverse'.format(subrxn)

                    split_reversible.append((sub_r1_id, sub_r2_id))
                    if dir == 1:
                        for_sub_rxn_list.append(sub_r1_id)
                        rev_sub_rxn_list.append(sub_r2_id)
                        mm_irrev.database.set_reaction(sub_r1_id, sub_r1)
                        mm_irrev.add_reaction(sub_r1_id)
                        mm_irrev.database.set_reaction(sub_r2_id, sub_r2)
                        mm_irrev.add_reaction(sub_r2_id)
                    elif dir == -1:
                        for_sub_rxn_list.append(sub_r2_id)
                        rev_sub_rxn_list.append(sub_r1_id)
                        mm_irrev.database.set_reaction(sub_r1_id, sub_r1)
                        mm_irrev.add_reaction(sub_r1_id)
                        mm_irrev.database.set_reaction(sub_r2_id, sub_r2)
                        mm_irrev.add_reaction(sub_r2_id)
                    mm_irrev.limits[sub_r1_id].lower = 0
                    mm_irrev.limits[sub_r1_id].upper = 100
                    mm_irrev.limits[sub_r2_id].lower = 0
                    mm_irrev.limits[sub_r2_id].upper = 100
                new_lump_rxn_dict[r_id] = for_sub_rxn_list
                new_lump_rxn_dict[r2_id] = rev_sub_rxn_list

    for rxn in lumped_rxns:
        mm_irrev.remove_reaction(rxn)

    return mm_irrev, split_reversible, new_lump_rxn_dict


def add_reaction_constraints(self, problem, mm, exclude_lumps, exclude_unknown, exclude_lumps_unknown, dgr_dict,
                             lump_rxn_list, split_rxns, transport_parameters, testing_list, scaled_compounds, water, hin, hout, temp, err_est=False, hamilton=False):

    dgf_scaling = {}
    if scaled_compounds is not None:
        scaled_compounds.seek(0)
        for row in csv.reader(scaled_compounds, delimiter=str('\t')):
            dgf_scaling[row[0]] = Decimal(row[1])

    R = Decimal(8.3144621 / 1000) # kJ/mol
    # R = Decimal(1.9858775 / 1000) # kcal/mol

    # T = Decimal(303.15)
    # T = Decimal(293.15) # 20 C
    # T = Decimal(288.15) # 15 C
    # T = Decimal(277.15)  # 4 C
    T = Decimal(temp) + Decimal(273.15)
    k = 500
    epsilon = 0.000001
    # epsilon = 0
    # h_e = problem.prob.var(str('h[e]'))
    # h_e = problem.prob.var(str('cpd_h[e]'))

    # problem._prob.add_linear_constraints(h_e <= 11)
    # problem._prob.add_linear_constraints(h_e >= 4)

    # h_p = problem.prob.var(str('C00080[e]'))
    h_p = self._xij(str(hout))

    problem.add_linear_constraints(h_p <= 11)
    problem.add_linear_constraints(h_p >= 4)
    # problem._prob.add_linear_constraints(h_e == 7.4)
    # h_c = problem.prob.var(str('h[c]'))

    # h_c = problem.prob.var(str('C00080[c]'))
    h_c = self._xij(str(hin))

    problem.add_linear_constraints(h_c >= 4)
    problem.add_linear_constraints(h_c <= 11)
    delta_ph = (h_c - h_p)

    F = Decimal(0.02306)

    excluded_cpd_list = ['cpd_h2o[e]', 'cpd_h2o[c]', 'cpd_h[c]', 'cpd_h[e]', 'cpd_h[p]', 'cpd_h2o[p]',
                         'C00080[e]', 'C00080[c]', 'C00001[e]', 'C00001[c]']
    excluded_cpd_list.append(hin)
    excluded_cpd_list.append(hout)
    for cpt in mm.compartments:
        excluded_cpd_list.append('{}[{}]'.format(water, cpt))
    logger.info('Excluded compounds: {}'.format(','.join(excluded_cpd_list)))
    logger.info('Temperature: {}'.format(T))
    logger.info('using h in {}'.format(hin))
    logger.info('using h out {}'.format(hout))
    logger.info('using water {}'.format(water))
    split_list = []
    for (f, r) in split_rxns:
        split_list.append(f)
        split_list.append(r)

    # dgri_var_dict = {}
    # for reaction in mm.reactions:
    # 	if reaction not in exclude_unknown:
    # 		problem.prob.define('dgri_{}'.format(reaction), types=lp.VariableType.Continuous, lower=-1000, upper=1000)
    # 		dgri = problem.prob.var('dgri_{}'.format(reaction))
    # 		dgri_var_dict[reaction] = dgri

    for (f,r) in split_rxns:
        if f not in exclude_unknown:
            dgrif = self._dgri(f)
            dgrir = self._dgri(r)
            problem.add_linear_constraints(dgrif == dgrir * -1)

    new_excluded_reactions = []
    for reaction in mm.reactions:
        if reaction not in exclude_unknown:
            rxn = mm.get_reaction(reaction)
            rhs_check = 0
            lhs_check = 0
            for (cpd, stoich) in rxn.compounds:
                if stoich < 0:
                    if str(cpd) not in excluded_cpd_list:
                        lhs_check += 1
                if stoich > 0:
                    if str(cpd) not in excluded_cpd_list:
                        rhs_check += 1
            if rhs_check == 0 or lhs_check == 0:
                new_excluded_reactions.append(reaction)
            # define variables for vmax, dgri, zi, yi, and vi

            # problem.prob.define('zi_{}'.format(reaction), types=lp.VariableType.Binary, lower=int(0), upper=int(1))
            zi = self._zi(reaction)
            dgri = self._dgri(reaction)
            vi = self._v(reaction)
            vmax = mm.limits[reaction].upper
            if reaction in testing_list:
                if reaction in transport_parameters.keys():
                    (c, h) = transport_parameters[reaction]
                    # dph = -math.log(math.exp(problem.prob.var('h[e]'))) - - math.log(math.exp(problem.prob.var('h[c]')))
                    # dph = Decimal(0.4)
                    ddph = Decimal(-2.3)*Decimal(h)*R*T*delta_ph
                    # dpsi = (33.33*dph)-143.33
                    dpsi = Decimal(33.33) * delta_ph - Decimal(143.33)
                    ddpsi = dpsi * Decimal(c) * Decimal(F)
                    dgr_trans = ddph + ddpsi
                else:
                    dgr_trans = 0

                # print('Reaction DGRTRANS {}: {}'.format(reaction, dgr_trans))
                (dgr0, err) = dgr_dict[reaction]

                # print('Reaction DGR0 dict lookup value: {}\t{}'.format(reaction, dgr_dict[reaction]))
                ssxi = 0
                # If no error then use this line

                if err_est:
                    problem.prob.define('dgr_err_{}'.format(reaction), types=lp.VariableType.Continuous, lower=-1000,
                                        upper=1000)
                    dgr_err = problem.prob.var('dgr_err_{}'.format(reaction))
                    problem.add_linear_constraints(dgr_err <= 2*err)
                    problem.add_linear_constraints(dgr_err >= -2*err)
                else:
                    dgr_err = 0
                # problem._prob.add_linear_constraints(dgr_err <= 0)
                # problem._prob.add_linear_constraints(dgr_err >= -20)

                # print(reaction, 'error', 2*err)
                for (cpd, stoich) in rxn.compounds:
                    if str(cpd) not in excluded_cpd_list:
                        scale = dgf_scaling.get(str(cpd), 1)
                        ssxi += self._xij(str(cpd)) * Decimal(stoich) * scale
                        # print('ssxi calc for {} compound {}\t{}'.format(reaction, problem.prob.var(str(cpd)), stoich))
                # print('Reaction dgri constraint calculation\t{}\t{}={}+({}*{}*({}))'.format(reaction, dgri, dgr0, R, T, ssxi))

                problem.add_linear_constraints(dgri == dgr0 + (R * T * (ssxi)) + dgr_err + dgr_trans)
                # print('Reaction dgri raw constraint calculation {}: '.format(reaction), (dgri == dgr0 + (R * T * (ssxi)) + dgr_err))
                if hamilton:
                    problem.add_linear_constraints(dgri <= 300-epsilon)
                    problem.add_linear_constraints(dgri >= -300+epsilon)

                if reaction not in exclude_lumps_unknown:
                    if rhs_check != 0 and lhs_check != 0:
                        problem.add_linear_constraints(dgri - k + (k * zi) <= -epsilon)
                        problem.add_linear_constraints(vi <= zi * vmax)
                        # problem._prob.add_linear_constraints(vi >= int(0))



                if reaction in lump_rxn_list.keys():
                    yi = problem.prob.var('yi_{}'.format(reaction))
                    if reaction not in new_excluded_reactions:
                        # print('TEST THIS IS A TEST', reaction)
                        vi = problem.get_flux_var(reaction)
                        # print('TEST VI TEST', vi)
                        yi = problem.prob.var('yi_{}'.format(reaction))
                        dgri = problem.prob.var('dgri_{}'.format(reaction))
                        problem._prob.add_linear_constraints(vi == 0)
                        # print('Lumped Reaction flux = 0 constraint\t{}=0'.format(vi))
                        # print('Lumped Reaction flux = 0 raw constraint {}=0'.format(reaction), (vi == 0))
                        problem._prob.add_linear_constraints(dgri - (k * yi) <= - epsilon)
                        # print('Lumped Reaction feasibility constraint\t{}\t{}-{}*{}<={}'.format(reaction, dgri, k, yi, epsilon))
                        # print('Lumped reaction feasibility raw constraint {}: '.format(reaction), (dgri - (k * yi) <= - epsilon))
                        sub_rxn_list = lump_rxn_list[reaction]
                        sszi = 0
                        for sub_rxn in sub_rxn_list:
                            sszi += problem.prob.var('zi_{}'.format(sub_rxn))
                        problem._prob.add_linear_constraints(yi + sszi <= len(sub_rxn_list))
                # print('Lump component constraints\t{}\t{}+{}<={}'.format(reaction, yi, sszi, sub_rxn_list))
                # print('Lumped component raw constraint {} :'.format(reaction), (yi + sszi <= len(sub_rxn_list)))

    for (forward, reverse) in split_rxns:
        problem.add_linear_constraints(
            self._zi(forward) + self._zi(reverse) <= int(1))
            # print('Split reaction Zi constraints\t{}\t{}\t{}+{}<=1'.format(forward, reverse, problem.prob.var('zi_{}'.format(forward)), problem.prob.var('zi_{}'.format(reverse))))
            # print('Split reaction zi raw constraint {} :'.format(forward), (problem.prob.var('zi_{}'.format(forward)) + problem.prob.var('zi_{}'.format(reverse)) <= 1))
    return problem








    # 	vmax = mm.limits[reaction].upper
    # 	problem.prob.define('yi_{}'.format(reaction), types=lp.VariableType.Binary)
    # 	yi = problem.prob.var('yi_{}'.format(reaction))
    # 	vi = problem.get_flux_var(reaction)
    #
    #
    #
    # 	# add flux constraint linking vi and zi for all reactions except lumps
    # 	if reaction not in exclude_lumps:
    # 		problem._prob.add_linear_constraints(vi <= zi * vmax)
    # 		problem._prob.add_linear_constraints(vi >= 0)
    # 		# print('Reaction Zi Vi constraint\t{}\t{}-{}*{}<=0'.format(reaction, vi, zi, vmax))
    # 		# print('Reaction Zi Vi constraint {}: '.format(reaction), (vi - zi * vmax <= 0))
    # 	# add thermodynamic feasibility constraint for all reactions where dgr0 is known except for lumps
    # 	if reaction not in exclude_lumps_unknown:
    # 		if reaction in testing_list:
    # 			if rhs_check != 0 and lhs_check != 0:
    # 				problem._prob.add_linear_constraints(dgri - k + (k * zi) <= 0 - epsilon)
    # 				# print('Reaction thermo feasibility constraint\t{}\t{}-{}+{}*{}<=-{}'.format(reaction, dgri, k, k, zi, epsilon))
    # 				# print('Reaction thermo feasibility raw constraint {}: '.format(reaction), (dgri - k + (k * zi) <= 0 - epsilon))
    # 	# add constraint to calculate dgri based on dgr0 and the concentrations of the metabolites
    # 	if reaction not in exclude_unknown:
    # 		if reaction in transport_parameters.keys():
    # 			(c, h) = transport_parameters[reaction]
    # 			# dph = -math.log(math.exp(problem.prob.var('h[e]'))) - - math.log(math.exp(problem.prob.var('h[c]')))
    # 			# dph = Decimal(0.4)
    # 			ddph = Decimal(-2.3)*Decimal(h)*R*T*delta_ph
    # 			# dpsi = (33.33*dph)-143.33
    # 			dpsi = Decimal(33.33) * delta_ph - Decimal(143.33)
    # 			ddpsi = dpsi * Decimal(c) * Decimal(F)
    # 			dgr_trans = ddph + ddpsi
    # 		else:
    # 			dgr_trans = 0
    #
    # 		# print('Reaction DGRTRANS {}: {}'.format(reaction, dgr_trans))
    # 		(dgr0, err) = dgr_dict[reaction]
    #
    # 		# print('Reaction DGR0 dict lookup value: {}\t{}'.format(reaction, dgr_dict[reaction]))
    # 		ssxi = 0
    # 		problem.prob.define('dgr_err_{}'.format(reaction))
    # 		# If no error then use this line
    #
    # 		# If you want to use the error estimates for dgr values then uncomment these lines
    # 		if err_est:
    # 			dgr_err = problem.prob.var('dgr_err_{}'.format(reaction))
    # 			problem._prob.add_linear_constraints(dgr_err <= 2*err)
    # 			problem._prob.add_linear_constraints(dgr_err >= -2*err)
    # 		else:
    # 			dgr_err = 0
    # 		# problem._prob.add_linear_constraints(dgr_err <= 0)
    # 		# problem._prob.add_linear_constraints(dgr_err >= -20)
    #
    # 		# print(reaction, 'error', 2*err)
    # 		for (cpd, stoich) in rxn.compounds:
    # 			if str(cpd) not in excluded_cpd_list:
    # 				scale = dgf_scaling.get(str(cpd), 1)
    # 				ssxi += problem.prob.var(str(cpd)) * Decimal(stoich) * scale
    # 				# print('ssxi calc for {} compound {}\t{}'.format(reaction, problem.prob.var(str(cpd)), stoich))
    # 		# print('Reaction dgri constraint calculation\t{}\t{}={}+({}*{}*({}))'.format(reaction, dgri, dgr0, R, T, ssxi))
    #
    # 		problem._prob.add_linear_constraints(dgri == dgr0 + (R * T * (ssxi)) + dgr_err + dgr_trans)
    # 		# print('Reaction dgri raw constraint calculation {}: '.format(reaction), (dgri == dgr0 + (R * T * (ssxi)) + dgr_err))
    # 		if hamilton:
    # 			problem._prob.add_linear_constraints(dgri <= 300-epsilon)
    # 			problem._prob.add_linear_constraints(dgri >= -300+epsilon)
    # # add constraints for thermodynamic feasibility of lump reactions and to constrain their constituent reactions
    # for reaction in mm.reactions:
    # 	if reaction in lump_rxn_list.keys():
    # 		if reaction not in new_excluded_reactions:
    # 			# print('TEST THIS IS A TEST', reaction)
    # 			vi = problem.get_flux_var(reaction)
    # 			# print('TEST VI TEST', vi)
    # 			yi = problem.prob.var('yi_{}'.format(reaction))
    # 			dgri = problem.prob.var('dgri_{}'.format(reaction))
    # 			problem._prob.add_linear_constraints(vi == 0)
    # 			# print('Lumped Reaction flux = 0 constraint\t{}=0'.format(vi))
    # 			# print('Lumped Reaction flux = 0 raw constraint {}=0'.format(reaction), (vi == 0))
    # 			problem._prob.add_linear_constraints(dgri - (k * yi) <= - epsilon)
    # 			# print('Lumped Reaction feasibility constraint\t{}\t{}-{}*{}<={}'.format(reaction, dgri, k, yi, epsilon))
    # 			# print('Lumped reaction feasibility raw constraint {}: '.format(reaction), (dgri - (k * yi) <= - epsilon))
    # 			sub_rxn_list = lump_rxn_list[reaction]
    # 			sszi = 0
    # 			for sub_rxn in sub_rxn_list:
    # 				sszi += problem.prob.var('zi_{}'.format(sub_rxn))
    # 			problem._prob.add_linear_constraints(yi + sszi <= len(sub_rxn_list))
    # 			# print('Lump component constraints\t{}\t{}+{}<={}'.format(reaction, yi, sszi, sub_rxn_list))
    # 			# print('Lumped component raw constraint {} :'.format(reaction), (yi + sszi <= len(sub_rxn_list)))
    # # Add linear constraint to disallow solutions that use both the forward and reverse reaction from a split reaction.
    # for (forward, reverse) in split_rxns:
    # 	problem._prob.add_linear_constraints(problem.prob.var('zi_{}'.format(forward)) + problem.prob.var('zi_{}'.format(reverse)) <= 1)
    # 	# print('Split reaction Zi constraints\t{}\t{}\t{}+{}<=1'.format(forward, reverse, problem.prob.var('zi_{}'.format(forward)), problem.prob.var('zi_{}'.format(reverse))))
    # 	# print('Split reaction zi raw constraint {} :'.format(forward), (problem.prob.var('zi_{}'.format(forward)) + problem.prob.var('zi_{}'.format(reverse)) <= 1))
    # return problem


def solve_objective(problem, objective, all_rxns=True):
    problem.prob.set_objective(problem.get_flux_var(objective))
    problem.prob.set_objective_sense(lp.ObjectiveSense.Maximize)
    problem.prob.solve()
    result = problem.prob.result
    # if all_rxns is True:
    # 	for rxn in sorted(problem._model.reactions):
    # 		try:
    # 			print('{}\t{}\t{}\t{}'.format(rxn, problem.get_flux(rxn), problem._model.get_reaction(rxn), problem.prob.result.get_value('dgri_{}'.format(rxn))))
    # 		except:
    # 			print('{}\t{}\t{}\t{}'.format(rxn, problem.get_flux(rxn), problem._model.get_reaction(rxn), 'No Err'))
    #
    # print('Biomass Objective: {}'.format(result.get_value(problem.get_flux_var(objective))))
    return result.get_value(problem.get_flux_var(objective))

def solve_objective_tmfa(problem, objective, metabolites):
    for metab in metabolites:
        try:
            problem.prob.set_objective(problem.prob.var(metab))
            problem.prob.set_objective_sense(lp.ObjectiveSense.Maximize)
            problem.prob.solve()
            result = problem.prob.result
            max = result.get_value(problem.prob.var(metab))
            problem.prob.set_objective_sense(lp.ObjectiveSense.Minimize)
            problem.prob.solve()
            result = problem.prob.result
            min = result.get_value(problem.prob.var(str(metab)))
            print('{}\t{}\t{}'.format(metab, min, max))
        except:
            print('{}\t{}\t{}'.format(metab, 'NA', 'NA'))

# if all_rxns is True:
    # 	for rxn in sorted(problem._model.reactions):
    # 		try:
    # 			print('{}\t{}\t{}\t{}'.format(rxn, problem.get_flux(rxn), problem._model.get_reaction(rxn), problem.prob.result.get_value('dgri_{}'.format(rxn))))
    # 		except:
    # 			print('{}\t{}\t{}\t{}'.format(rxn, problem.get_flux(rxn), problem._model.get_reaction(rxn), 'No Err'))
    #
    # print('Biomass Objective: {}'.format(result.get_value(problem.get_flux_var(objective))))

class TimeLimitExpired(Exception): pass

def timelimit(timeout, func, args=(), kwargs={}):
    import threading
    class FuncThread(threading.Thread):
        def __init__(self):
            threading.Thread.__init__(self)
            self.result = None

        def run(self):
            self.result = func(*args, **kwargs)

        def _stop(self):
            if self.isAlive():
                threading.Thread._Thread__stop(self)

    it = FuncThread()
    it.start()
    it.join(timeout)
    if it.isAlive():
        it._stop()
        raise(TimeLimitExpired())
    else:
        return it.result
