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
# Copyright 2017-2020  Keith Dufault-Thompson <keitht547@uri.edu>

from __future__ import unicode_literals
import logging
from ..command import SolverCommandMixin, MetabolicMixin, \
    Command, ObjectiveMixin, convert_to_unicode
from ..reaction import Compound
from ..lpsolver import lp
from ..datasource.reaction import parse_reaction
from decimal import Decimal
import csv
import math
import random
import yaml
import argparse
from psamm.datasource.entry import (DictReactionEntry as ReactionEntry)
from six import iteritems
logger = logging.getLogger(__name__)


class TMFACommand(MetabolicMixin, SolverCommandMixin, ObjectiveMixin, Command):
    """Run TMFA on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument('--config', type=argparse.FileType('r'),
                            help='Config file for TMFA settings')
        parser.add_argument('--threshold', default=None, type=Decimal,
                            help='value to fix biomass flux to during '
                                 'tmfa simulations (default = max biomass)')
        parser.add_argument('--temp', help='Temperature in Celsius',
                            default=25)
        parser.add_argument('--hamilton', action='store_true',
                            help='run model using Hamilton TMFA method')
        parser.add_argument('--err', action='store_true',
                            help='use error estimates when running TMFA')
        parser.add_argument('--phin', default='4,11',
                            help='Allowable range for internal compartment pH '
                            'format: Lower,Upper')
        parser.add_argument('--phout', default='4,11',
                            help='Allowable range for external compartment pH '
                            'format: Lower,Upper')

        subparsers = parser.add_subparsers(title='Functions', dest='which')
        parser_util = subparsers.add_parser('util', help='Utility functions')

        parser_util.add_argument('--generate-config', action='store_true')
        parser_util.add_argument('--random-addition', action='store_true',
                                 help='perform random reaction constraint '
                                      'addition in model')

        parser_sim = subparsers.add_parser('simulation',
                                           help='Simulation Functions')
        parser_sim.add_argument('--single-solution', type=str, default=None,
                                choices=['fba', 'l1min', 'random'],
                                help='get a single TMFA solution')

        super(TMFACommand, cls).init_parser(parser)

    def run(self):
        """Run TMFA command."""
        solver = self._get_solver()
        if solver._properties['name'] not in ['cplex', 'gurobi']:
            raise ImportError(
                'Unable to find an LP solver that satisfy the requirement of '
                'TMFA function (TMFA requires Cplex or Gurobi as LP solver)')

        mm = self._mm
        which_command = self._args.which

        if which_command == 'util':
            if self._args.generate_config:
                with open('./example-config.yaml', mode='w') as f:
                    f.write('deltaG: <path to deltaG input file>\n')
                    f.write('exclude: <path to excluded reactions file>\n')
                    f.write('transporters: <path to transporter parameters '
                            'file>\n')
                    f.write('concentrations: <path to set concentrations '
                            'file>\n')
                    f.write('proton-in: <id of the internal proton. '
                            'example cpd_h[c]>\n')
                    f.write('proton-out: <id of the external proton. '
                            'example cpd_h[p]>\n')
                    f.write('proton-other:\n')
                    f.write('  - <id of other proton compounds in other '
                            'compartments>\n')
                    f.write('water:\n')
                    f.write('  - <id of water compounds example cpd_h2o[c]. '
                            'list one per line\n')
                logger.info('Generated example config file '
                            'example-config.yaml')
                quit()

        if self._args.config is None:
            logger.warning('No configuration file was provided with'
                           ' --config option')
            quit()
        else:
            config_dict = yaml.safe_load(self._args.config)

        # Parse exclude file provided through command line.
        base_exclude_list = []
        if config_dict.get('exclude') is not None:
            for line in open(config_dict.get('exclude')).readlines():
                base_exclude_list.append(line.rstrip())

        # Parse set of reactions that will have their deltaG
        # constrained only by deltaG of the transport component.
        ph_difference_rxn = []
        if config_dict.get('rxn_conc_only') is not None:
            for line in open(config_dict.get('rxn_conc_only')).readlines():
                rxn = line.rstrip()
                ph_difference_rxn.append(rxn)
                ph_difference_rxn.append(u'{}_forward'.format(rxn))
                ph_difference_rxn.append(u'{}_reverse'.format(rxn))

        # Parse out the file of lumped reactions or set the lumps
        # dictionaries to be empty
        if config_dict.get('lumped_reactions') is not None:
            lumpid_to_rxn, rxn_lumpid, lump_to_rxnids, \
                lump_to_rxnids_dir = lump_parser(
                    open(config_dict.get('lumped_reactions')))
        else:
            logger.info('No lumped reactions were provided')
            lumpid_to_rxn = {}
            lump_to_rxnids = {}
            lump_to_rxnids_dir = {}

        # Add exchange reactions to the exclude list
        for reaction in mm.reactions:
            if mm.is_exchange(reaction):
                base_exclude_list.append(reaction)

        # these lines will add ph_difference list to the base exclude list.
        base_exclude_list = base_exclude_list + ph_difference_rxn

        # exclude lump_list is a list of all excluded reactions and
        # lump reactions. exclude_lump_unkown list is a list of all
        # excluded reactions, lump reactions and lumped reactions.
        exclude_lump_list = set(base_exclude_list)
        exclude_unknown_list = set(base_exclude_list)
        for lump_id, sub_rxn_list in iteritems(lump_to_rxnids):
            exclude_lump_list.add(lump_id)
            for sub_rxn in sub_rxn_list:
                exclude_unknown_list.add(sub_rxn)
        # Make list of all excluded, lumped, and unknown dgr.
        exclude_lump_unkown = exclude_unknown_list.union(exclude_lump_list)

        cpd_conc_dict = {}
        # Parse file containing set concentrations and ranges.
        if config_dict.get('concentrations') is not None:
            for row in csv.reader(
                open(config_dict.get('concentrations')),
                    delimiter=str('\t')):
                cpd, lower, upper = row
                cpd_conc_dict[convert_to_unicode(cpd)] = [lower, upper]

        # Add lump reactions to the mm_irreversible model
        for lump_id, rxn in iteritems(lumpid_to_rxn):
            reaction = parse_reaction(rxn)
            mm.database.set_reaction(lump_id, reaction)
            mm.add_reaction(lump_id)
            mm.limits[lump_id].lower = 0
            mm.limits[lump_id].upper = 0

        # make an irreversible version of the metabolic model:
        mm_irreversible, _, split_reversible, \
            reversible_lump_to_rxn_dict = mm.make_irreversible(
                gene_dict={}, exclude_list=exclude_lump_unkown,
                lumped_rxns=lump_to_rxnids_dir,
                all_reversible=self._args.hamilton)

        split_reactions = []
        for_rev_reactions = []
        for (i, j) in split_reversible:
            for_rev_reactions.append(i)
            for_rev_reactions.append(j)
            split_reactions.append(i[:-8])

        model = self._model

        for rx in [i for i in model.reactions]:
            if rx.id in split_reactions:
                f = ReactionEntry(
                    dict(id=u'{}_forward'.format(rx.id), genes=rx.genes))
                r = ReactionEntry(
                    dict(id=u'{}_reverse'.format(rx.id), genes=rx.genes))
                model.reactions.add_entry(f)
                model.reactions.add_entry(r)

        # update these lists for the new reversible lumps and constituent
        # reactions. exclude lump_list is a list of all excluded reactions
        # and lump reactions. exclude_lump_unkown list is a list of
        # all excluded reactions, lump reactions and lumped reactions.
        for lump_id, sub_rxn_list in iteritems(reversible_lump_to_rxn_dict):
            exclude_lump_list.add(lump_id)
            for sub_rxn in sub_rxn_list:
                exclude_unknown_list.add(sub_rxn)
        # Make list of all excluded, lumped, and unknown dgr.
        exclude_lump_unkown = exclude_unknown_list.union(exclude_lump_list)

        # Read DeltaG information from file.
        dgr_dict = {}
        if config_dict.get('deltaG') is not None:
            dgr_dict = parse_dgr_file(open(
                config_dict.get('deltaG')), mm_irreversible)
            print('using dgr file')
        elif config_dict.get('deltaGf') is not None:
            # Parse the deltaGf values for all metabolites
            # from a supplied file.
            dgf_dict = parse_dgf(mm_irreversible, open(
                config_dict.get('deltaGf')))
            # Calculate the deltaGr values for all reactions
            dgr_dict = calculate_dgr(
                mm_irreversible, dgf_dict, exclude_unknown_list,
                open(config_dict.get('transporters')),
                ph_difference_rxn, open(config_dict.get('scaled_compounds')))
            print('using dgf file')

        if config_dict.get('transporters') is not None:
            transport_parameters = parse_tparam_file(
                open(config_dict.get('transporters')))

        prob, v, zi, dgri, xij, cp_list = make_tmfa_problem(
            mm_irreversible, solver)

        # Parse pH settings from command line
        phin_s = self._args.phin
        phout_s = self._args.phout
        phinsplit = phin_s.split(',')
        phoutsplit = phout_s.split(',')
        phin = (float(phinsplit[0]), float(phinsplit[1]))
        phout = (float(phoutsplit[0]), float(phoutsplit[1]))
        ph_bounds = [phin, phout]

        # Run Utility based functions
        if which_command == 'util':
            if self._args.random_addition:
                testing_list = list(mm_irreversible.reactions)
                random.shuffle(testing_list)
                testing_list_tmp = []
                for rx in testing_list:
                    testing_list_iter = testing_list_tmp + [rx]
                    logger.info(u'current testing list: {}'.format(
                        testing_list_tmp))
                    logger.info(u'testing reaction: {}'.format(rx))

                    prob, v, zi, dgri, xij, cp_list = make_tmfa_problem(
                        mm_irreversible, solver)

                    prob, cpd_xij_dict = add_conc_constraints(
                        xij, prob, cpd_conc_dict, cp_list,
                        config_dict.get('water'), config_dict.get('proton-in'),
                        config_dict.get('proton-out'),
                        config_dict.get('proton-other'))

                    prob, excluded_compounds = add_reaction_constraints(
                        prob, v, zi, dgri, xij, mm_irreversible,
                        exclude_lump_list, exclude_unknown_list,
                        exclude_lump_unkown, dgr_dict,
                        reversible_lump_to_rxn_dict, split_reversible,
                        transport_parameters, testing_list_iter,
                        config_dict.get('scaled_compounds'),
                        config_dict.get('water'), config_dict.get('proton-in'),
                        config_dict.get('proton-out'),
                        config_dict.get('proton-other'), ph_bounds,
                        self._args.temp, self._args.err)
                    try:
                        biomass = get_var_bound(prob, v(self._get_objective()),
                                                lp.ObjectiveSense.Maximize)
                        logger.info(u'Current Biomass: {}'.format(biomass))
                        if biomass < self._args.threshold:
                            continue
                        else:
                            testing_list_tmp.append(rx)
                    except:
                        continue

                for rx in testing_list:
                    if rx not in testing_list_tmp:
                        print(u'{}\tBad Constraint'.format(rx))
                    else:
                        print(u'{}\tGood Constraint'.format(rx))
                quit()
        # Run simulation based functions
        elif which_command == 'simulation':
            prob, cpd_xij_dict = add_conc_constraints(
                xij, prob, cpd_conc_dict, cp_list,
                config_dict.get('water'), config_dict.get('proton-in'),
                config_dict.get('proton-out'), config_dict.get('proton-other'))
            testing_list_tmp = list(mm_irreversible.reactions)
            prob, excluded_compounds = add_reaction_constraints(
                prob, v, zi, dgri, xij, mm_irreversible, exclude_lump_list,
                exclude_unknown_list, exclude_lump_unkown, dgr_dict,
                reversible_lump_to_rxn_dict, split_reversible,
                transport_parameters, testing_list_tmp,
                config_dict.get('scaled_compounds'), config_dict.get('water'),
                config_dict.get('proton-in'), config_dict.get('proton-out'),
                config_dict.get('proton-other'), ph_bounds,
                self._args.temp, self._args.err)

            if self._args.threshold is not None:
                prob.add_linear_constraints(
                    v(self._get_objective()) == self._args.threshold)
                logger.info(u'Set biomass based on threshold to '
                            '{}'.format(self._args.threshold))
            else:
                max_biomass = get_var_bound(prob, v(self._get_objective()),
                                            lp.ObjectiveSense.Maximize)
                prob.add_linear_constraints(
                    v(self._get_objective()) == max_biomass)
                logger.info(u'Set biomass based on maxbiomass to {}'.format(
                    max_biomass))
            split_reversible_list = []
            for (f, r) in split_reversible:
                split_reversible_list.append(f)
                split_reversible_list.append(r)
            if self._args.single_solution is not None:
                print_fba(self._args.single_solution, prob,
                          self._get_objective(), v, zi, dgri, xij,
                          mm_irreversible, cp_list, excluded_compounds,
                          split_reversible_list)
            else:
                print_fva(prob, mm_irreversible, cp_list, exclude_unknown_list,
                          v, dgri, zi, xij, excluded_compounds,
                          split_reversible_list)
        quit()


def make_tmfa_problem(mm_irreversible, solver):
    """Make a flux balance problem with TMFA related variables

    The TMFA problem contains 4 classes of variables. Flux variables
    defined in the v namespace, binary reaction indicator variables
    defined in the zi namespace, Gibbs free energy variables defined
    in the dgri namespace, and concentration variables defined in the
    xij namespace. This function will add these variables with their
    default lower and upper bounds. The function will then return
    the flux problem, variable namespaces and a list of the compounds.

    Args:
        mm_irreversible: :class:`psamm.metabolicmodel.MetabolicModel`.
        solver: linear programming library to use.
    """
    prob = solver.create_problem()
    if solver._properties['name'] == 'gurobi':
        prob.integrality_tolerance.value = 1e-9
        logger.warning(
            'Gurobi supports minimum integrality tolerance of 1e-9. This may '
            'affect the results from this simulation')
    elif solver._properties['name'] == 'glpk':
        prob.integrality_tolerance.value = 1e-21
    else:
        prob.integrality_tolerance.value = 0
    v = prob.namespace(name='flux')
    zi = prob.namespace(name='zi')
    dgri = prob.namespace(name='dgri')
    xij = prob.namespace(name='xij')
    for reaction in mm_irreversible.reactions:
        lower, upper = mm_irreversible.limits[reaction]
        v.define([reaction], lower=lower, upper=upper,
                 types=lp.VariableType.Continuous)
        zi.define([reaction], lower=int(0), upper=int(1),
                  types=lp.VariableType.Binary)
        dgri.define([reaction], lower=-1000, upper=1000,
                    types=lp.VariableType.Continuous)

    massbalance_lhs = {compound: 0 for compound in mm_irreversible.compounds}
    for spec, value in iteritems(mm_irreversible.matrix):
        compound, reaction_id = spec
        massbalance_lhs[compound] += v(reaction_id) * value
    for compound, lhs in iteritems(massbalance_lhs):
        prob.add_linear_constraints(lhs == 0)

    cp_list = []
    for cp in mm_irreversible.compounds:
        cp_list.append(convert_to_unicode(str(cp)))
        xij.define([convert_to_unicode(str(cp))], lower=-50, upper=50,
                   types=lp.VariableType.Continuous)
    return prob, v, zi, dgri, xij, cp_list


def get_var_bound(prob, var, objective_sense):
    """Gets upper or lower bound of a variable in an LP problem.

    Args:
        prob: :class:`psamm.lpsolver.lp.Problem`.
        var: LP problem variable
        objective_sense: :class:`psamm.lpsolver.lp.ObjectiveSense`
    """
    prob.set_objective(var)
    result = prob.solve_unchecked(objective_sense)
    if not result.success:
        logger.error(u'Solution not optimal: {}'.format(result.status))
        quit()
    return result.get_value(var)


def print_fva(prob, mm_irreversible, cp_list, exclude_unknown_list,
              _v, _dgri, _zi, _xij, excluded_compounds, split_reversible_list):
    """Prints FVA like result from TMFA problem.

    This function will take a TMFA problem along with the associated
    variable namespaces and print out all of the upper and lower bounds
    of the variables.

    Args:
        prob: :class:`psamm.lpsolver.lp.Problem`.
        mm_irreversible: :class:`psamm.metabolicmodel.MetabolicModel`.
        cp_list: List of compounds in the metabolic model.
        exclude_unknown_list: List of reactions excluded from thermodynamic
            constraints.
        _v: variable namespace for flux variables.
        _dgri: variable namespace for gibbs free energy variables
        _zi: variables namespace for indicator variables
        _xij: variable namespace for concentrations
        excluded_compounds: list of compounds that are not constrained
        split_reversible_list: list of all split reactions

    """
    for step, reaction in enumerate(sorted(mm_irreversible.reactions)):
        # logger.info('Testing Reaction {}/{}'.format(step, len_rxns))
        min_flux = get_var_bound(prob, _v(reaction),
                                 lp.ObjectiveSense.Minimize)
        max_flux = get_var_bound(prob, _v(reaction),
                                 lp.ObjectiveSense.Maximize)
        print(u'Flux\t{}\t{}\t{}'.format(reaction, min_flux, max_flux))
        if reaction not in exclude_unknown_list:
            if reaction in split_reversible_list:
                if '_reverse' in reaction:
                    rxn_f = reaction
                    rxn_f = rxn_f.replace('_reverse', '')
                    rxn_f += '_forward'
                    max_dgr = -1*get_var_bound(prob, _dgri(rxn_f),
                                               lp.ObjectiveSense.Minimize)
                    min_dgr = -1*get_var_bound(prob, _dgri(rxn_f),
                                               lp.ObjectiveSense.Maximize)
                else:
                    min_dgr = get_var_bound(prob, _dgri(reaction),
                                            lp.ObjectiveSense.Minimize)
                    max_dgr = get_var_bound(prob, _dgri(reaction),
                                            lp.ObjectiveSense.Maximize)
            else:
                min_dgr = get_var_bound(prob, _dgri(reaction),
                                        lp.ObjectiveSense.Minimize)
                max_dgr = get_var_bound(prob, _dgri(reaction),
                                        lp.ObjectiveSense.Maximize)
            min_zi = get_var_bound(prob, _zi(reaction),
                                   lp.ObjectiveSense.Minimize)
            max_zi = get_var_bound(prob, _zi(reaction),
                                   lp.ObjectiveSense.Maximize)
            print(u'DGR\t{}\t{}\t{}'.format(reaction, min_dgr, max_dgr))
            print(u'Zi\t{}\t{}\t{}'.format(reaction, min_zi, max_zi))
        else:
            print(u'DGR\t{}\t{}\t{}'.format(reaction, 'NA', 'NA'))
            print(u'Zi\t{}\t{}\t{}'.format(reaction, 'NA', 'NA'))
    for step, cpd in enumerate(sorted(cp_list)):
        if cpd not in excluded_compounds:
            # logger.info('Testing Compound {}\{}'.format(step, len(cp_list)))
            min_cpd = get_var_bound(prob, _xij(cpd),
                                    lp.ObjectiveSense.Minimize)
            max_cpd = get_var_bound(prob, _xij(cpd),
                                    lp.ObjectiveSense.Maximize)
            print(u'CONC\t{}\t{}\t{}'.format(
                cpd, math.exp(min_cpd), math.exp(max_cpd)))


def print_fba(simulation, prob, objective, _v, _zi, _dgri, _xij, mm,
              cp_list, excluded_compounds, split_reversible_list):
    """Prints FBA like result from TMFA problem.

    This function will take a TMFA problem along with the associated
    variable namespaces and print out a single solution for the
    associated problem. The solution can be either a single FBA-like
    solution, a L1min solution, or a random solution from the
    boundary of the solution space.

    Args:
        simulation: type of solving for the problem. ['fba', 'l1min', 'random']
        prob: :class:`psamm.lpsolver.lp.Problem`.
        objective: Objective reaction ID.
        mm_irreversible: :class:`psamm.metabolicmodel.MetabolicModel`.
        cp_list: List of compounds in the metabolic model.
        exclude_unknown_list: List of reactions excluded
            from thermo constraints.
        _v: variable namespace for flux variables.
        _dgri: variable namespace for gibbs free energy variables
        _zi: variables namespace for indicator variables
        _xij: variable namespace for concentrations
        excluded_compounds: list of compounds that are not constrained
        split_reversible_list: list of all split reactions
    """
    prob.set_objective(_v(objective))
    result = prob.solve_unchecked(lp.ObjectiveSense.Maximize)
    if not result.success:
        logger.error(u'Solution not optimal: {}'.format(result.status))
        quit()
    if simulation == 'fba':
        for rxn in mm.reactions:
            print(u'Flux\t{}\t{}'.format(rxn, result.get_value(_v(rxn))))
            if rxn in split_reversible_list:
                if '_reverse' in rxn:
                    rxn_f = rxn
                    rxn_f = rxn_f.replace('_reverse', '')
                    rxn_f += '_forward'
                    print(u'DGR\t{}\t{}'.format(
                        rxn, -1*result.get_value(_dgri(rxn_f))))
                else:
                    print(u'DGR\t{}\t{}'.format(
                        rxn, result.get_value(_dgri(rxn))))
            else:
                print(u'DGR\t{}\t{}'.format(rxn, result.get_value(_dgri(rxn))))
            print(u'Zi\t{}\t{}'.format(rxn, result.get_value(_zi(rxn))))
        for cp in cp_list:
            if cp not in excluded_compounds:
                print(u'CPD\t{}\t{}'.format(cp, result.get_value(_xij(cp))))
    elif simulation == 'l1min':
        _z = prob.namespace(mm.reactions, lower=0)
        z = _z.set(mm.reactions)
        v = _v.set(mm.reactions)
        prob.add_linear_constraints(z >= v, v >= -z)

        objective = _z.expr(
            (reaction_id, -1)
            for reaction_id in mm.reactions)
        prob.set_objective(objective)
        result = prob.solve_unchecked(lp.ObjectiveSense.Maximize)
        if not result.success:
            logger.error(u'Solution not optimal: {}'.format(result.status))
            quit()
        for rxn in mm.reactions:
            print(u'Flux\t{}\t{}'.format(rxn, result.get_value(_v(rxn))))
            if rxn in split_reversible_list:
                if '_reverse' in rxn:
                    rxn_f = rxn
                    rxn_f = rxn_f.replace('_reverse', '')
                    rxn_f += '_forward'
                    print(u'DGR\t{}\t{}'.format(
                        rxn, -1*result.get_value(_dgri(rxn_f))))
                else:
                    print(u'DGR\t{}\t{}'.format(
                        rxn, result.get_value(_dgri(rxn))))
            else:
                print(u'DGR\t{}\t{}'.format(
                    rxn, result.get_value(_dgri(rxn))))
            print(u'Zi\t{}\t{}'.format(rxn, result.get_value(_zi(rxn))))
        for cp in cp_list:
            if cp not in excluded_compounds:
                print(u'CPD\t{}\t{}'.format(cp, result.get_value(_xij(cp))))
    elif simulation == 'random':
        optimize = None
        for rxn_id in mm.reactions:
            if optimize is None:
                optimize = _v(rxn_id) * random.random()
            else:
                optimize += _v(rxn_id) * random.random()
        prob.add_linear_constraints(
            _v(objective) == result.get_value(_v(objective)))

        prob.set_objective(optimize)

        result = prob.solve_unchecked()
        for rxn in mm.reactions:
            print(u'Flux\t{}\t{}'.format(rxn, result.get_value(_v(rxn))))
            if rxn in split_reversible_list:
                if '_reverse' in rxn:
                    rxn_f = rxn
                    rxn_f = rxn_f.replace('_reverse', '')
                    rxn_f += '_forward'
                    print(u'DGR\t{}\t{}'.format(
                        rxn, -1*result.get_value(_dgri(rxn_f))))
                else:
                    print(u'DGR\t{}\t{}'.format(
                        rxn, result.get_value(_dgri(rxn))))
            else:
                print(u'DGR\t{}\t{}'.format(rxn, result.get_value(_dgri(rxn))))
            print(u'Zi\t{}\t{}'.format(rxn, result.get_value(_zi(rxn))))
        for cp in cp_list:
            if cp not in excluded_compounds:
                print(u'CPD\t{}\t{}'.format(cp, result.get_value(_xij(cp))))


def lump_parser(lump_file):
    """Parses a supplied file and returns
    dictionaries containing lump information.

    The supplied file should be in a tab separated table in the format of
    lumpID  lump_deltaG rxn1:1,rxn2:-1,rxn3:1  lumpRXN
    Returns dictionaries storing this information, linking reactions to the
    lumps.
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
            subrxn, dir = convert_to_unicode(i).split(':')
            rx_dir_list.append((subrxn, dir))
            rx_list.append(subrxn)
            rxn_to_lump_id[i] = lump_id
        lump_to_rxn[lump_id] = lump_rxn
        lump_to_rxnids[lump_id] = rx_list
        lump_to_rxnids_dir[lump_id] = rx_dir_list
    return lump_to_rxn, rxn_to_lump_id, lump_to_rxnids, lump_to_rxnids_dir


def parse_tparam_file(file):
    """Parse a transport parameter file.

    This file contains reaction IDs, net charge transported into
    the cell, and net protons transported into the cell.

    """
    t_param = {}
    if file is not None:
        for row in csv.reader(file, delimiter=str('\t')):
            rxn, c, h = row
            rxn = convert_to_unicode(rxn)
            t_param[rxn] = (Decimal(c), Decimal(h))
            t_param[u'{}_forward'.format(rxn)] = (Decimal(c), Decimal(h))
            t_param[u'{}_reverse'.format(rxn)] = (-Decimal(c), -Decimal(h))
    return t_param


def parse_dgf(mm, dgf_file):
    """A function that will parse a supplied deltaG of formation file.

    Compound IDs in this file do not need to contain the compartments.
    compound deltaGf values should be in Kcal/mol

    Args:
        mm: a metabolic model object
        dgf_file: a file that containing 2 columns of
            compound ids and deltaGf values
    """

    cpd_dgf_dict = {}
    for row in csv.reader(dgf_file, delimiter=str('\t')):
        for cpt in mm.compartments:
            try:
                dg = Decimal(row[1])
                derr = Decimal(row[2])
                cpd_dgf_dict[Compound(row[0], cpt)] = (dg, derr)
            except:
                logger.info(u'Compound {} has an assigned detltGf value '
                            'of {}. This is not an number and will be treated '
                            'as a missing '
                            'value.'.format(Compound(row[0], cpt), row[1]))
    return cpd_dgf_dict


def add_conc_constraints(xij, problem, cpd_conc_dict,
                         cp_list, water, hin, hout, hother):
    """Add concentration constraints to TMFA problem
    based on parsed concentration dictionary.

    """
    # Water and hydrogen need to be excluded from
    # these concentration constraints.
    excluded_cpd_list = []
    excluded_cpd_list.append(hin)
    excluded_cpd_list.append(hout)
    excluded_cpd_list.append(hother)
    for wat in water:
        excluded_cpd_list.append(wat)
    cpdid_xij_dict = {}
    for cp in cp_list:
        # define concentration variable for compound.
        cpdid_xij_dict[cp] = xij(cp)
        var = xij(cp)
        # Define default constraints for anything not set in the conc file
        if cp not in cpd_conc_dict.keys():
            if cp not in excluded_cpd_list:
                # Add concentration constraints as the
                # ln of the concentration (M).
                problem.add_linear_constraints(var >= math.log(0.00001))
                problem.add_linear_constraints(var <= math.log(0.02))
        elif cp in cpd_conc_dict.keys():
            if cp not in excluded_cpd_list:
                conc_limits = cpd_conc_dict[cp]
                if conc_limits[0] > conc_limits[1]:
                    logger.error(u'lower bound for {} concentration higher '
                                 'than upper bound'.format(cp))
                    quit()
                if Decimal(conc_limits[0]) == Decimal(conc_limits[1]):
                    problem.add_linear_constraints(
                        var == math.log(Decimal(conc_limits[0])))
                else:
                    problem.add_linear_constraints(
                        var >= math.log(Decimal(conc_limits[0])))
                    problem.add_linear_constraints(
                        var <= math.log(Decimal(conc_limits[1])))
    return problem, cpdid_xij_dict


def parse_dgr_file(dgr_file, mm):
    """Parses DeltaG of reaction file and returns a dictionary of the values.

    """
    def is_number(val):
        try:
            float(val)
            return True
        except ValueError:
            return False
    dgr_dict = {}
    for row in csv.reader(dgr_file, delimiter=str('\t')):
        rxn, dgr, err = row
        rxn = convert_to_unicode(rxn)
        if is_number(dgr):
            if is_number(err):
                err = Decimal(err)
            else:
                err = Decimal(2)
            if rxn in mm.reactions:
                dgr_dict[rxn] = (Decimal(dgr), err)
            elif '{}_forward'.format(rxn) in mm.reactions:
                dgr_dict[u'{}_forward'.format(rxn)] = (Decimal(dgr), err)
                dgr_dict[u'{}_reverse'.format(rxn)] = (-Decimal(dgr), err)
        else:
            logger.info(
                u'Reaction {} was provided with dgr value of {}'.format(
                    rxn, dgr))
    return dgr_dict


def calculate_dgr(mm, dgf_dict, excluded_reactions,
                  transport_parameters, ph_difference_rxn, scaled_compounds):
    """Calculates DeltaG values from DeltaG of formation values of compounds.

    """
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
                for j in rxn.compounds:
                    if j[0] not in dgf_dict.keys():
                        print(j[0])
                if reaction not in ph_difference_rxn:
                    logger.error(u'Reaction {} contains at least 1 compound '
                                 'with an unknown deltaGf '
                                 'value'.format(reaction))
                    quit()
            else:
                dgr = 0
                # Make a variable dgf_sum that represents
                # the sum of sij *  (stoichiometry *
                # deltaGf for reaction j.
                for cpd in rxn.compounds:
                    cpd_0 = convert_to_unicode(str(cpd[0]))
                    if cpd_0 in dgf_scaling.keys():
                        dgscale = dgf_scaling[cpd_0]
                    else:
                        dgscale = 1
                    (dg, dge) = dgf_dict[cpd[0].name]
                    dgs = Decimal(dg) * (Decimal(cpd[1])*dgscale)
                    dgr += dgs
                    dgerr += Decimal(cpd[1])*Decimal(dgscale) * dge
            dgr_dict[reaction] = (dgr, dgerr)
    return dgr_dict


def add_reaction_constraints(
        problem, _v, _zi, _dgri, _xij, mm, exclude_lumps, exclude_unknown,
        exclude_lumps_unknown, dgr_dict, lump_rxn_list, split_rxns,
        transport_parameters, testing_list, scaled_compounds, water, hin,
        hout, hother, ph, temp, err_est=False, hamilton=False):
    """Adds reaction constraints to a TMFA problem

    This function will add in gibbs free energy constraints to a TMFA flux
    problem. These constraints will be added based on provided Gibbs free
    energy values and temperature.

    Args:
        problem: :class:`psamm.lpsolver.lp.Problem`.
        mm: :class:`psamm.metabolicmodel.MetabolicModel`.
        _v: variable namespace for flux variables.
        _dgri: variable namespace for gibbs free energy variables
        _zi: variables namespace for indicator variables
        _xij: variable namespace for concentrations
        exclude_unknown: List of reactions excluded from thermo constraints.
        exclude_lumps_unknown: List of excluded reactions and lumped reactions.
        dgr_dict: dictionary where keys are reaction ids and values
            are deltag values.
        lump_rxn_list: List of lump reaction IDs.
        split_rxns: List of tuples of (reaction_forward, reaction_reverse)
        transport_parameters: dictionary of reaction IDs to proton
            transport and charge transport values.
        testing_list: List of reactions to add deltaG constraints for.
        water: list of water compound IDs.
        hin: ID of proton compound from inside cell compartment
        hout: ID of proton compound from outside cell compartment
        hother: List of other proton IDs
        ph: list of two tuples containing pH bounds. [(in_lower, in_upper),
            (out_lower, out_upper)]
        temp: Temperature in Celsius
        err_est: True or False for using error estimates for deltaG values.
        hamilton: True or False for using Hamilton TMFA method.

    """
    dgf_scaling = {}
    if scaled_compounds is not None:
        scaled_compounds.seek(0)
        for row in csv.reader(scaled_compounds, delimiter=str('\t')):
            dgf_scaling[row[0]] = Decimal(row[1])

    # "idg" is the ideal gas constant in units of kJ/mol
    idg = Decimal(8.3144621 / 1000)  # kJ/mol
    # R = Decimal(1.987 / 1000) kcal/mol

    # "tkelvin" is the temperature on the Kelvin scale
    tkelvin = Decimal(temp) + Decimal(273.15)

    k = 500
    epsilon = 0.000001
    ph_in = ph[0]
    ph_out = ph[0]
    h_p = _xij(str(hout))

    problem.add_linear_constraints(h_p <= ph_out[1])
    problem.add_linear_constraints(h_p >= ph_out[0])

    h_c = _xij(str(hin))

    problem.add_linear_constraints(h_c >= ph_in[0])
    problem.add_linear_constraints(h_c <= ph_in[1])
    delta_ph = (h_c - h_p)

    # "fc" is the Faraday constant in units of kcal/(mV*mol)
    fc = Decimal(0.02306)

    excluded_cpd_list = []
    excluded_cpd_list.append(hin)
    excluded_cpd_list.append(hout)
    if hother is not None:
        excluded_cpd_list.append(hother)
    for wat in water:
        excluded_cpd_list.append(wat)
    logger.info(u'Excluded compounds: {}'.format(','.join(excluded_cpd_list)))
    logger.info(u'Temperature: {}'.format(tkelvin))
    logger.info(u'using h in {}'.format(hin))
    logger.info(u'using h out {}'.format(hout))
    logger.info(u'using water {}'.format(water))
    logger.info(u'using ph range of {} to {} for internal compartment'.format(
        ph_in[0], ph_in[1]))
    logger.info(u'using ph range of {} to {} for external compartment'.format(
        ph_out[0], ph_out[1]))
    # for (f, r) in split_rxns:
    #     split_list.append(f)
    #     split_list.append(r)
    #     if f not in exclude_unknown:
    #         dgrif = _dgri(f)
    #         dgrir = _dgri(r)
    #         problem.add_linear_constraints(dgrif == dgrir * -1)

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

            zi = _zi(reaction)
            split_rxns_l = []
            for (f, r) in split_rxns:
                split_rxns_l.append(f)
                split_rxns_l.append(r)
            if reaction in split_rxns_l:
                if '_reverse' in reaction:
                    f_rxn = reaction.replace('_reverse', '')
                    f_rxn = f_rxn+'_forward'
                    dgri = -1*_dgri(reaction)
                else:
                    dgri = _dgri(reaction)
            else:
                dgri = _dgri(reaction)
            vi = _v(reaction)
            vmax = mm.limits[reaction].upper
            if reaction in testing_list:
                if reaction in transport_parameters.keys():
                    (c, h) = transport_parameters[reaction]
                    ddph = Decimal(-2.3)*Decimal(h)*idg*tkelvin*delta_ph
                    dpsi = Decimal(33.33) * delta_ph - Decimal(143.33)
                    ddpsi = dpsi * Decimal(c) * Decimal(fc)
                    dgr_trans = ddph + ddpsi
                else:
                    dgr_trans = 0
                (dgr0, err) = dgr_dict[reaction]
                ssxi = 0

                if err_est:
                    problem.define(u'dgr_err_{}'.format(reaction),
                                   types=lp.VariableType.Continuous,
                                   lower=-1000, upper=1000)
                    dgr_err = problem.var(u'dgr_err_{}'.format(reaction))
                    problem.add_linear_constraints(dgr_err <= 2*err)
                    problem.add_linear_constraints(dgr_err >= -2*err)
                else:
                    dgr_err = 0

                for (cpd, stoich) in rxn.compounds:
                    cpd = convert_to_unicode(str(cpd))
                    if cpd not in excluded_cpd_list:
                        scale = dgf_scaling.get(cpd, 1)
                        ssxi += _xij(cpd) * Decimal(stoich) * scale

                problem.add_linear_constraints(
                    dgri == dgr0 + (idg * tkelvin * ssxi
                                    ) + dgr_err + dgr_trans)
                if hamilton:
                    problem.add_linear_constraints(dgri <= 300-epsilon)
                    problem.add_linear_constraints(dgri >= -300+epsilon)

                if reaction not in exclude_lumps_unknown:
                    if rhs_check != 0 and lhs_check != 0:
                        problem.add_linear_constraints(
                            dgri - k + (k * zi) <= -epsilon)
                        problem.add_linear_constraints(vi <= zi * vmax)

        if reaction in lump_rxn_list.keys():
            problem.define(u'yi_{}'.format(reaction), lower=int(0),
                           upper=int(1), types=lp.VariableType.Binary)
            if reaction not in new_excluded_reactions:
                vi = _v(reaction)
                yi = problem.var(u'yi_{}'.format(reaction))
                dgri = _dgri(reaction)
                problem.add_linear_constraints(vi == 0)
                problem.add_linear_constraints(dgri - (k * yi) <= - epsilon)
                sub_rxn_list = lump_rxn_list[reaction]
                sszi = 0
                for sub_rxn in sub_rxn_list:
                    sszi += _zi(sub_rxn)
                problem.add_linear_constraints(
                    yi + sszi <= len(sub_rxn_list))

    for (forward, reverse) in split_rxns:
        problem.add_linear_constraints(
            _zi(forward) + _zi(reverse) <= int(1))
    return problem, excluded_cpd_list
