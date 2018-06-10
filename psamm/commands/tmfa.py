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
# Copyright 2017-2017  Keith Dufault-Thompson <keitht547@uri.edu>

from __future__ import unicode_literals
import logging
from ..command import SolverCommandMixin, MetabolicMixin, Command, ObjectiveMixin
from ..reaction import Direction, Reaction, Compound
from .. import fluxanalysis
from ..lpsolver import lp
from ..datasource.reaction import parse_reaction

import csv
import math
from collections import Counter
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
            type=file)
        parser.add_argument('--trans-param', help='file containing transport c and h values.', type=file)
        parser.add_argument(
            '--set-concentrations', help='Tab seperated file with Reaction ID, [lower], [upper]', type=file)
        super(TMFACommand, cls).init_parser(parser)

    def run(self):
        """Run TMFA command."""
        solver = self._get_solver()
        objective = self._get_objective()
        mm = self._mm
        # Parse exclude file provided through command line.
        base_exclude_list = []
        for line in self._args.exclude.readlines():
            base_exclude_list.append(line.rstrip())

        # Parse out the file of lumped reactions or set the lumps dictionaries to be empty
        if self._args.lump_file is not None:
            lumpid_to_rxn, rxn_lumpid, lump_to_rxnids = lump_parser(self._args.lump_file)
        else:
            logger.info('No lumped reactions were provided')
            lumpid_to_rxn = {}
            rxn_lumpid = {}
            lump_to_rxnids = {}

        # Parse lump reaction equations and add them to the metabolic model
        for lump_id, rxn in lumpid_to_rxn.iteritems():
            reaction = parse_reaction(rxn)
            mm.database.set_reaction(lump_id, reaction)
            mm.add_reaction(lump_id)

        # Add any reactions that are in lumps to the exclude list
        base_exclude_lumped_list = list(base_exclude_list)
        for rxn in rxn_lumpid.iterkeys():
            base_exclude_lumped_list.append(rxn)

        base_exclude_lumped_exchange_list = list(base_exclude_lumped_list)
        for rxn in mm.reactions:
            if mm.is_exchange(rxn):
                base_exclude_lumped_exchange_list.append(rxn)

        cpd_dgf_dict = parse_dgf(mm, self._args.dgf_file)

        dg_r0_dict = calculate_dgr(mm, cpd_dgf_dict, base_exclude_lumped_exchange_list)



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
    for row in csv.reader(lump_file, delimiter=str('\t')):
        #print(row)
        lump_id, lump_rxn_list, lump_rxn = row
        for i in lump_rxn_list.split(','):
            rxn_to_lump_id[i] = lump_id
        lump_to_rxn[lump_id] = lump_rxn
        lump_to_rxnids[lump_id] = lump_rxn_list.split(',')
    return lump_to_rxn, rxn_to_lump_id, lump_to_rxnids


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
                dg = float(row[1])
                cpd_dgf_dict[Compound(row[0], cpt)] = dg
            except:
                logger.info('Compound {} has an assigned detltGf value of {}. This is not an number and will be treated\
                            as a missing value.'.format(Compound(row[0], cpt), row[1]))
    return cpd_dgf_dict


def add_conc_constraints(problem, conc_file):
    # Water needs to be excluded from these concentration constraints.
    excluded_compounds = ['h2o[c]', 'h2o[e]']
    cpdid_xij_dict = {}
    cpd_conc_dict = {}
    # Parse file containing set concentrations and ranges.
    if conc_file is not None:
        for row in csv.reader(conc_file, delimiter=str('\t')):
            if len(row) >= 3:
                cpd, lower, upper = row
                cpd_conc_dict[cpd] = [lower, upper]
            else:
                cpd, fixed = row
                cpd_conc_dict[cpd] = [fixed]
    for cp in problem._model.compounds:
        # define concentration variable for compound.
        problem.prob.define(cp)
        var = problem.prob.var(cp)
        cpdid_xij_dict[str(cp)] = var
        # Define default constraints for anything not set in the conc file
        if str(cp) not in cpd_conc_dict.keys():
            if str(cp) not in excluded_compounds:
                # Add concentration constraints as the ln of the concentration (M).
                problem.prob.add_linear_constraints(var >= math.log(0.00001))
                problem.prob.add_linear_constraints(var <= math.log(0.02))
        elif str(cp) in cpd_conc_dict.keys():
            conc_limits = cpd_conc_dict[str(cp)]
            if len(conc_limits) > 1:
                problem.prob.add_linear_constraints(var >= math.log(float(conc_limits[0])))
                problem.prob.add_linear_constraints(var <= math.log(float(conc_limits[1])))
                if float(conc_limits[0] >= conc_limits[1]):
                    print(conc_limits)
            else:
                problem.prob.add_linear_constraints(var == math.log(float(conc_limits[0])))
        else:
            logger.info('Compound {} had no concentration constraints added'.format(str(cp)))

    return problem, cpdid_xij_dict


def calculate_dgr(mm, dgf_dict, excluded_reactions):
    dgr_dict = {}
    for reaction in mm.reactions:
        if reaction not in excluded_reactions:
            rxn = mm.get_reaction(reaction)
            if any(dgf_dict.get(j[0]) is None for j in rxn.compounds):
                logger.error('Reaction {} contains at least 1 compound with an unknown deltaGf value'.format(reaction))
            else:
                dgf_sum = 0
                # Make a variable dgf_sum that represents the sum of sij *  (stoichiometry *
                # deltaGf for reaction j.
                for cpd in rxn.compounds:
                    dg = dgf_dict[cpd[0]]
                    dgs = dg * float(cpd[1])
                    dgf_sum += dgs
            dgr_dict[reaction] = dgf_sum
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