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
		parser.add_argument('--transport-parameters', help='file containing parameters for transport rxn dgr', type=file)
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

		for reaction in mm.reactions:
			if mm.is_exchange(reaction):
				base_exclude_list.append(reaction)

		# exclude lump_list is a list of all excluded reactions and lump reactions.
		# exclude_lump_unkown list is a list of all excluded reactions, lump reactions and lumped reactions.
		exclude_lump_list = list(base_exclude_list)
		exclude_unkown_list = list(base_exclude_list)
		for lump_id, sub_rxn_list in lump_to_rxnids.iteritems():
			exclude_lump_list.append(lump_id)
			for sub_rxn in sub_rxn_list:
				exclude_unkown_list.append(sub_rxn)
		# Make list of all excluded, lumped, and unknown dgr.
		exclude_lump_unkown = exclude_unkown_list + exclude_lump_list

		# make an irreversible version of the metabolic model:
		mm_irreversible, split_reversible = make_irreversible(mm, exclude_unkown_list, True)
		split_reactions = []
		for (i,j) in split_reversible:
			split_reactions.append(i[:-8])

		# Add lump reactions to the mm_irreversible model
		for lump_id, rxn in lumpid_to_rxn.iteritems():
			reaction = parse_reaction(rxn)
			mm_irreversible.database.set_reaction(lump_id, reaction)
			mm_irreversible.add_reaction(lump_id)

		dgf_dict = parse_dgf(mm_irreversible, self._args.dgf_file)
		dgr_dict = calculate_dgr(mm_irreversible, dgf_dict, exclude_unkown_list, self._args.transport_parameters)

		TMFA_Problem = fluxanalysis.FluxBalanceProblem(mm_irreversible, solver)
		TMFA_Problem, cpd_xij_dict = add_conc_constraints(TMFA_Problem, self._args.set_concentrations)

		TMFA_Problem.prob.set_objective(TMFA_Problem.get_flux_var(objective))
		TMFA_Problem.prob.set_objective_sense(lp.ObjectiveSense.Maximize)
		TMFA_Problem.prob.solve()
		result = TMFA_Problem.prob.result
		logger.info('TMFA Problem Status: {}'.format(result.get_value(TMFA_Problem.get_flux_var(objective))))

		TMFA_Problem = add_reaction_constraints(TMFA_Problem, mm_irreversible, exclude_lump_list, exclude_unkown_list,
												exclude_lump_unkown, dgr_dict, lump_to_rxnids, split_reversible)

		# TMFA_Problem.prob.set_objective(TMFA_Problem.get_flux_var(objective))
		# TMFA_Problem.prob.set_objective_sense(lp.ObjectiveSense.Maximize)
		# TMFA_Problem.prob.solve()
		# result = TMFA_Problem.prob.result
		# logger.info('TMFA Problem Status: {}'.format(result.get_value(TMFA_Problem.get_flux_var(objective))))
		# for reaction in mm_irreversible.reactions:
		# 	print(result.get_value(TMFA_Problem.get_flux_var(reaction)), result.get_value('dgri_{}'.format(reaction)))


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
		problem.prob.define(str(cp))
		var = problem.prob.var(str(cp))
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


def calculate_dgr(mm, dgf_dict, excluded_reactions, transport_parameters):
	F = 0.02306
	dpsi = -130
	R = 1.9858775 / 1000
	T = 298
	dph = 0.4
	t_param = {}
	for row in csv.reader(transport_parameters, delimiter=str('\t')):
		rxn, c, h = row
		t_param[rxn] = (c, h)

	dgr_dict = {}
	for reaction in mm.reactions:
		dgt = 0
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
			if reaction in t_param.keys():
				(c, h) = t_param[reaction]
				dgt = c * F * dpsi - 2.3 * h * R * T * dph
		dgr_dict[reaction] = dgf_sum + dgt

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


def make_irreversible(mm, exclude_list, reversible):
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
	for rxn in mm.reactions:
		if rxn not in exclude_list:
			reaction = mm_irrev.get_reaction(rxn)
			if reversible is False:
				if reaction.direction == Direction.Both:
					mm_irrev.remove_reaction(rxn)
					r = Reaction(Direction.Forward, reaction.left, reaction.right)
					r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
					r_id = str('{}_forward'.format(rxn))
					r2_id = str('{}_reverse'.format(rxn))
					mm_irrev.database.set_reaction(r_id, r)
					mm_irrev.database.set_reaction(r2_id, r2)
					mm_irrev.add_reaction(r_id)
					mm_irrev.add_reaction(r2_id)
					split_reversible.append((r_id, r2_id))
			elif reversible is True:
				mm_irrev.remove_reaction(rxn)
				r = Reaction(Direction.Forward, reaction.left, reaction.right)
				r2 = Reaction(Direction.Forward, reaction.right, reaction.left)
				r_id = str('{}_forward'.format(rxn))
				r2_id = str('{}_reverse'.format(rxn))
				mm_irrev.database.set_reaction(r_id, r)
				mm_irrev.database.set_reaction(r2_id, r2)
				mm_irrev.add_reaction(r_id)
				mm_irrev.add_reaction(r2_id)
				split_reversible.append((r_id, r2_id))
	return mm_irrev, split_reversible


def add_reaction_constraints(problem, mm, exclude_lumps, exclude_unknown, exclude_lumps_unknown, dgr_dict,
							 lump_rxn_list, split_rxns):
	R = 1.9858775 / 1000
	T = 298
	k = 1000000
	epsilon = 1e-6
	excluded_cpd_list = ['h2o[e]', 'h2o[c]', 'h[e]', 'h[c]']
	for reaction in mm.reactions:
		# define variables for vmax, dgri, zi, yi, and vi
		vmax = mm.limits[reaction].upper
		problem.prob.define('zi_{}'.format(reaction), types=lp.VariableType.Binary)
		problem.prob.define('yi_{}'.format(reaction), types=lp.VariableType.Binary)
		problem.prob.define('dgri_{}'.format(reaction))
		zi = problem.prob.var('zi_{}'.format(reaction))
		yi = problem.prob.var('yi_{}'.format(reaction))
		vi = problem.get_flux_var(reaction)
		dgri = problem.prob.var('dgri_{}'.format(reaction))
		rxn = mm.get_reaction(reaction)
		# add flux constraint linking vi and zi for all reactions except lumps
		if reaction not in exclude_lumps:
			problem.prob.add_linear_constraints(vi - zi * vmax <= 0)
		# add thermodynamic feasibility constraint for all reactions where dgr0 is known except for lumps
		if reaction not in exclude_lumps_unknown:
			problem.prob.add_linear_constraints(dgri - k + k * zi <= 0 - epsilon)
		# add constraint to calculate dgri based on dgr0 and the concentrations of the metabolites
		if reaction not in exclude_unknown:
			dgr0 = dgr_dict[reaction]
			ssxi = 0
			for (cpd, stoich) in rxn.compounds:
				if cpd not in excluded_cpd_list:
					ssxi += problem.prob.var(str(cpd)) * float(stoich)
			problem.prob.add_linear_constraints(dgri == dgr0 + R * T * ssxi)
		# add constraints for thermodynamic feasibility of lump reactions and to constrain their constituent reactions
	for reaction in mm.reactions:
		if reaction in lump_rxn_list.keys():
			problem.prob.add_linear_constraints(vi == 0)
			problem.prob.add_linear_constraints(dgri - k * yi <= 0 - epsilon)
			sub_rxn_list = lump_rxn_list[reaction]
			sszi = 0
			for sub_rxn in sub_rxn_list:
				sszi += problem.prob.var('zi_{}'.format(sub_rxn))
			problem.prob.add_linear_constraints(yi + sszi <= len(sub_rxn_list))
	# Add linear constraint to disallow solutions that use both the forward and reverse reaction from a split reaction.
	for (forward, reverse) in split_rxns:
		problem.prob.add_linear_constraints(problem.prob.var('zi_{}'.format(forward)) + problem.prob.var('zi_{}'.format(reverse)) <= 1)
	return problem
