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
# Copyright 2018-2020  Ke Zhang <kzhang@my.uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth1@my.uri.edu>

from __future__ import unicode_literals

import logging
import io
import os
import csv
import argparse
from collections import defaultdict, Counter
from six import text_type
from ..datasource.context import FilePathContext, FileMark
from ..command import MetabolicMixin, Command, FilePrefixAppendAction, \
	convert_to_unicode, SolverCommandMixin
from .. import graph
import sys
from ..formula import Atom, _AtomType
try:
	from graphviz import render, FORMATS
except ImportError:
	render = None
	FORMATS = None
from six import string_types
from psamm.expression import boolean
from psamm.importer import count_genes
from psamm.lpsolver import generic
from psamm import fluxanalysis
from psamm.lpsolver import glpk, lp
from psamm.graph import make_network_dict, write_network_dict
from psamm.datasource.reaction import Reaction, Compound, Direction
import plotly.express as px
from psamm.lpsolver import glpk, lp
import json
import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash import ctx
from dash.dependencies import Input, Output, State
import dash_cytoscape as cyto
import re
import pandas as pd
import webbrowser
from psamm.commands.chargecheck import charge_check
from psamm.commands.formulacheck import formula_check
from psamm.datasource.native import ModelWriter, NativeModel
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,DictReactionEntry as ReactionEntry, DictCompartmentEntry as CompartmentEntry)
from decimal import *
import yaml
from psamm.datasource.native import parse_reaction_equation_string
from psamm.datasource.reaction import parse_reaction
from pkg_resources import resource_filename
import webbrowser


# GLOBALS ---------------------------------------------------------------------

logger = logging.getLogger(__name__)
REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'
EPSILON = 1e-5

nodes = []
edges = []
pathway_default = "Central carbon metabolism, Glycolysis or Gluconeogenesis"


# FUNCTIONS -------------------------------------------------------------------

# Gets model and network information from Model files
def read_model(model, mm, el):
	if isinstance(el, list) or el is None:
		el = "C"
	excluded_reactions = []
	for rxn in model.reactions:
		for cpd, v in rxn.equation.compounds:
			if (float(v)).is_integer() is False or float(v) > 10:
				excluded_reactions.append(rxn.id)
	network = make_network_dict(model, mm, subset=None, method='fpp',
		element=el, excluded_reactions=excluded_reactions, reaction_dict={}, analysis=None)
	return model, network

# Builds a reaction set from a model based on a pathway of interest
def get_pathway_list(nm, pathway):
	pathway_list = set()
	rxn_set = set()
	if isinstance(pathway, str):
		pathway = [pathway]
	for path in pathway:
		if path == "None":
			rxn_set = set()
			for i in nm.reactions:
				if 'pathways' in i.properties:
					if isinstance(i.properties['pathways'], list):
						for j in [i.properties['pathways']]:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['pathways'])
				elif 'subsystem' in i.properties:
					if isinstance(i.properties['subsystem'], list):
						for j in i.properties['subsystem']:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['subsystem'])
		elif path != "All":
			for i in nm.reactions:
				if 'pathways' in i.properties:
					if isinstance(i.properties['pathways'], list):
						for j in [i.properties['pathways']]:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['pathways'])
					if path in i.properties['pathways']:
						rxn_set.add(i.id)
				elif 'subsystem' in i.properties:
					if isinstance(i.properties['subsystem'], list):
						for j in i.properties['subsystem']:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['subsystem'])
					if path in i.properties['subsystem']:
						rxn_set.add(i.id)

		else:
			for i in nm.reactions:
				if 'pathways' in i.properties:
					if isinstance(i.properties['pathways'], list):
						for j in i.properties['pathways']:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['pathways'])
				elif 'subsystem' in i.properties:
					if isinstance(i.properties['subsystem'], list):
						for j in i.properties['subsystem']:
							pathway_list.add(j)
					else:
						pathway_list.add(i.properties['subsystem'])
				rxn_set.add(i.id)
	pathway_list = ["None"] + ["All"] + list(pathway_list)
	return pathway_list, rxn_set

# Retreive list of compounds from model 
def get_compounds_list(nm):
	compounds_list = []
	for i in nm.compounds:
		compounds_list.append(i.id)
	return compounds_list

# Breadth first search function for model
def bfs_compounds(start, end, network, rxn_list, rxn_set, middle2, middle3):
	# initialize useful variables
	generic = ['h2o', 'h', 'accoa', 'coa', 'nad', 'nadh', 'nadp', 'nadph', \
		'pi', 'co2', 'ppi', 'q8h2', 'q8', 'atp', 'gtp', 'utp', 'ttp']
	depth_start = {(start, "none"): 0}
	depth_end = {(end, "none"): 0}
	i = 1
	frontier_start = [(start, "none")]
	frontier_end = [(end, "none")]
	parent_start = {}
	parent_start[(start, "none")] = "none"
	parent_end = {}
	parent_end[(end, "none")] = "none"
	frontier_check = "False"
	move_start = {}
	move_end = {}
	final_list = {}
	while frontier_check == "False":
		next = []
		for u in frontier_start:
			adj = {}
			for rxn in network[0]:
				left = []
				for le in rxn._properties['equation'].__dict__['_left']:
					left.append(le[0].__dict__['_name'])
				right = []
				for ri in rxn._properties['equation'].__dict__['_right']:
					right.append(ri[0].__dict__['_name'])
				for cpd in network[0][rxn][0]:
					for j in cpd:
						if j.name == u[0]:
							for k in cpd:
								if j.name in left and k.name in right and (str(rxn._properties['equation']. __dict__['_direction']) == "Direction.Forward" \
									or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
										adj[k.name] = rxn.id
								elif j.name in right and k.name in left and \
								(str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse" \
									or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
										adj[k.name] = rxn.id
			for v, r in adj.items():
				if (v, r) not in depth_start and v not in generic and \
						(v, r) != middle2 and (v, r) != middle3:
						depth_start[(v, r)] = i
						parent_start[(v, r)] = u
						next.append((v, r))
						move_start[(v, r)] = adj[v]
				if (v, r) in depth_end and v not in generic and (v, r) \
						!= middle2 and (v, r) != middle3:
						frontier_check = "True"
						middle = (v, r)
		frontier_start = next
		next = []
		for u in frontier_end:
			adj = {}
			for rxn in network[0]:
				left = []
				for le in rxn._properties['equation'].__dict__['_left']:
					left.append(le[0].__dict__['_name'])
				right = []
				for ri in rxn._properties['equation'].__dict__['_right']:
					right.append(ri[0].__dict__['_name'])
				for cpd in network[0][rxn][0]:
					for j in cpd:
						if j.name == u[0]:
							for k in cpd:
								if j.name in left and k.name in right and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Forward" \
									or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
										adj[k.name] = rxn.id
								elif j.name in right and k.name in left and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse" \
									or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
										adj[k.name] = rxn.id
			for v, r in adj.items():
				if (v, r) not in depth_end and v not in generic and \
						(v, r) != middle2 and (v, r) != middle3:
						depth_end[(v, r)] = i
						parent_end[(v, r)] = u
						next.append((v, r))
						move_end[(v, r)] = adj[v]
				if (v, r) in depth_start and frontier_check != "True" and \
						v not in generic and (v, r) != \
						middle2 and (v, r) != middle3:
						frontier_check = "True"
						middle = (v, r)
		frontier_end = next
		i += 1
		if i > 50:
			return([], [])
	# collect the rxns from the start frontier into a final list of rxns
	i = depth_start[middle]
	j = 1
	final_list[i] = move_start[middle]
	parent = parent_start[middle]
	while parent != (start, "none"):
		i -= 1
		j += 1
		final_list[i] = move_start[parent]
		parent = parent_start[parent]
	# Checks to make sure the solution isnt a 0 path solution
	if depth_end[middle] > 0:
		# Collects the moves from the end frontier into a final list of moves
		j += 1
		final_list[j] = move_end[middle]
		parent = parent_end[middle]
		while parent != (end, "none"):
			j += 1
			final_list[j] = move_end[parent]
			parent = parent_end[parent]
	sorted_list = []
	# Sorts the moves by their index and store them in a final list
	for k in range(1, len(final_list) + 1):
		sorted_list.append(final_list[k])
	return(middle, sorted_list)

# Builds the subset network, updates nodes and edges for network graph
def build_network(self, nm, mm, rxn_set, network, fba_dropdown, nodes, edges):
	# setting up storage data structures for compound info
	name = {}
	formula = {}
	charge = {}
	for i in nm.compounds:
		if i.name:
			name[i.id] = i.name
		else:
			name[i.id] = i.id
		formula[i.id] = i.formula
		charge[i.id] = i.charge
	obj = 0
	if not isinstance(fba_dropdown, list):
		objective = str(fba_dropdown)
		# Set up the solver
		solver = self._get_solver()
		# Define the flux balance
		flux_carrying = defaultdict(lambda: 0)
		try:
			problem = fluxanalysis.FluxBalanceProblem(mm, solver)
			problem.maximize(fba_dropdown)
			obj = problem.get_flux(fba_dropdown)
			for i in mm.reactions:
				flux_carrying[i] = problem.get_flux(i)
		except:
			flux_carrying = {}
			for i in mm.reactions:
				flux_carrying[i] = 0
			logger.info("Solution Not Optimal")
	else:
		flux_carrying = {}
		for i in mm.reactions:
			flux_carrying[i] = 0
	comp = []
	for i in nm.compartments:
		comp.append(i.id)
	# adding compounds and reactions to the network 
	for rxn in network[0]:
		if rxn.id in rxn_set:
			if rxn.name:
				rname = rxn.name
			else:
				rname = rxn.id
			rxn_num = 0
			visited = set()
			# looking through each compound in the reaction
			for cpd in network[0][rxn][0]:
				if cpd[0] not in visited:
					visited.add(cpd[0])
					visited.add(cpd[1])
					# set color based on compartment
					color = "#102A5F" # c = blue
					if not cpd[0].compartment == "c":
						color = "#105f41" # e = green
					# adding nodes for compounds
					nodes.append({'data': {'id': str(cpd[0]), 'label': name[str(cpd[0])[0:(str(cpd[0]).find('['))]],
						'formula': formula[str(cpd[0])[0:(str(cpd[0]).find('['))]], 'compartment': cpd[0].compartment,
						'col': color, 'charge': charge[str(cpd[0])[0:(str(cpd[0]).find('['))]],
						'orig_id': str(cpd[0])[0:(str(cpd[0]).find('['))], 'type': 'compound'}})
					nodes.append({'data': {'id': str(cpd[1]), 'label': name[str(cpd[1])[0:(str(cpd[1]).find('['))]],
						'formula': formula[str(cpd[1])[0:(str(cpd[1]).find('['))]],'compartment': cpd[0].compartment,
						'col': color,'charge': charge[str(cpd[0])[0:(str(cpd[0]).find('['))]],
						'orig_id': str(cpd[0])[0:(str(cpd[0]).find('['))],'type': 'compound'}})
					if 'pathways' in rxn.properties:
						path = rxn.properties['pathways']
					elif 'subsystem' in rxn.properties:
						path = rxn.properties['subsystem']
					else:
						path = ['No pathway exists']
					if "{}_node".format(rxn.id) not in nodes:
						if "name" not in rxn.properties:
							rxn.properties['name'] = "None"
						if "genes" not in rxn.properties:
							rxn.properties['genes'] = "None"
						if "pathways" in rxn.properties:
							path_prop = rxn.properties['pathways']
						elif "subsystem" in rxn.properties:
							path_prop = rxn.properties['subsystem']
						else:
							path_prop = "None"
						if "EC" not in rxn.properties:
							rxn.properties['EC'] = "None"
						# adding nodes for reactions (purple squares)
						nodes.append({'data': {'id': "{}_node".format(rxn.id),'label': rxn.id,'type': 'reaction',
							'equation': str(rxn.properties["equation"]),'pathways': path,'orig_id': rxn.id,'name': rxn.properties["name"],
							'genes': rxn.properties["genes"],'pathways': path_prop,'EC': rxn.properties["EC"],'flux_node': flux_carrying[rxn.id]},
							'style': {'shape': 'rectangle', "background-color": "#41105f"}})
						# adding edges (connections) between nodes 
						# arrow direction on edge depends on the direction of the reaction
						if str(rxn._properties['equation'].__dict__
							   ['_direction']) == "Direction.Forward":
							flux = flux_carrying[rxn.id]
							edges.append({'data': {'id': "".join([rxn.id, "_",str(cpd[1]), "for"]),
								'orig_id': rxn.id,
								'source': "{}_node".format(rxn.id),
								'target': str(cpd[1]),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux}})
							rxn_num += 1
							edges.append({'data': {
								'id': "".join([rxn.id, "_",str(cpd[0]), "for"]),
								'orig_id': rxn.id,
								'source': str(cpd[0]),
								'target': "{}_node".format(rxn.id),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux}})
							rxn_num += 1
						elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse":
							flux = 0 - float(flux_carrying[rxn.id])
							edges.append({'data': {
								'id': "".join([rxn.id, "_",str(cpd[0]), "rev"]),
								'orig_id': rxn.id,
								'source': "{}_node".format(rxn.id),
								'target': str(cpd[0]),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux}})
							rxn_num += 1
							edges.append({'data': {
								'id': "".join([rxn.id, "_",str(cpd[1]), "rev"]),
								'orig_id': rxn.id,
								'source': str(cpd[1]),
								'target': "{}_node".format(rxn.id),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux}})
							rxn_num += 1
						elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both":
							flux = flux_carrying[rxn.id]
							edges.append({'data': {
								'id': "".join([rxn.id, "_", str(cpd[1]),"for"]),
								'orig_id': rxn.id,
								'source': "{}_node".format(rxn.id),
								'target': str(cpd[1]),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux,
								'classes': 'bidirectional'}})
							rxn_num += 1
							flux = 0 - float(flux_carrying[rxn.id])
							edges.append({'data': {
								'id': "".join([rxn.id, "_", str(cpd[0]),"for"]),
								'orig_id': rxn.id,
								'source': "{}_node".format(rxn.id),
								'target': str(cpd[0]),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux,
								'classes': 'bidirectional'}})
							rxn_num += 1
							flux = flux_carrying[rxn.id]
							edges.append({'data': {
								'id': "".join([rxn.id, "_", str(cpd[0]),"rev"]),
								'orig_id': rxn.id,
								'source': str(cpd[0]),
								'target': "{}_node".format(rxn.id),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux,
								'classes': 'bidirectional'}})
							rxn_num += 1
							flux = 0 - float(flux_carrying[rxn.id])
							edges.append({'data': {
								'id': "".join([rxn.id, "_", str(cpd[1]),"rev"]),
								'orig_id': rxn.id,
								'source': str(cpd[1]),
								'target': "{}_node".format(rxn.id),
								'rid': rxn.id,
								'label': "".join([rxn.name]),
								'equation': str(rxn.properties["equation"]),
								'pathways': path,
								'flux': flux,
								'classes': 'bidirectional'}})
							rxn_num += 1
	return obj

# builds the network model for individual reactions (charge check and formula check)
def rxn_network(nm, rxn_eq, reaction_name):
	nodes = []
	edges = []
	# storing compound ids to track compounds involved in the reaction
	rxn_sides = []
	cpds_left = []
	cpds_right = []
	clean_cpds_left = []
	cc_left_comp = []
	clean_cpds_right = []
	cc_right_comp = []
	arrow = ""
	dir_both = False
	# parsing out compounds involved 
	if "<=>" in str(rxn_eq):
		dir_both = True
		arrow = "double"
		rxn_sides = str(rxn_eq).strip().split("<=>")
		cpds_left.extend(rxn_sides[0].strip().split(" + "))
		cpds_right.extend(rxn_sides[1].strip().split(" + "))
		for cpd in cpds_left:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_left_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_left_comp.append("e")
			else:
				append("none")
			clean_cpds_left.append(temp2.strip())
		for cpd in cpds_right:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_right_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_right_comp.append("e")
			else:
				append("none")
			clean_cpds_right.append(temp2.strip())
	elif "<=" in str(rxn_eq):
		arrow = "back"
		rxn_sides = str(rxn_eq).strip().split("<=")
		cpds_left.extend(rxn_sides[0].strip().split(" + "))
		cpds_right.extend(rxn_sides[1].strip().split(" + "))
		for cpd in cpds_left:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_left_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_left_comp.append("e")
			else:
				append("none")
			clean_cpds_left.append(temp2.strip())
		for cpd in cpds_right:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_right_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_right_comp.append("e")
			else:
				append("none")
			clean_cpds_right.append(temp2.strip())
	elif "=>" in str(rxn_eq):
		arrow = "forward"
		rxn_sides = str(rxn_eq).strip().split("=>")
		cpds_left.extend(rxn_sides[0].strip().split(" + "))
		cpds_right.extend(rxn_sides[1].strip().split(" + "))
		for cpd in cpds_left:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_left_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_left_comp.append("e")
			else:
				append("none")
			clean_cpds_left.append(temp2.strip())
		for cpd in cpds_right:
			split_index = cpd.find(")")
			temp = cpd[split_index+1:].strip()
			temp2 = ""
			if "[c]" in temp:
				temp2 = temp.replace("[c]", "")
				cc_right_comp.append("c")
			elif "[e]" in temp:
				temp2 = temp.replace("[e]", "")
				cc_right_comp.append("e")
			else:
				append("none")
			clean_cpds_right.append(temp2.strip())
	# creating dictonaries to store relevant compound info
	names = dict()
	formulas = dict()
	charges = dict()
	for i in nm.compounds:
		if i.name:
			names[i.id] = i.name
		else:
			names[i.id] = i.id
		formulas[i.id] = i.formula
		charges[i.id] = i.charge
	# creating a node for the reaction
	nodes.append({'data': {'id': 'rxn', 'label': str(reaction_name),'type': 'reaction',
		'equation': str(rxn_eq),}, 'style': {'shape': 'rectangle', "background-color":"#41105f"}})
	count = 0
	# creating nodes for all compounds involved in the reaction
	# creating edges to connect all compounds to the reaction
	for cpd in clean_cpds_left:
		color = "#102A5F"
		if cc_left_comp[count] == "e":
			color = "#105f41"
		nodes.append({'data': {'id': str(cpd), 'label': names[str(cpd)], 'formula': formulas[str(cpd)], 
			'charge': charges[str(cpd)], 'type': 'compound'}})
		if dir_both:
			edges.append({'data': {'source': str(cpd), 'target': 'rxn', 'classes': 'bidirectional'}})
		else:
			edges.append({'data': {'source': str(cpd), 'target': 'rxn'}})
		count += 1
	count = 0
	for cpd in clean_cpds_right:
		color = "#102A5F"
		if cc_right_comp[count] == "e":
			color = "#105f41"
		nodes.append({'data': {'id': str(cpd), 'label': names[str(cpd)], 'formula': formulas[str(cpd)], 
			'charge': charges[str(cpd)], 'type': 'compound'}})
		if dir_both:
			edges.append({'data': {'source': 'rxn', 'target': str(cpd), 'classes': 'bidirectional'}})
		else:
			edges.append({'data': {'source': 'rxn', 'target': str(cpd)}})
		count += 1
	return nodes, edges

# returns the initial setup for the app 
def build_app(self):
	'''
	Generates the initial required data for the visualization and
	launches the app on a local port.
	'''
	# styles for the home/simulate/curate/documentation sidebar
	sidebar_style = {
		"position": "fixed",
		"top": 0,
		"left": 0,
		"bottom": 0,
		"width": "16rem",
		"padding": "2rem 1rem",
		"background-color": "#E6EBEE",
	}
	# creating the site navigation sidebar 
	sidebar = html.Div(
		# PSAMM title and description
		[html.H2("PSAMM", className="display-4", style={"color": "#102A5F"}),
		html.Hr(),
		html.P("Interactive Metabolic Modelling", className="lead", style={"color": "#102A5F"},),
		# home, curate, simulate, and documentation buttons 
		dbc.Nav([
			dbc.NavLink("Home", href="/", active="exact", style={"background-color": "#102A5F", "border-radius": "8px", "color": "#E6EBEE"}),
			dbc.NavLink("Simulate", href="/simulate", active="exact", style={"background-color": "#102A5F", "border-radius": "8px",
				"margin-top": "0.25rem", "color": "#E6EBEE"}),
			dbc.NavLink("Curate", href="/curate", active="exact", style={"background-color": "#102A5F", "border-radius": "8px",
				"margin-top": "0.25rem", "color": "#E6EBEE"}),
			dbc.NavLink("Documentation", href="https://psamm.readthedocs.io", active="exact", style={"background-color": "#102A5F",
			 "border-radius": "8px","margin-top": "0.25rem", "color": "#E6EBEE"}
			 )],
			vertical=True,
			pills=True,),
		# open and close menu buttons 
		# only visible on the simulate tab (display element in style is used to do this)
		# CALLBACKS ON LINE 739
		dbc.Nav([
			dbc.Button("Close Menu",id="rsb_btn", style={"margin-top": "5rem", "margin-left":"3rem",
				"background-color": "#5d80c9", "border-radius": "8px", "display": "block"}),
			dbc.Button("Open Menu",id="open_menu_btn", style={"margin-top": "0.25rem", "margin-left":"3rem",
				"background-color": "#5d80c9", "border-radius": "8px", "display": "block"}
				),],
			vertical=True,
			),
		],
		style=sidebar_style,
	)
	# create a container for the page content
	content = html.Div(id="page-content", style={"margin-left": "18rem", "margin-top": "1rem"})
	# return the layout for the app, including the sidebar and page content
	return html.Div([dcc.Location(id="url"), sidebar, content])


# CLASSES ---------------------------------------------------------------------

# HiddenPrints suppresses print statements to avoid cluttering the terminal
# from running backend commands (ex: charge, formulacheck)
class HiddenPrints:
	def __enter__(self):
		self._original_stdout = sys.stdout
		sys.stdout = open(os.devnull, 'w')
	def __exit__(self, exc_type, exc_val, exc_tb):
		sys.stdout.close()
		sys.stdout = self._original_stdout


# Launches the app for PSAMM
class InteractiveCommand(MetabolicMixin, SolverCommandMixin,
						 Command, FilePrefixAppendAction):
	"""Launches an interactive interface for the analysis """

	# INTIALIZERS --------------------------------------------------------------
	
	@classmethod
	def init_parser(cls, parser):
		parser.add_argument(
			'--exclude', metavar='reaction', type=convert_to_unicode,
			default=[], action=FilePrefixAppendAction,
			help='Reaction(s) to exclude from metabolite pair prediction')
		super(InteractiveCommand, cls).init_parser(parser)

	def run(self):
		# sets the bootstrap theme for the app 
		bs_theme = "https://cdn.jsdelivr.net/npm/bootswatch@5.2.3/dist/lux/bootstrap.min.css"
		# launch app 
		self._app = dash.Dash(
			__name__, external_stylesheets=[
				bs_theme], 
				suppress_callback_exceptions=True,
				prevent_initial_callbacks="initial_duplicate",
			)
		self._app.title = "PSAMM"
		server = self._app.server
		if self._app is not None and hasattr(self, "callbacks"):
			self.callbacks(self._app)
		application = build_app(self)
		self._app.layout = application
		# server link 
		webbrowser.open_new("http://localhost:{}".format("8050"))
		self._app.run_server(debug=False)


	# CALLBACKS ----------------------------------------------------------------

	def callbacks(self, _app):

		# OPERATES THE PSAMM SIDEBAR TO NAVIGATE THE APP ----------------------- 
		@_app.callback(
			Output("page-content", "children"), 
			Output("rsb_btn", "style"),
			Output("open_menu_btn", "style"),
			[Input("url", "pathname")],
			Input("rsb_btn", "n_clicks"),
			Input("open_menu_btn", "n_clicks"))
		def render_page_content(pathname, rsb_clicks, open_clicks):
			# HOME PAGE SETUP --------------------------------------------------
			# background image
			image = html.Img(src='https://docs.google.com/drawings/d/e/2PACX-1vQvQsjym3p5llDbQYGRW15GDw-TFiBitfkDye9ZN4BnZqgEJW06XC18gQ-9DZ66F8meeOw1CTyY9CRy/pub?w=1152&h=432',
				style={"width": "100%", "height": "100%"})
			# help page content (explanation of features)
			# CALLBACK TO OPEN HELP PAGE ON LINE 929
			help_content = html.Div(id='help_dd', children=[
				dbc.Row([
					dbc.Col([
						html.Div([
							html.H4("Choose a page:", className="text-center", 
								style={"color": "#102A5F", "text-transform": "none",
								"margin-left": "12rem", "margin-top": "1rem"}),
							], style={"text-align": "right"}),
						], width="auto"),
					dbc.Col([
						html.Div([
							dcc.Dropdown(id="help_dropdown",options=["simulate", "curate"],
								placeholder="", style={"margin-top": "0.5rem", "color": "#102A5F"},),
							]),
						], width=5),
					]),
				])
			# home page layout (contains the image, the help title, and the help content)
			home = html.Div(id="home_page", style={'position':'absolute', 'z-index':'1000', 
				"top": "0px", "left": "16rem", "width": "80%"}, children=[
				image,
				html.H1("PSAMM HELP", id="help_title", className="text-center", style={"color": "#102A5F", "margin-top": "1rem"}),
				help_content,
				html.Div(id="help_space", children=[])
				])

			# SIMULATE MENU SETUP ----------------------------------------------
			# set styles for the menu
			sidebar_style = {
				"width": "22rem",
				"padding": "1rem 1rem",
				"background-color": "#E6EBEE",
			}
			# menu setup (contains the buttons for each simulate tab and a space for them to open up)
			# CALLBACKKS FOR SIMULATE MENU BUTTONS ON LINE 1035
			menu = html.Div(id="menu", children=
				[dbc.Nav([
					dbc.Row([
						html.H3("Model Menu", className="text-center", style={"color": "#102A5F"})
						]),
					dbc.Row([
						html.Div([
								dbc.Button("controls", id='ctrl',color="primary", style={"background-color": "#102A5F", 
									"margin-left": "1.25rem", "border-radius": "8px", "margin-top": "0.5rem"}),
								dbc.Button("constraints", id='cnst',color="primary", style={"margin-top": "0.5rem",
									"margin-left": "0.5rem", "background-color": "#102A5F", "border-radius": "8px"}),
							]),
						]),
					dbc.Row([
						html.Div([
								dbc.Button("exchange", id='exg',color="primary", style={"margin-bottom": "1.25rem",
									"margin-left": "1.25rem", "background-color": "#102A5F", "border-radius": "8px", "margin-top": "0.5rem"},),
								dbc.Button("model stats", id='ms',color="primary", style={"margin-left": "0.5rem", "margin-bottom": "1.25rem",
									"background-color": "#102A5F", "border-radius": "8px", "margin-top": "0.5rem"}),
							]),
						]),
					dbc.Row([
						html.Div(id="menu_content", className="menu", children=[])
						]),
					],
					vertical=True,
					pills=True,
				),], style=sidebar_style
			)

			# SIMULATE MODEL SETUP ---------------------------------------------
			# create the network model using cytoscape
			s_n = cyto.Cytoscape(id='net', layout={'name': 'cose','animate': True,
				'animationDuration': 1000}, style={'width': '100%', 'height': '37.5rem'},
				elements=nodes + edges, stylesheet=[{'selector': 'node', 'style': 
				{'background-color': 'data(col)', 'label': 'data(label)'}},
				{"selector": "edge", "style": {"width": 1, "curve-style": "straight", "target-arrow-shape":"triangle"}},
				{"selector": "[flux > 1e-10]", "style": {"line-color": "blue", "target-arrow-color": "blue", "width": 1,
				"curve-style": "bezier","target-arrow-shape":"triangle","line-dash-pattern": [6, 3],"line-dash-offset": 24}},
				{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
				minZoom=0.2, maxZoom=3, responsive=True)

			# SIMULATE PAGE SETUP ----------------------------------------------
			# return the model, the menu, and the flux analysis report
			sim = dbc.Container([
				dbc.Row([
					html.H5("Metabolic Network Visualized here:", style={"margin-top": "0.5rem", 
						"color": "#102A5F", "margin-left": "1rem"}),
					]),
				dbc.Row([
					dbc.Col([s_n], width=7),
					dbc.Col([menu], width=1),
					]),
				dbc.Row([
					html.Div(id="fba-data", className="results", children=[], style={"margin-top": "1.25rem", 
						"background-color": "#d2def7", "border-radius": "8px", "padding": "1.25rem", "color": "#102A5F"}),
					]),
				], fluid=True)

			# CLOSED MENU MODE -------------------------------------------------
			# check if the button clicked was the 'close menu' button
			if dash.callback_context.triggered[0]['prop_id'] == 'rsb_btn.n_clicks':
				# return the layout for the simulate page with the model and node analysis
				return dbc.Container([
					dbc.Row([s_n]),
					dbc.Row([
						dbc.Col([
							html.Div(id="node-data", className="text-center", children=[
								html.H6("Click on a node to see its details", style={"margin-top": "0.5rem", "color": "#102A5F"})
								], style={"background-color": "#d2def7", "border-radius": "8px", "padding": "0.5rem"}),
							]),
						]),
					]), {"margin-top": "7.25rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}, \
						{"margin-top": "0.25rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}

			# OPENED MENU MODE -------------------------------------------------
			# if 'open menu' button clicked, return the simulate layout
			if dash.callback_context.triggered[0]['prop_id'] == \
					'open_menu_btn.n_clicks':
					return sim, {"margin-top": "7.5rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}, \
					{"margin-top": "0.25rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}

			# CURATE PAGE SETUP ------------------------------------------------
			# styles for the topbar to navigate the curate page
			cur_topbar_style = {"position": "static", "top": 0, "left": "16rem",
				"height": "5rem", "width": "98%", "padding": "1rem", "background-color": "#102A5F",
				"border-radius": "8px",}
			# creating the topbar for the curate page
			# CALLBACK FOR ADD/EDIT/CHARGECHECK/FORMULACHECK/MODELWIDESCREEN TAB BUTTON ON LINE 1925
			cur_menu = html.Div(
				dbc.Row([
					html.Div([
						dbc.Button("add", id='add_btn', color="secondary", style={"border-radius": "8px", "color": "#102A5F",
							"margin-left": "3.45rem"},),
						dbc.Button("edit", id='edit_btn', color="secondary",style={"border-radius": "8px", "color": "#102A5F",
							"margin-left": "2.5rem"},),
						dbc.Button("charge check", id='cc_btn', color="secondary",style={"border-radius": "8px", "color": "#102A5F",
							"margin-left": "2.5rem"},),
						dbc.Button("formula check", id='fc_btn', color="secondary",style={"border-radius": "8px", "color": "#102A5F",
							"margin-left": "2.5rem"},),
						dbc.Button("model widescreen", id='mws_btn', color="secondary",style={"border-radius": "8px", "color": "#102A5F",
							"margin-left": "2.5rem"},),
						], style={"width": "100%"}),
					], style={"width": "110%"}),
				style=cur_topbar_style,
			)
			cur_content = html.Div(id="cur_content", children=[
				html.H5("Start curating the metabolic model with the menu above", className="text-center", 
					style={"margin-top": "0.5rem", "color": "#102A5F"}),
				s_n
				])

			# SIDEBAR PAGE BUTTON FUNCTIONALITY --------------------------------
			# returns the proper page layout based on the button clicked on the sidebar
			if pathname == "/":
				# home page returned, along with style components to hide
				# the open/close menu buttons (displayed on the simulate page only)
				return home, {"display": "none"}, {"display": "none"}
			elif pathname == "/simulate":
				# simulate page returned, along with style components to show
				# the open/close menu buttons (displayed on the simulate page only)
				return sim, {"margin-top": "7.5rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}, \
				{"margin-top": "0.25rem", "background-color": "#5d80c9", "border-radius": "8px", "display": "block"}
			elif pathname == "/curate":
				# curate page returned, along with style components to hide
				# the open/close menu buttons (displayed on the simulate page only)
				return html.Div([cur_menu, cur_content]), {"display": "none"}, {"display": "none"}
			# if the user tries to reach a different page, return a 404 message
			return html.Div([
					html.H1("404: Not found", className="text-danger"),
					html.Hr(),
					html.P(f"The pathname {pathname} was not recognised..."),
					],className="p-3 bg-light rounded-3",)
		# END OF SIDEBAR OPERATION CALLBACK ------------------------------------

		# OPERATES HELP PAGE MENU ----------------------------------------------
		@_app.callback(
			Output("help_space", "children"),
			Input("help_dropdown", "value"),
		)
		def open_help(page):
			# opens different help messages based on the page selected in the dropdown
			if page == "simulate":
				# simulate page: explains each dropdown in the control menu
				return html.Div([
					dbc.Row([
					  dbc.Col([
						html.Div([
						  html.H4("Pathways", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F"}),
						  html.H5("Use the dropdown to display pathways in the metabolic model. Multiple can be selected.",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  dbc.Col([
						html.Div([
						  html.H4("Reactions", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F"}),
						  html.H5("Use the dropdown to display reactions in the metabolic model. Multiple can be selected.",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  dbc.Col([
						html.Div([
						  html.H4("Gene Delete", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F"}),
						  html.H5("""Use the dropdown to delete genes from the model. If all genes associated with a reaction are
							deleted, the reaction will be removed from analysis.""",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  ]),
					dbc.Row([
					  dbc.Col([
						html.Div([
						  html.H4("Element Transfer Networks", className="text-center", style={"margin-top": "1rem", "color": "#102A5F"}),
						  html.H5("Use the dropdown to choose a chemical element for the network (default Carbon).",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  dbc.Col([
						html.Div([
						  html.H4("Path Search", className="text-center", style={"margin-top": "1rem", "color": "#102A5F"}),
						  html.H5("Use the path search to conduct a bidirectional breadth first search to show the shortest path between two compounds.",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  dbc.Col([
						html.Div([
						  html.H4("Compound Search", className="text-center", style={"margin-top": "1rem", "color": "#102A5F"}),
						  html.H5("Use the compound search to show all reactions containing this compound. Select ALL for pathways to see all reactions.",
							className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						  ]),
						], width=4),
					  ]),
					dbc.Row([
					  html.Div([
						html.H4("Flux Analysis", className="text-center", style={"margin-top": "1rem", "color": "#102A5F"}),
						html.H5("Use the dropdown to choose a reaction to optimize for flux balance analysis. Flux is visualized in blue.",
						  className="text-center", style={"color": "#102A5F", "text-transform": "none"}),
						]),
					  ]),
					])
			elif page == "curate":
				# curate page: explains what each tab in the contorl topbar is used for
				return html.Div([
					dbc.Row([
					  dbc.Col([
						html.Div([
						  html.H4("Add", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", "margin-left": "6.25rem"}),
						  html.H5("""To add a new reaction to the model: fill in id, name, equation, and pathway
							To add a new compound to the model: fill in id, name, charge, and formula""",
							className="text-center", style={"color": "#102A5F", "text-transform": "none", "margin-left": "6.25rem"}),
						  ]),
						]),
					  dbc.Col([
						html.Div([
						  html.H4("Edit", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", "margin-left": "0.5rem"}),
						  html.H5("""To edit reactions in the model: fill in direction and details for compounds on the left and right.
							To edit compounds in the model: fill in charge and formula.""",
							className="text-center", style={"color": "#102A5F", "text-transform": "none", "margin-left": "0.5rem"}),
						  ]),
						]),
					  ]),
					dbc.Row([
					  dbc.Col([
						html.Div([
						  html.H4("Chargecheck / Formulacheck", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", "margin-left": "6.25rem"}),
						  html.H5("""Select an unbalanced reaction from the dropdown to view relevant data. To download
							all unbalanced reactions and associated data, click the 'Download List & Data' button. Edit current
							compounds involved in the reaction or add new compounds to the reaction by filling in the id, stoich,
							formula, charge, and component fields.""",
							className="text-center", style={"color": "#102A5F", "text-transform": "none", "margin-left": "6.25rem"}),
						  ]),
						]),
					  dbc.Col([
						html.Div([
						  html.H4("Model Widescreen", className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", "margin-left": "0.5rem"}),
						  html.H5("View your changes to the model here.",
							className="text-center", style={"color": "#102A5F", "text-transform": "none", "margin-left": "0.5rem"}),
						  ]),
						])
					  ]),
					])
		# END OF HELP PAGE MENU CALLBACK ---------------------------------------

		# OPERATES SIMULATE PAGE MENU ------------------------------------------
		@_app.callback(
			Output('menu_content', 'children'),
			Input('ctrl', 'n_clicks'),
			Input("cnst", "n_clicks"),
			Input("exg", "n_clicks"),
			Input("ms", "n_clicks"),
			prevent_initial_call=True
		)
		# REDIRECTS APP TO OPEN CORRECT MENU -----------------------------------
		def control_panel_director(ctrl_clicks, cnst_clicks, exg_clicks, ms_clicks):
			triggered_id = ctx.triggered_id
			if triggered_id == "ctrl":
				return open_ctrl(ctrl_clicks)
			elif triggered_id == "cnst":
				return open_cnst(cnst_clicks)
			elif triggered_id == "exg":
				return open_exg(exg_clicks)
			elif triggered_id == "ms":
				return open_ms(ms_clicks)
		# RETURNS THE CONTROL TAB TO THE MENU ----------------------------------
		def open_ctrl(nclicks):
			# DEFINITIONS ------------------------------------------------------
			# enable svg export
			cyto.load_extra_layouts()
			# Start with default of carbon tranfer network
			el = 'C'
			# generate the initial network
			pathway = 'All'
			pathway_list, rxn_set = get_pathway_list(self._model, "All")
			rxns = list(rxn_set)
			# get a list of all compound in the model
			compounds_list = get_compounds_list(self._model)
			# generate information about number of genes in the model
			count = 0
			rxns_with_genes = 0
			gene_content = set()
			rxns_full = []
			for reaction in self._model.reactions:
				rxns_full.append(reaction.id)
				count += 1
				if reaction.genes is None or reaction.genes == '':
					continue
				elif isinstance(reaction.genes, string_types):
					rxns_with_genes += 1
					assoc = boolean.Expression(reaction.genes)
				else:
					rxns_with_genes += 1
					variables = [boolean.Variable(g) for g in reaction.genes]
					assoc = boolean.Expression(boolean.And(*variables))
				gene_content.update(v.symbol for v in assoc.variables)
			# CREATING CONTROL TAB LAYOUT --------------------------------------
			if dash.callback_context.triggered[0]['prop_id'] == \
					'ctrl.n_clicks':
				return html.Div(id='control-tab', className='control-tab', children=[
					# contains dropdowns for each model option
					# dropdowns are loaded with options from above definitions and model
					dbc.Row([
						html.Div([
							html.H5("Pathways", className="text-center",
								style={"margin-top": "0.5rem", "color": "#102A5F"}),
							dcc.Dropdown(id="pathways_dropdown",options=[
								{"label": i,"value": i, } for i in list(pathway_list)],
								value=pathway_list[0], multi=True, placeholder="",
								optionHeight=80, style={"margin-top": "0.5rem", "color": "#102A5F"}),
							]),
						]),
					dbc.Row([
						html.Div([
							html.H5("Reactions", className="text-center",
								style={"margin-top": "0.5rem", "color": "#102A5F"}),
							dcc.Dropdown(id="reaction_dropdown", options=[
								{"label": i.id,"value": i.id, } for i in list(self._model.reactions)],
								value=pathway_list[0],multi=True,placeholder="",
								style={"margin-top": "0.5rem", "color": "#102A5F"}),
							])
						]),
					dbc.Row([
						html.Div([
							html.H5("Element Transfer Networks", className="text-center", style={"margin-top": "0.75rem", "color": "#102A5F"}),
							dcc.Dropdown(id="element_dropdown",options=[
								{"label": i,"value": i, } for i in ["C", "N", "S", "P"]],
								value=["C", "N", "S", "P"],multi=False),
							])
						]),
					dbc.Row([
						html.Div([
							html.H5("Compound Search", className="text-center", style={"margin-top": "0.75rem", "color": "#102A5F"}),
							dcc.Dropdown(id="compounds_dropdown",options=[
								{"label": i,"value": i, } for i in list(compounds_list)],
								value=compounds_list,multi=False,style={"width": "100%", "color": "#102A5F"}),
							])
						]),
					dbc.Row(html.H5("Path Search", className="text-center", style={"margin-top": "0.75rem", "color": "#102A5F"}),),
					dbc.Row([
						dbc.Col([
							html.Div([
								dcc.Dropdown(id="filter1_dropdown",options=[
									{"label": i,"value": i, } for i in list(compounds_list)],
									value=compounds_list, multi=False),
								]),
							]),
						dbc.Col([
							html.Div([
								dcc.Dropdown(id="filter2_dropdown", options=[
									{"label": i,"value": i, } for i in list(compounds_list)],
									value=compounds_list, multi=False),
								])
							]),
						]),
					dbc.Row([
						html.Div([
							html.H5("Gene Delete", className="text-center", style={"margin-top": "0.75rem", "color": "#102A5F"}),
							dcc.Dropdown(id="delete_dropdown",options=[
								{"label": i,"value": i, } for i in list(gene_content)],
								value=None,multi=True,style={"width": "100%", "color": "#102A5F"}),
							]),
						]),
					dbc.Row([
						html.Div([
							html.H5("Flux Analysis", className="text-center", style={"margin-top": "0.75rem", "color": "#102A5F"}),
							dcc.Dropdown(id="fba_dropdown",options=[{
								"label": i,"value": i, } for i in list(rxns_full)],
								value=rxns,multi=False,style={"width": "100%", "color": "#102A5F"},),
							]),
						]),
					# buttons to submit and download the current model
					# CALLBACK FOR SUBMIT BUTTON ON LINE 1310
					# CALLBACK FOR DOWNLOAD CSV BUTTON ON LINE 1418
					# CALLBACK FOR DOWNLOAD PNG BUTTON ON LINE 1474
					dbc.Row([
						html.Div([
							dbc.Button("Submit",id="btn_sub", size="sm", style={"margin-top": "1.5rem", "margin-left": "0.25rem",
								"background-color": "#5d80c9", "border-radius": "8px"}),
							dbc.Button("Download CSV",id="btn_tsv", size="sm", style={"margin-top": "1.5rem",
								"background-color": "#5d80c9","border-radius": "8px", "margin-left": "0.25rem"}),
							dcc.Download(id="download_csv"),
							dbc.Button("Download png", id="btn-get-png", size="sm", style={"margin-top": "1.5rem",
								"background-color": "#5d80c9", "border-radius": "8px", "margin-left": "0.25rem"}, n_clicks=None),
							])
						]),
					])
		# RETURNS THE CONTRAINTS LAYOUT TO THE MENU ----------------------------
		def open_cnst(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == 'cnst.n_clicks':
				return html.Div([
							# text areas for column headers and user warnings / confirmations
							dbc.Row([
								html.Div([
									html.H5(id="constraints", className="text-center", children="",
									style={"color":"#102A5F", "margin-top":"0.5rem"}),
									html.H5(id="constraint-confirm", className="text-center", children="",
										style={"color":"#102A5F", "margin-top":"0.5rem", "text-transform": "none"}),
									html.H5(id="constraint-save-confirm", className="text-center", children="",
										style={"color":"#102A5F", "margin-top":"0.5rem", "text-transform": "none"}),
									]),
								]),
							# buttons to load and update constraints
							# CALLBACK FOR LOAD CONSTRAINTS BUTTON ON LINE 1569
							# CALLBACK FOR UPDATE CONSTRAINTS BUTTON ON LINE 1528
							dbc.Row([
								html.Div([
									dbc.Button("Load constraints", id="constraint-load", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-left":"1.1rem", "margin-top": "0.5rem"}),
									dbc.Button("Update constraints", id="constraint-modify", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-left": "0.5rem", "margin-top": "0.5rem"}),
									]),
								]),
							# button to save constraints
							# CALLBACK FOR SAVE CONSTRAINTS BUTTON ON LINE 1489
							dbc.Row([
								html.Div([
									dbc.Button("Save constraints", id="constraint-save", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-left": "6rem", "margin-top": "0.5rem"})
									]),
								]),
						])
		# RETURNS THE EXCHANGE LAYOUT TO THE MENU ------------------------------
		def open_exg(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == 'exg.n_clicks':
				return html.Div([
						# text area for column headers
						html.H6(id="exchange", children="",
							className="text-center", style={"color": "#102A5F", "margin-top": "0.5rem"}),
						html.Div([
							# buttons to load, update, and save exchange 
							# CALLBACK FOR LOAD EXCHANGE BUTTON ON LINE 1758
							# CALLBACK FOR UPDATE EXCHANGE BUTTON ON LINE 1706
							# CALLBACK FOR SAVE EXCHANGE BUTTON ON LINE 1670
							dbc.Row([
								html.Div([
									dbc.Button("Load exchange", id="exchange-load", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-top": "0.5rem", "margin-left": "2rem"}),
									dbc.Button("Update exchange", id="exchange-modify", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-top": "0.5rem", "margin-left": "0.5rem"}),
									]),
								html.Div([
									dbc.Button("Save exchange", id="exchange-save", size="sm", style={"background-color": "#5d80c9",
										"border-radius": "8px", "margin-top": "0.5rem", "margin-left": "6.5rem"}),
									]),
								])
							]),
						# text areas for user warnings / confirmations
						html.H6(id="exchange-confirm", children="",
							className="text-center", style={"color": "#102A5F", "margin-top": "0.5rem", "text-transform": "none"}),
						html.H6(id="exchange-save-confirm", children="",
							className="text-center", style={"color": "#102A5F", "margin-top": "0.5rem", "text-transform": "none"}),
					])
		# RETURNS MODEL STATISTICS LAYOUT TO THE MENU --------------------------
		def open_ms(nclicks):
			# DEFINITIONS ------------------------------------------------------
			# enable svg export
			cyto.load_extra_layouts()
			# Start with default of carbon tranfer network
			el = 'C'
			# generate the initial network
			pathway = 'All'
			pathway_list, rxn_set = get_pathway_list(self._model, "All")
			rxns = list(rxn_set)
			# generate information about number of genes in the model
			count = 0
			rxns_with_genes = 0
			# stores the initial model stats in content_model and generate the
			# initial stats
			content_model = []
			content_model.append(html.H5("General Model Statistics:", id="gen_title", className="text-center", 
				style={"color": "#102A5F"}))
			content_model.append(html.P("Reactions in model: " + str(count), id="stats_rxns", className="text-center", 
				style={"color": "#5d80c9"}))
			count = 0
			for i in self._model.compounds:
				count += 1
			content_model.append(html.P("Compounds in model: " + str(count), id="stats_cmpd", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(
				html.P("Pathways in mode: " + str(len(pathway_list)), id="stats_paths", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(html.P("Reactions with gene assocations: " + str(rxns_with_genes),
				id="stats_rxngenes", className="text-center", style={"color": "#5d80c9"}))
			try:
				content_model.append(html.P("Total genes: " + str(count_genes(self._model)), 
					id="stats_genes", className="text-center", style={"color": "#5d80c9"}))
			except BaseException:
				content_model.append(html.P("Total genes: Unable to count genes", id="nogenes", className="text-center", 
				style={"color": "#5d80c9"}))
			with HiddenPrints():
				unbalanced, count, unchecked, exclude, unbalanced_list, \
					unbalance_dict = charge_check(self, 1e-6)
			content_model.append(html.H5("Chargecheck results:", id="cc_title", className="text-center", 
				style={"color": "#102A5F"}))
			content_model.append(html.P("Unblanced reactions: " + str(unbalanced), id="ubchargecheck", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked), id="bchargecheck", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude)), id="exclrxns", className="text-center", 
				style={"color": "#5d80c9"}))
			with HiddenPrints():
				unbalanced, count, unchecked, exclude, form_list, \
					form_dict = formula_check(self)
			content_model.append(html.H5("Formulacheck results:", id="fc_title", className="text-center", 
				style={"color": "#102A5F"}))
			content_model.append(html.P("Unbalanced reactions: " + str(unbalanced), id="ubrxns", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked), id="ucrxns", className="text-center", 
				style={"color": "#5d80c9"}))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude)), id="exclrxns2", className="text-center", 
				style={"color": "#5d80c9"}))
			# CREATING MODEL STATS TAB CONTENT ---------------------------------
			if dash.callback_context.triggered[0]['prop_id'] == 'ms.n_clicks':
				return html.Div([
					# returns information from above
					html.H4("Model Statistics", id="model-stuff", className="text-center", 
						style={"color": "#102A5F", "margin-top": "0.5rem"}),
					html.Div(content_model),
					])
		# END OF OPERATION OF SIMULATE MENU CALLBACK ---------------------------


		# CALLBACKS TO OPERATE FUNCTIONS ON SIMULATE PAGE ----------------------

		# FUNCTIONALITY OF CONTROL TAB SUBMIT BUTTON ---------------------------
		@_app.callback(
			[Output("net", "elements"),
			Output("fba-data", "children"),],
			[Input("btn_sub", "n_clicks"), ],
			[State("pathways_dropdown", "value"),
			 State("element_dropdown", "value"),
			 State("compounds_dropdown", "value"),
			 State("fba_dropdown", "value"),
			 State("filter1_dropdown", "value"),
			 State("filter2_dropdown", "value"),
			 State("delete_dropdown", "value"),
			 State("reaction_dropdown", "value")],
			prevent_initial_call=True, )
		# takes input from all dropdowns on control menu and creates nodes and edges for model
		def filter_nodes(n_clicks, pathways_dropdown, element_dropdown,
			compounds_dropdown, fba_dropdown,filter1_dropdown, 
			filter2_dropdown, delete_dropdown, r_dropdown):
			# check that a dropdown has a value to display in the model
			if isinstance(pathways_dropdown, list) \
					or isinstance(delete_dropdown, list) \
					or isinstance(element_dropdown, str) \
					or isinstance(compounds_dropdown, str) \
					or isinstance(fba_dropdown, str) \
					or (isinstance(filter1_dropdown, str)
						and isinstance(filter1_dropdown, str)) \
					or isinstance(r_dropdown, list):
				model = self._mm.copy()
				# add genes and reactions to model 
				if delete_dropdown is not None:
					genes = set()
					gene_assoc = {}
					deleted_reactions = set()
					for reaction in self._model.reactions:
						assoc = None
						if reaction.genes is None or reaction.genes == '':
							continue
						elif isinstance(reaction.genes, string_types):
							assoc = boolean.Expression(reaction.genes)
						else:
							variables = [boolean.Variable(g) for g in reaction.genes]
							assoc = boolean.Expression(boolean.And(*variables))
						genes.update(v.symbol for v in assoc.variables)
						gene_assoc[reaction.id] = assoc
						reactions = set(model.reactions)
					# take out reactions associated with delete genes dropdown value
					for reaction in reactions:
						if reaction not in gene_assoc:
							continue
						assoc = gene_assoc[reaction]
						if any(boolean.Variable(gene) in assoc.variables
								for gene in delete_dropdown):
							new_assoc = assoc.substitute(
								lambda v: v if v.symbol not in
								delete_dropdown else False)
							if new_assoc.has_value() and not new_assoc.value:
								logger.info('Deleting reaction {}...'.format(reaction))
								deleted_reactions.add(reaction)
					for r in deleted_reactions:
						model.remove_reaction(r)
				# create the model network
				nm, network = read_model(self._model, model, element_dropdown)
				pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)
				if r_dropdown is not None:
					for i in r_dropdown:
						rxn_set.add(i)
				if isinstance(filter1_dropdown, str) and isinstance(filter2_dropdown, str):
					rxn_list = set()
					cpd_list = [filter1_dropdown, filter2_dropdown]
					middle2 = []
					middle3 = []
					middle2, rxn_list = bfs_compounds(filter1_dropdown, filter2_dropdown,
						network, rxn_list,rxn_set, middle2, middle3)
					if len(rxn_list) == 0:
						return []
				elif isinstance(compounds_dropdown, str):
					rxn_list = []
					for rxn in network[0]:
						for cpd in network[0][rxn][0]:
							for i in cpd:
								if i.name == compounds_dropdown and rxn.id in rxn_set:
									rxn_list.append(rxn.id)
				else:
					rxn_list = rxn_set
				if fba_dropdown is None:
					fba_dropdown = []
				nodes.clear()
				edges.clear()
				# add nodes and edges to the network graph 
				obj_flux = build_network(self, nm, model,rxn_list, network,
					fba_dropdown, nodes, edges)
				elements = nodes + edges
				# report object flux 
				if obj_flux != 0:
					contents = []
					contents.append(html.H5("{} flux ""is {}".format(fba_dropdown, obj_flux), 
						style={"text-transform": "none", "color": "#102A5F"}))
				else:
					contents = []
					contents.append(html.H5("No flux through the network", 
						style={"text-transform": "none", "color": "#102A5F"}))
				# return the nodes/edges and the flux reaction text to their elements
				return elements, contents
			else:
				return tuple(), list()

		# FUNCTIONALITY OF THE DOWNLOAD CSV BUTTON ON CONTROL TAB --------------
		@_app.callback(
			Output("download_csv", "data"),
			Input("btn_tsv", "n_clicks"),
			prevent_initial_call=True
		)
		# creates a csv based on data currently in nodes and edges global variables
		def create_csv(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "btn_tsv.n_clicks":
				node_info = dict()
				# collecting node info
				for node in nodes: 
					if node['data']['type'] == "reaction":
						continue
					temp_dict = dict()
					temp_dict["id"] = node['data']['id']
					temp_dict["name"] = node['data']['label']
					if "formula" in node['data'].keys():
						temp_dict["formula"] = node['data']['formula']
					temp_dict["equations"] = set()
					temp_dict["pathways"] = set()
					temp_dict["rxns_involved"] = set()
					node_info[node['data']['id']] = temp_dict.copy()
				# collecting edge info 
				for edge in edges:
					source = edge['data']['source']
					if source in node_info.keys():
						node_info[source]["equations"].add(edge['data']['equation'])
						node_info[source]["pathways"].add(edge['data']['pathways'])
						target = edge['data']['target']
						if target[0] == "R":
							node_info[source]["rxns_involved"].add(target.replace("_node", ""))
					target = edge['data']['target']
					if target in node_info.keys():
						node_info[target]["equations"].add(edge['data']['equation'])
						node_info[target]["pathways"].add(edge['data']['pathways'])
						source = edge['data']['source']
						if source[0] == "R":
							node_info[target]["rxns_involved"].add(source.replace("_node", ""))
				# creating the csv dataframe object 
				fields = ["Name", "ID", "Formula", "Equations Involved in", "Pathways Involved in", "Reactions Involved in"]
				translations = {"Name": "name", "ID": "id", "Formula": "formula", "Equations Involved in": "equations", \
					"Pathways Involved in": "pathways", "Reactions Involved in": "rxns_involved"}
				all_data = []
				for field in fields:
					temp_list = []
					for cpd in node_info.keys():
						if type(node_info[cpd][translations[field]]) == set:
							temp_list.append("; ".join(node_info[cpd][translations[field]]))
						else:
							temp_list.append(node_info[cpd][translations[field]])
					all_data.append(temp_list)
				dataframe = pd.DataFrame({"Name": all_data[0], "ID": all_data[1], "Formula": all_data[2], \
					"Equations Involved in": all_data[3], "Pathways Involved in": all_data[4], "Reactions Involved in": all_data[5]})
				return dcc.send_data_frame(dataframe.to_csv, "ModelNodeData.csv")

		# FUNCTIONALITY OF THE DOWNLOAD PNG BUTTON ON THE CONTROL TAB ----------
		@_app.callback(
			Output("net", "generateImage"),
			Input("btn-get-png", "n_clicks"),
			prevent_initial_call=True,
		)
		# downloads the cytoscape graph
		def create_png(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "btn-get-png.n_clicks" and \
				dash.callback_context.triggered[0]['value'] is not None:
				return {"type": "png", "action": "download"}
			else:
				return {}

		# FUNCTIONALITY OF THE SAVE CONSTRAINT BUTTON --------------------------
		@_app.callback(
			Output("constraint-save-confirm", "children"),
			Input("constraint-save", "n_clicks"),
			State("constraints", "children")
		)
		# saves to YAML and sends message to constraint menu page
		def save_file(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == 'constraint-save.n_clicks':
				# helps parse through the datalist 
				# the datalist stores all the elements on the screen and their values
				def rec_props(data, out):
					for i in data:
						if isinstance(i, dict):
							if 'children' in i['props']:
								rec_props(i['props']['children'], out)
							elif 'value' in i['props']:
								out.append(i['props']['value'])
							elif 'placeholder' in i['props']:
								out.append(i['props']['placeholder'])
					return out
				# find all changes the user inputted on the page
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)
				# add changes to YAML file 
				path = FilePathContext(self._args.model)
				with open('{}/limits_curated.yaml'.format(path), "w") as f:
					for e in ex:
						if e[0] != '':
							f.write(
								"{}\t{}\t{}\n".format(
									e[0], e[1], e[2]))
				return("Constraints saved!")

		# FUNCTIONALITY OF THE UPDATE CONSTRAINTS BUTTON -----------------------
		@_app.callback(
			Output("constraint-confirm", "children"),
			Input("constraint-modify", "n_clicks"),
			State("constraints", "children")
		)
		# adds constraint to model and sends message to constraints menu
		def modify_constraint(nclicks, datalist):
			# helps parse through the datalist 
			# the datalist stores all the elements on the screen and their values
			def rec_props(data, out):
				for i in data:
					if isinstance(i, dict):
						if 'children' in i['props']:
							rec_props(i['props']['children'], out)
						elif 'value' in i['props']:
							out.append(i['props']['value'])
						elif 'placeholder' in i['props']:
							out.append(i['props']['placeholder'])
				return out
			if dash.callback_context.triggered[0]['prop_id'] == 'constraint-modify.n_clicks':
				# find all changes the user inputted on the page
				const = []
				for i in datalist:
					if isinstance(i, dict):
						const_temp = rec_props(i['props']['children'], [])
						if len(const_temp) > 0:
							const.append(const_temp)
				# add changes to model
				rxns = []
				for i in self._model.reactions:
					rxns.append(i.id)
				for c in const:
					if c[0] in rxns:
						self._model._limits[c[0]] = (c[0], Decimal(c[1]), Decimal(c[2]))
						self._mm._limits_lower[c[0]] = Decimal(c[1])
						self._mm._limits_upper[c[0]] = Decimal(c[2])
					else:
						return("Added constraints for reaction not in model")

		# FUNCTIONALITY OF THE LOAD CONSTRAINTS BUTTON
		@_app.callback(
			Output("constraints", "children"),
			Input("constraint-load", "n_clicks")
		)
		# puts the input boxes and current constraints on constaints menu page
		def load_constraints(nclicks):
			rxn_set = set()
			for i in self._model.reactions:
				rxn_set.add(i.id)
			rxn_list = list(rxn_set)
			contents = []
			contents.append(
				# column titles 
				dbc.Row([
					dbc.Col([
						html.P("ID")], width=4),
					dbc.Col([
						html.P("Lower Bound")], width=4),
					dbc.Col([
						html.P("Upper Bound")], width=4),
				]),)
			if dash.callback_context.triggered[0]['prop_id'] == 'constraint-load.n_clicks':
				# get constraint data from model 
				constraint = set()
				for c in self._model._limits:
					lower = self._model._limits[c][1]
					upper = self._model._limits[c][2]
					if lower is not None:
						constraint.add(c)
					if upper is not None:
						constraint.add(c)
				rlist = set()
				for i in self._model.reactions:
					rlist.add(i.id)
				for c in constraint:
					if self._model._limits[c][1] is None:
						lower = -float(self._model.default_flux_limit)
					else:
						lower = float(self._model._limits[c][1])
					if self._model._limits[c][2] is None:
						upper = float(self._model.default_flux_limit)
					else:
						upper = float(self._model._limits[c][2])
					contents.append(
						# adding current constraints to the page
						dbc.Row([
							dbc.Col([
								dbc.Input(id="id{}".format(str(c)),
									type="text",
									size="sm",
									placeholder=str(c),
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
							dbc.Col([
								dbc.Input(id="lower{}".format(str(c)),
									type="number",
									size="sm",
									placeholder=lower,
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
							dbc.Col([
								dbc.Input(id="upper{}".format(str(c)),
									type="number",
									size="sm",
									placeholder=upper,
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
						]),)
				contents.append(html.P("Add New Constraint:", style={"margin-top": "0.5rem"}))
				contents.append(
					# input boxes to add new constraints to the model
					dbc.Row([
						dbc.Col([
							dbc.Input(id="id{}".format("new"),
								type="text",
								size="sm",
								placeholder=str(''),
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
						dbc.Col([
							dbc.Input(id="lower{}".format("new"),
								type="number",
								size="sm",
								placeholder=-1000,
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
						dbc.Col([
							dbc.Input(id="upper{}".format(str(c)),
								type="number",
								size="sm",
								placeholder=1000,
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=4),
					]),
				)
			return(contents)

		# FUNCTIONALITY OF THE SAVE EXCHANGE BUTTON ----------------------------
		@_app.callback(
			Output("exchange-save-confirm", "children"),
			Input("exchange-save", "n_clicks"),
			State("exchange", "children")
		)
		# saves the changes to a YAML and sends a message to exchange menu
		def save_file(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == 'exchange-save.n_clicks':
				# helps parse through the datalist 
				# the datalist stores all the elements on the screen and their values
				def rec_props(data, out):
					for i in data:
						if isinstance(i, dict):
							if 'children' in i['props']:
								rec_props(i['props']['children'], out)
							elif 'value' in i['props']:
								out.append(i['props']['value'])
							elif 'placeholder' in i['props']:
								out.append(i['props']['placeholder'])
					return out
				# getting all user inputs from the page 
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)
				# saving user inputs back to the model YAML files 
				path = FilePathContext(self._args.model)
				with open('{}/exchange_curated.yaml'.format(path), "w") as f:
					for e in ex:
						if e[0] != '':
							f.write("{}\t{}\t{}\t{}\n".format(e[0], e[1], e[2], e[3]))
				return("Exchange saved!")

		# FUNCTIONALITY OF THE UPDATE EXCHANGE BUTTON --------------------------
		@_app.callback(
			Output("exchange-confirm", "children"),
			Input("exchange-modify", "n_clicks"),
			State("exchange", "children")
		)
		# modifies the model and sends a message to the exchange menu 
		def modify_exchange(nclicks, datalist):
			# helps parse through the datalist 
			# the datalist stores all the elements on the screen and their values
			def rec_props(data, out):
				for i in data:
					if isinstance(i, dict):
						if 'children' in i['props']:
							rec_props(i['props']['children'], out)
						elif 'value' in i['props']:
							out.append(i['props']['value'])
						elif 'placeholder' in i['props']:
							out.append(i['props']['placeholder'])
				return out
			if dash.callback_context.triggered[0]['prop_id'] == 'exchange-modify.n_clicks':
				# getting the user inputs from the page 
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)
				# saving user inputs back to the model instance 
				c = []
				for i in self._model.compounds:
					c.append(i.id)
				rxns = []
				for i in self._model.reactions:
					rxns.append(i.id)
				for e in ex:
					if e[0] in c:
						tup = (Compound(e[0]).in_compartment(e[1]), None, e[2], e[3])
						self._model._exchange[Compound(e[0]).in_compartment(e[1])] = tup
						self._mm._limits_lower["EX_{}[{}]".format(e[0], e[1])] = float(e[2])
						self._mm._limits_upper["EX_{}[{}]".format(e[0], e[1])] = float(e[3])
						if self._mm._database.has_reaction("EX_{}[{}]".format(e[0], e[1])) is False:
							properties = {"id": "EX_{}[{}]".format(e[0], e[1]), 
								"equation": Reaction(Direction.Forward, [(Compound(e[0]).in_compartment(e[1]), -1)])}
							reaction = ReactionEntry(properties, filemark=None)
							self._model.reactions.add_entry(reaction)
							self._mm._reaction_set.add(reaction.id)
							self._mm._database.set_reaction(reaction.id, Reaction(Direction.Forward, [
								(Compound(e[0]).in_compartment(e[1]), -1)]))
					else:
						return("Added exchange compound not in model")

		# FUNCTIONALITY OF THE LOAD EXCHANGE BUTTON ----------------------------
		@_app.callback(
			Output("exchange", "children"),
			Input("exchange-load", "n_clicks")
		)
		# creates the input rows to view current exchanges and add new exchanges to model
		def load_exchange(nclicks):
			# setting compound list 
			comp_set = set()
			for i in self._model.reactions:
				for j in i.equation.compounds:
					comp_set.add(j[0].compartment)
			comp_list = list(comp_set)
			contents = []
			contents.append(
				# column titles 
				dbc.Row([
					dbc.Col([
						html.P("ID")], width=3, className="text-center"),
					dbc.Col([
						html.P("COMPRT")], width=3, className="text-center"),
					dbc.Col([
						html.P("Lower")], width=3, className="text-center"),
					dbc.Col([
						html.P("Upper")], width=3, className="text-center"),
				]),)
			if dash.callback_context.triggered[0]['prop_id'] == 'exchange-load.n_clicks':
				# getting exchange data from model 
				exchange = set()
				for e in self._model._exchange:
					lower = self._model._exchange[e][2]
					upper = self._model._exchange[e][3]
					if lower is None:
						lower = -self._model.default_flux_limit
					if upper is None:
						upper = self._model.default_flux_limit
					if float(lower) != 0 and float(lower) != -1000:
						exchange.add(e)
					if float(upper) != 0 and float(upper) != 1000:
						exchange.add(e)
				for e in exchange:
					contents.append(
						# creating input rows to show current exchange in model
						dbc.Row([
							dbc.Col([
								dbc.Input(id="id{}".format(str(e.name)),
									type="text",
									size="sm",
									placeholder=str(e.name),
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
							dbc.Col([
								dcc.Dropdown(
									id="compartment{}".format(str(e.name)),
									options=[{"label": x,
											  "value": x, }
											 for x in list(comp_list)],
									placeholder=e.compartment,
									multi=False,
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"}
								)
							], width=3),
							dbc.Col([
								dbc.Input(id="lower{}".format(str(e.name)),
									size="sm",
									type="number",
									placeholder=float(self._model._exchange[e][2]),
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
							dbc.Col([
								dbc.Input(id="upper{}".format(str(e.name)),
									size="sm",
									type="number",
									placeholder=float(self._model._exchange[e][3]),
									style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
						]),)
				contents.append(html.H5("Add New Exchange Constraint:", className="text-center", style={"color": "#102A5F",
					"margin-bottom": "0.5rem", "margin-top": "1.1rem"}))
				contents.append(
					# input row to add new exchange to the model 
					dbc.Row([
						dbc.Col([
							dbc.Input(id="id{}".format("new"),
								size="sm",
								type="text",
								placeholder=str(''),
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
						dbc.Col([
							dcc.Dropdown(
								id="compartment{}".format("new"),
								options=[{"label": x,"value": x, } for x in list(comp_list)],
								placeholder='',
								multi=False,
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"}
							)
						], width=3),
						dbc.Col([
							dbc.Input(id="lower{}".format("new"),
								size="sm",
							   	type="number",
								placeholder=-1000,
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
						dbc.Col([
							dbc.Input(id="upper{}".format(str(e.name)),
								type="number",
								size="sm",
								placeholder=1000,
								style={"color": "#102A5F",
									"margin-bottom": "0.25rem"})], width=3),
					]),
				)
			return(contents)

		# FUNCTIONALITY TO CLICK NODES AND SEE DATA (SIMULATE) -----------------
		@_app.callback(
			Output("node-data", "children"),
			[Input("net", "selectedNodeData")]
		)
		# displays data on nodes based on user clicks
		# viewed when menu is closed only 
		def display_nodedata(datalist):
			if datalist is not None:
				if len(datalist) > 0:
					contents = []
					data = datalist[-1]
					# adding node data to the simulate page 
					if "formula" in data:
						if len(datalist) > 0:
							data = datalist[-1]
							contents.append(html.H5("ID: " + data["id"].title(), style={"text-transform": "none", "color": "#102A5F"}))
							contents.append(html.P("Name: " + data["label"], style={"text-transform": "none", "color": "#345591"}))
							contents.append(html.P("Formula: " + str(data["formula"]), style={"text-transform": "none", "color": "#345591"}))
							contents.append(html.P("Charge: " + str(data["charge"]), style={"text-transform": "none", "color": "#345591"}))
					elif "equation" in data:
						if isinstance(data["pathways"], str):
							path = [data["pathways"]]
						else:
							path = data["pathways"]
						contents.append(html.H5("ID: " + data["orig_id"].title(), style={"text-transform": "none", "color": "#102A5F"}))
						contents.append(html.P("Name: " + data["label"], style={"text-transform": "none", "color": "#345591"}))
						contents.append(html.P("Equation: " + str(data["equation"]), style={"text-transform": "none", "color": "#345591"}))
						if isinstance(data["pathways"], str):
							contents.append(html.P("Pathways: " + data["pathways"], style={"text-transform": "none", "color": "#345591"}))
						elif isinstance(data["pathways"], list):
							contents.append(html.P("Pathways: " + ";".join(data["pathways"]), style={"text-transform": "none", "color": "#345591"}))
						contents.append(html.P("Flux: " + str(data["flux_node"]), style={"text-transform": "none", "color": "#345591"}))
					return html.Div(contents)
			else:
				# in the case where no data is found on the node, have the user click another node
				return html.Div([
					html.H5("Select a node to see data", style={"text-transform": "none", "color": "#102A5F"})
					])

		# END OF CALLBACKS TO OPERATE FUNCTIONS ON SIMULATE PAGE ---------------


		# CALLBACKS TO OPERATE THE CURATE PAGE ---------------------------------

		# OPERATES THE CURATE TOPBAR MENU --------------------------------------
		@_app.callback(
			Output('cur_content', 'children'),
			Input('add_btn', 'n_clicks'),
			Input("edit_btn", "n_clicks"),
			Input("cc_btn", "n_clicks"),
			Input("fc_btn", "n_clicks"),
			Input("mws_btn", "n_clicks"),
			prevent_initial_call=True
		)
		# redirects the app to the correct function based on which button was pressed
		def cur_topbar_director(add_clicks, edit_clicks, cc_clicks, fc_clicks, mws_clicks):
			triggered_id = ctx.triggered_id
			if triggered_id == "add_btn":
				return open_add(add_clicks)
			elif triggered_id == "edit_btn":
				return open_edit(edit_clicks)
			elif triggered_id == "cc_btn":
				return open_cc(cc_clicks)
			elif triggered_id == "fc_btn":
				return open_fc(fc_clicks)
			elif triggered_id == "mws_btn":
				return open_mws(mws_clicks)
		# opens add tab
		def open_add(nclicks):
			# styles for the add and edit menus
			ae_style = {
				"width": "22rem",
				"padding": "1rem 1rem",
				"background-color": "#E6EBEE",
				"margin-top": "1.25rem",
			}
			# DEFINITIONS ------------------------------------------------------
			# generate the initial network
			pathway = 'All'
			pathway_list, rxn_set = get_pathway_list(self._model, "All")
			# add menu: inputs and buttons to add new reactions and compounds to the model 
			# CALLBACK FOR ADD REACTION SAVE BUTTON ON LINE 2267
			# CALLBACK FOR ADD COMPOUND SAVE BUTTON ON LINE 2301
			# CALLBACK FOR SAVE AND VIEW MODEL BUTTON ON LINE 2332
			return dbc.Row([dbc.Col([
						html.Div(id="add_model_space", children=[
							html.H4("Choose a pathway and save to view visualization", id="add_model_notice", className="text-center",
								style={"color": "#102A5F", "margin-top": "1.25rem"})
							],)
						], width=7),
					dbc.Col([
						dbc.Nav([
						dbc.Row([
							html.H3("Add to Model", className="text-center", style={"color": "#102A5F"})
							]),
						dbc.Row([
							html.H5("Add a new reaction", id="add_rxn_descript", className="text-center", 
								style={"color": "#102A5F", "text-transform": "none", "margin-bottom": "0.5rem"}),
							]),
						dbc.Row([
							dbc.Col([
								dbc.Input(id="add_rxn_id_input", type="text", placeholder="ID",
									className="text-center",style={"color": "#102A5F", "margin-bottom": "0.25rem",
									"border-radius": "8px"}),
								]),
							dbc.Col([
								dbc.Input(id="add_rxn_name_input", type="text", placeholder="NAME",
									className="text-center",style={"color": "#102A5F", "margin-bottom": "0.25rem",
									"border-radius": "8px"}),
								]),
							]),
						dbc.Row([
							dbc.Col([
								dbc.Input(id="add_rxn_equation_input", type="text", placeholder="EQUATION",
									className="text-center",style={"color": "#102A5F", "margin-bottom": "0.25rem",
									"border-radius": "8px"}),
								]),
							]),
						dbc.Row([
							dbc.Col([
								dbc.Input(id="add_rxn_pathway_input", type="text", placeholder="PATHWAY",
									className="text-center",style={"color": "#102A5F", "margin-bottom": "0.5rem",
									"border-radius": "8px"}),
								]),
							]),
						dbc.Row([
							html.Div([dbc.Button("save reaction", id="add_rxn_btn", style={"background-color": "#102A5F",
								"margin-left": "5.25rem", "margin-bottom": "1.25rem", "border-radius": "8px"}),]),
							]),
						dbc.Row([
							html.H5("Add a new compound", id="add_cpd_descript", className="text-center", 
								style={"color": "#102A5F", "text-transform": "none", "margin-bottom": "0.5rem"}),
							]),
						dbc.Row([
							dbc.Col([
								dbc.Input(id="add_cpd_id_input", type="text", placeholder="ID",
									className="text-center", style={"color": "#102A5F", "margin-bottom": "0.25rem",
									"border-radius": "8px"}),
								]),
							dbc.Col([
								dbc.Input(id="add_cpd_name_input", type="text", placeholder="NAME",
									className="text-center", style={"color": "#102A5F", "margin-bottom": "0.25rem",
									"border-radius": "8px"}),
								]),
							]),
						dbc.Row([
							dbc.Col([
								dbc.Input(id="add_cpd_charge_input", type="text", placeholder="CHARGE",
									className="text-center", style={"color": "#102A5F", "margin-bottom": "0.5rem",
									"border-radius": "8px"}),
								]),
							dbc.Col([
								dbc.Input(id="add_cpd_formula_input", type="text", placeholder="FORMULA",
									className="text-center", style={"color": "#102A5F", "margin-bottom": "0.5rem",
									"border-radius": "8px"}),
								]),
							]),
						dbc.Row([
							html.Div([dbc.Button("save compound", id="add_cpd_btn", style={"background-color": "#102A5F",
								"margin-left": "5rem", "border-radius": "8px"})]),
							]),
						dbc.Row([
							html.Div([
								html.H5("Select a pathway to view model...", className="text-center",
									style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}),
								dcc.Dropdown(id="add_path_dropdown",options=[
									{"label": i,"value": i, } for i in list(pathway_list)],
									value=pathway_list[0], multi=True, placeholder="",
									optionHeight=80, style={"margin-top": "0.5rem", "color": "#102A5F"}),
								dbc.Button("save & view model", id="add_save_btn", style={"background-color": "#102A5F",
								"margin-left": "4.25rem", "border-radius": "8px", "margin-top": "0.25rem"})]),
							]),
						dbc.Row([
							html.Div(id="add_rxn_confirm", children=[]),
							html.Div(id="add_cpd_confirm", children=[])]),
					], style=ae_style, vertical=True, pills=True,
				)],),
			])
		# opens edit tab
		def open_edit(nclicks):
			# styles for the add and edit menus
			ae_style = {
				"width": "22rem",
				"padding": "1rem 1rem",
				"background-color": "#E6EBEE",
				"margin-top": "1.25rem",
			}
			# DEFINITIONS ------------------------------------------------------
			# generate the initial network
			pathway = 'All'
			pathway_list, rxn_set = get_pathway_list(self._model, "All")
			rxns_full = []
			for reaction in self._model.reactions:
				rxns_full.append(reaction.id)
			compounds_list = get_compounds_list(self._model)
			# edit menu: dropdowns and submit buttons to confirm which reaction 
			# or compound the user wants to edit 
			# CALLBACK FOR CONFIRM REACTION EDIT BUTTON ON LINE 2416
			# CALLBACK FOR ADD REACTION TO MODEL BUTTON ON LINE 2629
			# CALLBACK FOR CONFIRM COMPOUND EDIT BUTTON ON LINE 2569
			# CALLBACK FOR ADD COMPOUND TO MODEL BUTTON ON LINE 2727
			# CALLBACK FOR SAVE AND VIEW MODEL BUTTON ON LINE 2761
			return dbc.Row([dbc.Col([
						html.Div(id="edit_model_space", children=[
							html.H4("Choose a pathway and save to view visualization", id="edit_model_notice", className="text-center",
								style={"color": "#102A5F", "margin-top": "1.25rem"})
							],)
						], width=7),
					dbc.Col([
						dbc.Nav([
						dbc.Row([
							html.H3("Edit Model", className="text-center", style={"color": "#102A5F"})
							]),
						dbc.Row([
							html.Div([
								html.H5("Reaction to Edit", id="edit_rxn_descript", className="text-center", 
									style={"color": "#102A5F", "text-transform": "none", "margin-bottom": "0.5rem"}),
								dbc.Row([
									dbc.Col([
										dcc.Dropdown(id="edit_rxn_dropdown",options=rxns_full,
											value='', multi=True, placeholder="",
											optionHeight=80, style={"margin-top": "1.1rem", "color": "#102A5F"}),
										]),
									dbc.Col([
										dbc.Button("confirm edit", id="cfrm_edit_btn", style={"background-color": "#102A5F",
											"border-radius": "8px", "margin-top": "0.5rem"}),
										]),
									]),
								html.Div(id="edit_rxn_space", children=[]),
								dbc.Row([
									html.Div([
										dbc.Button("Add to Model", id="submit_rxn_btn", style={"display": "none"}),
										]),
									]),
								]),
							]),
						dbc.Row([
							html.Div([
								html.H5("Compound to Edit", id="edit_cpd_descript", className="text-center", 
									style={"color": "#102A5F", "text-transform": "none", "margin-bottom": "0.5rem",
									"margin-top": "0.75rem"}),
								dbc.Row([
									dbc.Col([
										dcc.Dropdown(id="edit_cpd_dropdown",options=[{"label": i,
											"value": i} for i in list(compounds_list)], value='', multi=True, placeholder="",
											optionHeight=80, style={"margin-top": "1.1rem", "color": "#102A5F"}),
										]),
									dbc.Col([
										dbc.Button("confirm edit", id="cfrm_edit_btn2", style={"background-color": "#102A5F",
											"border-radius": "8px", "margin-top": "0.5rem"}),
										]),
									]),
								html.Div(id="edit_cpd_space", children=[]),
								dbc.Row([
									html.Div([
										dbc.Button("Add to Model", id="submit_rxn_btn2", style={"display": "none"}),
										]),
									]),
								]),
							]),
						dbc.Row([
							html.Div([
								html.H5("Select a pathway to view model...", className="text-center",
									style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}),
								dcc.Dropdown(id="edit_path_dropdown",options=[
									{"label": i,"value": i, } for i in list(pathway_list)],
									value=pathway_list[0], multi=True, placeholder="",
									optionHeight=80, style={"margin-top": "0.5rem", "color": "#102A5F"}),
								dbc.Button("save & view model", id="edit_save_btn", style={"background-color": "#102A5F",
								"margin-left": "4.25rem", "border-radius": "8px", "margin-top": "0.25rem"})]),
							]),
					], style=ae_style, vertical=True, pills=True,
				)],),
			])
		# opens the charge check tab
		def open_cc(nclicks):
			# getting information on which reactions have a charge imbalance
			with HiddenPrints(): 
				unbalanced, count, unchecked, exclude, unbalanced_list, unbalance_dict = charge_check(self, 1e-6)
			contents = []
			# dislaying a menu that shows the unbalanced reactions
			# includes a button to download data on the unbalanced reactions
			# CALLBACK FOR UNBALANCED REACTIONS DROPDOWN ON LINE 2849
			# CALLBACK FOR DOWNLOAD LIST AND DATA BUTTON ON LINE 3150
			contents.append(
				dbc.Row([
					dbc.Col([
						html.Div([
							html.H4("Unbalanced Reactions:", style={"margin-top": "1.5rem", 
								"color": "#102A5F", "margin-bottom": "1.5rem"}),
							], style={"text-align": "left"})
						], width="auto"),
					dbc.Col([
						html.Div([
							dcc.Dropdown(id="unbalanced_rxns",options=[{"label": i,
								"value": i, } for i in list(unbalanced_list)],
								multi=False, placeholder="", style={"color": "#102A5F", "margin-top": "17px"}),
							], style={"text-align": "left", "width": "19rem"}),
						], width="auto"),
					dbc.Col([
						html.Div([
							dbc.Button("Download List & Data", id="download_unbal", style={
							"background-color": "#102A5F", "border-radius": "8px", "margin-top": "0.5rem",
							"margin-bottom": "0.5rem"}),
							dcc.Download(id="down_list_cc")
							], style={"text-align": "right"}),
						], width="auto"),
					], style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
							"margin-top":"0.5rem","margin-left": "1rem"}),
				)
			contents.append(html.Div(id="cc_space", children=[]))
			return contents
		# opens the formula check tab 
		def open_fc(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "fc_btn.n_clicks":
				# getting information on which reactions have a stiochiometric imbalance
				with HiddenPrints():
					unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(self)
				contents = []
				# dislaying a menu that shows the unbalanced reactions
				# includes a button to download data on the unbalanced reactions
				# CALLBACK FOR UNBALANCED REACTIONS DROPDOWN ON LINE 3391
				# CALLBACK FOR DOWNLOAD LIST AND DATA BUTTON ON LINE 3709
				contents.append(
					dbc.Row([
						dbc.Col([
							html.Div([
								html.H4("Unbalanced Reactions:", style={"margin-top": "1.5rem", 
									"color": "#102A5F", "margin-bottom": "1.5rem"}),
								], style={"text-align": "left"})
							], width="auto"),
						dbc.Col([
							html.Div([
								dcc.Dropdown(id="unbalanced_rxns_fc",options=[{"label": i,
									"value": i, } for i in list(form_list)],
									multi=False, placeholder="", style={"color": "#102A5F", "margin-top": "17px"}),
								], style={"text-align": "left", "width": "19rem"}),
							], width="auto"),
						dbc.Col([
							html.Div([
								dbc.Button("Download List & Data", id="download_unbal_fc", style={
								"background-color": "#102A5F", "border-radius": "8px", "margin-top": "0.5rem",
								"margin-bottom": "0.5rem"}),
								dcc.Download(id="down_list_fc")
								], style={"text-align": "right"}),
							], width="auto"),
						], style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
								"margin-top":"0.5rem","margin-left": "1rem"}),
					)
				contents.append(html.Div(id="fc_space", children=[]))
				return contents
		# opens the model widescreen tab
		def open_mws(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "mws_btn.n_clicks":
				# creating the default model 
				# default modle is automatically displayed on the model widescreen page
				model = self._mm
				nm, network = read_model(self._model, model, "C")
				pathway_list, rxn_set = get_pathway_list(nm, pathway_default)
				if nodes is None or edges is None or len(nodes) == 0 or len(edges) == 0:
					rxn_list = rxn_set
					fba_dropdown = []
					obj_flux = build_network(self, nm, model, rxn_list, network, fba_dropdown, nodes, edges)
				# cytoscape for default pathway
				c_n = cyto.Cytoscape(id='net_curate', layout={'name': 'cose', 'animate': True,
					'animationDuration': 1000}, style={'width': '100%', 'height': '37.5rem'}, 
					elements=nodes + edges, stylesheet=[{'selector': 'node', 'style': {'background-color':
						'data(col)','label': 'data(label)'}}, {"selector": "edge", "style": {"width": 1,
						"curve-style": "straight","target-arrow-shape":"triangle"}},{"selector": "[flux > 0]",
					"style": {"line-color": "blue","target-arrow-color":"blue","width": 1,"curve-style": "bezier",
						"target-arrow-shape":"triangle","line-dash-pattern":[6, 3], "line-dash-offset": 24}}],
					minZoom=0.2, maxZoom=3, responsive=True)
				# menu for selecting a pathway
				# CALLBACK FOR PATHWAY DROPDOWN ON LINE 3962
				return [
					dbc.Row([
						dbc.Col([
							html.Div([
								dcc.Dropdown(id="ws_path_dropdown",options=[{"label": i,
									"value": i,} for i in list(pathway_list)],
									multi=True, placeholder="Pathway...", optionHeight=80, style={
									"color": "#102A5F", "margin-top": "0.75rem","margin-left": "0.5rem"}),
								]),
							], width=5),
						dbc.Col([
							html.Div([
								dbc.Button("Submit", id="ws_sub_btn", size="sm", style={
								"background-color": "#102A5F", "border-radius": "8px", "margin-top": "0.5rem"}),
								]),
							], width=3),
						]),
					dbc.Row([html.Div(id="ws_model_space", children=[c_n])]),
					dbc.Row([html.Div(id="ws_node_data", children=[
						html.H5("Click on a node to see more information", className="text-center",
							style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}),
							], style={"background-color": "#d2def7", "border-radius": "8px", "padding": "0.5rem"})
						]),
					]

		# FUNCTIONALITY OF THE ADD TAB SUBMIT REACTION BUTTON ------------------ 
		@_app.callback(
			Output("add_rxn_confirm", "children"),
			Input("add_rxn_btn", "n_clicks"),
			[State("add_rxn_id_input", "value"),
			 State("add_rxn_name_input", "value"),
			 State("add_rxn_equation_input", "value"),
			 State("add_rxn_pathway_input", "value")]
		)
		# saves the new reaction to the model 
		def save_rxn(nclicks, id_, name, eq, path):
			if dash.callback_context.triggered[0]['prop_id'] == 'add_rxn_btn.n_clicks':
				# check if all fields are filled out
				if id_ is not None and name is not None and eq is not None and path is not None:
					# add reaction info from user to the model 
					equation = parse_reaction_equation_string(eq, "c")
					properties = {"id": id_,"name": name,"equation": equation,"pathways": [path]}
					reaction = ReactionEntry(properties, filemark=None)
					self._model.reactions.add_entry(reaction)
					self._mm._reaction_set.add(reaction.id)
					self._mm._database.set_reaction(reaction.id, equation)
					return([
						html.H5("Added {} to the model".format(id_), id="add_save",
							className="text-center", style={"color": "#5d80c9", "margin-top": "0.5rem",
							"text-transform": "none"}),
						])
				else:
					# error message if not all fields were filled out 
					return([
						html.H5("Error: Please enter a reaction to add to the model", id="add_save",
							className="text-center", style={"color": "#5d80c9", "margin-top": "0.5rem",
							"text-transform": "none"}),
						])

		# FUNCTIONALITY OF THE ADD TAB SUBMIT COMPOUND BUTTON ------------------
		@_app.callback(
			Output("add_cpd_confirm", "children"),
			Input("add_cpd_btn", "n_clicks"),
			[State("add_cpd_id_input", "value"),
			 State("add_cpd_name_input", "value"),
			 State("add_cpd_charge_input", "value"),
			 State("add_cpd_formula_input", "value")]
		)
		# saves a new compound to the model 
		def save_rxn(nclicks, id_, name, charge, form):
			if dash.callback_context.triggered[0]['prop_id'] == 'add_cpd_btn.n_clicks':
				# check if all fields are filled out 
				if id_ is not None and name is not None and charge is not None and form is not None:
					# add compound info from user to the model 
					properties = {"id": id_,"name": name,"charge": charge,"formula": form}
					compound = CompoundEntry(properties, filemark=None)
					self._model.compounds.add_entry(compound)
					return([
						html.H5("Added {} to the model".format(id_), id="add_save",
							className="text-center", style={"color": "#5d80c9", "margin-top": "0.5rem",
							"text-transform": "none"}),
						])
				else:
					# error message if not all fields were filled out
					return([
						html.H5("Error: Please enter a compound to add to the model", id="add_save",
							className="text-center", style={"color": "#5d80c9", "margin-top": "0.5rem",
							"text-transform": "none"}),
						])

		# FUNCTIONALITY OF THE ADD SAVE AND VIEW MODEL BUTTON ------------------
		@_app.callback(
			[Output("add_model_space", "children")],
			[Input("add_save_btn", "n_clicks"),
			Input("add_path_dropdown", "value")],
			prevent_initial_call=True,)
		# saves user input to the model YAML files and displays the network graph
		def save_model(clicks, add_path_dropdown):
			if dash.callback_context.triggered[0]['prop_id'] == 'add_save_btn.n_clicks':
				# set up the model to be added to 
				out_mw = ModelWriter()
				path = FilePathContext(self._args.model)
				with open('{}/model.yaml'.format(path), 'r') as f:
					modelfile = yaml.safe_load(f)
				exported = []
				# add reaction entries to the model and YAML files 
				for i in modelfile['reactions']:
					out_nm = NativeModel()
					with open('{}/{}'.format(path, i['include']), 'r') as f:
						reactionfile = yaml.safe_load(f)
						for j in reactionfile:
							exported.append(j['id'])
							if j['id'] in self._model.reactions:
								out_nm.reactions.add_entry(self._model.reactions[j['id']])
						with open('{}/{}'.format(path, i['include']), "w") as f:
							out_mw.write_reactions(f, out_nm.reactions)
				out_nm = NativeModel()
				for j in self._model.reactions:
					if j.id not in exported:
						out_nm.reactions.add_entry(self._model.reactions[j.id])
				with open('{}/{}'.format(path,
						  "reactions_curated.yaml"), "a+") as f:
					out_mw.write_compounds(f, out_nm.reactions)
				# add compound entries to the model and YAML files 
				exported = []
				for i in modelfile['compounds']:
					out_nm = NativeModel()
					with open('{}/{}'.format(path, i['include']), 'r') \
							as f:
						compoundfile = yaml.safe_load(f)
						for j in compoundfile:
							exported.append(j['id'])
							if j['id'] in self._model.compounds:
								out_nm.compounds.add_entry(
									self._model.compounds[j['id']])
						with open('{}/{}'.format(path,
												 i['include']), "w") as f:
							out_mw.write_compounds(f, out_nm.compounds)
				out_nm = NativeModel()
				for j in self._model.compounds:
					if j.id not in exported:
						out_nm.compounds.add_entry(self._model.compounds[j.id])
				with open('{}/{}'.format(path, "compounds_curated.yaml"),
						  "a+") as f:
					out_mw.write_compounds(f, out_nm.compounds)
				# creating the cytoscape visualization to show network
				model = self._mm
				nm, network = read_model(self._model, model, "C")
				pathway_list, rxn_set = get_pathway_list(nm, add_path_dropdown)
				rxn_list = rxn_set
				fba_dropdown = []
				nodes.clear()
				edges.clear()
				obj_flux = build_network(self, nm, model, rxn_list, network, fba_dropdown, nodes, edges)
				c_n = cyto.Cytoscape(id='net_curate', layout={'name': 'cose','animate': True,
					'animationDuration': 1000},style={'width': '100%', 'height': '37.5rem'},
					elements=nodes + edges,stylesheet=[{'selector': 'node',
						'style': {'background-color':'data(col)','label': 'data(label)'}},
						{"selector": "edge","style": {"width": 1,"curve-style": "straight",
							"target-arrow-shape":"triangle"}},
						{"selector": "[flux > 0]","style": {"line-color": "blue",
							"target-arrow-color":"blue","width": 1,"curve-style": "bezier",
							"target-arrow-shape":"triangle","line-dash-pattern":[6, 3],
							"line-dash-offset": 24}},
						{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
					minZoom=0.2, maxZoom=3, responsive=True)
				return [c_n]
			else:
				# error message because user did not select a pathway to display
				return ["Please select a pathway"]

		# FUNCTIONALITY OF THE CONFIRM REACTION EDIT BUTTON --------------------
		@_app.callback(
			[Output("edit_rxn_space", "children"),
			Output("submit_rxn_btn", "style")],
			[Input("cfrm_edit_btn", "n_clicks"),
			Input("edit_rxn_dropdown", "value")],
		)
		# opens the info panel to edit a specific reaction, selected from the dropdown
		def open_ers(clicks, edit_rxn_dropdown):
			if dash.callback_context.triggered[0]['prop_id'] == "cfrm_edit_btn.n_clicks":
				# check if a reaction was inputted, send error message if not
				if edit_rxn_dropdown == "":
					return html.Div(id="ret_edit", children=[
						html.H5("Error: Please select a reaction to edit", id="rxn_error",
							className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
						]), {"display":"none"}
				# find the reaction in the model corresponding with the reaction 
				# chosen in the dropdown 
				model = self._model
				r = edit_rxn_dropdown
				if r is not None and len(r) > 0:
					for i in model.reactions:
						if isinstance(r, list):
							if i.id == r[0]:
								reaction = i
						else:
							if i.id == r:
								reaction = i
					# find info for the compounds in the selected reaction
					c_list = []
					formula_dict = {}
					charge_dict = {}
					comp_set = set()
					for i in model.compounds:
						c_list.append(i.id)
						formula_dict[i.id] = i.formula
						charge_dict[i.id] = i.charge
					for i in model.reactions:
						for j in i.equation.compounds:
							comp_set.add(j[0].compartment)
					comp_list = list(comp_set)
					# print out info for the reaction 
					contents = []
					contents.append(html.H5("ID: {}".format(reaction.id), id="ers_id_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "0.5rem"}),)
					contents.append(html.H5("Name: {}".format(reaction.name), id="ers_name_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "0.5rem"}),)
					# dropdown for the user to select a direction for the reaction
					contents.append(dbc.Row(id="dir", children=[
						dbc.Col([
							html.H5("Direction:", id="ers_dir_text", className="text-center",
								style={"color": "#102A5F", "text-transform": "none", "margin-top": "0.5rem"}),
							]),
						dbc.Col([
							dcc.Dropdown(id="direction", options=[{"label": x,
								"value": x, } for x in ["forward", "reverse", 
								"both"]], multi=False, style={"color": "#102A5F", "border-radius": "8px"})
							]),
						]),)
					# inputs and dropdowns for the user to change info about each compound
					# in the reaction (on left and right)
					contents.append(html.H5("Current Left Side Compounds", id="ers_lsc_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "1.1rem"}),)
					dropid = 0
					for j in reaction.equation.left:
						contents.append(
							dbc.Row([
								dbc.Col([
									dbc.Input(id="stoich{}".format(str(dropid)),type="number", size="sm",
										placeholder=str(j[1]), style={"color": "#102A5F", "border-radius": "8px",
										"margin-top": "0.25rem"})
									]),
								dbc.Col([
									dcc.Dropdown(id="drop{}".format(str(dropid)), options=[{"label": x, 
										"value": x, } for x in list(c_list)], placeholder=j[0].name, multi=False,
										style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
									]),
								dbc.Col([
									dcc.Dropdown(id="compartment{}".format(str(dropid)), options=[{"label": x, 
										"value": x, } for x in list(comp_list)], placeholder=j[0].compartment, multi=False,
										style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
									]), 
								], id="left"),)
						dropid += 1
					# user can also add a new compound to each side 
					contents.append(html.H5("Add a new Compound to the Left Side", id="ers_lsc_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "1.1rem"}),)
					contents.append(
						dbc.Row([
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(dropid)),type="number", size="sm",
									placeholder=0, style={"color": "#102A5F", "border-radius": "8px",
									"margin-top": "0.25rem"}),
								]),
							dbc.Col([
								dcc.Dropdown(id="drop{}".format(str(dropid)), options=[{"label": x, 
									"value": x, } for x in list(c_list)], placeholder="", multi=False,
									style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
							]),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(dropid)), options=[{"label": x, 
									"value": x, } for x in list(comp_list)], placeholder="", multi=False,
									style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
								]), 
							], id="left_new"))
					dropid += 1
					contents.append(html.H5("Current Right Side Compounds", id="ers_lsc_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "1.1rem"}),)
					for j in reaction.equation.right:
						contents.append(
							dbc.Row([
								dbc.Col([
									dbc.Input(id="stoich{}".format(str(dropid)),type="number", size="sm",
										placeholder=str(j[1]), style={"color": "#102A5F", "border-radius": "8px",
										"margin-top": "0.25rem"}),
									]),
								dbc.Col([
									dcc.Dropdown(id="drop{}".format(str(dropid)), options=[{"label": x, 
										"value": x, } for x in list(c_list)], placeholder=j[0].name, multi=False,
										style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
									]),
								dbc.Col([
									dcc.Dropdown(id="compartment{}".format(str(dropid)), options=[{"label": x, 
										"value": x, } for x in list(comp_list)], placeholder=j[0].compartment, multi=False,
										style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
									]), 
								], id="right"), )
						dropid += 1
					contents.append(html.H5("Add a new Compound to the Right Side", id="ers_lsc_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "1.1rem"}),)
					contents.append(
						dbc.Row([
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(dropid)),type="number", size="sm",
									placeholder=0, style={"color": "#102A5F", "border-radius": "8px",
									"margin-top": "0.25rem"}),
								]),
							dbc.Col([
								dcc.Dropdown(id="drop{}".format(str(dropid)), options=[{"label": x, 
									"value": x, } for x in list(c_list)], placeholder="", multi=False,
									style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
								]),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(dropid)), options=[{"label": x, 
									"value": x, } for x in list(comp_list)], placeholder="", multi=False,
									style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"})
								]), ], id="right_new"))
					dropid += 1
					return html.Div(id="ret_edit", children=contents), {"background-color": "#102A5F",
						"margin-left": "5.25rem", "margin-bottom": "1.25rem", "border-radius": "8px", "margin-top": "0.5rem"}
			else:
				return html.Div(id="ret_edit", children=[]), {"display":"none"}

		# FUNCTIONALITY OF THE CONFIRM COMPOUND EDIT BUTTON --------------------
		@_app.callback(
			[Output("edit_cpd_space", "children"),
			Output("submit_rxn_btn2", "style")],
			[Input("cfrm_edit_btn2", "n_clicks"),
			Input("edit_cpd_dropdown", "value")],
		)
		# opens the info panel to edit a specific compound, selected from the dropdown
		def open_ecs(clicks, edit_cpd_dropdown):
			if dash.callback_context.triggered[0]['prop_id'] == "cfrm_edit_btn2.n_clicks":
				# check if a reaction was inputted, send error message if not
				if edit_cpd_dropdown == "":
					return html.Div(id="ret_edit", children=[
						html.H5("Error: Please select a compound to edit", id="cpd_error",
							className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
						]), {"display":"none"}
				# find the compound in the model corresponding with the compound 
				# chosen in the dropdown 
				model = self._model
				c = edit_cpd_dropdown[0]
				for i in model.compounds:
					if i.id == c:
						compound = i
				contents = []
				# printing compound info 
				contents.append(
					html.H5("ID: " + compound.id, id="ecs_id_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "0.5rem"})
					)
				contents.append(
					html.H5("Name: " + compound.name, id="ecs_id_text", className="text-center",
						style={"color": "#102A5F", "text-transform": "none", "margin-top": "0.5rem"})
					)
				# dropdowns and inputs to change details for this compound
				contents.append(
					dbc.Row([
						dbc.Col([
							html.Div([
								html.H6("Charge", id="charge_text", className="text-center", style=
									{"color": "#102A5F", "margin-top": "0.25rem"}),
								dbc.Input(id="charge" ,type="number", size="sm", placeholder=compound.charge, 
									style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"}),
								]),
							]),
						dbc.Col([
							html.H6("Formula", id="charge_text", className="text-center", style=
									{"color": "#102A5F", "margin-top": "0.25rem"}),
							dbc.Input(id="formula" ,type="text", size="sm", placeholder=compound.formula, 
								style={"color": "#102A5F", "border-radius": "8px", "margin-top": "0.25rem"}),
							]),
						], id='compounds'), )
				return html.Div(id="ret_edit2", children=contents), {"background-color": "#102A5F",
						"margin-left": "5.25rem", "margin-bottom": "1.25rem", "border-radius": "8px", "margin-top": "0.5rem"}
			else:
				return html.Div(id="ret_edit2", children=[]), {"display":"none"}

		# FUNCTIONALITY OF SAVE REACTION BUTTON ON EDIT TAB --------------------
		@_app.callback(
			[Output("ret_edit", "children"),
			Output("submit_rxn_btn", "style", allow_duplicate=True)],
			[Input("submit_rxn_btn", "n_clicks"), ],
			[State("ret_edit", "children"), ],
			prevent_initial_call=True,
		)
		# makes changes to a reaction in this instance of the model from user input
		def update_model_rxn(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == 'submit_rxn_btn.n_clicks':
				# storing data inputted to the page by the user
				direction = []
				left = []
				right = []
				# datalist is used to find all data inputted by the user 
				for i in datalist:
					if isinstance(i, dict):
						if 'id' in i['props']:
							if i['props']['id'] == 'dir':
								id_temp = i['props']['children'][1]['props']['children'][0]['props']
								if id_temp['id'] == "direction":
									if 'value' in id_temp:
										direction.append(id_temp['value'])
							elif i['props']['id'] == 'left':
								left_list = []
								for j in range(3):
									temp = i['props']['children'][j]['props']['children'][0]['props']
									if 'value' in temp:
										left_list.append(temp['value'])
									elif 'placeholder' in temp:
										left_list.append(temp['placeholder'])
								left.append(left_list)
							elif i['props']['id'] == 'left_new':
								left_list = []
								for j in range(3):
									temp = i['props']['children'][j]['props']['children'][0]['props']
									if 'value' in temp:
										left_list.append(temp['value'])
								add = True
								if left_list == []:
									add = False
								for l in left_list:
									if l == 0 or l == "":
										add = False
								if add:
									left.append(left_list)
							elif i['props']['id'] == 'right':
								right_list = []
								for j in range(3):
									temp = i['props']['children'][j]['props']['children'][0]['props']
									if 'value' in temp:
										right_list.append(temp['value'])
									elif 'placeholder' in temp:
										right_list.append(temp['placeholder'])
								right.append(right_list)
							elif i['props']['id'] == 'right_new':
								right_list = []
								for j in range(3):
									temp = i['props']['children'][j]['props']['children'][0]['props']
									if 'value' in temp:
										right_list.append(temp['value'])
								add = True
								if right_list == []:
									add = False
								for r in right_list:
									if r == 0 or r == "":
										add = False
								if add:
									right.append(right_list)
				# figuring out which reaction is being edited by the user
				curr_r = datalist[0]['props']['children'].strip("ID:").strip()
				for r in self._model.reactions:
					if r.id.upper() == curr_r.upper():
						reaction = r
				# adding changes to this instance of the model 
				eq_dict = {}
				for le in left:
					eq_dict[Compound(le[1]).in_compartment(le[2])] = int(le[0]) * - 1
				for r in right:
					eq_dict[Compound(r[1]).in_compartment(r[2])] = int(r[0])
				if direction == []:
					# error message if no direction is given 
					return [html.H5("Error: Please give direction. Click confirm edit to try again.", id="error_edit", 
						className="text-center", style={"color":"#5d80c9", "text-transform":"none", "margin-top": "0.5rem"})],\
					{"display":"none"}
				if direction[0] == "forward":
					new_eq = Reaction(Direction.Forward, eq_dict)
				elif direction[0] == "reverse":
					new_eq = Reaction(Direction.Reverse, eq_dict)
				elif direction[0] == "both":
					new_eq = Reaction(Direction.Both, eq_dict)
				else:
					quit("Nonvalid direction")
				reaction.equation = new_eq
				self._model.reactions.add_entry(reaction)
			# on success, send a message to the user to confirm their edits
			return [html.H5("Edit Saved! Select confirm edit to make another.", id="error_edit", className="text-center", style={"color":"#102A5F", "text-transform":"none", 
				"margin-top": "0.5rem"})], {"display":"none"}

		# FUNCTIONALITY OF SAVE COMPOUND BUTTON ON EDIT TAB --------------------
		@_app.callback(
			[Output("ret_edit2", "children"),
			Output("submit_rxn_btn2", "style", allow_duplicate=True)],
			[Input("submit_rxn_btn2", "n_clicks"), ],
			[State("ret_edit2", "children"), ],
			prevent_initial_call=True,
		)
		# makes changes to a compound in this instance of the model from user input
		def update_model_cpd(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == 'submit_rxn_btn2.n_clicks':
				# collecting charge and formula inputted by the user 
				charge_info = datalist[-1]['props']['children'][0]['props']['children'][0]['props']['children'][1]['props']
				formula_info = datalist[-1]['props']['children'][1]['props']['children'][1]['props']
				charge = ""
				formula = ""
				if 'value' in charge_info:
				  charge = charge_info['value']
				elif 'placeholder' in charge_info:
				  charge = charge_info['placeholder']
				if 'value' in formula_info:
				  formula = formula_info['value']
				elif 'placeholder' in formula_info:
				  formula = formula_info['placeholder']
				cid = datalist[0]['props']['children'].strip("ID: ").strip()
				# making changes to this compound in the model instance 
				for c in self._model.compounds:
					if c.id.upper() == cid.upper():
						compound = c
				compound.charge = charge
				compound.formula = formula
				self._model.compounds.add_entry(compound)
			# send a confirmation message to the user on success
			return [html.H5("Edit Saved! Select confirm edit to make another.", id="error_edit", className="text-center", style={"color":"#102A5F", "text-transform":"none", 
				"margin-top": "0.5rem"})], {"display":"none"}

		# FUNCTIONALITY OF THE SAVE AND VIEW MODEL BUTTON ON THE EDIT TAB ------
		@_app.callback(
			[Output("edit_model_space", "children")],
			[Input("edit_save_btn", "n_clicks"),
			Input("edit_path_dropdown", "value")],
			prevent_initial_call=True,)
		# saves the user's reaction/compound changes and displays a network graph
		def save_model(clicks, add_path_dropdown):
			if dash.callback_context.triggered[0]['prop_id'] == 'edit_save_btn.n_clicks':
				# setting up the model files for writing
				out_mw = ModelWriter()
				path = FilePathContext(self._args.model)
				with open('{}/model.yaml'.format(path), 'r') as f:
					modelfile = yaml.safe_load(f)
				exported = []
				# saving all reactions currently in the model
				for i in modelfile['reactions']:
					out_nm = NativeModel()
					with open('{}/{}'.format(path, i['include']), 'r') as f:
						reactionfile = yaml.safe_load(f)
						for j in reactionfile:
							exported.append(j['id'])
							if j['id'] in self._model.reactions:
								out_nm.reactions.add_entry(
									self._model.reactions[j['id']])
						with open('{}/{}'.format(path, i['include']), "w") \
								as f:
							out_mw.write_reactions(f, out_nm.reactions)
				out_nm = NativeModel()
				# writing all reactions to the YAML files 
				for j in self._model.reactions:
					if j.id not in exported:
						out_nm.reactions.add_entry(self._model.reactions[j.id])
				with open('{}/{}'.format(path,
						  "reactions_curated.yaml"), "a+") as f:
					out_mw.write_compounds(f, out_nm.reactions)
				# saving all compounds currently in the model 
				exported = []
				for i in modelfile['compounds']:
					out_nm = NativeModel()
					with open('{}/{}'.format(path, i['include']), 'r') \
							as f:
						compoundfile = yaml.safe_load(f)
						for j in compoundfile:
							exported.append(j['id'])
							if j['id'] in self._model.compounds:
								out_nm.compounds.add_entry(
									self._model.compounds[j['id']])
						with open('{}/{}'.format(path,
												 i['include']), "w") as f:
							out_mw.write_compounds(f, out_nm.compounds)
				out_nm = NativeModel()
				# writing all compounds to the YAML files 
				for j in self._model.compounds:
					if j.id not in exported:
						out_nm.compounds.add_entry(self._model.compounds[j.id])
				with open('{}/{}'.format(path, "compounds_curated.yaml"),
						  "a+") as f:
					out_mw.write_compounds(f, out_nm.compounds)
				# creating the cytoscape visualization
				model = self._mm
				nm, network = read_model(self._model, model, "C")
				pathway_list, rxn_set = get_pathway_list(nm, add_path_dropdown)
				rxn_list = rxn_set
				fba_dropdown = []
				nodes.clear()
				edges.clear()
				obj_flux = build_network(self, nm, model, rxn_list, network, fba_dropdown, nodes, edges)
				c_n = cyto.Cytoscape(id='net_curate',layout={'name': 'cose',
					'animate': True,'animationDuration': 1000},style={'width': '100%', 
					'height': '37.5rem'},elements=nodes + edges,stylesheet=[{'selector': 'node',
					'style': {'background-color':'data(col)','label': 'data(label)'}},
					{"selector": "edge","style": {"width": 1,"curve-style": "straight",
					"target-arrow-shape":"triangle"}},{"selector": "[flux > 0]",
					"style": {"line-color": "blue","target-arrow-color":"blue",
					"width": 1,"curve-style": "bezier","target-arrow-shape":"triangle",
					"line-dash-pattern":[6, 3],"line-dash-offset": 24}},{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
					minZoom=0.2, maxZoom=3, responsive=True)
				return [c_n]
			else:
				# error message if no pathway is selected
				return ["Please select a pathway"]

		# FUNCTIONALITY OF THE UNBALANCED REACTIONS DROPDOWN (CHARGECHECK) -----
		@_app.callback(
			Output("cc_space", "children"),
			Input("unbalanced_rxns", "value")
		)
		# opens the charge check page when an unbalanced reaction is selected from the dropdown
		def open_cc_details(reaction):
			# function to get compound name from model 
			def compound_name(id):
				if id not in self._model.compounds:
					return id
				return self._model.compounds[id].properties.get('name', id)
			# if no reaction is provided, skip this callback
			if reaction == "" or reaction == [] or reaction == None:
				return []
			else:
				# store compound and reaction info from the model
				model = self._model
				c_list = []
				formula_dict = {}
				charge_dict = {}
				comp_set = set()
				for i in model.compounds:
					c_list.append(i.id)
					formula_dict[i.id] = i.formula
					charge_dict[i.id] = i.charge
				for i in model.reactions:
					for j in i.equation.compounds:
						comp_set.add(j[0].compartment)
				comp_list = list(comp_set)
				with HiddenPrints(): 
					unbalanced, count, unchecked, exclude, unbalanced_list, unbalance_dict = charge_check(self, 1e-6)
				contents = []
				# identify which unbalanced reaction was selected in the dropdown
				if isinstance(reaction, str):
					for i in self._model.reactions:
						if i.id == reaction:
							r = i
					# create a cytoscape network graph for the unbalanced reaction
					mmodel = self._mm
					nm, network = read_model(self._model, mmodel, "C")
					temp_nodes, temp_edges = rxn_network(nm, r.equation, reaction)
					c_n = cyto.Cytoscape(id='rxn_cc',layout={'name': 'cose','animate': True,'animationDuration': 1000},
						style={'width': '100%', 'height': '37.5rem'}, elements=temp_nodes + temp_edges,
						stylesheet=[{'selector': 'node',
						'style': {'background-color': '#102A5F', 'label': 'data(label)'}}, {"selector": "edge",
						"style": {"width": 1, "curve-style": "straight"}},
						{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
						minZoom=0.2, maxZoom=3, responsive=True)
					# display reaction, equation, and charge balance for the selected unbalanced reaction
					contents.append(html.H4("REACTION: {}".format(r.id), className="text-center",
						style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H5("EQUATION: ", className="text-center",
						style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H6("{}".format(r.equation), className="text-center",
						style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none",
						"padding":"0.5rem"}))
					contents.append(html.H5("TRANSLATED EQUATION:", className="text-center",
						style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H6("{}".format(r.equation.translated_compounds(compound_name)),
						className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
						"text-transform": "none", "padding":"0.5rem"}))
					try:
						contents.append(html.H5("Charge off by {}".format(unbalance_dict[r.id]),
							className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", 
							"text-transform": "none"}))
					except KeyError:
						contents.append(html.H5("Charge successfully balanced",
							className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", 
							"text-transform": "none"}))
				# display all compounds on the left side of the equation
				# allow the user to edit the ID, stoich, formula, and component of each compound
				left_contents = []
				left_contents.append(html.H4("Left:", className="text-center",
					style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
				left_contents.append(
					dbc.Row([
						dbc.Col([
							html.H6("ID", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Stoich", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Formula", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Q", style={
								"color": "#102A5F", "text-transform": "none"})
							], width=1),
						dbc.Col([
							html.H6("Comp", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=2),
					]))
				cur_id = 0
				for i in r.equation.left:
					left_contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Dropdown(
									id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
									placeholder=i[0].name, multi=False, style={"color": "#102A5F", "height": "2.25rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"},
									)], width=3),
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder=str(i[1]),
									style={'width': '80%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								dbc.Input(id="formula{}".format(str(cur_id)), type="text", placeholder=str(formula_dict[i[0].name]),
									style={'width': '90%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								html.H5("{}".format(str(charge_dict[i[0].name])), id="charge{}".format(str(cur_id)),
									className="text-center", style={"color": "#102A5F", "border-radius": "8px", 
									"margin-bottom": "0.25rem", "margin-top": "0.5rem"})
								], width=1),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
									"value": x, } for x in list(comp_list)], placeholder=i[0].compartment,
									multi=False, style={"color": "#102A5F", "height": "2.25rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"})
								], width=2),
						], justify="center"),
					)
					cur_id += 1
				# display dropdowns and inputs for the user to add their own compound to the left
				# side of the reaction
				left_contents.append(html.H4("Add to Left:", className="text-center",
					style={"margin-top": "1.1rem", "color": "#102A5F", "text-transform": "none"}))
				left_contents.append(
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Div([
									dcc.Dropdown(
										id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
										placeholder="COMPOUND ID", multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"},)
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder="STOICH",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							]),
						dbc.Row([
							dbc.Col([
								html.Div([
									dbc.Input(id="formula{}".format(str(cur_id)), type="text", 
										placeholder="FORMULA", style={'width': '90%', "color": "#102A5F", 
										"border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="charge{}".format(str(cur_id)), type="number", placeholder="CHARGE",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
										"value": x, } for x in list(comp_list)], placeholder="COMPARTMENT",
										multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"})
									]),
								]),
							]),
						]))
				right_contents = []
				# display all compounds on the right side of the equation
				# allow the user to edit the ID, stoich, formula, and component of each compound
				right_contents.append(html.H4("Right:", className="text-center",
					style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
				right_contents.append(
					dbc.Row([
						dbc.Col([
							html.H6("ID", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Stoich", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Formula", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Q", style={
								"color": "#102A5F", "text-transform": "none"})
							], width=1),
						dbc.Col([
							html.H6("Comp", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=2),
					]))
				cur_id = 0
				for i in r.equation.right:
					right_contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Dropdown(
									id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
									placeholder=i[0].name, multi=False, style={"color": "#102A5F", "height": "2.25rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"},
									)], width=3),
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder=str(i[1]),
									style={'width': '80%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								dbc.Input(id="formula{}".format(str(cur_id)), type="text", placeholder=str(formula_dict[i[0].name]),
									style={'width': '90%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								html.H5("{}".format(str(charge_dict[i[0].name])), id="charge{}".format(str(cur_id)),
									className="text-center", style={"color": "#102A5F", "border-radius": "8px", 
									"margin-bottom": "0.25rem", "margin-top": "0.5rem"})
								], width=1),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
									"value": x, } for x in list(comp_list)], placeholder=i[0].compartment,
									multi=False, style={"color": "#102A5F", "height": "2.25rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"})
								], width=2),
						], justify="center"),
					)
					cur_id += 1
				# display inputs and dropdowns to allow the user to add their own compound
				# to the right side of the reaction 
				right_contents.append(html.H4("Add to Right:", className="text-center",
					style={"margin-top": "1.1rem", "color": "#102A5F", "text-transform": "none"}))
				right_contents.append(
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Div([
									dcc.Dropdown(
										id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
										placeholder="COMPOUND ID", multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"},)
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder="STOICH",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							]),
						dbc.Row([
							dbc.Col([
								html.Div([
									dbc.Input(id="formula{}".format(str(cur_id)), type="text", 
										placeholder="FORMULA", style={'width': '90%', "color": "#102A5F", 
										"border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="charge{}".format(str(cur_id)), type="number", placeholder="CHARGE",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
										"value": x, } for x in list(comp_list)], placeholder="COMPARTMENT",
										multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"})
									]),
								]),
							]),
						]))
				# return all of the contents for the chargecheck page
				# CALLBACK FOR SAVE CHANGES TO MODEL BUTTON ON LINE 3224
				return html.Div([
					html.Div(id="cc_info", children=contents, style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
					html.Div(id="cc_model", children=[c_n], style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
					dbc.Row([
						html.Div(id="cc_left", children=left_contents, style={"background-color": "#E6EBEE", "width": "46%", "border-radius": "8px",
							"margin-top":"0.5rem", "padding": "0.5rem", "margin-left": "0.5rem"}),
						html.Div(id="cc_right", children=right_contents, style={"background-color": "#E6EBEE", "width": "46%", "border-radius": "8px",
							"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
						], style={"margin-left": "1rem"}),
					html.Div(id="cc_save", children=[
						dbc.Button("Save Changes to Model", id="cc_save_btn", style={
							"background-color": "#102A5F", "border-radius": "8px", "margin-left": "8rem",
							"padding": "1.25rem"})
						], style={"background-color": "#E6EBEE", "width": "50%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "15rem", "padding": "0.5rem"}),
					])

		# FUNCTIONALITY OF THE DOWNLOAD LIST AND DATA BUTTON (CHARGECHECK) -----
		@_app.callback(
			Output("down_list_cc", "data"),
			Input("download_unbal", "n_clicks"),
		)
		# creates a csv file containing the unbalanced reactions and compound info
		def create_csv(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "download_unbal.n_clicks":
				# function to get compound name from model 
				def compound_name(id):
					if id not in self._model.compounds:
						return id
					return self._model.compounds[id].properties.get('name', id)
				# storing necessary info on compounds and reactions from model
				with HiddenPrints(): 
					unbalanced, count, unchecked, exclude, unbalanced_list, unbalance_dict = charge_check(self, 1e-6)
				model = self._model
				c_list = []
				formula_dict = {}
				charge_dict = {}
				comp_set = set()
				for i in model.compounds:
					c_list.append(i.id)
					formula_dict[i.id] = i.formula
					charge_dict[i.id] = i.charge
				for i in model.reactions:
					for j in i.equation.compounds:
						comp_set.add(j[0].compartment)
				comp_list = list(comp_set)
				data = {"Compound ID": [], "Stoich": [], "Formula": [], "Charge": [], "Component": [], \
					"Reaction Involved with":[], "Equation": [], "Translated Equation": [], "Reaction Charge off by": []}
				# filling the data dictionary with info on each compound
				# stores ID, stoich, formula, charge, component, reactions involved with, equation,
				# translated equation, and reaction charge imbalance 
				for rxn in unbalanced_list:
					if isinstance(rxn, str):
						for i in self._model.reactions:
							if i.id == rxn:
								r = i
					# compounds on the left side of the reaction
					for cpd in r.equation.left:
						data["Compound ID"].append(cpd[0].name)
						data["Stoich"].append(str(cpd[1]))
						data["Formula"].append(str(formula_dict[cpd[0].name]))
						data["Charge"].append(str(charge_dict[cpd[0].name]))
						data["Component"].append(cpd[0].compartment)
						data["Reaction Involved with"].append(r.id)
						data["Equation"].append(r.equation)
						data["Translated Equation"].append(r.equation.translated_compounds(compound_name))
						if r.id in unbalance_dict.keys():
							data["Reaction Charge off by"].append(unbalance_dict[r.id])
						else:
							data["Reaction Charge off by"].append("balanced charge")
					# compounds on the right side of the reaction 
					for cpd in r.equation.right:
						data["Compound ID"].append(cpd[0].name)
						data["Stoich"].append(str(cpd[1]))
						data["Formula"].append(str(formula_dict[cpd[0].name]))
						data["Charge"].append(str(charge_dict[cpd[0].name]))
						data["Component"].append(cpd[0].compartment)
						data["Reaction Involved with"].append(r.id)
						data["Equation"].append(r.equation)
						data["Translated Equation"].append(r.equation.translated_compounds(compound_name))
						if r.id in unbalance_dict.keys():
							data["Reaction Charge off by"].append(unbalance_dict[r.id])
						else:
							data["Reaction Charge off by"].append("balanced charge")
				# creating a dataframe object with the compound dictionary 
				dataframe = pd.DataFrame(data)
				# returning a downloadable object to the user 
				return dcc.send_data_frame(dataframe.to_csv, "ChargeCheckData.csv")

		# FUNCTIONALITY OF THE CHARGECHECK SAVE BUTTON -------------------------
		@_app.callback(
			Output("cc_save", "children"),
			Input("cc_save_btn", "n_clicks"),
			State("cc_left", "children"),
			State("cc_right", "children"),
			)
		# parses through datalist to add compound/reaction info to the model instance
		def save_cc(nclicks, left_data, right_data):
			if dash.callback_context.triggered[0]['prop_id'] == "cc_save_btn.n_clicks":
				eq_dict = dict()
				for i in range(2, len(left_data)-2):
					# getting id data from dropdowns
					if 'value' in (left_data[i]['props']['children'][0]['props']['children'][0]['props']).keys():
						id_ = left_data[i]['props']['children'][0]['props']['children'][0]['props']['value']
					else:
						id_ = left_data[i]['props']['children'][0]['props']['children'][0]['props']['placeholder']
					# getting stoich from input 
					if 'value' in (left_data[i]['props']['children'][1]['props']['children'][0]['props']).keys():
						stoich = float(left_data[i]['props']['children'][1]['props']['children'][0]['props']['value'])
					else:
						stoich = float(left_data[i]['props']['children'][1]['props']['children'][0]['props']['placeholder'])
					# getting formula from input
					if 'value' in (left_data[i]['props']['children'][2]['props']['children'][0]['props']).keys():
						formula = left_data[i]['props']['children'][2]['props']['children'][0]['props']['value']
					else:
						formula = left_data[i]['props']['children'][2]['props']['children'][0]['props']['placeholder']
					# getting charge from text
					charge = int(left_data[i]['props']['children'][3]['props']['children'][0]['props']['children'])
					# getting component from dropdown
					if 'value' in (left_data[i]['props']['children'][4]['props']['children'][0]['props']).keys():
						comp = left_data[i]['props']['children'][4]['props']['children'][0]['props']['value']
					else:
						comp = left_data[i]['props']['children'][4]['props']['children'][0]['props']['placeholder'] 
					# adding compound to database and reaction information dictionary
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				last_row = left_data[-1]
				all_found = True
				# getting id data from dropdowns
				if 'value' in (last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					id_ = last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting stoich from input 
				if 'value' in (last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					stoich = float(last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting formula from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					formula = last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting charge from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					charge = int(last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting component from dropdown
				if 'value' in (last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]).keys():
					comp = last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# adding compound to database and reaction information dictionary
				if all_found:
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				else:
					return html.H5("Error: Please fill in all fields to add a compound to the left or right", 
						id="cc_save_confirm", className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
				# right side 
				for i in range(2, len(right_data)-2):
					# getting id data from dropdowns
					if 'value' in (right_data[i]['props']['children'][0]['props']['children'][0]['props']).keys():
						id_ = right_data[i]['props']['children'][0]['props']['children'][0]['props']['value']
					else:
						id_ = right_data[i]['props']['children'][0]['props']['children'][0]['props']['placeholder']
					# getting stoich from input 
					if 'value' in (right_data[i]['props']['children'][1]['props']['children'][0]['props']).keys():
						stoich = float(right_data[i]['props']['children'][1]['props']['children'][0]['props']['value'])
					else:
						stoich = float(right_data[i]['props']['children'][1]['props']['children'][0]['props']['placeholder'])
					# getting formula from input
					if 'value' in (right_data[i]['props']['children'][2]['props']['children'][0]['props']).keys():
						formula = right_data[i]['props']['children'][2]['props']['children'][0]['props']['value']
					else:
						formula = right_data[i]['props']['children'][2]['props']['children'][0]['props']['placeholder']
					# getting charge from text
					charge = int(right_data[i]['props']['children'][3]['props']['children'][0]['props']['children'])
					# getting component from dropdown
					if 'value' in (right_data[i]['props']['children'][4]['props']['children'][0]['props']).keys():
						comp = right_data[i]['props']['children'][4]['props']['children'][0]['props']['value']
					else:
						comp = right_data[i]['props']['children'][4]['props']['children'][0]['props']['placeholder'] 
					# adding compound to database and reaction information dictionary
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				last_row = right_data[-1]
				all_found = True
				# getting id data from dropdowns
				if 'value' in (last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					id_ = last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting stoich from input 
				if 'value' in (last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					stoich = float(last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting formula from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					formula = last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting charge from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					charge = int(last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting component from dropdown
				if 'value' in (last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]).keys():
					comp = last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# adding compound to database and reaction information dictionary
				if all_found:
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				else:
					return html.H5("Error: Please fill in all fields to add a compound to the left or right", 
						id="cc_save_confirm", className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
				# adding to reaction dictionary
				for r in self._model.reactions:
					if r.id.upper() == rxn.upper():
						reaction = r
						new_eq = Reaction(reaction.equation.__dict__['_direction'], eq_dict)
						reaction.equation = new_eq
						self._model.reactions.add_entry(reaction)
				# return a save message
				return html.H5("Model saved!", id="cc_save_confirm",
					className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
			else:
				return dbc.Button("Save Changes to Model", id="cc_save_btn", className="text-center", style={
					"background-color": "#102A5F", "border-radius": "8px", "margin-left": "7rem",
					"padding": "1.25rem"})

		# FUNCTIONALITY OF THE UNBALANCED REACTIONS DROPDOWN (FORMULACHECK) ----
		@_app.callback(
			Output("fc_space", "children"),
			Input("unbalanced_rxns_fc", "value")
		)
		# opens the formula check tab when an unbalanced reaction is selected
		def open_fc_details(reaction):
			# function to get the compound name from the model
			def compound_name(id):
				if id not in self._model.compounds:
					return id
				return self._model.compounds[id].properties.get('name', id)
			# if no reaction is provided, skip this callback 
			if reaction == "" or reaction == [] or reaction == None:
				return []
			else:
				# store compound and reaction info from the model
				model = self._model
				c_list = []
				formula_dict = {}
				charge_dict = {}
				comp_set = set()
				for i in model.compounds:
					c_list.append(i.id)
					formula_dict[i.id] = i.formula
					charge_dict[i.id] = i.charge
				for i in model.reactions:
					for j in i.equation.compounds:
						comp_set.add(j[0].compartment)
				comp_list = list(comp_set)
				with HiddenPrints():
					unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(self)
				contents = []
				# identify which unbalanced reaction was selected in the dropdown
				if isinstance(reaction, str):
					for i in self._model.reactions:
						if i.id == reaction:
							r = i
					# create a cytoscape network graph for the unbalanced reaction
					mmodel = self._mm
					nm, network = read_model(self._model, mmodel, "C")
					temp_nodes, temp_edges = rxn_network(nm, r.equation, reaction)
					c_n = cyto.Cytoscape(id='rxn_cc',layout={'name': 'cose','animate': True,'animationDuration': 1000},
						style={'width': '100%', 'height': '37.5rem'}, elements=temp_nodes + temp_edges,
						stylesheet=[{'selector': 'node',
						'style': {'background-color': '#102A5F', 'label': 'data(label)'}}, {"selector": "edge",
						"style": {"width": 1, "curve-style": "straight"}},
						{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
						minZoom=0.2, maxZoom=3, responsive=True)
					# display reaction, equation, and stioch balance for reaction selected
					contents.append(html.H4("REACTION: {}".format(r.id), className="text-center",
						style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H5("EQUATION: ", className="text-center",
						style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H6("{}".format(r.equation), className="text-center",
						style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none",
						"padding":"0.5rem"}))
					contents.append(html.H5("TRANSLATED EQUATION:", className="text-center",
						style={"margin-top": "1.25rem", "color": "#102A5F", "text-transform": "none"}))
					contents.append(html.H6("{}".format(r.equation.translated_compounds(compound_name)),
						className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
						"text-transform": "none", "padding":"0.5rem"}))
					try:
						contents.append(
							dbc.Row([
								dbc.Col([
									html.H5("Left side total: {}".format(form_dict[r.id][0]), id="fc_lst_txt", 
										className="text-center", style={"color": "#5d80c9", "text-transform": "none"}),
									html.H5("Left side off by: {}".format(form_dict[r.id][1]), id="fc_lst_txt", 
										className="text-center", style={"color": "#5d80c9", "text-transform": "none"}),
								]),
								dbc.Col([
									html.H5("Right side total: {}".format(form_dict[r.id][2]), id="fc_lst_txt", 
										className="text-center", style={"color": "#5d80c9", "text-transform": "none"}),
									html.H5("Right side off by: {}".format(form_dict[r.id][3]), id="fc_lst_txt", 
										className="text-center", style={"color": "#5d80c9", "text-transform": "none"}),
								]),
							]),

						)
					except KeyError:
						contents.append(html.H5("Reaction successfully balanced",
							className="text-center", style={"margin-top": "1.25rem", "color": "#102A5F", 
							"text-transform": "none"}))
				# display all compounds on the left side of the equation 
				# allow the user to edit the ID, stoich, formula, and component of each compound
				left_contents = []
				left_contents.append(html.H4("Left:", className="text-center",
					style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
				left_contents.append(
					dbc.Row([
						dbc.Col([
							html.H6("ID", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Stoich", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Formula", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Q", style={
								"color": "#102A5F", "text-transform": "none"})
							], width=1),
						dbc.Col([
							html.H6("Comp", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=2),
					]))
				cur_id = 0
				for i in r.equation.left:
					left_contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Dropdown(
									id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
									placeholder=i[0].name, multi=False, style={"color": "#102A5F", "height": "3rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"},
									)], width=3),
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder=str(i[1]),
									style={'width': '80%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								dbc.Input(id="formula{}".format(str(cur_id)), type="text", placeholder=str(formula_dict[i[0].name]),
									style={'width': '90%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								html.H5("{}".format(str(charge_dict[i[0].name])), id="charge{}".format(str(cur_id)),
									className="text-center", style={"color": "#102A5F", "border-radius": "8px", 
									"margin-bottom": "0.25rem", "margin-top": "0.5rem"})
								], width=1),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
									"value": x, } for x in list(comp_list)], placeholder=i[0].compartment,
									multi=False, style={"color": "#102A5F", "height": "3rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"})
								], width=2),
						], justify="center"),
					)
					cur_id += 1
				# display dropdowns and inputs for the user to add their own compound to
				# the left side of the reaction
				left_contents.append(html.H4("Add to Left:", className="text-center",
					style={"margin-top": "1.1rem", "color": "#102A5F", "text-transform": "none"}))
				left_contents.append(
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Div([
									dcc.Dropdown(
										id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
										placeholder="COMPOUND ID", multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"},)
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder="STOICH",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							]),
						dbc.Row([
							dbc.Col([
								html.Div([
									dbc.Input(id="formula{}".format(str(cur_id)), type="text", 
										placeholder="FORMULA", style={'width': '90%', "color": "#102A5F", 
										"border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="charge{}".format(str(cur_id)), type="number", placeholder="CHARGE",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
										"value": x, } for x in list(comp_list)], placeholder="COMPARTMENT",
										multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"})
									]),
								]),
							]),
						]))
				right_contents = []
				# display all compounds on the right side of the equation
				# allow the user to edit the ID, stoich, formula, and component of each compound
				right_contents.append(html.H4("Right:", className="text-center",
					style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}))
				right_contents.append(
					dbc.Row([
						dbc.Col([
							html.H6("ID", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Stoich", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Formula", style={"margin-left":"1rem",
								"color": "#102A5F", "text-transform": "none"})
							], width=3),
						dbc.Col([
							html.H6("Q", style={
								"color": "#102A5F", "text-transform": "none"})
							], width=1),
						dbc.Col([
							html.H6("Comp", className="text-center", style={ 
								"color": "#102A5F", "text-transform": "none"})
							], width=2),
					]))
				cur_id = 0
				for i in r.equation.right:
					right_contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Dropdown(
									id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
									placeholder=i[0].name, multi=False, style={"color": "#102A5F", "height": "3rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"},
									)], width=3),
							dbc.Col([
								dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder=str(i[1]),
									style={'width': '80%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								dbc.Input(id="formula{}".format(str(cur_id)), type="text", placeholder=str(formula_dict[i[0].name]),
									style={'width': '90%', "color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
								], width=3),
							dbc.Col([
								html.H5("{}".format(str(charge_dict[i[0].name])), id="charge{}".format(str(cur_id)),
									className="text-center", style={"color": "#102A5F", "border-radius": "8px", 
									"margin-bottom": "0.25rem", "margin-top": "0.5rem"})
								], width=1),
							dbc.Col([
								dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
									"value": x, } for x in list(comp_list)], placeholder=i[0].compartment,
									multi=False, style={"color": "#102A5F", "height": "3rem",
									"border-radius": "8px", "margin-bottom": "0.5rem"})
								], width=2),
						], justify="center"),
					)
					cur_id += 1
				# display inputs and dropdowns to allow the user to add their own compound
				# to the right side of the reaction 
				right_contents.append(html.H4("Add to Right:", className="text-center",
					style={"margin-top": "1.1rem", "color": "#102A5F", "text-transform": "none"}))
				right_contents.append(
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Div([
									dcc.Dropdown(
										id="drop{}".format(str(cur_id)), options=[{"label": x,"value": x, } for x in list(c_list)],
										placeholder="COMPOUND ID", multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"},)
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="stoich{}".format(str(cur_id)), type="number", placeholder="STOICH",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							]),
						dbc.Row([
							dbc.Col([
								html.Div([
									dbc.Input(id="formula{}".format(str(cur_id)), type="text", 
										placeholder="FORMULA", style={'width': '90%', "color": "#102A5F", 
										"border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dbc.Input(id="charge{}".format(str(cur_id)), type="number", placeholder="CHARGE",
										style={"color": "#102A5F", "border-radius": "8px", "margin-bottom": "0.25rem"})
									]),
								]),
							dbc.Col([
								html.Div([
									dcc.Dropdown(id="compartment{}".format(str(cur_id)),options=[{"label": x,
										"value": x, } for x in list(comp_list)], placeholder="COMPARTMENT",
										multi=False, style={"color": "#102A5F", "height": "3rem",
										"border-radius": "8px", "margin-bottom": "0.5rem"})
									]),
								]),
							]),
						]))
				# return all contents for the formula check page 
				# CALLBACK FOR SAVE CHANGES TO MODEL BUTTON ON LINE 3794
				return html.Div([
					html.Div(id="fc_info", children=contents, style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
					html.Div(id="fc_model", children=[c_n], style={"background-color": "#E6EBEE", "width": "95%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
					dbc.Row([
						html.Div(id="fc_left", children=left_contents, style={"background-color": "#E6EBEE", "width": "46%", "border-radius": "8px",
							"margin-top":"0.5rem", "padding": "0.5rem", "margin-left": "0.5rem"}),
						html.Div(id="fc_right", children=right_contents, style={"background-color": "#E6EBEE", "width": "46%", "border-radius": "8px",
							"margin-top":"0.5rem","margin-left": "1rem", "padding": "0.5rem"}),
						], style={"margin-left": "1rem"}),
					html.Div(id="fc_save", children=[
						dbc.Button("Save Changes to Model", id="fc_save_btn", className="text-center", style={
							"background-color": "#102A5F", "border-radius": "8px", "margin-left": "7rem",
							"padding": "1.25rem"})
						], style={"background-color": "#E6EBEE", "width": "50%", "border-radius": "8px",
						"margin-top":"0.5rem","margin-left": "15rem", "padding": "0.5rem"}),
					])

		# FUNCTIONALITY OF THE DOWNLOAD LITS AND DATA BUTTON (FORMULA CHECK) ---
		@_app.callback(
			Output("down_list_fc", "data"),
			Input("download_unbal_fc", "n_clicks"),
		)
		# creates a csv file containing the unbalanced reaction and compound info
		def create_csv(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "download_unbal_fc.n_clicks":
				# function to get compound name from the model 
				def compound_name(id):
					if id not in self._model.compounds:
						return id
					return self._model.compounds[id].properties.get('name', id)
				# storing necessary info on the compounds and reactions from the model
				with HiddenPrints():
					unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(self)
				model = self._model
				c_list = []
				formula_dict = {}
				charge_dict = {}
				comp_set = set()
				for i in model.compounds:
					c_list.append(i.id)
					formula_dict[i.id] = i.formula
					charge_dict[i.id] = i.charge
				for i in model.reactions:
					for j in i.equation.compounds:
						comp_set.add(j[0].compartment)
				comp_list = list(comp_set)
				data = {"Compound ID": [], "Stoich": [], "Formula": [], "Charge": [], "Component": [], \
					"Reaction Involved with":[], "Equation": [], "Translated Equation": [], \
					"Left Total": [], "Left Side off by": [], "Right Total": [], "Right Side off by": []}
				# filling the data fictionary with info on each compound
				# stores ID, stoich, formula, charge, component, reactions involved in, 
				# equations, translated equations, and stoich imbalance for each side of the equation
				for rxn in form_list:
					if isinstance(rxn, str):
						for i in self._model.reactions:
							if i.id == rxn:
								r = i
					# compounds on the left side of the reaction 
					for cpd in r.equation.left:
						data["Compound ID"].append(cpd[0].name)
						data["Stoich"].append(str(cpd[1]))
						data["Formula"].append(str(formula_dict[cpd[0].name]))
						data["Charge"].append(str(charge_dict[cpd[0].name]))
						data["Component"].append(cpd[0].compartment)
						data["Reaction Involved with"].append(r.id)
						data["Equation"].append(r.equation)
						data["Translated Equation"].append(r.equation.translated_compounds(compound_name))
						if r.id in form_dict.keys():
							data["Left Total"].append(form_dict[r.id][0])
							data["Left Side off by"].append(form_dict[r.id][1])
							data["Right Total"].append(form_dict[r.id][2])
							data["Right Side off by"].append(form_dict[r.id][3])
						else:
							data["Left Total"].append("balanced")
							data["Left Side off by"].append("balanced")
							data["Right Total"].append("balanced")
							data["Right side off by"].append("balanced")
					# compounds on the right side of the reaction
					for cpd in r.equation.right:
						data["Compound ID"].append(cpd[0].name)
						data["Stoich"].append(str(cpd[1]))
						data["Formula"].append(str(formula_dict[cpd[0].name]))
						data["Charge"].append(str(charge_dict[cpd[0].name]))
						data["Component"].append(cpd[0].compartment)
						data["Reaction Involved with"].append(r.id)
						data["Equation"].append(r.equation)
						data["Translated Equation"].append(r.equation.translated_compounds(compound_name))
						if r.id in form_dict.keys():
							data["Left Total"].append(form_dict[r.id][0])
							data["Left Side off by"].append(form_dict[r.id][1])
							data["Right Total"].append(form_dict[r.id][2])
							data["Right Side off by"].append(form_dict[r.id][3])
						else:
							data["Left Total"].append("balanced")
							data["Left Side off by"].append("balanced")
							data["Right Total"].append("balanced")
							data["Right side off by"].append("balanced")
				# creating a dataframe object with the compound dictionary
				dataframe = pd.DataFrame(data)
				# returning a downloadable object to the user
				return dcc.send_data_frame(dataframe.to_csv, "FormulaCheckData.csv")

		# FUNCTIONALITY OF THE FORMULA CHECK SAVE BUTTON -----------------------
		@_app.callback(
			Output("fc_save", "children"),
			Input("fc_save_btn", "n_clicks"),
			State("fc_left", "children"),
			State("fc_right", "children"),
			)
		# parses through datalist to add compound/reaction info to the model instance
		def save_cc(nclicks, left_data, right_data):
			if dash.callback_context.triggered[0]['prop_id'] == "fc_save_btn.n_clicks":
				eq_dict = dict()
				for i in range(2, len(left_data)-2):
					# getting id data from dropdowns
					if 'value' in (left_data[i]['props']['children'][0]['props']['children'][0]['props']).keys():
						id_ = left_data[i]['props']['children'][0]['props']['children'][0]['props']['value']
					else:
						id_ = left_data[i]['props']['children'][0]['props']['children'][0]['props']['placeholder']
					# getting stoich from input 
					if 'value' in (left_data[i]['props']['children'][1]['props']['children'][0]['props']).keys():
						stoich = float(left_data[i]['props']['children'][1]['props']['children'][0]['props']['value'])
					else:
						stoich = float(left_data[i]['props']['children'][1]['props']['children'][0]['props']['placeholder'])
					# getting formula from input
					if 'value' in (left_data[i]['props']['children'][2]['props']['children'][0]['props']).keys():
						formula = left_data[i]['props']['children'][2]['props']['children'][0]['props']['value']
					else:
						formula = left_data[i]['props']['children'][2]['props']['children'][0]['props']['placeholder']
					# getting charge from text
					charge = int(left_data[i]['props']['children'][3]['props']['children'][0]['props']['children'])
					# getting component from dropdown
					if 'value' in (left_data[i]['props']['children'][4]['props']['children'][0]['props']).keys():
						comp = left_data[i]['props']['children'][4]['props']['children'][0]['props']['value']
					else:
						comp = left_data[i]['props']['children'][4]['props']['children'][0]['props']['placeholder'] 
					# adding compound to database and reaction information dictionary
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				last_row = left_data[-1]
				all_found = True
				# getting id data from dropdowns
				if 'value' in (last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					id_ = last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting stoich from input 
				if 'value' in (last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					stoich = float(last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting formula from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					formula = last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting charge from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					charge = int(last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting component from dropdown
				if 'value' in (last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]).keys():
					comp = last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# adding compound to database and reaction information dictionary
				if all_found:
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				else:
					return html.H5("Error: Please fill in all fields to add a compound to the left or right", 
						id="cc_save_confirm", className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
				# right side 
				for i in range(2, len(right_data)-2):
					# getting id data from dropdowns
					if 'value' in (right_data[i]['props']['children'][0]['props']['children'][0]['props']).keys():
						id_ = right_data[i]['props']['children'][0]['props']['children'][0]['props']['value']
					else:
						id_ = right_data[i]['props']['children'][0]['props']['children'][0]['props']['placeholder']
					# getting stoich from input 
					if 'value' in (right_data[i]['props']['children'][1]['props']['children'][0]['props']).keys():
						stoich = float(right_data[i]['props']['children'][1]['props']['children'][0]['props']['value'])
					else:
						stoich = float(right_data[i]['props']['children'][1]['props']['children'][0]['props']['placeholder'])
					# getting formula from input
					if 'value' in (right_data[i]['props']['children'][2]['props']['children'][0]['props']).keys():
						formula = right_data[i]['props']['children'][2]['props']['children'][0]['props']['value']
					else:
						formula = right_data[i]['props']['children'][2]['props']['children'][0]['props']['placeholder']
					# getting charge from text
					charge = int(right_data[i]['props']['children'][3]['props']['children'][0]['props']['children'])
					# getting component from dropdown
					if 'value' in (right_data[i]['props']['children'][4]['props']['children'][0]['props']).keys():
						comp = right_data[i]['props']['children'][4]['props']['children'][0]['props']['value']
					else:
						comp = right_data[i]['props']['children'][4]['props']['children'][0]['props']['placeholder'] 
					# adding compound to database and reaction information dictionary
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				last_row = right_data[-1]
				all_found = True
				# getting id data from dropdowns
				if 'value' in (last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					id_ = last_row['props']['children'][0]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting stoich from input 
				if 'value' in (last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					stoich = float(last_row['props']['children'][0]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting formula from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]).keys():
					formula = last_row['props']['children'][1]['props']['children'][0]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# getting charge from input
				if 'value' in (last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]).keys():
					charge = int(last_row['props']['children'][1]['props']['children'][1]['props']['children'][0]['props']['children'][0]['value'])
				else:
					all_found = False
				# getting component from dropdown
				if 'value' in (last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]).keys():
					comp = last_row['props']['children'][1]['props']['children'][2]['props']['children'][0]['props']['children'][0]['value']
				else:
					all_found = False
				# adding compound to database and reaction information dictionary
				if all_found:
					for c in self._model.compounds:
						if c.id.upper() == id_.upper():
							compound = c
							compound.charge = charge
							compound.formula = formula
							self._model.compounds.add_entry(compound)
							eq_dict[Compound(id_).in_compartment(comp)] = stoich * -1
				else:
					return html.H5("Error: Please fill in all fields to add a compound to the left or right", 
						id="cc_save_confirm", className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
				# adding to reaction dictionary
				for r in self._model.reactions:
					if r.id.upper() == rxn.upper():
						reaction = r
						new_eq = Reaction(reaction.equation.__dict__['_direction'], eq_dict)
						reaction.equation = new_eq
						self._model.reactions.add_entry(reaction)
				# return a save message
				return html.H5("Model saved!", id="fc_save_confirm",
					className="text-center", style={"color": "#5d80c9", "text-transform": "none"})
			else:
				return dbc.Button("Save Changes to Model", id="fc_save_btn", className="text-center", style={
					"background-color": "#102A5F", "border-radius": "8px", "margin-left": "7rem",
					"padding": "1.25rem"})
		
		# FUNCTIONALITY OF MODEL WIDESCREEN DROPDOWN ---------------------------
		@_app.callback(
			Output("ws_model_space", "children"),
			Input("ws_path_dropdown", "value"),
			Input("ws_sub_btn", "n_clicks"),
			prevent_initial_call=True,
		)
		# display cytoscape network graph when a pathway is selected 
		def display_model(pathway, nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == "ws_sub_btn.n_clicks":
				if pathway == "":
					# creating a default cytoscape
					c_n = cyto.Cytoscape(id='net_curate', layout={'name': 'cose',
						'animate': True,'animationDuration': 1000}, style={'width': '100%', 
						'height': '37.5rem'},elements=nodes + edges, stylesheet=[{'selector': 'node',
						'style': {'background-color':'data(col)','label': 'data(label)'}},
						{"selector": "edge","style": {"width": 1,"curve-style": "straight",
						"target-arrow-shape":"triangle"}},{"selector": "[flux > 0]",
						"style": {"line-color": "blue","target-arrow-color":"blue",
						"width": 1,"curve-style": "bezier","target-arrow-shape":"triangle",
						"line-dash-pattern":[6, 3],"line-dash-offset": 24}},
						{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
						minZoom=0.2, maxZoom=3, responsive=True)
					return [c_n]
				else:
					# if a pathway was selected in the dropdown, create a cytoscape to display it
					temp_nodes = []
					temp_edges = []
					contents = []
					model = self._mm
					nm, network = read_model(self._model, model, "C")
					pathway_list, rxn_set = get_pathway_list(nm, pathway)
					rxn_list = rxn_set
					fba_dropdown = []
					build_network(self, nm, model, rxn_list, network, fba_dropdown, temp_nodes, temp_edges)
					c_n = cyto.Cytoscape(id='net_curate', layout={'name': 'cose', 'animate': True,
						'animationDuration': 1000}, style={'width': '100%', 'height': '37.5rem'},
						elements=temp_nodes + temp_edges, stylesheet=[{'selector': 'node',
						'style': {'background-color': 'data(col)', 'label': 'data(label)'}},
						{"selector": "edge", "style": {"width": 1, "curve-style": "straight",
						"target-arrow-shape": "triangle"}}, {"selector": "[flux > 0]", "style": 
						{"line-color": "blue","target-arrow-color":"blue","width": 1,"curve-style": "bezier",
						"target-arrow-shape":"triangle","line-dash-pattern":[6, 3],"line-dash-offset": 24}},
						{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
						minZoom=0.2, maxZoom=3, responsive=True)
					return [c_n]
			else:
				# default cytoscape if something goes wrong
				c_n = cyto.Cytoscape(id='net_curate',layout={'name': 'cose','animate': True,
					'animationDuration': 1000},style={'width': '100%', 'height': '37.5rem'},
					elements=nodes + edges,stylesheet=[{'selector': 'node','style': 
					{'background-color':'data(col)','label': 'data(label)'}},{"selector": "edge",
					"style": {"width": 1,"curve-style": "straight","target-arrow-shape":"triangle"}},
					{"selector": "[flux > 0]", "style": {"line-color": "blue","target-arrow-color":"blue",
					"width": 1,"curve-style": "bezier","target-arrow-shape":"triangle","line-dash-pattern":[6, 3],
					"line-dash-offset": 24}},{"selector": ".bidirectional", "style":{'source-arrow-shape': 'triangle-backcurve', 'target-arrow-shape': 'triangle-backcurve'}}],
					minZoom=0.2, maxZoom=3, responsive=True)
				return [c_n]

		# FUNCTIONALITY OF SELECT NODE TO SEE DATA ON MODEL WIDESCREEN TAB ----- 
		@_app.callback(
			Output("ws_node_data", "children"),
			[Input("net_curate", "selectedNodeData"), ],
			prevent_initial_call=True,
		)
		# shows data for node when clicked on the model widescreen tab
		def updated_selected(node):
			if dash.callback_context.triggered[0]['prop_id'] == 'net_curate.selectedNodeData' and node is not None:
				# if node has no data, skip
				if len(node) <= 0:
					return []
				# display data for node 
				contents = []
				# whether the node is a reaction or compound
				contents.append(
					dbc.Row([
						html.H4("{}".format(node[0]['type']), className="text-center", 
							style={"margin-top": "0.5rem", "color": "#102A5F", 
							"text-transform": "uppercase"}),
						]))
				# the ID and name for the node
				contents.append(
					dbc.Row([
						html.H5("ID: {} & Name: {}".format(node[0]['id'], node[0]['label']), 
							className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
							"text-transform": "none"}),
						]),
					)
				# info displayed if the node is a compound: formula and charge
				if node[0]['type'] == 'compound': 
					contents.append(
						dbc.Row([
							html.H5("Formula: {}".format(node[0]['formula']), 
								className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
								"text-transform": "none"}),
							]),
						)
					contents.append(
						dbc.Row([
							html.H5("Charge: {}".format(node[0]['charge']), 
								className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
								"text-transform": "none"}),
							]),
						)
				# info displayed if the node is a reaction: equation, pathways and genes
				elif node[0]['type'] == "reaction":
					contents.append(
						dbc.Row([
							html.H5("Equation: {}".format(node[0]['equation']), 
								className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
								"text-transform": "none"}),
							]),
						)
					contents.append(
						dbc.Row([
							html.H5("Pathways: {}".format(node[0]['pathways']), 
								className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
								"text-transform": "none"}),
							]),
						)
					contents.append(
						dbc.Row([
							html.H5("Genes: {}".format(node[0]['genes']), 
								className="text-center", style={"margin-top": "0.5rem", "color": "#102A5F", 
								"text-transform": "none"}),
							]),
						)
				return contents
			else:
				return html.H5("Click on a node to see more information", className="text-center",
					style={"margin-top": "0.5rem", "color": "#102A5F", "text-transform": "none"}),
	
		# END OF CALLBACKS TO OPERATE FUNCTIONS ON CURATE PAGE -----------------
