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
import dash_core_components as dcc
import dash_html_components as html
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


# FUNCTIONS -------------------------------------------------------------------

# Gets network from PSAMM
def read_model(model, mm, el):
	if isinstance(el, list) or el is None:
		el = "C"
	excluded_reactions = []
	for rxn in model.reactions:
		for cpd, v in rxn.equation.compounds:
			if (float(v)).is_integer() is False or float(v) > 10:
				excluded_reactions.append(rxn.id)
	network = make_network_dict(model, mm, subset=None, method='fpp',
								element=el,
								excluded_reactions=excluded_reactions,
								reaction_dict={}, analysis=None)
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
	generic = ['h2o', 'h', 'accoa', 'coa', 'nad', 'nadh', 'nadp', 'nadph',
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

# Builds the subset network, updates nodes and edges for graph
def build_network(self, nm, mm, rxn_set, network, fba_dropdown, nodes, edges):
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
	color = ["#6FB2D2", "#85C88A", "#EBD671", "#EEEEEE"]
	col_dict = defaultdict(lambda: "c")
	for col in range(0, len(comp)):
		col_dict[comp[col]] = color[col]
	for rxn in network[0]:
		if rxn.id in rxn_set:
			if rxn.name:
				rname = rxn.name
			else:
				rname = rxn.id
			rxn_num = 0
			visited = set()
			for cpd in network[0][rxn][0]:
				if cpd[0] not in visited:
					visited.add(cpd[0])
					visited.add(cpd[1])
					nodes.append({'data': {'id': str(cpd[0]), 'label': name[str(cpd[0])[0:(str(cpd[0]).find('['))]],
						'formula': formula[str(cpd[0])[0:(str(cpd[0]).find('['))]], 'compartment': cpd[0].compartment,
						'col': col_dict[cpd[0].compartment], 'charge': charge[str(cpd[0])[0:(str(cpd[0]).find('['))]],
						'orig_id': str(cpd[0])[0:(str(cpd[0]).find('['))], 'type': 'compound'}})
					nodes.append({'data': {'id': str(cpd[1]), 'label': name[str(cpd[1])[0:(str(cpd[1]).find('['))]],
						'formula': formula[str(cpd[1])[0:(str(cpd[1]).find('['))]],'compartment': cpd[0].compartment,
						'col': col_dict[cpd[0].compartment],'charge': charge[str(cpd[0])[0:(str(cpd[0]).find('['))]],
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
						nodes.append({'data': {'id': "{}_node".format(rxn.id),'label': rxn.id,'type': 'reaction',
							'equation': str(rxn.properties["equation"]),'pathways': path,'orig_id': rxn.id,'name': rxn.properties["name"],
							'genes': rxn.properties["genes"],'pathways': path_prop,'EC': rxn.properties["EC"],'flux_node': flux_carrying[rxn.id]},
							'style': {'shape': 'rectangle'}})
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
								'flux': flux}})
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
								'flux': flux}})
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
								'flux': flux}})
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
								'flux': flux
							}})
							rxn_num += 1
	return obj


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
	# INTIALIZERS --------------------------------------------------------------
	@classmethod
	def init_parser(cls, parser):
		parser.add_argument(
			'--exclude', metavar='reaction', type=convert_to_unicode,
			default=[], action=FilePrefixAppendAction,
			help='Reaction(s) to exclude from metabolite pair prediction')
		super(InteractiveCommand, cls).init_parser(parser)

	def run(self):
		# bootstrap theme CDN
		bs_theme = "https://cdn.jsdelivr.net/npm/bootswatch@5.2.3/dist/lux/bootstrap.min.css"
		self._app = dash.Dash(
			__name__, external_stylesheets=[
				bs_theme], 
				#use_pages=True
				suppress_callback_exceptions=True,
				prevent_initial_callbacks="initial_duplicate",
			)
		server = self._app.server
		if self._app is not None and hasattr(self, "callbacks"):
			self.callbacks(self._app)
		application = build_app(self)
		self._app.layout = application
		webbrowser.open_new("http://localhost:{}".format("8050"))
		self._app.run_server(debug=True)

	
	def callbacks(self, _app):


		@_app.callback(
			Output("page-content", "children"), 
			Output("rsb_btn", "style"),
			Output("open_menu_btn", "style"),
			[Input("url", "pathname")],
			Input("rsb_btn", "n_clicks"),
			Input("open_menu_btn", "n_clicks"))

		def render_page_content(pathname, rsb_clicks, open_clicks):

			# HOME PAGE --------------------------------------------------------------------------

			home = html.Img(src='')

			"""

			# DEFINITIONS -------------------------------------------------------------------------

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

			# stores the initial model stats in content_model and generate the
			# initial stats
			content_model = []
			content_model.append(html.H5("General Model Statistics:"))
			content_model.append(html.P("Reactions in model: " + str(count)))
			count = 0
			for i in self._model.compounds:
				count += 1
			content_model.append(html.P("Compounds in model: " + str(count)))
			content_model.append(
				html.P("Pathways in mode: " + str(len(pathway_list))))
			content_model.append(
				html.P(
					"Reactions with gene assocations: " +
					str(rxns_with_genes)))
			try:
				content_model.append(
					html.P("Total genes: " + str(count_genes(self._model))))
			except BaseException:
				content_model.append(html.P("Total genes: Unable to count genes"))
			with HiddenPrints():
				unbalanced, count, unchecked, exclude, unbalanced_list, \
					unbalance_dict = charge_check(self, 1e-6)
			content_model.append(html.H5("Chargecheck results:"))
			content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude))))

			with HiddenPrints():
				unbalanced, count, unchecked, exclude, form_list, \
					form_dict = formula_check(self)
			content_model.append(html.H5("Formulacheck results:"))
			content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude))))

			# initialize an empty list. the full is messy
			"""
			#nodes, edges = [], []

			sidebar_style = {
				#"position": "fixed",
				#"top": 0,
				"vertical-align": "top",
				"horizontal-align": "right",
				#"right": 0,
				#"bottom": 0,
				"width": "42rem",
				"padding": "2rem 1rem",
				"background-color": "#E6EBEE",
			}

			CONTENT_STYLE = {"margin-left": "18rem", "margin-right": "2rem", "padding": "2rem 1rem"}

			menu = html.Div(id="menu", children=
				[dbc.Row([
					html.H2("Model Menu", className="text-center",),
					]), 
					dbc.Nav(
						[
							dbc.Row([
								html.Div([
										dbc.Button("controls", id='ctrl',color="primary", style={"margin-bottom": "10px",
											"background-color": "#102A5F", "margin-left": "3px", "border-radius": "8px", "margin-top": "10px"}),
										dbc.Button("constraints", id='cnst',color="primary", style={"margin-top": "10px",
											"margin-left": "4px", "margin-bottom": "10px", "background-color": "#102A5F", "border-radius": "8px"}),
										dbc.Button("exchange", id='exg',color="primary", style={"margin-bottom": "10px",
											"margin-left": "4px", "background-color": "#102A5F", "border-radius": "8px", "margin-top": "10px"},),
										dbc.Button("model stats", id='ms',color="primary", style={"margin-left": "4px", "margin-bottom": "10px",
											"background-color": "#102A5F", "border-radius": "8px", "margin-top": "10px"}),
										dbc.Button("help", id="help", color="primary", style={"margin-left": "4px", "margin-bottom": "10px",
											"background-color": "#102A5F", "border-radius": "8px", "margin-top": "10px"},),
									])
								]),
							dbc.Row([
								html.Div(id="menu_content", className="menu", children=[])
								])
						],
						vertical=True,
						pills=True,
					),
				],
				style=sidebar_style
			)

			'''
			# TABS ---------------------------------------------------------------------
			s_c = "jfdhf"

			simulate_stats = dcc.Tab(
				label='Model Stats',
				value='Statistics',
				children=html.Div(
					className='about-tab',
					children=[
						html.H4(
						className='Statistics',
						children='Model Statistics:'),
							html.P("Below shows some basic information for the model."),
							dbc.Alert(
								id="model-stats",
								children=content_model,
								color="white",
								style={"width": "75%"}),
					]),
			)

			s_con = dcc.Tab(label='Constraints', value='constraints',
				children=html.Div(className='about-tab', children=[
					dbc.Row([html.P("""
							Select Load constraints to load the reaction
							constraints. To modify, select update constraints.
							To export changes as a new limits file, select save
							constraints.
							"""),
					dbc.Alert(id="constraints", children="display exchange here", color="secondary",),
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Button("Load constraints", id="constraint-load")
								]),
							dbc.Col([
								html.Button("Update constraints", id="constraint-modify")
								]),
							dbc.Col([
								html.Button("Save constraints", id="constraint-save")
								])]), ]),
					dbc.Alert(id="constraint-confirm", children="This will confirm that the constraint was modified", color="secondary",),
					dbc.Alert(id="constraint-save-confirm", children="confirm that constraint is saved", color="secondary",),
					]),
					]),
				)
			s_ex = dcc.Tab(label='Exchange', value='what-is',
				children=html.Div(className='about-tab', children=[
					dbc.Row([html.P("""
							Select Load exchange to load the exchange
							constraints. To modify, select update exchange.
							To export changes as a new media, select save
							exchange."""
							),
					dbc.Alert(id="exchange", children="display exchange here", color="secondary",),
					html.Div([
						dbc.Row([
							dbc.Col([
								html.Button("Load exchange", id="exchange-load")
								]), 
							dbc.Col([
								html.Button("Update exchange", id="exchange-modify")
								]),
							dbc.Col([
								html.Button("Save exchange", id="exchange-save")
								])]), ]),
					dbc.Alert(id="exchange-confirm", children="This will confirm that the exchange was modified", color="secondary",),
					dbc.Alert(id="exchange-save-confirm", children="confirm that exchange is saved", color="secondary",),
					]), ]),
				)

				'''

			s_n = cyto.Cytoscape(id='net', layout={'name': 'cose','animate': True,
				'animationDuration': 1000}, style={'width': '100%', 'height': '600px'},
				elements=nodes + edges, stylesheet=[{'selector': 'node', 'style': 
				{'background-color': 'data(col)', 'label': 'data(label)'}},
				{"selector": "edge", "style": {"width": 1, "curve-style": "bezier", "target-arrow-shape":"triangle"}},
				{"selector": "[flux > 1e-10]", "style": {"line-color": "blue", "target-arrow-color": "blue", "width": 1,"curve-style": "bezier","target-arrow-shape":"triangle","line-dash-pattern": [6, 3],"line-dash-offset": 24}}, ],
				minZoom=0.06)

			"""
			dbc.Row([html.Div(id='model-control-tabs', className='control-tabs',
			style={'width': '69%', 'font-size': '50%','height': '50%'}, 
			children=[simulate_tabs]), ]),"""
			if dash.callback_context.triggered[0]['prop_id'] == \
					'rsb_btn.n_clicks':
				return dbc.Container([
					dbc.Row([
						s_n,
						dbc.Row([
							html.Div(id="node-data", className="node_info", children=[
								html.H6("Click on a node to see its details", style={"margin-top": "10px"})
								], style={"margin-top": "20px", "background-color": "#CCEAEB"}),
							]),
						]),
					dbc.Row([dcc.Markdown('''\\* Data analysis carried out for demonstration of data 
						visualisation purposes only.''')], style={"fontSize": 11,"color": "gray"},), ]), \
				{"margin-top": "120px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}, \
				{"margin-top": "5px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}


			sim = dbc.Container([
				dbc.Row([
					dbc.Col([s_n]),
					dbc.Col([menu]),
					]),
				dbc.Row([
					html.Div(id="fba-data", className="results", children=[], style={"margin-top": "20px", "background-color": "#CCEAEB"}),
					]),
				])

			if dash.callback_context.triggered[0]['prop_id'] == \
					'open_menu_btn.n_clicks':
					return sim, {"margin-top": "120px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}, {"margin-top": "5px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}

			if pathname == "/":
				return home, {"display": "none"}, {"display": "none"}
			elif pathname == "/simulate":
				return sim, {"margin-top": "120px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}, {"margin-top": "5px", "background-color": "#102A5F", "border-radius": "8px", "display": "block"}
			elif pathname == "/curate":
				return html.P("Oh cool, this is page 2!"), {"display": "none"}, {"display": "none"}
			# If the user tries to reach a different page, return a 404 message
			return html.Div(
				[
					html.H1("404: Not found", className="text-danger"),
					html.Hr(),
					html.P(f"The pathname {pathname} was not recognised..."),
				],
				className="p-3 bg-light rounded-3",
			   	)


		@_app.callback(
			Output('menu_content', 'children'),
			Input('ctrl', 'n_clicks'),
			Input('help', 'n_clicks'),
			Input("cnst", "n_clicks"),
			prevent_initial_call=True

			# add the other callbacks (for buttons like control)
			# here --> consolidates code
		)
		def control_panel_director(ctrl_clicks, help_clicks, cnst_clicks):
			triggered_id = ctx.triggered_id
			if triggered_id == "ctrl":
				return open_ctrl(ctrl_clicks)
			elif triggered_id == "help":
				return open_help(help_clicks)
			elif triggered_id == "cnst":
				return open_cnst(cnst_clicks)


		def open_ctrl(nclicks):
			# DEFINITIONS -------------------------------------------------------------------------
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

			# stores the initial model stats in content_model and generate the
			# initial stats
			content_model = []
			content_model.append(html.H5("General Model Statistics:"))
			content_model.append(html.P("Reactions in model: " + str(count)))
			count = 0
			for i in self._model.compounds:
				count += 1
			content_model.append(html.P("Compounds in model: " + str(count)))
			content_model.append(
				html.P("Pathways in mode: " + str(len(pathway_list))))
			content_model.append(
				html.P(
					"Reactions with gene assocations: " +
					str(rxns_with_genes)))
			try:
				content_model.append(
					html.P("Total genes: " + str(count_genes(self._model))))
			except BaseException:
				content_model.append(html.P("Total genes: Unable to count genes"))
			with HiddenPrints():
				unbalanced, count, unchecked, exclude, unbalanced_list, \
					unbalance_dict = charge_check(self, 1e-6)
			content_model.append(html.H5("Chargecheck results:"))
			content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude))))

			with HiddenPrints():
				unbalanced, count, unchecked, exclude, form_list, \
					form_dict = formula_check(self)
			content_model.append(html.H5("Formulacheck results:"))
			content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
			content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
			content_model.append(
				html.P("Excluded reactions: " + str(len(exclude))))

			# initialize an empty list. the full is messy
			nodes, edges = [], []

			if dash.callback_context.triggered[0]['prop_id'] == \
					'ctrl.n_clicks':
				return html.Div(id='control-tab', className='control-tab', children=[
					dbc.Row([
						dbc.Col([
							html.Div([
								html.H5("Pathways", className="text-center",
									style={"margin-top": "10px"}),
								dcc.Dropdown(id="pathways_dropdown",options=[
									{"label": i,"value": i, } for i in list(pathway_list)],
									value=pathway_list[0], multi=True, placeholder="",
									style={"margin-top": "10px"}),
								])
							]),
						dbc.Col([
							html.Div([
								html.H5("Reactions", className="text-center",
									style={"margin-top": "10px"}),
								dcc.Dropdown(id="reaction_dropdown", options=[
									{"label": i.id,"value": i.id, } for i in list(self._model.reactions)],
									value=pathway_list[0],multi=True,placeholder="",
									style={"margin-top": "10px"}),
								]),
							]),
						]),
					dbc.Row([
						html.Div([
							html.H5("Element Transfer Networks", className="text-center", style={"margin-top": "15px"}),
							dcc.Dropdown(id="element_dropdown",options=[
								{"label": i,"value": i, } for i in ["C", "N", "S", "P"]],
								value=["C", "N", "S", "P"],multi=False),
							])
						]),
					dbc.Row([
						html.Div([
							html.H5("Compound Search", className="text-center", style={"margin-top": "15px"}),
							dcc.Dropdown(id="compounds_dropdown",options=[
								{"label": i,"value": i, } for i in list(compounds_list)],
								value=compounds_list,multi=False,style={"width": "100%"}),
							])
						]),
					dbc.Row(html.H5("Path Search", className="text-center", style={"margin-top": "15px"}),),
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
						dbc.Col([
							html.Div([
								html.H5("Gene Delete", className="text-center", style={"margin-top": "15px"}),
								dcc.Dropdown(id="delete_dropdown",options=[
									{"label": i,"value": i, } for i in list(gene_content)],
									value=None,multi=True,style={"width": "100%"}),
								]),
							]),
						dbc.Col([
							html.Div([
								html.H5("Flux Analysis", className="text-center", style={"margin-top": "15px"}),
								dcc.Dropdown(id="fba_dropdown",options=[{
									"label": i,"value": i, } for i in list(rxns_full)],
									value=rxns,multi=False,style={"width": "100%"},),
								]),
							]),
						]),
					dbc.Row([
						html.Div([
							dbc.Button("Submit",id="btn_sub", style={"margin-top": "20px", "margin-left": "60px"}, outline=True, n_clicks=0),
							dbc.Button("Download CSV",id="btn_tsv", style={"margin-top": "20px", "margin-left": "60px"}, outline=True),
							dcc.Download(id="download_r"),
							dcc.Download(id="download_c"),
							dbc.Button("Download png", id="btn-get-png", style={"margin-top": "20px", "margin-left": "60px"}, outline=True),
							dcc.Download(id="downloadpng"),
							])
						]),
					])
		def open_cnst(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == 'cnst.n_clicks':
				return html.Div([
							dbc.Row([
								dbc.Col([
									html.Button("Load constraints", id="constraint-load")
									]),
								dbc.Col([
									html.Button("Update constraints", id="constraint-modify")
									]),
								dbc.Col([
									html.Button("Save constraints", id="constraint-save")
									])
								]), 
						]),
		def open_help(nclicks):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'help.n_clicks':
				return html.Div([
					dbc.Row([
						html.H4("Welcome to the help page", className="text-center", style={"margin-top": "20px", "color": "#102A5F"}),
						]),
					dbc.Row([
						html.Div([
							dbc.Button("pathways & reactions", id='pr_help',color="secondary", style={"margin-bottom": "10px", "margin-left": "70px",
								"margin-top": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							dbc.Button("element transfer network", id='etn_help',color="secondary", outline=True, style={
								"margin-left": "15px", "margin-bottom": "10px", "margin-top": "10px", "border-radius": "8px", "color": "#102A5F"}),
							]) 
						]),
					dbc.Row([
						html.Div([
							dbc.Button("gene delete", id='gd_help',color="secondary", outline=True, style={"margin-bottom": "10px",
								"margin-left": "32px", "border-radius": "8px", "color": "#102A5F"},),
							dbc.Button("path & compound search", id='pcs_help',color="secondary", style={"margin-left": "15px", 
								"margin-bottom": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							dbc.Button("flux analysis", id='fa_help',color="secondary", style={"margin-left": "15px", 
								"margin-bottom": "10px", "border-radius": "8px", "color": "#102A5F"}, outline=True),
							])
						]),
					dbc.Row([html.Div(id="help_area", className="help", children=[])]),
					])

		@_app.callback(
			Output("help_area", "children"),
			Input("pr_help", "n_clicks"),
			Input("etn_help", "n_clicks"),
			Input("gd_help", "n_clicks"),
			Input("pcs_help", "n_clicks"),
			Input("fa_help", "n_clicks"),
		)
		def det_help(pr_clicks, etn_clicks, gd_clicks, pcs_clicks, fa_clicks):
			triggered_id = ctx.triggered_id
			if triggered_id == "pr_help":
				return open_pr_help(pr_clicks)
			elif triggered_id == "etn_help":
				return open_etn_help(etn_clicks)
			elif triggered_id == "gd_help":
				return open_gd_help(gd_clicks)
			elif triggered_id == "pcs_help":
				return open_pcs_help(pcs_clicks)
			elif triggered_id == "fa_help":
				return open_fa_help(fa_clicks)
		def open_pr_help(n_clicks):
			return html.Div([
				dbc.Row([
					html.Div([
						html.H5("Pathways", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the dropdown to display pathways in the metabolic model. Multiple can be selected.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				dbc.Row([
					html.Div([
						html.H5("Reactions", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the dropdown to display reactions in the metabolic model. Multiple can be selected.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				])
		def open_etn_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H5("Element Transfer Networks", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the dropdown to choose a chemical element for the network (default Carbon).",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				])
		def open_gd_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H5("Gene Delete", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("""Use the dropdown to delete genes from the model. If all genes associated with a reaction are
							deleted, the reaction will be removed from analysis.""",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				])
		def open_pcs_help(nclicks):
			return html.Div([
				dbc.Row([
					html.Div([
						html.H5("Path Search", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the path search to conduct a bidirectional breadth first search to show the shortest path between two compounds.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				dbc.Row([
					html.Div([
						html.H5("Compound Search", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the compound search to show all reactions containing this compound. Select ALL for pathways to see all reactions.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				])
		def open_fa_help(nclicks):
			return html.Div([
				dbc.Row([
					dbc.Col([
						html.H5("Flux Analysis", className="text-center", style={"margin-top": "30px", "color": "#102A5F"}),
						html.H6("Use the dropdown to choose a reaction to optimize for flux balance analysis. Flux is visualized in blue.",
							className="text-center", style={"margin-left": "60px", "margin-right": "60px", "color": "#102A5F"}),
						]),
					]),
				])


		@_app.callback(
			Output("rxn-save-confirm", "children"),
			Input("btn_rxn", "n_clicks"),
			[State("r_id_input", "value"),
			 State("r_name_input", "value"),
			 State("r_eq_input", "value"),
			 State("r_pathway_input", "value")]
		)
		def save_rxn(nclicks, id, name, eq, path):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'btn_rxn.n_clicks':
				equation = parse_reaction_equation_string(eq, "c")
				properties = {
					"id": id,
					"name": name,
					"equation": equation,
					"pathways": [path]}
				reaction = ReactionEntry(properties, filemark=None)
				self._model.reactions.add_entry(reaction)
				self._mm._reaction_set.add(reaction.id)
				self._mm._database.set_reaction(reaction.id, equation)
				return("Added {} to the model".format(id))

		@_app.callback(
			Output("cpd-save-confirm", "children"),
			Input("btn_cpd", "n_clicks"),
			[State("c_id_input", "value"),
			 State("c_name_input", "value"),
			 State("c_charge_input", "value"),
			 State("c_form_input", "value")]
		)
		def save_rxn(nclicks, id, name, charge, form):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'btn_cpd.n_clicks':
				properties = {
					"id": id,
					"name": name,
					"charge": charge,
					"formula": form}
				compound = CompoundEntry(properties, filemark=None)
				self._model.compounds.add_entry(compound)
				# reaction = ReactionEntry(properties, filemark = None)
				# self._model.reactions.add_entry(reaction)
				# self._mm._reaction_set.add(reaction.id)
				# self._mm._database.set_reaction(reaction.id, equation)
				return("Added {} to the model".format(id))

		### Constraint Start
		@_app.callback(
			Output("constraint-save-confirm", "children"),
			Input("constraint-save", "n_clicks"),
			State("constraints", "children")
		)
		def save_file(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'constraint-save.n_clicks':
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
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)

				path = FilePathContext(self._args.model)
				with open('{}/limits_curated.yaml'.format(path), "w") as f:
					for e in ex:
						if e[0] != '':
							f.write(
								"{}\t{}\t{}\n".format(
									e[0], e[1], e[2]))
				return("Constraints saved!")

		@_app.callback(
			Output("constraint-confirm", "children"),
			Input("constraint-modify", "n_clicks"),
			State("constraints", "children")
		)
		def modify_constraint(nclicks, datalist):
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
			if dash.callback_context.triggered[0]['prop_id'] == \
					'constraint-modify.n_clicks':
				const = []
				for i in datalist:
					if isinstance(i, dict):
						const_temp = rec_props(i['props']['children'], [])
						if len(const_temp) > 0:
							const.append(const_temp)
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

		@_app.callback(
			Output("constraints", "children"),
			Input("constraint-load", "n_clicks")
		)
		def load_constraints(nclicks):
			rxn_set = set()
			for i in self._model.reactions:
				rxn_set.add(i.id)
			rxn_list = list(rxn_set)
			contents = []
			contents.append(
				dbc.Row([
					dbc.Col([
						html.P("ID")], width=4),
					dbc.Col([
						html.P("Lower Bound")], width=4),
					dbc.Col([
						html.P("Upper Bound")], width=4),
				]),)
			if dash.callback_context.triggered[0]['prop_id'] == \
					'constraint-load.n_clicks':
				print(self._model._limits)
				constraint = set()
				for c in self._model._limits:
					print(c)
					lower = self._model._limits[c][1]
					upper = self._model._limits[c][2]
					print(self._model._limits[c], lower, upper)
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
						dbc.Row([
							dbc.Col([
								dcc.Input(id="id{}".format(str(c)),
									type="text",
									placeholder=str(c),
									style={'width': '75%'})], width=4),
							dbc.Col([
								dcc.Input(id="lower{}".format(str(c)),
									type="number",
									placeholder=lower,
									style={'width': '85%'})], width=4),
							dbc.Col([
								dcc.Input(id="upper{}".format(str(c)),
									type="number",
									placeholder=upper,
									style={'width': '85%'})], width=4),
						]),)
				contents.append(html.P("Add New Constraint:"))
				contents.append(
					dbc.Row([
						dbc.Col([
							dcc.Input(id="id{}".format("new"),
								type="text",
								placeholder=str(''),
								style={'width': '75%'})], width=4),
						dbc.Col([
							dcc.Input(id="lower{}".format("new"),
								type="number",
								placeholder=-1000,
								style={'width': '85%'})], width=4),
						dbc.Col([
							dcc.Input(id="upper{}".format(str(c)),
								type="number",
								placeholder=1000,
								style={'width': '85%'})], width=4),
					]),
				)
			return(contents)

			### Constraint End

		@_app.callback(
			Output("exchange-save-confirm", "children"),
			Input("exchange-save", "n_clicks"),
			State("exchange", "children")
		)
		def save_file(nclicks, datalist):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'exchange-save.n_clicks':
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
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)

				path = FilePathContext(self._args.model)
				with open('{}/exchange_curated.yaml'.format(path), "w") as f:
					for e in ex:
						if e[0] != '':
							f.write("{}\t{}\t{}\t{}\n".format(e[0], e[1], e[2], e[3]))
				return("Exchange saved!")

		@_app.callback(
			Output("exchange-confirm", "children"),
			Input("exchange-modify", "n_clicks"),
			State("exchange", "children")
		)
		def modify_exchange(nclicks, datalist):
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
			if dash.callback_context.triggered[0]['prop_id'] == \
					'exchange-modify.n_clicks':
				ex = []
				for i in datalist:
					if isinstance(i, dict):
						ex_temp = rec_props(i['props']['children'], [])
						if len(ex_temp) > 0:
							ex.append(ex_temp)
				c = []
				for i in self._model.compounds:
					c.append(i.id)
				rxns = []
				for i in self._model.reactions:
					rxns.append(i.id)
				for e in ex:
					if e[0] in c:
						tup = (Compound(e[0]).in_compartment(e[1]),
							   None, e[2], e[3])
						self._model._exchange[Compound(
							e[0]).in_compartment(e[1])] = tup
						self._mm._limits_lower["EX_{}[{}]".format(
							e[0], e[1])] = float(e[2])
						self._mm._limits_upper["EX_{}[{}]".format(
							e[0], e[1])] = float(e[3])
						if self._mm._database.has_reaction(
								"EX_{}[{}]".format(e[0], e[1])) is False:
							properties = {
								"id": "EX_{}[{}]".format(e[0], e[1]), "equation": Reaction(
									Direction.Forward, [
										(Compound(
											e[0]).in_compartment(
											e[1]), -1)])}
							reaction = ReactionEntry(properties, filemark=None)
							self._model.reactions.add_entry(reaction)
							self._mm._reaction_set.add(reaction.id)
							self._mm._database.set_reaction(
								reaction.id, Reaction(
									Direction.Forward, [
										(Compound(
											e[0]).in_compartment(
											e[1]), -1)]))
					else:
						return("Added exchange compound not in model")

		@_app.callback(
			Output("exchange", "children"),
			Input("exchange-load", "n_clicks")
		)
		def load_exchange(nclicks):
			comp_set = set()
			for i in self._model.reactions:
				for j in i.equation.compounds:
					comp_set.add(j[0].compartment)
			comp_list = list(comp_set)
			contents = []
			contents.append(
				dbc.Row([
					dbc.Col([
						html.P("ID")], width=3),
					dbc.Col([
						html.P("compartment")
					], width=3),
					dbc.Col([
						html.P("Lower")], width=3),
					dbc.Col([
						html.P("Upper")], width=3),
				]),)
			if dash.callback_context.triggered[0]['prop_id'] == \
					'exchange-load.n_clicks':
				exchange = set()
				for e in self._model._exchange:
					lower = self._model._exchange[e][2]
					upper = self._model._exchange[e][3]
					if float(lower) != 0 and float(lower) != -1000:
						exchange.add(e)
					if float(upper) != 0 and float(upper) != 1000:
						exchange.add(e)
				for e in exchange:
					contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Input(id="id{}".format(str(e.name)),
									type="text",
									placeholder=str(e.name),
									style={'width': '75%'})], width=3),
							dbc.Col([
								dcc.Dropdown(
									id="compartment{}".format(str(e.name)),
									options=[{"label": x,
											  "value": x, }
											 for x in list(comp_list)],
									placeholder=e.compartment,
									multi=False
								)
							], width=3),
							dbc.Col([
								dcc.Input(id="lower{}".format(str(e.name)),
									type="number",
									placeholder=float(self._model._exchange[e][2]),
									style={'width': '85%'})], width=3),
							dbc.Col([
								dcc.Input(id="upper{}".format(str(e.name)),
									type="number",
									placeholder=float(self._model._exchange[e][3]),
									style={'width': '85%'})], width=3),
						]),)
				contents.append(html.P("Add New Exchange Constraint:"))
				contents.append(
					dbc.Row([
						dbc.Col([
							dcc.Input(id="id{}".format("new"),
								type="text",
								placeholder=str(''),
								style={'width': '75%'})], width=3),
						dbc.Col([
							dcc.Dropdown(
								id="compartment{}".format("new"),
								options=[{"label": x,"value": x, } for x in list(comp_list)],
								placeholder='',
								multi=False
							)
						], width=3),
						dbc.Col([
							dcc.Input(id="lower{}".format("new"),
							   	type="number",
								placeholder=-1000,
								style={'width': '85%'})], width=3),
						dbc.Col([
							dcc.Input(id="upper{}".format(str(e.name)),
								type="number",
								placeholder=1000,
								style={'width': '85%'})], width=3),
					]),
				)
			return(contents)
		@_app.callback(
			Output("node-data", "children"),
			[Input("net", "selectedNodeData")]
		)
		def display_nodedata(datalist):
			if datalist is not None:
				if len(datalist) > 0:
					contents = []
					data = datalist[-1]
					if "formula" in data:
						if len(datalist) > 0:
							data = datalist[-1]
							contents.append(
								html.H5(
									"ID: " + data["id"].title()
								)
							)
							contents.append(
								html.P(
									"Name: " + data["label"]
								)
							)
							contents.append(
								html.P(
									"Formula: "
									+ str(data["formula"])
								)
							)
							contents.append(
								html.P(
									"Charge: "
									+ str(data["charge"])
								)
							)
					elif "equation" in data:
						if isinstance(data["pathways"], str):
							path = [data["pathways"]]
						else:
							path = data["pathways"]
						contents.append(html.H5("ID: " +
												data["orig_id"].title()))
						contents.append(html.P("Name: " + data["label"]))
						contents.append(html.P("Equation: " +
											   str(data["equation"])))
						if isinstance(data["pathways"], str):
							contents.append(html.P("Pathways: " +
												   data["pathways"]))
						elif isinstance(data["pathways"], list):
							contents.append(html.P("Pathways: " +
												   ";".join(data["pathways"])))
						contents.append(html.P("Flux: " +
											   str(data["flux_node"])))
					return html.Div(contents)
			else:
				return html.Div([
					html.H5("Select a node to see data")
					])

		@_app.callback(
			Output("net_curate", "elements"),
			Input("btn_update_cur_net", "n_clicks"),
			State("pathways_dropdown_curate", "value"),
			prevent_initial_call=True,
		)
		def update_net_cur(nclicks, pathway):
			if dash.callback_context.triggered[0]['prop_id'] == \
					'btn_update_cur_net.n_clicks':
				if pathway != '':
					model = self._mm
					nm, network = read_model(self._model, model, "C")
					pathway_list, rxn_set = get_pathway_list(nm, pathway)
					# nodes, edges, 
					obj_flux = build_network(
						self, nm, model, rxn_set, network, [], nodes, edges)
					elements = nodes + edges
				else:
					elements = []
				return(elements)

		@_app.callback(
			[Output("reaction_dropdown_curate", "value"),
			 Output("compound_dropdown_curate", "value"),
			 ],
			[Input("net_curate", "selectedNodeData"), ],
			prevent_initial_call=True,
		)
		def updated_selected(node):
			# if dash.callback_context.triggered[0]['prop_id'] == \
			#		 'net_curate.selectedEdgeData' and len(edge) > 0:
			#	 return(edge[0]['rid'], [])
			if dash.callback_context.triggered[0]['prop_id'] == \
					'net_curate.selectedNodeData' and len(node) > 0:
				if node[0]['type'] == "reaction":
					return([node[0]['label']], [])
				else:
					return([], node[0]['id'][0:(str(node[0]['id']).find('['))])
			else:
				return([], [])

		@_app.callback(
			[Output("edit-confirm", "children"),
			 ],
			[Input("btn_save", "n_clicks"), ],
			[State("reaction_dropdown_curate", "value"),
			 State("compound_dropdown_curate", "value"),
			 State("edit-data", "children"), ],
			prevent_initial_call=True,
		)
		def update_model(nclicks, r_dropdown_curate,
						 c_dropdown_curate, datalist):

			changed_id = [p['prop_id'] for p in
						  dash.callback_context.triggered][0]
			contents = []

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
			if dash.callback_context.triggered[0]['prop_id'] == \
					'btn_save.n_clicks':
				if isinstance(r_dropdown_curate, str) and \
						len(r_dropdown_curate) > 0:
					direction = []
					for i in datalist:
						if isinstance(i, dict):
							if 'id' in i['props']:
								if i['props']['id'] == 'direction':
									direction_temp = rec_props(i['props']
															   ['children'],
															   [])
									direction.append(direction_temp)
					left = []
					for i in datalist:
						if isinstance(i, dict):
							if 'id' in i['props']:
								if i['props']['id'] == 'left':
									left_temp = rec_props(i['props']
														  ['children'], [])
									left.append(left_temp)
								elif i['props']['id'] == 'left_new':
									left_temp = rec_props(i['props']
														  ['children'], [])
									if left_temp[1] != '':
										left.append(left_temp)
					right = []
					for i in datalist:
						if isinstance(i, dict):
							if 'id' in i['props']:
								if i['props']['id'] == 'right':
									right_temp = rec_props(i['props']
														   ['children'], [])
									right.append(right_temp)
								elif i['props']['id'] == 'right_new':
									right_temp = rec_props(i['props']
														   ['children'], [])
									if right_temp[1] != '':
										right.append(right_temp)

					rid = re.split(" ", datalist[0]['props']['children'])[1]

					for r in self._model.reactions:
						if r.id.upper() == rid.upper():
							reaction = r

					eq_dict = {}
					for le in left:
						eq_dict[Compound(le[1]).in_compartment(le[2])] = \
							int(le[0]) * - 1
					for r in right:
						eq_dict[Compound(r[1]).in_compartment(r[2])] = \
							int(r[0])
					dir = direction[0][0]

					if dir == "Direction.Forward":
						new_eq = Reaction(Direction.Forward, eq_dict)
					elif dir == "Direction.Reverse":
						new_eq = Reaction(Direction.Reverse, eq_dict)
					elif dir == "Direction.Both":
						new_eq = Reaction(Direction.Both, eq_dict)
					else:
						quit("Nonvalid direction")

					reaction.equation = new_eq
					self._model.reactions.add_entry(reaction)

					contents.append(html.H5(
									"{} has been updated".format(rid)))
				elif isinstance(c_dropdown_curate, str) and \
						len(c_dropdown_curate) > 0:
					comp = rec_props(datalist[-1]['props']['children'], [])

					cid = re.split(" ", datalist[0]['props']['children'])[1]

					for c in self._model.compounds:
						if c.id.upper() == cid.upper():
							compound = c

					compound.charge = comp[0]
					compound.formula = comp[1]
					self._model.compounds.add_entry(compound)
					contents.append(html.H5(
									"{} has been updated".format(cid)))
				return contents

		@_app.callback([Output("edit-data", "children"),
						Output("edit-data", "is_open")],
					   [Input("btn_update", "n_clicks"),
						Input("reaction_dropdown_curate", "value"),
						Input("compound_dropdown_curate", "value")])
		def display_nodedata(clicks, r, c):
			print(r)
			changed_id = [p['prop_id'] for p in
						  dash.callback_context.triggered][0]
			contents = "Click on an edge to see its details here"
			model = self._model
			c_list = []
			comp_set = set()
			for i in model.compounds:
				c_list.append(i.id)
			for i in model.reactions:
				for j in i.equation.compounds:
					comp_set.add(j[0].compartment)
			comp_list = list(comp_set)
			if dash.callback_context.triggered[0]['prop_id'] == \
					'btn_update.n_clicks':
				if r is not None and len(r) > 0:
					for i in model.reactions:
						if isinstance(r, list):
							if i.id == r[0]:
								reaction = i
						else:
							if i.id == r:
								reaction = i
					print(reaction)
					contents = []
					contents.append(html.H5(
									"ID: " + reaction.id))
					contents.append(html.H5(
									"Name: " + reaction.name))
					contents.append(
						dbc.Row([
							dbc.Col([
								html.H5("Direction: ")
							]),
							dbc.Col([
								dcc.Dropdown(
									id="direction",
									options=[{"label": x,
											  "value": x, }
											 for x in ["Direction.Forward",
													   "Direction.Reverse",
													   "Direction.Both"]],
									placeholder=str(
										reaction.equation.direction),
									multi=False
								)
							]), ], id='direction'),)
					contents.append(html.H5("Left Equation:"))
					dropid = 0
					for j in reaction.equation.left:
						contents.append(
							dbc.Row([
								dbc.Col([
									dcc.Input(id="stoich{}".format(
										str(dropid)),
											  type="number",
											  placeholder=str(j[1]))]),
								dbc.Col([
									dcc.Dropdown(
										id="drop{}".format(str(dropid)),
										options=[{"label": x,
												  "value": x, }
												 for x in list(c_list)],
										placeholder=j[0].name,
										multi=False
									)
								]),
								dbc.Col([
									dcc.Dropdown(
										id="compartment{}".format(str(dropid)),
										options=[{"label": x,
												  "value": x, }
												 for x in list(comp_list)],
										placeholder=j[0].compartment,
										multi=False
									)
								]), ], id="left"),)
						dropid += 1
					contents.append(html.P("Add new compound to left side"))
					contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Input(id="stoich{}".format(str(dropid)),
										  type="number",
										  placeholder=0)]),
							dbc.Col([
								dcc.Dropdown(
									id="drop{}".format(str(dropid)),
									options=[{"label": x,
											  "value": x, }
											 for x in list(c_list)],
									placeholder='',
									multi=False
								)
							]),
							dbc.Col([
								dcc.Dropdown(
									id="compartment{}".format(str(dropid)),
									options=[{"label": x,
											  "value": x, }
											 for x in list(comp_list)],
									placeholder='',
									multi=False
								)
							]), ], id="left_new"))
					dropid += 1
					contents.append(html.H5("Right Equation:"))
					for j in reaction.equation.right:
						contents.append(
							dbc.Row([
								dbc.Col([
									dcc.Input(id="stoich{}".format(
										str(dropid)),
											  type="number",
											  placeholder=str(j[1]))]),
								dbc.Col([
									dcc.Dropdown(
										id="drop{}".format(str(dropid)),
										options=[{"label": x,
												  "value": x, }
												 for x in list(c_list)],
										placeholder=j[0].name,
										multi=False
									)
								]),
								dbc.Col([
									dcc.Dropdown(
										id="compartment{}".format(str(dropid)),
										options=[{"label": x,
												  "value": x, }
												 for x in list(comp_list)],
										placeholder=j[0].compartment,
										multi=False)
								]), ], id="right"), )
						dropid += 1
					contents.append(html.P("Add new compound to right side"))
					contents.append(
						dbc.Row(
							[dbc.Col([
								dcc.Input(id="stoich{}".format(str(dropid)),
									 type="number",
									 placeholder=0)]),
								dbc.Col([
									dcc.Dropdown(id="drop{}".format(
										str(dropid)),
												 options=[{"label": x,
														   "value": x}
														  for x in
														  list(c_list)],
												 placeholder='',
												 multi=False)]),
								dbc.Col([
									dcc.Dropdown(id="compart{}".
												 format(str(dropid)),
												 options=[{"label": x,
														   "value": x}
														  for x in
														  list(comp_list)],
												 placeholder=j[0].
												 compartment,
												 multi=False)])],
							id="right_new"))
					dropid += 1
				elif c is not None and len(c) > 0:
					for i in model.compounds:
						if i.id == c:
							compound = i
					contents = []
					contents.append(html.H5(
									"ID: " + compound.id))
					contents.append(html.H5(
									"Name: " + compound.name))
					contents.append(
						dbc.Row([
							dbc.Col([
								dcc.Input(id="charge",
									type="number",
									placeholder=compound.charge)]),
								dbc.Col([
									dcc.Input(id="formula", type='text',
											  placeholder=compound.formula)]),
								], id='compounds'), )

				return contents, True
			else:
				return [], False

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
		def filter_nodes(n_clicks, pathways_dropdown, element_dropdown,
						 compounds_dropdown, fba_dropdown,
						 filter1_dropdown, filter2_dropdown, delete_dropdown,
						 r_dropdown):
			if isinstance(pathways_dropdown, list) \
					or isinstance(delete_dropdown, list) \
					or isinstance(element_dropdown, str) \
					or isinstance(compounds_dropdown, str) \
					or isinstance(fba_dropdown, str) \
					or (isinstance(filter1_dropdown, str)
						and isinstance(filter1_dropdown, str)) \
					or isinstance(r_dropdown, list):
				model = self._mm.copy()
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
							variables = [boolean.Variable(g)
										 for g in reaction.genes]
							assoc = boolean.Expression(boolean.And(*variables))
						genes.update(v.symbol for v in assoc.variables)
						gene_assoc[reaction.id] = assoc
						reactions = set(model.reactions)
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
								logger.info('Deleting reaction {}...'.
											format(reaction))
								deleted_reactions.add(reaction)
					for r in deleted_reactions:
						model.remove_reaction(r)

				nm, network = read_model(self._model, model, element_dropdown)
				pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)

				if r_dropdown is not None:
					for i in r_dropdown:
						rxn_set.add(i)

				if isinstance(filter1_dropdown, str) and \
						isinstance(filter2_dropdown, str):
					rxn_list = set()
					cpd_list = [filter1_dropdown, filter2_dropdown]
					middle2 = []
					middle3 = []
					middle2, rxn_list = bfs_compounds(filter1_dropdown,
													  filter2_dropdown,
													  network, rxn_list,
													  rxn_set, middle2,
													  middle3)
					if len(rxn_list) == 0:
						return []
				elif isinstance(compounds_dropdown, str):
					rxn_list = []
					for rxn in network[0]:
						for cpd in network[0][rxn][0]:
							for i in cpd:
								if i.name == compounds_dropdown and \
										rxn.id in rxn_set:
									rxn_list.append(rxn.id)
				else:
					rxn_list = rxn_set

				if fba_dropdown is None:
					fba_dropdown = []
				# nodes, edges, 
				obj_flux = build_network(self, nm, model,
													   rxn_list, network,
													   fba_dropdown, nodes, edges)
				elements = nodes + edges

				if obj_flux != 0:
					contents = []
					contents.append(html.H5("{} flux "
									"is {}".format(fba_dropdown, obj_flux)))
				else:
					contents = []
					contents.append(html.H5("No flux through the network"))

				return elements, contents
			else:
				return tuple(), list()



def build_app(self):
	'''
	Generates the initial required data for the visualization and
	launches the app on a local port.
	'''

	# Initialize the app and launch the server
	# self._app = dash.Dash(
	#	 __name__, external_stylesheets=[
	#		 dbc.themes.BOOTSTRAP])
	# self._app = app
	# server = self._app.server
	# if self._app is not None and hasattr(self, "callbacks"):
	#	 self.callbacks(self._app)

	# Set styles for the app
	styles = {'output': {'border': 'thin lightgrey solid','overflowX': 'scroll'}}
	col_swatch = px.colors.qualitative.Dark24
	default_stylesheet = [{'selector': 'node','style': {'background-color': '#BFD7B5',
	'label': 'data(label)'}}, {"selector": "edge", "style": {"width": 1, "curve-style": "bezier"}}]

	sidebar_style = {
		"position": "fixed",
		"top": 0,
		"left": 0,
		"bottom": 0,
		"width": "16rem",
		"padding": "2rem 1rem",
		"background-color": "#E6EBEE",
	}

	sidebar = html.Div(
		[
			html.H2("PSAMM", className="display-4", style={"color": "#102A5F"}),
			html.Hr(),
			html.P(
				"Interactive Metabolic Modelling", className="lead", style={"color": "#102A5F"},
			),
			dbc.Nav(
				[
					dbc.NavLink("Home", href="/", active="exact", style={"background-color": "#102A5F", "border-radius": "8px", "color": "#E6EBEE"}),
					dbc.NavLink("Simulate", href="/simulate", active="exact", style={"background-color": "#102A5F", "border-radius": "8px",
						"margin-top": "5px", "color": "#E6EBEE"}),
					dbc.NavLink("Curate", href="/curate", active="exact", style={"background-color": "#102A5F", "border-radius": "8px",
						"margin-top": "5px", "color": "#E6EBEE"}),
					dbc.NavLink("Documentation", href="https://psamm.readthedocs.io", active="exact", style={"background-color": "#102A5F",
					 "border-radius": "8px","margin-top": "5px", "color": "#E6EBEE"})
				],
				vertical=True,
				pills=True,
			),
			dbc.Nav(
				[
					dbc.Button("Close Menu",id="rsb_btn", style={"margin-top": "80px", "margin-left":"50px",
						"background-color": "#102A5F", "border-radius": "8px", "display": "block"}),
					dbc.Button("Open Menu",id="open_menu_btn", style={"margin-top": "5px", "margin-left":"60px",
						"background-color": "#102A5F", "border-radius": "8px", "display": "block"}),
				],
				vertical=True,
			),

		],
		style=sidebar_style,
	)

	content = html.Div(id="page-content", style={"margin-left": "18rem", "margin-top": "1rem"})

	return html.Div([dcc.Location(id="url"), sidebar, content])