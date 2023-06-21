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
from psamm.datasource.context import FilePathContext, FileMark
from psamm.command import MetabolicMixin, Command, FilePrefixAppendAction, \
	convert_to_unicode, SolverCommandMixin
from psamm import graph
import sys
from psamm.formula import Atom, _AtomType
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
from dash import callback
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import dash_cytoscape as cyto
import re
import pandas as pd
import webbrowser
from psamm.commands.chargecheck import charge_check
from psamm.commands.formulacheck import formula_check
from psamm.datasource.native import ModelWriter, NativeModel
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry, DictReactionEntry as ReactionEntry, DictCompartmentEntry as CompartmentEntry)
from decimal import *
import yaml
from psamm.datasource.native import parse_reaction_equation_string
from psamm.datasource.reaction import parse_reaction
from pkg_resources import resource_filename

logger = logging.getLogger(__name__)
dash.register_page(__name__, path='/')

REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'
EPSILON = 1e-5







def build_app(self):
    '''
    Generates the initial required data for the visualization and
    launches the app on a local port.
    '''

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

    # Initialize the app and launch the server
    # self._app = dash.Dash(
    #     __name__, external_stylesheets=[
    #         dbc.themes.BOOTSTRAP])
    # self._app = app
    # server = self._app.server
    # if self._app is not None and hasattr(self, "callbacks"):
    #     self.callbacks(self._app)

    # Set styles for the app
    styles = {
        'output': {
            'border': 'thin lightgrey solid',
            'overflowX': 'scroll'},
    }
    col_swatch = px.colors.qualitative.Dark24
    default_stylesheet = [
        {
            'selector': 'node',
            'style': {'background-color': '#BFD7B5',
                      'label': 'data(label)'}},
        {"selector": "edge", "style": {"width": 1,
                                       "curve-style": "bezier"}}]


    # Set the navigation bar for the app
    navbar = dbc.Navbar(
        dbc.Container(
            [
                html.A(
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.NavbarBrand(
                                    "Interactive Metabolic Modelling "
                                    "with PSAMM",
                                    className="ms-2"), 
                                    align="center"
                                ),
                            dbc.Col(
                                html.Div(
                                    [
                                        dcc.Dropdown(
                                            options=["Simulate", "Curate", "Documentation"],
                                            id="startmenu"), 
                                    ],
                                    style={"width": "200%"},
                                ),
                            ),
                        ],
                        align="center",
                        className="g-0",
                    ),
                    style={"textDecoration": "none"},
                    #href="https://psamm.readthedocs.io",
                ),
            ], fluid=True,
        ),
        color="dark",
        dark=True, 
    )


    # control_tabs is an object that stores the control tabs under the
    # s_c - simulate tab
    s_c = dcc.Tab(label='Controls', id='Controls',
                  children=html.Div(id='control-tab', className='control-tab', children=[
                      html.H5(
                          "Pathways:"),
                      dcc.ConfirmDialog(id="pathway_dialog",
                                        message="Select the pathways you "
                                        "want to display from the "
                                        "metabolic model. Multiple can be "
                                        "selected. Reactions can be added "
                                        "to this display from the "
                                        "Reactions dropdown."),
                      dbc.Row([
                              dbc.Col(
                                  dcc.Dropdown(
                                      id="pathways_dropdown",
                                      options=[
                                          {"label": i,
                                           "value": i, }
                                          for i in list(pathway_list)],
                                      value=pathway_list[0],
                                      multi=True,
                                      placeholder="",
                                      style={"width": "100%"}), width=11),
                              dbc.Col(dbc.Button("i", id='pathway_help',
                                                 outline=True,
                                                 color="primary",
                                                 size="sm"),
                                      width=1), ]),


                      html.H5("Reactions:"),
                      dcc.ConfirmDialog(id="reaction_dialog",
                                        message="Select the reactions you "
                                        "want to display from the "
                                        "metabolic model. Multiple can "
                                        "be selected."),
                      dbc.Row([
                              dbc.Col(
                                  dcc.Dropdown(
                                      id="reaction_dropdown",
                                      options=[
                                          {"label": i.id,
                                           "value": i.id, }
                                          for i in
                                          list(self._model.reactions)],
                                      value=pathway_list[0],
                                      multi=True,
                                      placeholder="",
                                      style={"width": "100%"},), width=11),
                              dbc.Col(dbc.Button("i", id='reaction_help',
                                                 outline=True,
                                                 color="primary",
                                                 size="sm"), width=1), ]),
                      html.H5("Element Transfer Networks:"),
                      dcc.ConfirmDialog(id="element_dialog",
                                        message="Choose the chemical "
                                        "element for the network. Default "
                                        "is Carbon (C). Options are "
                                        "Carbon (C), Nitrogen (N), Sulfur "
                                        "(S), or Phosphorous (P)"),
                      dbc.Row([
                              dbc.Col(
                                  dcc.Dropdown(
                                      id="element_dropdown",
                                      options=[
                                          {"label": i,
                                           "value": i, }
                                          for i in ["C", "N", "S", "P"]],
                                      value=[
                                          "C", "N", "S", "P"],
                                      multi=False),
                                  width=11),
                              dbc.Col(dbc.Button("i", id='element_help',
                                                 outline=True,
                                                 color="primary",
                                                 size="sm"),
                                      width=1), ]),
                      html.H5("Compound Search:"),
                      dcc.ConfirmDialog(id="compound_dialog",
                                        message="Select a compound to "
                                        "show all reactions containing "
                                        "that compound. Note that this "
                                        "combines with the pathways and "
                                        "reactions options above. If you "
                                        "want to see all reactions with "
                                        "that comound, select All for "
                                        "pathways."),
                      dbc.Row([
                              dbc.Col(
                                  dcc.Dropdown(
                                      id="compounds_dropdown",
                                      options=[
                                          {"label": i,
                                           "value": i, }
                                          for i in list(compounds_list)],
                                      value=compounds_list,
                                      multi=False,
                                      style={"width": "100%"}, ),
                                  width=11),
                              dbc.Col(
                                  dbc.Col(dbc.Button("i",
                                                     id='compound_help',
                                                     outline=True,
                                                     color="primary",
                                                     size="sm"), width=1),
                              ), ]),
                      html.H5(
                          "Path Search:"),
                      dcc.ConfirmDialog(id="bfs_dialog",
                                        message="Select two compounds to "
                                        "conduct a bidirectional breadth "
                                        "first search to show the "
                                        "shortest path."),
                      dbc.Row([
                              dbc.Col([
                                  html.P(
                                      "Compound 1:"),
                                  dcc.Dropdown(
                                      id="filter1_dropdown",
                                      options=[{"label": i,
                                                "value": i, }
                                               for i in
                                               list(compounds_list)],
                                      value=compounds_list,
                                      multi=False,
                                      style={"width": "100%"},), ],
                                      width=6),
                              dbc.Col([
                                  html.P(
                                      "Compound 2:"),
                                  dcc.Dropdown(
                                      id="filter2_dropdown",
                                      options=[{"label": i,
                                                "value": i, }
                                               for i in
                                               list(compounds_list)],
                                      value=compounds_list,
                                      multi=False,
                                      style={"width": "100%"},)], width=5),
                              dbc.Col(dbc.Button("i", id='bfs_help',
                                                 outline=True,
                                                 color="primary",
                                                 size="sm"), width=1)]),
                      html.H5(
                          "Gene Delete:"),
                      dcc.ConfirmDialog(id="gene_dialog",
                                        message="Select genes to delete "
                                        "the gene from the model. If all "
                                        "of the genes selected constitute "
                                        "all of the genes associated with "
                                        "a reaction. That reaction will "
                                        "be removed from the analysis."),
                      dbc.Row([dbc.Col(
                              dcc.Dropdown(

                                  id="delete_dropdown",
                                  options=[{"label": i,
                                            "value": i, }
                                           for i in list(gene_content)],
                                  value=None,
                                  multi=True,
                                  style={"width": "100%"}), width=11),
                          dbc.Col(dbc.Button("i", id='gene_help',
                                             outline=True, color="primary",
                                             size="sm"), width=1)
                      ]),
                      html.H5(
                          "Flux Analysis:"),
                      dcc.ConfirmDialog(id="fba_dialog",
                                        message="Choose a reaction to "
                                        "optimize for flux balance "
                                        "analysis. Reactions carrying "
                                        "flux will be visualized with "
                                        "a blue color."),
                      dbc.Row([
                              dbc.Col(dcc.Dropdown(
                                  id="fba_dropdown",
                                  options=[{"label": i,
                                            "value": i, }
                                           for i in list(rxns_full)],
                                  value=rxns,
                                  multi=False,
                                  style={"width": "100%"},), width=11),
                              dbc.Col(dbc.Button("i", id='fba_help',
                                                 outline=True,
                                                 color="primary",
                                                 size="sm"), width=1)]),
                      dbc.Alert(id="fba-data",
                                children="Select a reaction "
                                "to see the flux here",
                                color="secondary",
                                style={"width": "75%"}),
                      dbc.Col([
                              html.Div([
                                  html.Button("Submit",
                                              id="btn_sub")]),
                              html.Div([
                                  html.Button("Download CSV",
                                              id="btn_tsv"),
                                  dcc.Download(
                                      id="download_r"),
                                 dcc.Download(
                                     id="download_c"),
                              ]),
                              html.Div([
                                  html.Button(
                                      "Download png", id="btn-get-png"),
                                  dcc.Download(
                                      id="downloadpng"),
                              ]),
                              ]),
                  ]),
                  )
    simulate_stats = dcc.Tab(
        label='Model Stats',
        value='Statistics',
        children=html.Div(
            className='about-tab',
            children=[
                html.H4(
                    className='Statistics',
                    children='Model Statistics:'),
                html.P("""
                    Below shows some basic information for the
                    model.
                    """),
                dbc.Alert(
                    id="model-stats",
                    children=content_model,
                    color="white",
                    style={
                        "width": "75%"}),
            ]),
    )
    s_con = dcc.Tab(label='Constraints', value='constraints',
                   children=html.Div(className='about-tab', children=[
                         dbc.Row([html.P(
                             """
                    Select Load constraints to load the reaction
                    constraints. To modify, select update constraints.
                    To export changes as a new limits file, select save
                    constraints.
                    """
                         ),
                             dbc.Alert(id="constraints",
                                       children="display exchange here",
                                       color="secondary",),
                             html.Div([
                                      dbc.Row([
                                          dbc.Col([
                                              html.Button(
                                                  "Load constraints",
                                                  id="constraint-load")
                                          ]),
                                          dbc.Col([
                                              html.Button(
                                                  "Update constraints",
                                                  id="constraint-modify")
                                          ]),
                                          dbc.Col([
                                              html.Button(
                                                  "Save constraints",
                                                  id="constraint-save")
                                          ])]), ]),
                             dbc.Alert(id="constraint-confirm",
                                       children="This will confirm that "
                                       "the constraint was modified",
                                       color="secondary",),
                             dbc.Alert(id="constraint-save-confirm",
                                       children="confirm that constraint "
                                       "is saved",
                                       color="secondary",),
                         ]),
                   ]),
                   )
    s_ex = dcc.Tab(label='Exchange', value='what-is',
                   children=html.Div(className='about-tab', children=[
                         dbc.Row([html.P(
                             """
                    Select Load exchange to load the exchange
                    constraints. To modify, select update exchange.
                    To export changes as a new media, select save
                    exchange.
                    """
                         ),
                             dbc.Alert(id="exchange",
                                       children="display exchange here",
                                       color="secondary",),
                             html.Div([
                                      dbc.Row([
                                          dbc.Col([
                                              html.Button(
                                                  "Load exchange",
                                                  id="exchange-load")
                                          ]),
                                          dbc.Col([
                                              html.Button(
                                                  "Update exchange",
                                                  id="exchange-modify")
                                          ]),
                                          dbc.Col([
                                              html.Button(
                                                  "Save exchange",
                                                  id="exchange-save")
                                          ])]), ]),
                             dbc.Alert(id="exchange-confirm",
                                       children="This will confirm that "
                                       "the exchange was modified",
                                       color="secondary",),
                             dbc.Alert(id="exchange-save-confirm",
                                       children="confirm that exchange "
                                       "is saved",
                                       color="secondary",),
                         ]),
                   ]),
                   )

    simulate_tabs = dcc.Tabs(id="model-tabs", value="model", children=[
        s_c, s_ex, s_con, simulate_stats,
    ], style={'display': 'inline-block', 'width': '100%'})

    s_n = cyto.Cytoscape(id='net',
                         layout={'name': 'cose',
                                 'animate': True,
                                 'animationDuration': 1000},
                         style={'width': '100%',
                                'height': '600px'},
                         elements=nodes + edges,
                         stylesheet=[{'selector': 'node',
                                      'style': {'background-color':
                                                'data(col)',
                                                'label': 'data(label)'}},
                                     {"selector": "edge",
                                      "style": {"width": 1,
                                                "curve-style": "bezier",
                                                "target-arrow-shape":
                                                "triangle"}},
                                     {"selector": "[flux > 1e-10]",
                                     "style": {"line-color": "blue",
                                               "target-arrow-color":
                                               "blue",
                                               "width": 1,
                                               "curve-style": "bezier",
                                               "target-arrow-shape":
                                               "triangle",
                                               "line-dash-pattern": [6, 3],
                                               "line-dash-offset": 24}}, ],
                            minZoom=0.06)

    
    # Define the main body of the app
    body_layout = dbc.Container(
        [dbc.Row([
            dbc.Row([
                dbc.Col(
                    html.Div(id='model-control-tabs',className='control-tabs',
                        style={'width': '69%', 'font-size': '50%', 'height': '50%'}, 
                        children=[simulate_tabs])
                    ),
                dbc.Col([s_n,
                    dbc.Row([dbc.Alert(id="node-data",children="Click on a node to see its details", color="secondary"),
                        ]), 
                    ]),
                ]),
            dbc.Row([dcc.Markdown('''\\* Data analysis carried out for 
                demonstration of data visualisation purposes only.''')], style={"fontSize": 11, "color": "gray"},), ],
            style={"marginTop": 20},)],

        # take up the whole screen
        fluid = True, style={"height": "100vh"} 
        )


    layout = html.Div([navbar, body_layout])
    return layout

# layout definition


