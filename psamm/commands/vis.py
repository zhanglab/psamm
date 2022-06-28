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
from dash.dependencies import Input, Output, State
import dash_cytoscape as cyto
import re
import pandas as pd
import webbrowser
from psamm.commands.chargecheck import charge_check
from psamm.commands.formulacheck import formula_check
from psamm.datasource.native import ModelWriter, NativeModel
from psamm.datasource.entry import (DictCompoundEntry as CompoundEntry,
                                    DictReactionEntry as ReactionEntry,
                                    DictCompartmentEntry as CompartmentEntry)
from decimal import *
import yaml
from psamm.datasource.native import parse_reaction_equation_string
from psamm.datasource.reaction import parse_reaction

logger = logging.getLogger(__name__)

REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'

EPSILON = 1e-5


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class InteractiveCommand(MetabolicMixin, SolverCommandMixin,
                         Command, FilePrefixAppendAction):
    """Generate graphical representations of the model"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method', type=str, default='fpp',
            help='Compound pair prediction method',
            choices=['fpp', 'no-fpp'])
        parser.add_argument(
            '--exclude', metavar='reaction', type=convert_to_unicode,
            default=[], action=FilePrefixAppendAction,
            help='Reaction(s) to exclude from metabolite pair prediction')
        parser.add_argument(
            '--subset', type=argparse.FileType('rU'), default=None,
            help='File containing a list of reactions IDs or compound IDs '
                 '(with compartment). This file determines which reactions '
                 'will be visualized')
        parser.add_argument(
            '--compartment', action='store_true',
            help='Include compartments in final visualization')
        parser.add_argument(
            '--output', type=str,
            help='Full path for visualization outputs, including prefix of '
                 'output files. '
                 'e.g. "--output <output-directory-path>/output-prefix"')
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '--fba', type=argparse.FileType('rU'), default=None,
            help='File containing fba reaction flux')
        group.add_argument(
            '--fva', type=argparse.FileType('rU'), default=None,
            help='File containing fva reaction flux')
        super(InteractiveCommand, cls).init_parser(parser)

    def run(self):

        # enable svg export
        cyto.load_extra_layouts()

        el = 'C'
        pathway = 'All'

        pathway_list, rxn_set = get_pathway_list(self._model, "All")
        compounds_list = get_compounds_list(self._model)
        rxns = list(rxn_set)
        rxns_full = []

        content_model = []
        count = 0
        rxns_with_genes = 0
        gene_content = set()
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
            unbalanced, count, unchecked, exclude, unbalanced_list, unbalance_dict = charge_check(
                self, 1e-6)
        content_model.append(html.H5("Chargecheck results:"))
        content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
        content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
        content_model.append(
            html.P("Excluded reactions: " + str(len(exclude))))

        with HiddenPrints():
            unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(
                self)
        content_model.append(html.H5("Formulacheck results:"))
        content_model.append(html.P("Unblanced reactions: " + str(unbalanced)))
        content_model.append(html.P("Unchecked reactions: " + str(unchecked)))
        content_model.append(
            html.P("Excluded reactions: " + str(len(exclude))))

        # nodes, edges = build_network(nm, rxn_set, network)
        # initialize an empty list. the full is messy
        nodes, edges = [], []

        # Initialize the app
        self._app = dash.Dash(
            __name__, external_stylesheets=[
                dbc.themes.BOOTSTRAP])
        server = self._app.server
        if self._app is not None and hasattr(self, "callbacks"):
            self.callbacks(self._app)

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
        navbar = dbc.NavbarSimple(
            children=[
                dbc.NavItem(
                    dbc.NavLink(
                        "Source Code",
                        href="https://github.com/cpowers11060/metavis",
                    )
                ),
                dbc.NavItem(
                    dbc.NavLink(
                        "Psamm documentation",
                        href="https://psamm.readthedocs.io",
                    )
                ),
                dbc.NavItem(
                    dbc.NavLink(
                        "Psamm publication",
                        href="https://journals.plos.org/ploscompbiol/"
                             "article?id=10.1371/journal.pcbi.1004732",
                    )
                ), ],

            brand="Psamm web-based visualization of metabolic models",
            brand_href="#",
            color="dark",
            dark=True,
        )

        # Define the body of the app
        simtab = dbc.Row([
                    dbc.Col([
                        html.Div(id='model-control-tabs',
                                 className='control-tabs',
                                 children=[
                                     dcc.Tabs(id="model-tabs", value="what is", children=[
                                         dcc.Tab(label='Controls', value='Controls',
                                                 children=html.Div(className='control-tab', children=[
                                                     html.H5(
                                                         "Pathways:"),
                                                     dcc.ConfirmDialog(id="pathway_dialog",
                                                                       message="Select the pathways you want to display "
                                                                       "from the metabolic model. Multiple can be "
                                                                       "selected. Reactions can be added to this display "
                                                                       "from the Reactions dropdown."),
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
                                                         dbc.Col(dbc.Button("i", id='pathway_help', outline=True,
                                                                            color="primary", size="sm"), width=1), ]),


                                                     html.H5(
                                                         "Reactions:"),
                                                     dcc.ConfirmDialog(id="reaction_dialog",
                                                                       message="Select the reactions you want "
                                                                       "to display from the metabolic model. "
                                                                       "Multiple can be selected."),
                                                     dbc.Row([
                                                     dbc.Col(
                                                         dcc.Dropdown(
                                                             id="reaction_dropdown",
                                                             options=[
                                                                 {"label": i.id,
                                                                  "value": i.id, }
                                                                 for i in list(self._model.reactions)],
                                                             value=pathway_list[0],
                                                             multi=True,
                                                             placeholder="",
                                                             style={"width": "100%"},), width=11),
                                                         dbc.Col(dbc.Button("i", id='reaction_help',
                                                                            outline=True, color="primary",
                                                                            size="sm"), width=1),
                                                     ]),
                                                     html.H5(
                                                         "Element Transfer Networks:"),
                                                     dcc.ConfirmDialog(id="element_dialog",
                                                                       message="Choose the chemical element "
                                                                       "for the network to be based on. Default "
                                                                       "is Carbon (C). Options are Carbon (C), "
                                                                       "Nitrogen (N), Sulfur (S), or Phosphorous (P)"),
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
                                                         dbc.Col(dbc.Button("i", id='element_help', outline=True,
                                                                            color="primary", size="sm"), width=1), ]),
                                                     html.H5(
                                                         "Compound Search:"),
                                                     dcc.ConfirmDialog(id="compound_dialog",
                                                                       message="Select a compound to show all "
                                                                       "reactions containing that compound. Note "
                                                                       "that this combines with the pathways and "
                                                                       "reactions options above. If you want to "
                                                                       "see all reactions with that comound, "
                                                                       "select All for pathways."),
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
                                                         style={"width": "100%"}, ), width=11),
                                                         dbc.Col(
                                                         dbc.Col(dbc.Button("i", id='compound_help',
                                                                            outline=True, color="primary",
                                                                            size="sm"), width=1),
                                                         ), ]),
                                                     html.H5(
                                                         "Path Search:"),
                                                     dcc.ConfirmDialog(id="bfs_dialog",
                                                                       message="Select two compounds to conduct "
                                                                       "a bidirectional breadth first search to "
                                                                       "show the shortest path."),
                                                     dbc.Row([
                                                         dbc.Col([
                                                             html.P(
                                                                 "Compound 1:"),
                                                             dcc.Dropdown(
                                                                 id="filter1_dropdown",
                                                                 options=[{"label": i,
                                                                           "value": i, }
                                                                      for i in list(compounds_list)],
                                                                 value=compounds_list,
                                                                 multi=False,
                                                                 style={"width": "100%"},),], width=6),
                                                         dbc.Col([
                                                             html.P(
                                                                 "Compound 2:"),
                                                             dcc.Dropdown(
                                                                 id="filter2_dropdown",
                                                                 options=[{"label": i,
                                                                           "value": i, }
                                                                          for i in list(compounds_list)],
                                                                 value=compounds_list,
                                                                 multi=False,
                                                                 style={"width": "100%"},)], width=5),
                                                         dbc.Col(dbc.Button("i", id='bfs_help',
                                                                            outline=True, color="primary",
                                                                            size="sm"), width=1)]),
                                                     html.H5(
                                                         "Gene Delete:"),
                                                     dcc.ConfirmDialog(id="gene_dialog",
                                                                       message="Select genes to delete the "
                                                                       "gene from the model. If all of the "
                                                                       "genes selected constitute all of the "
                                                                       "genes associated with a reaction. "
                                                                       "That reaction will be removed from "
                                                                       "the analysis."),
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
                                                                       message="Choose a reaction to optimize "
                                                                       "for flux balance analysis. Reactions "
                                                                       "carrying flux will be visualized with "
                                                                       "a blue color."),
                                                     dbc.Row([
                                                         dbc.Col(dcc.Dropdown(
                                                         id="fba_dropdown",
                                                         options=[{"label": i,
                                                                   "value": i, }
                                                              for i in list(rxns_full)],
                                                         value=rxns,
                                                         multi=False,
                                                         style={"width": "100%"},), width = 11),
                                                         dbc.Col(dbc.Button("i", id='fba_help',
                                                                            outline=True, color="primary",
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
                                                                 id="download"),
                                                         ]),
                                                         html.Div([
                                                             html.Button(
                                                                 "Download png", id="btn-get-png"),
                                                             dcc.Download(
                                                                 id="downloadpng"),
                                                         ]),
                                                     ]),
                                                 ]),
                                                 ),
                                         dcc.Tab(label='Model Stats', value='Statistics',
                                                 children=html.Div(className='about-tab', children=[
                                                     html.H4(
                                                         className='Statistics',
                                                         children='Model Statistics:'
                                                     ),

                                                     html.P(
                                                         """
                            Below shows some basic information for the
                            model.
                            """
                                                     ),
                                                     dbc.Alert(id="model-stats",
                                                               children=content_model,
                                                               color="white",
                                                               style={"width": "75%"}),
                                                 ]),
                                                 ),
                                         dcc.Tab(label='Exchange', value='what-is',
                                                 children=html.Div(className='about-tab', children=[
                                                     dbc.Row([html.P(
                                                         """
                            Select Load exchange to load the exchange constraints.
                            To modify, select update exchange.
                            To export changes as a new media, select save exchange.
                            """
                                                     ),
                                                         dbc.Alert(id="exchange",
                                                                   children="display exchange here",
                                                                   color="secondary",),
                                                         html.Div([
                                                             dbc.Row([
                                                                 dbc.Col([
                                                                     html.Button(
                                                                         "Load exchange", id="exchange-load")
                                                                 ]),
                                                                 dbc.Col([
                                                                     html.Button(
                                                                         "Update exchange", id="exchange-modify")
                                                                 ]),
                                                                 dbc.Col([
                                                                     html.Button(
                                                                         "Save exchange", id="exchange-save")
                                                                 ])]), ]),
                                                         dbc.Alert(id="exchange-confirm",
                                                                   children="This will confirm that the exchange was modified",
                                                                   color="secondary",),
                                                         dbc.Alert(id="exchange-save-confirm",
                                                                   children="confirm that exchange is saved",
                                                                   color="secondary",),
                                                     ]),
                                                 ]),
                                                 ),

                                     ], style={'display': 'inline-block', 'width': '100%'}),
                                     # sm=12,
                                     # md=4,
                                 ]), ]),
                    dbc.Col([
                        cyto.Cytoscape(id='net',
                                       layout={'name': 'cose',
                                               'animate': True,
                                               'animationDuration': 1000},
                                       style={
                                           'width': '100%', 'height': '600px'},
                            elements=nodes + edges,
                            stylesheet=[{
                                'selector': 'node',
                                'style': {
                                    'background-color': 'data(col)',
                                    'label': 'data(label)'}},
                                           {
                                "selector": "edge",
                                "style": {
                                    "width": 1,
                                    "curve-style": "bezier",
                                    "target-arrow-shape": "triangle"}},
                                           {
                                "selector": "[flux > 1e-10]",
                                "style": {
                                    "line-color": "blue",
                                    "target-arrow-color": "blue",
                                    "width": 1,
                                    "curve-style": "bezier",
                                    "target-arrow-shape": "triangle",
                                    "line-dash-pattern": [6, 3],
                                    "line-dash-offset": 24
                                }
                            },

                                       ],
                            minZoom=0.06
                        ),
                        dbc.Row([
                            dbc.Alert(
                                id="node-data",
                                children="Click on a node to see its details here",
                                color="secondary",
                            ),
                        ]), ]),
                    dbc.Row([dcc.Markdown(
                        """
                \\* Data analysis carried out for demonstration of data visualisation purposes only.
                """
                    )],
                        style={"fontSize": 11, "color": "gray"},),

                ], style={"marginTop": 20},)
        body_layout = dbc.Container(
            [dbc.Row([
                dcc.Markdown(
                    """
                -----
                ### Filter / Explore metabolic models
                Use these filters to highlight reactions and compounds associated with different reactions.

                -----
                """
                ),
                dbc.Row([
                    html.Div(id='modelling-tabs', className='control-tabs', children=[
                        dcc.Tabs(id="parent-tabs", value="what is", children=[
                            dcc.Tab(label='Simulate', value='sim', children=[
                                simtab ]),
                            dcc.Tab(label='Curate', value='cur', children=[
                                html.Div(id='curate-control-tabs', className='curate-tabs', children=[
                                    dcc.Tabs(id="curate-subtabs", value="what is", children=[
                                            dcc.Tab(label='Add', value='add',
                                                    children=html.Div(className='add-tab', children=[
                                                        dbc.Row([
                                                            dbc.Col([
                                                                html.H5(
                                                                    "Add a new reaction"),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "ID:"),
                                                                        dcc.Input(
                                                                            id="r_id_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Name:"),
                                                                        dcc.Input(
                                                                            id="r_name_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Equation:"),
                                                                        dcc.Input(
                                                                            id="r_eq_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Pathway:"),
                                                                        dcc.Input(
                                                                            id="r_pathway_input"),
                                                                    ], width=11),
                                                                ]),
                                                                html.Div([
                                                                    html.Button(
                                                                        "Save reaction", id="btn_rxn")
                                                                ]),
                                                                dbc.Alert(
                                                                    id="rxn-save-confirm")
                                                            ]),
                                                            dbc.Col([
                                                                html.H5(
                                                                    "Add a new compound"),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "ID:"),
                                                                        dcc.Input(
                                                                            id="c_id_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Name:"),
                                                                        dcc.Input(
                                                                            id="c_name_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Charge:"),
                                                                        dcc.Input(
                                                                            id="c_charge_input"),
                                                                    ], width=11),
                                                                ]),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        html.P(
                                                                            "Formula:"),
                                                                        dcc.Input(
                                                                            id="c_form_input"),
                                                                    ], width=11),
                                                                ]),
                                                                html.Div([
                                                                    html.Button(
                                                                        "Save compound", id="btn_cpd")
                                                                ]),
                                                                dbc.Alert(
                                                                    id="cpd-save-confirm")
                                                            ]),
                                                        ]),
                                                    ]),),
                                            dcc.Tab(label='Edit', value='edit',
                                                    children=html.Div(className='control-tab', children=[
                                                        dbc.Row([
                                                            dbc.Col([
                                                                html.H5(
                                                                    "Reaction to edit:"),
                                                                dbc.Row(
                                                                    [dcc.Dropdown(
                                                                        id="reaction_dropdown_curate",
                                                                        options=[
                                                                            {"label": i,
                                                                             "value": i, }
                                                                            for i in list(rxns_full)],
                                                                        value='',
                                                                        multi=False,
                                                                        placeholder="",
                                                                        style={"width": "100%"},), ]),
                                                                html.H5(
                                                                    "Compound to edit:"),
                                                                dbc.Row(
                                                                    [dcc.Dropdown(
                                                                        id="compound_dropdown_curate",
                                                                        options=[
                                                                            {"label": i,
                                                                             "value": i, }
                                                                            for i in list(compounds_list)],
                                                                        value='',
                                                                        multi=False,
                                                                        placeholder="",
                                                                        style={"width": "100%"},), ]),
                                                                html.Div([
                                                                    html.Button(
                                                                        "Confirm Edit", id="btn_update")
                                                                ]),
                                                                dbc.Alert(
                                                                    id="edit-data",
                                                                    children="Click on a node to see its details here",
                                                                    color="secondary",
                                                                ),
                                                                html.Div([
                                                                    html.Button(
                                                                        "Save", id="btn_save")
                                                                ]),
                                                                dbc.Alert(
                                                                    id="edit-confirm",
                                                                    children="Click save to confirm changes to model",
                                                                    color="secondary",
                                                                ),
                                                            ]),
                                                            dbc.Col([
                                                                html.H5(
                                                                    "Pathways:"),
                                                                dbc.Row(
                                                                    [dcc.Dropdown(
                                                                        id="pathways_dropdown_curate",
                                                                        options=[
                                                                            {"label": i,
                                                                             "value": i, }
                                                                            for i in list(pathway_list)],
                                                                        value=pathway_list[0],
                                                                        multi=True,
                                                                        placeholder="",
                                                                        style={"width": "100%"},), ]),
                                                                html.Div([
                                                                    html.Button(
                                                                        "Submit", id="btn_update_cur_net")
                                                                ]),
                                                                cyto.Cytoscape(id='net_curate',
                                                                               layout={'name': 'cose',
                                                                                       'animate': True,
                                                                                       'animationDuration': 1000},
                                                                               style={
                                                                                   'width': '100%', 'height': '600px'},
                                                                               elements=nodes + edges,
                                                                               stylesheet=[{
                                                                                   'selector': 'node',
                                                                                   'style': {
                                                                                       'background-color': 'data(col)',
                                                                                       'label': 'data(label)'}},
                                                                                   {
                                                                                   "selector": "edge",
                                                                                   "style": {
                                                                                       "width": 1,
                                                                                       "curve-style": "bezier",
                                                                                       "target-arrow-shape": "triangle"}},
                                                                                   {
                                                                                   "selector": "[flux > 0]",
                                                                                   "style": {
                                                                                       "line-color": "blue",
                                                                                       "target-arrow-color": "blue",
                                                                                       "width": 1,
                                                                                       "curve-style": "bezier",
                                                                                       "target-arrow-shape": "triangle",
                                                                                       "line-dash-pattern": [6, 3],
                                                                                       "line-dash-offset": 24
                                                                                   }
                                                                               },

                                                                               ],
                                                                               minZoom=0.06
                                                                               ),
                                                            ]), ]),
                                                    ]),),
                                            dcc.Tab(label='Charge Check', value='chargeCheck',
                                                    children=html.Div(className='control-tab', children=[
                                                        dbc.Row(
                                                            [html.H5("Unbalanced Reactions:"),
                                                             dcc.Dropdown(
                                                                id="unbalanced_rxns",
                                                                options=[
                                                                    {"label": i,
                                                                     "value": i, }
                                                                    for i in list(unbalanced_list)],
                                                                value=pathway_list[0],
                                                                multi=False,
                                                                placeholder="",),
                                                                dbc.Alert(
                                                                id="chargeCheckCuration",
                                                                children="Select a reaction above to curate the charge balance",
                                                                color="secondary",
                                                            ),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        dbc.Alert(
                                                                            id="chargeLeftCuration",
                                                                            children="The left side of the reaction will display here",
                                                                            color="secondary",
                                                                        ),
                                                                        html.Div([
                                                                            html.Button(
                                                                                "Save Changes", id="btn_save_charge")
                                                                        ]),
                                                                    ]),
                                                                    dbc.Col([
                                                                        dbc.Alert(
                                                                            id="chargeRightCuration",
                                                                            children="The right side of the reaction will display here",
                                                                            color="secondary",
                                                                        ),
                                                                    ]),
                                                                ]),
                                                            ]),

                                                    ]),),
                                            dcc.Tab(label='Formula Check', value='formulaCheck',
                                                    children=html.Div(className='control-tab', children=[
                                                        dbc.Row(
                                                            [html.H5("Unbalanced Reactions:"),
                                                             dcc.Dropdown(
                                                                id="unbalanced_rxns_formula",
                                                                options=[
                                                                    {"label": i,
                                                                     "value": i, }
                                                                    for i in list(form_list)],
                                                                value=pathway_list[0],
                                                                multi=False,
                                                                placeholder="",),
                                                                dbc.Alert(
                                                                id="formCheckCuration",
                                                                children="Select a reaction above to curate the formula balance",
                                                                color="secondary",
                                                            ),
                                                                dbc.Row([
                                                                    dbc.Col([
                                                                        dbc.Alert(
                                                                            id="formLeftCuration",
                                                                            children="The left side of the reaction will display here",
                                                                            color="secondary",
                                                                        ),
                                                                        html.Div([
                                                                            html.Button(
                                                                                "Save Changes", id="btn_save_form")
                                                                        ]),
                                                                    ]),
                                                                    dbc.Col([
                                                                        dbc.Alert(
                                                                            id="formRightCuration",
                                                                            children="The right side of the reaction will display here",
                                                                            color="secondary",
                                                                        ),
                                                                    ]),
                                                                ]),
                                                            ]),
                                                    ]),),
                                    ]),

                                ]),
                                dbc.Row([
                                    dbc.Alert(
                                        id="save_confirmation",
                                        children="Press Save Model to save the model",
                                        color="seconday",
                                    ),
                                    dbc.Row([
                                        dbc.Col([
                                            html.Button(
                                                "Save Model", id="btn_save_model"),
                                        ]), ]), ],
                                    style={"fontSize": 16, "color": "gray"},),



                            ]),
                        ]),
                    ]), ]),
            ]),
            ])

        #
        #
        #         ],
        #         style={"marginTop": 20},
        # )

        self._app.layout = html.Div([navbar, body_layout])
        webbrowser.open_new("http://localhost:{}".format("8050"))
        self._app.run_server(debug=True)

    def callbacks(self, _app):
        @_app.callback(
            Output('pathway_dialog', 'displayed'),
            Input('pathway_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'pathway_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('reaction_dialog', 'displayed'),
            Input('reaction_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'reaction_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('gene_dialog', 'displayed'),
            Input('gene_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'gene_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('fba_dialog', 'displayed'),
            Input('fba_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'fba_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('bfs_dialog', 'displayed'),
            Input('bfs_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'bfs_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('compound_dialog', 'displayed'),
            Input('compound_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'compound_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('element_dialog', 'displayed'),
            Input('element_help', 'n_clicks')
        )
        def pathway_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'element_help.n_clicks':
                return True
            else:
                return False

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
                            f.write(
                                "{}\t{}\t{}\t{}\n".format(
                                    e[0], e[1], e[2], e[3]))
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
                                "EX_{}[{}]".format(e[0], e[1])) == False:
                            properties = {
                                "id": "EX_{}[{}]".format(
                                    e[0], e[1]), "equation": Reaction(
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
            if dash.callback_context.triggered[0]['prop_id'] == 'exchange-load.n_clicks':
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
                                          placeholder=float(self._model._exchange[e][2]), style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="upper{}".format(str(e.name)),
                                          type="number",
                                          placeholder=float(self._model._exchange[e][3]),
                                          style={'width': '85%'})], width=3),
                        ]),)
                contents.append(html.P("Add new constraint:"))
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
                                options=[{"label": x,
                                          "value": x, }
                                         for x in list(comp_list)],
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
            [Output("save_confirmation", "children"), ],
            [Input("btn_save_model", "n_clicks")],
            prevent_initial_call=True,)
        def save_model(clicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'btn_save_model.n_clicks':
                out_mw = ModelWriter()
                path = FilePathContext(self._args.model)
                with open('{}/model.yaml'.format(path), 'r') as f:
                    modelfile = yaml.safe_load(f)

                exported = []
                for i in modelfile['reactions']:
                    out_nm = NativeModel()
                    with open('{}/{}'.format(path, i['include']), 'r') as f:
                        reactionfile = yaml.safe_load(f)
                        for j in reactionfile:
                            exported.append(j['id'])
                            if j['id'] in self._model.reactions:
                                out_nm.reactions.add_entry(
                                    self._model.reactions[j['id']])
                        with open('{}/{}'.format(path, i['include']), "w") as f:
                            out_mw.write_reactions(f, out_nm.reactions)
                out_nm = NativeModel()
                for j in self._model.reactions:
                    if j.id not in exported:
                        out_nm.reactions.add_entry(self._model.reactions[j.id])
                with open('{}/{}'.format(path,
                          "reactions_curated.yaml"), "a+") as f:
                    out_mw.write_compounds(f, out_nm.reactions)

                exported = []
                for i in modelfile['compounds']:
                    out_nm = NativeModel()
                    with open('{}/{}'.format(path, i['include']), 'r') as f:
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
                return(["Model successfully saved"])
            else:
                return(["Press Save Model to save changes"])

        @_app.callback(
            [Output("formCheckCuration", "children"),
             Output("formLeftCuration", "children"),
             Output("formRightCuration", "children"),
             Output("unbalanced_rxns_formula", "options")],
            [Input("unbalanced_rxns_formula", "value"),
             Input("btn_save_form", "n_clicks")],
            [State("formLeftCuration", "children"),
             State("formRightCuration", "children")],
            prevent_initial_call=True,
        )
        def select_formulaCheck(rxn, formClick, leftdata, rightdata):
            changed_id = [p['prop_id']
                          for p in dash.callback_context.triggered][0]

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

            def compound_name(id):
                if id not in self._model.compounds:
                    return id
                return self._model.compounds[id].properties.get('name', id)

            if isinstance(rxn, str):
                if dash.callback_context.triggered[0]['prop_id'] == \
                        'btn_save_form.n_clicks':
                    left_out = []
                    eq_dict = {}
                    for i in leftdata:
                        if isinstance(i, dict):
                            left_out_temp = rec_props(
                                i['props']['children'], [])
                            if len(left_out_temp) > 0:
                                if left_out_temp[0] != '':
                                    left_out.append(left_out_temp)
                    for l in left_out:
                        for c in self._model.compounds:
                            if c.id.upper() == l[0].upper():
                                compound = c
                        compound.charge = int(l[3])
                        compound.formula = l[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(l[0]).in_compartment(
                            l[4])] = int(l[1]) * -1

                    right_out = []
                    for i in rightdata:
                        if isinstance(i, dict):
                            right_out_temp = rec_props(
                                i['props']['children'], [])
                            if len(right_out_temp) > 0:
                                if right_out_temp[0] != '':
                                    right_out.append(right_out_temp)
                    for r in right_out:
                        for c in self._model.compounds:
                            if c.id.upper() == r[0].upper():
                                compound = c
                        compound.charge = int(r[3])
                        compound.formula = r[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(r[0]).in_compartment(
                            r[4])] = int(r[1])

                    for r in self._model.reactions:
                        if r.id.upper() == rxn.upper():
                            reaction = r
                    new_eq = Reaction(
                        reaction.equation.__dict__['_direction'], eq_dict)
                    reaction.equation = new_eq
                    self._model.reactions.add_entry(reaction)

                with HiddenPrints():
                    unbalanced, count, unchecked, exclude, form_list, form_dict = formula_check(
                        self)
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
                if isinstance(rxn, str):
                    for i in self._model.reactions:
                        if i.id == rxn:
                            r = i
                    contents = []
                    contents.append(
                        html.H5("Reaction: {}".format(r.id))
                    )
                    contents.append(
                        html.H5("Equation: ")
                    )
                    contents.append(
                        html.P("{}".format(r.equation))
                    )

                    contents.append(
                        html.H5("Translated Equation:")
                    )
                    contents.append(
                        html.P(
                            "{}".format(
                                r.equation.translated_compounds(compound_name))))
                    try:
                        contents.append(
                            dbc.Row([
                                dbc.Col([
                                    html.P("Left side total: {}".format(form_dict[r.id][0])),
                                    html.P("Left side off by: {}".format(form_dict[r.id][1])),
                                ]),
                                dbc.Col([
                                    html.P("Right side total: {}".format(form_dict[r.id][2])),
                                    html.P("Right side off by: {}".format(form_dict[r.id][3])),
                                ]),
                            ]),

                        )
                    except KeyError:
                        contents = [
                            "{} was successfully balanced! Select a reaction above to curate the charge balance".format(
                                r.id)]
                        left_contents = [
                            "The left side of the reaction will display here"]
                        right_contents = [
                            "The right side of the reaction will display here"]
                        return contents, left_contents, right_contents, [
                            {"label": i, "value": i, } for i in list(form_list)]

                    left_contents = []
                    # for i in r.equation.left:
                    #     left_contents.append(html.H5("test"))

                    left_contents.append(html.H5("Left:"))
                    left_contents.append(
                        dbc.Row([
                            dbc.Col([
                                html.P("ID")
                            ], width=3),
                            dbc.Col([
                                html.P("Stoich.")
                            ], width=2),
                            dbc.Col([
                                html.P("Formula")
                            ], width=3),
                            dbc.Col([
                                html.P("Charge")
                            ], width=2),
                            dbc.Col([
                                html.P("Comp")
                            ], width=2),
                        ]))
                    cur_id = 0
                    for i in r.equation.left:
                        left_contents.append(
                            dbc.Row([
                                dbc.Col([
                                    dcc.Dropdown(
                                        id="drop{}".format(str(cur_id)),
                                        options=[{"label": x,
                                                  "value": x, }
                                                 for x in list(c_list)],
                                        placeholder=i[0].name,
                                        multi=False
                                    )], width=3),
                                dbc.Col([
                                    dcc.Input(id="stoich{}".format(str(cur_id)),
                                              type="number",
                                              placeholder=str(i[1]), style={'width': '75%'})], width=2),
                                dbc.Col([
                                    dcc.Input(id="formula{}".format(str(cur_id)),
                                              type="text",
                                              placeholder=str(formula_dict[i[0].name]), style={'width': '85%'})], width=3),
                                dbc.Col([
                                    dcc.Input(id="charge{}".format(str(cur_id)),
                                              type="number",
                                              placeholder=str(charge_dict[i[0].name]),
                                              style={'width': '75%'})], width=2),
                                dbc.Col([
                                    dcc.Dropdown(
                                        id="compartment{}".format(str(cur_id)),
                                        options=[{"label": x,
                                                  "value": x, }
                                                 for x in list(comp_list)],
                                        placeholder=i[0].compartment,
                                        multi=False
                                    )
                                ], width=2),
                            ]),
                        )
                        cur_id += 1
                    left_contents.append(html.P("Add to left:"))
                    left_contents.append(
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    id="drop{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(c_list)],
                                    placeholder='',
                                    multi=False
                                )], width=3),
                            dbc.Col([
                                dcc.Input(id="stoich{}".format(str(cur_id)),
                                          type="number",
                                          placeholder="0", style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder='', style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder="0",
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Dropdown(
                                    id="compartment{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(comp_list)],
                                    placeholder='',
                                    multi=False
                                )], width=2),

                        ]),
                    )

                    right_contents = []

                    right_contents.append(html.H5("Right:"))
                    right_contents.append(
                        dbc.Row([
                            dbc.Col([
                                html.P("ID")
                            ], width=3),
                            dbc.Col([
                                html.P("Stoich.")
                            ], width=2),
                            dbc.Col([
                                html.P("Formula")
                            ], width=3),
                            dbc.Col([
                                html.P("Charge")
                            ], width=2),
                            dbc.Col([
                                html.P("Comp")
                            ], width=2),
                        ]))
                    cur_id = 0
                    for i in r.equation.right:
                        right_contents.append(
                            dbc.Row([
                                dbc.Col([
                                    dcc.Dropdown(
                                        id="drop{}".format(str(cur_id)),
                                        options=[{"label": x,
                                                  "value": x, }
                                                 for x in list(c_list)],
                                        placeholder=i[0].name,
                                        multi=False
                                    )], width=3),
                                dbc.Col([
                                    dcc.Input(id="stoich{}".format(str(cur_id)),
                                              type="number",
                                              placeholder=str(i[1]), style={'width': '75%'})], width=2),
                                dbc.Col([
                                    dcc.Input(id="formula{}".format(str(cur_id)),
                                              type="text",
                                              placeholder=str(formula_dict[i[0].name]), style={'width': '85%'})], width=3),
                                dbc.Col([
                                    dcc.Input(id="charge{}".format(str(cur_id)),
                                              type="number",
                                              placeholder=str(charge_dict[i[0].name]),
                                              style={'width': '75%'})], width=2),
                                dbc.Col([
                                    dcc.Dropdown(
                                        id="compartment{}".format(str(cur_id)),
                                        options=[{"label": x,
                                                  "value": x, }
                                                 for x in list(comp_list)],
                                        placeholder=i[0].compartment,
                                        multi=False
                                    )
                                ], width=2),
                            ]),
                        )
                        cur_id += 1
                    right_contents.append(html.P("Add to right:"))
                    right_contents.append(
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    id="drop{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(c_list)],
                                    placeholder='',
                                    multi=False
                                )], width=3),
                            dbc.Col([
                                dcc.Input(id="stoich{}".format(str(cur_id)),
                                          type="number",
                                          placeholder="0", style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder='', style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder="0",
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Dropdown(
                                    id="compartment{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(comp_list)],
                                    placeholder='',
                                    multi=False
                                )], width=2),
                        ]),
                    )
                    return contents, left_contents, right_contents, [
                        {"label": i, "value": i, } for i in list(form_list)]
                else:
                    contents = [
                        "Select a reaction above to curate the formula balance"]
                    left_contents = [
                        "The left side of the reaction will display here"]
                    right_contents = [
                        "The right side of the reaction will display here"]

        @_app.callback(
            [Output("chargeCheckCuration", "children"),
             Output("chargeLeftCuration", "children"),
             Output("chargeRightCuration", "children"),
             Output("unbalanced_rxns", "options")],
            [Input("unbalanced_rxns", "value"),
             Input("btn_save_charge", "n_clicks")],
            [State("chargeLeftCuration", "children"),
             State("chargeRightCuration", "children")],
            prevent_initial_call=True,
        )
        def select_chargeCheck(rxn, chargeClick, leftdata, rightdata):
            changed_id = [p['prop_id']
                          for p in dash.callback_context.triggered][0]

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

            def compound_name(id):
                if id not in self._model.compounds:
                    return id
                return self._model.compounds[id].properties.get('name', id)
            if isinstance(rxn, str):
                if dash.callback_context.triggered[0]['prop_id'] == 'btn_save_charge.n_clicks':
                    left_out = []
                    eq_dict = {}
                    for i in leftdata:
                        if isinstance(i, dict):
                            left_out_temp = rec_props(
                                i['props']['children'], [])
                            if len(left_out_temp) > 0:
                                if left_out_temp[0] != '':
                                    left_out.append(left_out_temp)
                    for l in left_out:
                        for c in self._model.compounds:
                            if c.id.upper() == l[0].upper():
                                compound = c
                        compound.charge = int(l[3])
                        compound.formula = l[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(l[0]).in_compartment(
                            l[4])] = int(l[1]) * -1

                    right_out = []
                    for i in rightdata:
                        if isinstance(i, dict):
                            right_out_temp = rec_props(
                                i['props']['children'], [])
                            if len(right_out_temp) > 0:
                                if right_out_temp[0] != '':
                                    right_out.append(right_out_temp)
                    for r in right_out:
                        for c in self._model.compounds:
                            if c.id.upper() == r[0].upper():
                                compound = c
                        compound.charge = int(r[3])
                        compound.formula = r[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(r[0]).in_compartment(
                            r[4])] = int(r[1])

                    for r in self._model.reactions:
                        if r.id.upper() == rxn.upper():
                            reaction = r
                    new_eq = Reaction(
                        reaction.equation.__dict__['_direction'], eq_dict)
                    reaction.equation = new_eq
                    self._model.reactions.add_entry(reaction)

            with HiddenPrints():
                unbalanced, count, unchecked, exclude, unbalanced_list, unbalance_dict = charge_check(
                    self, 1e-6)
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
            if isinstance(rxn, str):
                for i in self._model.reactions:
                    if i.id == rxn:
                        r = i
                contents = []
                contents.append(
                    html.H5("Reaction: {}".format(r.id))
                )
                contents.append(
                    html.H5("Equation: ")
                )
                contents.append(
                    html.P("{}".format(r.equation))
                )

                contents.append(
                    html.H5("Translated Equation:")
                )
                contents.append(
                    html.P(
                        "{}".format(
                            r.equation.translated_compounds(compound_name))))
                try:
                    contents.append(
                        html.P("Charge off by {}".format(unbalance_dict[r.id]))

                    )
                except KeyError:
                    contents = [
                        "{} was successfully balanced! Select a reaction above to curate the charge balance".format(
                            r.id)]
                    left_contents = [
                        "The left side of the reaction will display here"]
                    right_contents = [
                        "The right side of the reaction will display here"]
                    return contents, left_contents, right_contents, [
                        {"label": i, "value": i, } for i in list(unbalanced_list)]

                left_contents = []
                # for i in r.equation.left:
                #     left_contents.append(html.H5("test"))

                left_contents.append(html.H5("Left:"))
                left_contents.append(
                    dbc.Row([
                        dbc.Col([
                            html.P("ID")
                        ], width=3),
                        dbc.Col([
                            html.P("Stoich.")
                        ], width=2),
                        dbc.Col([
                            html.P("Formula")
                        ], width=3),
                        dbc.Col([
                            html.P("Charge")
                        ], width=2),
                        dbc.Col([
                            html.P("Comp")
                        ], width=2),
                    ]))
                cur_id = 0
                for i in r.equation.left:
                    left_contents.append(
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    id="drop{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(c_list)],
                                    placeholder=i[0].name,
                                    multi=False
                                )], width=3),
                            dbc.Col([
                                dcc.Input(id="stoich{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(i[1]), style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder=str(formula_dict[i[0].name]), style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(charge_dict[i[0].name]),
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Dropdown(
                                    id="compartment{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(comp_list)],
                                    placeholder=i[0].compartment,
                                    multi=False
                                )
                            ], width=2),
                        ]),
                    )
                    cur_id += 1
                left_contents.append(html.P("Add to left:"))
                left_contents.append(
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                id="drop{}".format(str(cur_id)),
                                options=[{"label": x,
                                          "value": x, }
                                         for x in list(c_list)],
                                placeholder='',
                                multi=False
                            )], width=3),
                        dbc.Col([
                            dcc.Input(id="stoich{}".format(str(cur_id)),
                                      type="number",
                                      placeholder="0", style={'width': '75%'})], width=2),
                        dbc.Col([
                            dcc.Input(id="formula{}".format(str(cur_id)),
                                      type="text",
                                      placeholder='', style={'width': '85%'})], width=3),
                        dbc.Col([
                            dcc.Input(id="charge{}".format(str(cur_id)),
                                      type="number",
                                      placeholder="0",
                                      style={'width': '75%'})], width=2),
                        dbc.Col([
                            dcc.Dropdown(
                                id="compartment{}".format(str(cur_id)),
                                options=[{"label": x,
                                          "value": x, }
                                         for x in list(comp_list)],
                                placeholder='',
                                multi=False
                            )], width=2),

                    ]),
                )

                right_contents = []

                right_contents.append(html.H5("Right:"))
                right_contents.append(
                    dbc.Row([
                        dbc.Col([
                            html.P("ID")
                        ], width=3),
                        dbc.Col([
                            html.P("Stoich.")
                        ], width=2),
                        dbc.Col([
                            html.P("Formula")
                        ], width=3),
                        dbc.Col([
                            html.P("Charge")
                        ], width=2),
                        dbc.Col([
                            html.P("Comp")
                        ], width=2),
                    ]))
                cur_id = 0
                for i in r.equation.right:
                    right_contents.append(
                        dbc.Row([
                            dbc.Col([
                                dcc.Dropdown(
                                    id="drop{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(c_list)],
                                    placeholder=i[0].name,
                                    multi=False
                                )], width=3),
                            dbc.Col([
                                dcc.Input(id="stoich{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(i[1]), style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder=str(formula_dict[i[0].name]), style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(charge_dict[i[0].name]),
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Dropdown(
                                    id="compartment{}".format(str(cur_id)),
                                    options=[{"label": x,
                                              "value": x, }
                                             for x in list(comp_list)],
                                    placeholder=i[0].compartment,
                                    multi=False
                                )
                            ], width=2),
                        ]),
                    )
                    cur_id += 1
                right_contents.append(html.P("Add to right:"))
                right_contents.append(
                    dbc.Row([
                        dbc.Col([
                            dcc.Dropdown(
                                id="drop{}".format(str(cur_id)),
                                options=[{"label": x,
                                          "value": x, }
                                         for x in list(c_list)],
                                placeholder='',
                                multi=False
                            )], width=3),
                        dbc.Col([
                            dcc.Input(id="stoich{}".format(str(cur_id)),
                                      type="number",
                                      placeholder="0", style={'width': '75%'})], width=2),
                        dbc.Col([
                            dcc.Input(id="formula{}".format(str(cur_id)),
                                      type="text",
                                      placeholder='', style={'width': '85%'})], width=3),
                        dbc.Col([
                            dcc.Input(id="charge{}".format(str(cur_id)),
                                      type="number",
                                      placeholder="0",
                                      style={'width': '75%'})], width=2),
                        dbc.Col([
                            dcc.Dropdown(
                                id="compartment{}".format(str(cur_id)),
                                options=[{"label": x,
                                          "value": x, }
                                         for x in list(comp_list)],
                                placeholder='',
                                multi=False
                            )], width=2),
                    ]),
                )
                return contents, left_contents, right_contents, [
                    {"label": i, "value": i, } for i in list(unbalanced_list)]
            else:
                contents = [
                    "Select a reaction above to curate the charge balance"]
                left_contents = [
                    "The left side of the reaction will display here"]
                right_contents = [
                    "The right side of the reaction will display here"]

        @_app.callback(
            Output("node-data", "children"),
            [Input("net", "selectedNodeData")]
        )
        def display_nodedata(datalist):
            contents = "Click on a node to see its details here"
            if datalist is not None:
                if len(datalist) > 0:
                    data = datalist[-1]
                    if "formula" in data:
                        if len(datalist) > 0:
                            data = datalist[-1]
                            contents = []
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
                        contents = []
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
            return contents

        @_app.callback(
            Output("net_curate", "elements"),
            Input("btn_update_cur_net", "n_clicks"),
            State("pathways_dropdown_curate", "value"),
            prevent_initial_call=True,
        )
        def update_net_cur(nclicks, pathway):
            if dash.callback_context.triggered[0]['prop_id'] == 'btn_update_cur_net.n_clicks':
                if pathway != '':
                    model = self._mm
                    nm, network = read_model(self._model, model, "C")
                    pathway_list, rxn_set = get_pathway_list(nm, pathway)
                    nodes, edges, obj_flux = build_network(
                        self, nm, model, rxn_set, network, [])
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
            #         'net_curate.selectedEdgeData' and len(edge) > 0:
            #     return(edge[0]['rid'], [])
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
                    for l in left:
                        eq_dict[Compound(l[1]).in_compartment(l[2])] = \
                            int(l[0]) * - 1
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
                        if type(r) == list:
                            if i.id == r[0]:
                                reaction = i
                        else:
                            if i.id == r:
                                reaction = i
                    contents = []
                    contents.append(html.H5(
                                    "ID: " + i.id))
                    contents.append(html.H5(
                                    "Name: " + i.name))
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
                                    placeholder=str(reaction.equation.direction),
                                    multi=False
                                )
                            ]), ], id='direction'),)
                    contents.append(html.H5("Left Equation:"))
                    dropid = 0
                    for j in reaction.equation.left:
                        contents.append(
                            dbc.Row([
                                dbc.Col([
                                    dcc.Input(id="stoich{}".format(str(dropid)),
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
                                    dcc.Input(id="stoich{}".format(str(dropid)),
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
                                    dcc.Dropdown(id="drop{}".format(str(dropid)),
                                                 options=[{"label": x, "value": x}
                                                          for x in list(c_list)],
                                                 placeholder='',
                                                 multi=False)]),
                                dbc.Col([
                                    dcc.Dropdown(id="compart{}".
                                                 format(str(dropid)),
                                                 options=[{"label": x, "value": x}
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
             Output("fba-data", "children")],
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
                nodes, edges, obj_flux = build_network(self, nm, model,
                                                       rxn_list, network,
                                                       fba_dropdown)
                elements = nodes + edges

                if obj_flux != 0:
                    contents = []
                    contents.append(html.H5("{} flux "
                                    "is {}".format(fba_dropdown, obj_flux)))
                else:
                    contents = []
                    contents.append(html.H5("No flux through the network"))

                return elements, contents

        @_app.callback(Output("download", "data"),
                       [Input("btn_tsv", "n_clicks"),
                        Input("pathways_dropdown", "value"),
                        Input("element_dropdown", "value"),
                        Input("compounds_dropdown", "value"),
                        Input("fba_dropdown", "value"),
                        Input("filter1_dropdown", "value"),
                        Input("filter2_dropdown", "value"), ],
                       prevent_initial_call=True,)
        def func(n_clicks, pathways_dropdown, element_dropdown,
                 compounds_dropdown, fba_dropdown,
                 filter1_dropdown, filter2_dropdown):
            if n_clicks is not None:
                nm, network = read_model(self._model, self._mm,
                                         element_dropdown)
                pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)

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
                    if rxn_list is str:
                        raise PreventUpdate
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
                id = []
                name = []
                equation = []
                for rxn in nm.reactions:
                    if rxn.id in rxn_list:
                        id.append(rxn.id)
                        name.append(rxn.name)
                        equation.append(str(rxn.properties["equation"]))
                df = pd.DataFrame({"id": id, "name": name,
                                   "equation": equation})
                return(dcc.send_data_frame(df.to_csv, "exported_rxns.csv"))

        @_app.callback(Output("net", "generateImage"),
                       [Input('btn-get-png', 'n_clicks')],
                       prevent_initial_call=True,)
        def get_image(n_clicks):
            if n_clicks is not None:
                ftype = "png"
                action = "download"
                return {'type': ftype, 'action': action}


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


def get_compounds_list(nm):
    compounds_list = []
    for i in nm.compounds:
        compounds_list.append(i.id)
    return compounds_list


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
                                if j.name in left and k.name in right \
                                    and (str(rxn._properties['equation'].
                                             __dict__['_direction']) ==
                                         "Direction.Forward" or
                                         str(rxn._properties['equation'].
                                             __dict__['_direction']) ==
                                         "Direction.Both"):
                                    adj[k.name] = rxn.id
                                elif j.name in right and k.name in left and \
                                    (str(rxn._properties['equation'].
                                         __dict__['_direction']) ==
                                     "Direction.Reverse" or
                                     str(rxn._properties['equation'].
                                         __dict__['_direction']) ==
                                     "Direction.Both"):
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
                                if j.name in left and k.name in right and \
                                    (str(rxn._properties['equation'].
                                     __dict__['_direction']) ==
                                     "Direction.Forward" or
                                     str(rxn._properties['equation'].
                                         __dict__['_direction']) ==
                                        "Direction.Both"):
                                    adj[k.name] = rxn.id
                                elif j.name in right and k.name in left and \
                                    (str(rxn._properties['equation'].
                                     __dict__['_direction']) ==
                                     "Direction.Reverse" or
                                     str(rxn._properties['equation'].
                                         __dict__['_direction']) ==
                                        "Direction.Both"):
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


# Useful function to build the subset network. Returns nodes and edges
# from the network associated with the rxn_set of interest
def build_network(self, nm, mm, rxn_set, network, fba_dropdown):
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
    nodes = []
    edges = []

    obj = 0
    if not isinstance(fba_dropdown, list):
        objective = str(fba_dropdown)

        # Set up the solver
        solver = self._get_solver()
        # Define the flux balance
        problem = fluxanalysis.FluxBalanceProblem(mm, solver)
        problem.maximize(fba_dropdown)
        obj = problem.get_flux(fba_dropdown)
        flux_carrying = defaultdict(lambda: 0)
        for i in mm.reactions:
            flux_carrying[i] = problem.get_flux(i)
    else:
        flux_carrying = {}
        for i in mm.reactions:
            flux_carrying[i] = 0
    comp = []
    for i in nm.compartments:
        comp.append(i.id)
    color = ["#6FB2D2", "#85C88A", "#EBD671", "#EEEEEE"]
    col_dict = {}
    for col in range(0, len(comp)):
        col_dict[comp[col]] = color[col]
    for rxn in network[0]:
        if rxn.id in rxn_set:
            if rxn.name:
                rname = rxn.name
            else:
                rname - rxn.id
            rxn_num = 0
            visited = set()
            for cpd in network[0][rxn][0]:
                if cpd[0] not in visited:
                    visited.add(cpd[0])
                    visited.add(cpd[1])
                    nodes.append({'data': {
                                  'id': str(cpd[0]),
                                  'label': name[str(cpd[0])
                                                [0:(str(cpd[0]).
                                                 find('['))]],
                                  'formula': formula[str(cpd[0])
                                                     [0:(str(cpd[0]).
                                                      find('['))]],
                                  'compartment': cpd[0].compartment,
                                  'col': col_dict[cpd[0].compartment],
                                  'charge': charge[str(cpd[0])
                                                   [0:(str(cpd[0]).
                                                    find('['))]],
                                  'type': 'compound'
                                  }})
                    nodes.append({'data': {
                                  'id': str(cpd[1]),
                                  'label': name[str(cpd[1])
                                                [0:(str(cpd[1]).
                                                 find('['))]],
                                  'formula': formula[str(cpd[1])
                                                     [0:(str(cpd[1]).
                                                      find('['))]],
                                  'compartment': cpd[0].compartment,
                                  'col': col_dict[cpd[0].compartment],
                                  'charge': charge[str(cpd[0])
                                                   [0:(str(cpd[0]).
                                                    find('['))]],
                                  'type': 'compound'
                                  }})
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
                        nodes.append({'data': {
                                      'id': "{}_node".format(rxn.id),
                                      'label': rxn.id,
                                      'type': 'reaction',
                                      'equation': str(rxn.properties
                                                      ["equation"]),
                                      'pathways': path,
                                      'orig_id': rxn.id,
                                      'name': rxn.properties["name"],
                                      'genes': rxn.properties["genes"],
                                      'pathways': path_prop,
                                      'EC': rxn.properties["EC"],
                                      'flux_node': flux_carrying[rxn.id]
                                      },
                                      'style': {
                                      'shape': 'rectangle'
                        }})
                        if str(rxn._properties['equation'].__dict__
                               ['_direction']) == "Direction.Forward":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_",
                                               str(cpd[1]), "for"]),
                                'orig_id': rxn.id,
                                'source': "{}_node".format(rxn.id),
                                'target': str(cpd[1]),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_",
                                               str(cpd[0]), "for"]),
                                'orig_id': rxn.id,
                                'source': str(cpd[0]),
                                'target': "{}_node".format(rxn.id),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                        elif str(rxn._properties['equation'].__dict__
                                 ['_direction']) == "Direction.Reverse":
                            flux = 0 - float(flux_carrying[rxn.id])
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_",
                                              str(cpd[0]), "rev"]),
                                'orig_id': rxn.id,
                                'source': "{}_node".format(rxn.id),
                                'target': str(cpd[0]),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_",
                                               str(cpd[1]), "rev"]),
                                'orig_id': rxn.id,
                                'source': str(cpd[1]),
                                'target': "{}_node".format(rxn.id),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                        elif str(rxn._properties['equation'].
                                 __dict__['_direction']) \
                                == "Direction.Both":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_", str(cpd[1]),
                                              "for"]),
                                'orig_id': rxn.id,
                                'source': "{}_node".format(rxn.id),
                                'target': str(cpd[1]),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                            flux = 0 - float(flux_carrying[rxn.id])
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_", str(cpd[0]),
                                               "for"]),
                                'orig_id': rxn.id,
                                'source': "{}_node".format(rxn.id),
                                'target': str(cpd[0]),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                            flux = flux_carrying[rxn.id]
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_", str(cpd[0]),
                                               "rev"]),
                                'orig_id': rxn.id,
                                'source': str(cpd[0]),
                                'target': "{}_node".format(rxn.id),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.properties
                                                ["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
                            flux = 0 - float(flux_carrying[rxn.id])
                            edges.append({'data': {
                                'id': "".join([rxn.id, "_", str(cpd[1]),
                                              "rev"]),
                                'orig_id': rxn.id,
                                'source': str(cpd[1]),
                                'target': "{}_node".format(rxn.id),
                                'rid': rxn.id,
                                'label': "".join([rxn.name]),
                                'equation': str(rxn.
                                                properties["equation"]),
                                'pathways': path,
                                'flux': flux
                            }})
                            rxn_num += 1
    return nodes, edges, obj


class VisualizationCommand(MetabolicMixin,
                           Command, FilePrefixAppendAction):
    """Generate graphical representations of the model"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method', type=str, default='fpp',
            help='Compound pair prediction method',
            choices=['fpp', 'no-fpp'])
        parser.add_argument(
            '--exclude', metavar='reaction', type=convert_to_unicode,
            default=[], action=FilePrefixAppendAction,
            help='Reaction(s) to exclude from metabolite pair prediction')
        parser.add_argument(
            '--element', type=str, default='C',
            help='Element transfer to show on the graph (C, N, O, S).')
        parser.add_argument(
            '--cpd-detail', type=str, default=None, action='append',
            nargs='+', help='List of properties to include on '
                            'compound nodes.')
        parser.add_argument(
            '--rxn-detail', type=str, default=None, action='append',
            nargs='+', help='List of properties to include on reaction '
                            'nodes.')
        parser.add_argument(
            '--subset', type=argparse.FileType('rU'), default=None,
            help='File containing a list of reactions IDs or compound IDs '
                 '(with compartment). This file determines which reactions '
                 'will be visualized')
        parser.add_argument(
            '--color', type=argparse.FileType('rU'), default=None, nargs='+',
            help='File containing node color mappings')
        parser.add_argument(
            '--image', metavar='Image Format', type=str, choices=FORMATS,
            required='--image-size' in sys.argv,
            help='Image output file format '
                 '(e.g. pdf, png, eps)')
        parser.add_argument(
            '--hide-edges', type=argparse.FileType('rU'), default=None,
            help='File containing compound edges to exclude '
                 'from visualization')
        parser.add_argument(
            '--combine', metavar='Combine Level', type=int, choices=range(3),
            default=0, help='Combined reaction nodes in three '
                            'different levels.')
        parser.add_argument(
            '--compartment', action='store_true',
            help='Include compartments in final visualization')
        parser.add_argument(
            '--output', type=str,
            help='Full path for visualization outputs, including prefix of '
                 'output files. '
                 'e.g. "--output <output-directory-path>/output-prefix"')
        parser.add_argument(
            '--image-size', metavar=('Width', 'Height'),
            default=('None', 'None'), nargs=2, type=float,
            help='Set the width and height of the graph image. '
                 '(width height)(inches)')
        parser.add_argument(
            '--array', type=int, default=None,
            help='The number of columns to use in the final '
                 'network image')
        group = parser.add_mutually_exclusive_group()
        group.add_argument(
            '--fba', type=argparse.FileType('rU'), default=None,
            help='File containing fba reaction flux')
        group.add_argument(
            '--fva', type=argparse.FileType('rU'), default=None,
            help='File containing fva reaction flux')

        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""
        for reaction in self._model.reactions:
            if reaction.equation is None:
                logger.warning(
                    'Reaction {} was excluded from visualization due to '
                    'missing reaction equation'.format(reaction.id))

        compound_formula = graph.get_compound_dict(self._model)

        vis_rxns = rxnset_for_vis(
            self._mm, self._args.subset, self._args.exclude)

        if self._args.array is not None and int(self._args.array) <= 0:
            logger.error(
                "'--array' should be followed by a positive integer, number "
                "'{}' is invalid. Visualization has stopped, please fix the "
                "number first".format(self._args.array))
            quit()

        if self._args.element.lower() == 'all':
            self._args.element = None
        else:
            if self._args.element not in _AtomType._ELEMENTS:
                logger.error(
                    "Given element '{}' doesn't represent any chemical element"
                    ", visualization has terminated. Please check your "
                    "--element parameter".format(self._args.element))
                quit()

        self.analysis = None
        reaction_dict = {}
        if self._args.fba is not None:
            self.analysis = 'fba'
            for row in csv.reader(self._args.fba, delimiter=str('\t')):
                row[0] = convert_to_unicode(row[0])
                if row[0] in vis_rxns:
                    try:
                        if abs(float(row[1])) <= EPSILON:
                            row[1] = 0
                    except ValueError:
                        logger.warning(
                            'Reaction {} has an invalid flux value of {} '
                            'and will be considered as a flux of 0'.format(
                                row[0], row[1]))
                        row[1] = 0
                    reaction_dict[row[0]] = (float(row[1]), 1)
                else:
                    logger.warning(
                        'Reaction {} in input fba file was excluded from '
                        'visualization'.format(row[0]))

        if self._args.fva is not None:
            self.analysis = 'fva'
            for row in csv.reader(self._args.fva, delimiter=str('\t')):
                row[0] = convert_to_unicode(row[0])
                if row[0] in vis_rxns:
                    try:
                        if abs(float(row[1])) <= EPSILON:
                            row[1] = 0
                        if abs(float(row[2])) <= EPSILON:
                            row[2] = 0
                    except ValueError:
                        logger.warning(
                            'Reaction {} has an invalid flux value of {},{} '
                            'and will be considered as a flux of 0,0'.format(
                                row[0], row[1], row[2]))
                        row[1], row[2] = 0, 0
                    reaction_dict[row[0]] = (float(row[1]), float(row[2]))
                else:
                    logger.warning(
                        'Reaction {} in input fva file was excluded from '
                        'visualization due to not being defined in the '
                        'model'.format(row[0]))

        exclude_for_fpp = [self._model.biomass_reaction] + self._args.exclude
        filter_dict, style_flux_dict = graph.make_network_dict(
            self._model, self._mm, vis_rxns, self._args.method,
            self._args.element, exclude_for_fpp, reaction_dict, self.analysis)

        model_compound_entries, model_reaction_entries = {}, {}
        for c in self._model.compounds:
            model_compound_entries[c.id] = c
        for r in self._model.reactions:
            model_reaction_entries[r.id] = r

        hide_edges = []
        if self._args.hide_edges is not None:
            for row in csv.reader(self._args.hide_edges, delimiter=str('\t')):
                hide_edges.append((row[0], row[1]))
                hide_edges.append((row[1], row[0]))

        cpair_dict, new_id_mapping, new_style_flux_dict = \
            graph.make_cpair_dict(
                filter_dict, self._args.method, self._args.combine,
                style_flux_dict, hide_edges)

        g = graph.make_bipartite_graph_object(
            cpair_dict, new_id_mapping, self._args.method, self._args.combine,
            model_compound_entries, new_style_flux_dict, self.analysis)

        exchange_cpds = set()
        for rid in vis_rxns:
            if self._mm.is_exchange(rid) and rid != \
                    self._model.biomass_reaction:
                exchange_rxn = self._mm.get_reaction(rid)
                for c, _ in exchange_rxn.compounds:
                    if c not in g.nodes_id_dict:
                        g = add_ex_cpd(g, c, model_compound_entries[c.name],
                                       compound_formula, self._args.element)
                    exchange_cpds.add(c)
                g = add_exchange_rxns(
                    g, rid, exchange_rxn, style_flux_dict)

        recolor_dict = {}
        if self._args.color is not None:
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str('\t')):
                    recolor_dict[row[0]] = row[1]
        g = add_node_props(g, recolor_dict)
        g = add_node_label(g, self._args.cpd_detail, self._args.rxn_detail)
        bio_cpds_sub = set()
        bio_cpds_pro = set()

        if self._model.biomass_reaction in vis_rxns:
            if self._model.biomass_reaction in model_reaction_entries:
                nm_biomass_rxn = \
                    model_reaction_entries[self._model.biomass_reaction]
                g = add_biomass_rxns(g, nm_biomass_rxn)
                for cpd, _ in nm_biomass_rxn.equation.left:
                    bio_cpds_sub.add(text_type(cpd))
                for cpd, _ in nm_biomass_rxn.equation.right:
                    bio_cpds_pro.add(text_type(cpd))
            else:
                logger.warning(
                    'Biomass reaction {} was excluded from visualization '
                    'due to missing reaction entry'.format(
                        self._model.biomass_reaction))

        for node in g.nodes:
            if node.props['id'] in bio_cpds_sub:
                node.props['type'] = 'cpd_biomass_substrate'
            elif node.props['id'] in bio_cpds_pro:
                node.props['type'] = 'cpd_biomass_product'
            elif node.props['id'] in exchange_cpds:
                node.props['type'] = 'cpd_exchange'
            else:
                continue

        if self._args.method == 'no-fpp' and self._args.combine != 0:
            logger.warning(
                '--combine option is not compatible with no-fpp method')

        hw = self._args.image_size
        width = hw[0]
        height = hw[1]

        if self._args.output is not None:
            output = self._args.output
        else:
            output = 'reactions'

        if self._args.compartment:
            boundaries, extracellular = get_cpt_boundaries(self._model)
            boundary_tree, extracellular = make_cpt_tree(boundaries,
                                                         extracellular)
            with io.open('{}.dot'.format(output), 'w', encoding='utf-8',
                         errors='backslashreplace') as f:
                g.write_graphviz_compartmentalized(
                    f, boundary_tree, extracellular, width, height)
        else:
            with io.open('{}.dot'.format(output), 'w', encoding='utf-8',
                         errors='backslashreplace') as f:
                g.write_graphviz(f, width, height)
        with io.open('{}.nodes.tsv'.format(output), 'w', encoding='utf-8',
                     errors='backslashreplace') as f:
            g.write_nodes_tables(f)
        with io.open('{}.edges.tsv'.format(output), 'w', encoding='utf-8',
                     errors='backslashreplace') as f:
            g.write_edges_tables(f)

        if self._args.image is not None:
            if render is None:
                raise ImportError(
                    'Making an image directly requires the '
                    'graphviz python bindings and the graphviz program '
                    'to be installed ("pip install graphviz")')
            else:
                if len(g.nodes_id_dict) > 700:
                    logger.info(
                        'This graph contains {} reactions, '
                        'graphs of this size may take a long time to '
                        'create'.format(len(filter_dict.keys())))
                if self._args.array is None:
                    render('dot', self._args.image, '{}.dot'.format(output))
                else:
                    out = '{}.dot'.format(output)
                    format = self._args.image
                    image = '{}.dot.{}'.format(output, format)
                    os.system("ccomps -x {} |dot |gvpack -array{} |neato "
                              "-T{} -n2 -o {}".format(out, self._args.array,
                                                      format, image))
        else:
            if self._args.array is not None:
                out = '{}.dot'.format(output)
                os.system("ccomps -x {} |dot |gvpack -array{} "
                          ">{}_new.dot".format(out, self._args.array, output))
                os.system("mv {}_new.dot {}.dot".format(output, output))


def rxnset_for_vis(mm, subset_file, exclude):
    """ Create a collection of reaction IDs that need to be visualized.

    Generate a set that contains IDs of all reactions that need to be
    visualized. This set can be generated from a file containing
    reaction or compound IDs.

    Args:
        mm: Metabolic model, class 'psamm.metabolicmodel.MetabolicModel'.
        subset_file: None or an open file containing a list of reactions
                     or compound ids.
        exclude: rxns to exclude
    """
    all_cpds = set()
    for cpd in mm.compounds:
        all_cpds.add(text_type(cpd))
    if subset_file is None:
        if len(exclude) == 0:
            final_rxn_set = set(mm.reactions)
        else:
            final_rxn_set = set(
                [rxn for rxn in mm.reactions if rxn not in exclude])
    else:
        final_rxn_set = set()
        cpd_set = set()
        rxn_set = set()
        for line in subset_file.readlines():
            if not convert_to_unicode(line).startswith('#'):
                value = convert_to_unicode(line).strip()
                if value in all_cpds:
                    cpd_set.add(value)
                elif mm.has_reaction(value):
                    rxn_set.add(value)
                else:
                    raise ValueError(u'{} is in subset file but is not a '
                                     'compound or reaction ID'.format(value))

        if all(i > 0 for i in [len(cpd_set), len(rxn_set)]):
            raise ValueError('Subset file is a mixture of reactions and '
                             'compounds.')
        else:
            if len(cpd_set) > 0:
                for rx in mm.reactions:
                    rxn = mm.get_reaction(rx)
                    if any(text_type(c) in cpd_set
                           for (c, _) in rxn.compounds):
                        final_rxn_set.add(rx)
            elif len(rxn_set) > 0:
                final_rxn_set = rxn_set

    return final_rxn_set


def add_biomass_rxns(g, nm_bio_reaction):
    """ Adds biomass reaction nodes and edges to a graph object.

    This function is used to add nodes and edges of biomass reaction to
    the graph object. All the following properties are defined
    for added nodes: shape, style(filled), fill color, label.

    Args:
        g: Graph object.
        nm_bio_reaction: Biomass reaction DictReactionEntry.
    """
    bio_reaction = nm_bio_reaction.equation
    biomass_rxn_id = nm_bio_reaction.id
    direction = graph.dir_value(bio_reaction.direction)
    reactant_list, product_list = [], []
    for c, _ in bio_reaction.left:
        reactant_list.append(c)
    for c, _ in bio_reaction.right:
        product_list.append(c)
    bio_pair = Counter()
    for c, _ in bio_reaction.compounds:
        if text_type(c) in g.nodes_id_dict:
            bio_pair[biomass_rxn_id] += 1
            node_bio = graph.Node({
                'id': u'{}_{}'.format(biomass_rxn_id,
                                      bio_pair[biomass_rxn_id]),
                'entry': [nm_bio_reaction],
                'shape': 'box',
                'style': 'filled',
                'label': biomass_rxn_id,
                'type': 'bio_rxn',
                'fillcolor': ALT_COLOR,
                'compartment': c.compartment})
            g.add_node(node_bio)

            if c in reactant_list:
                g.add_edge(graph.Edge(g.get_node(text_type(c)), node_bio,
                                      {'dir': direction}))
            if c in product_list:
                g.add_edge(graph.Edge(node_bio, g.get_node(text_type(c)),
                                      {'dir': direction}))
    return g


def add_ex_cpd(g, mm_cpd, nm_cpd, compound_formula, element):
    node_dict = {'id': text_type(mm_cpd),
                 'entry': [nm_cpd],
                 'compartment': mm_cpd.compartment,
                 'type': 'cpd'}
    if element is not None:
        if mm_cpd.name not in compound_formula:
            logger.error(u'Compound formulas are required for fpp or specific '
                         'element visualizations, but compound {} does not '
                         'have valid formula, add its formula, or try '
                         '--element all to visualize all pathways without '
                         'compound formula input.'.format(mm_cpd.name))
        else:
            if Atom(element) in compound_formula[mm_cpd.name]:
                node = graph.Node(node_dict)
                g.add_node(node)
    else:
        node = graph.Node(node_dict)
        g.add_node(node)
    return g


def add_node_props(g, recolor_dict):
    """ Update node color in Graph object based on a mapping dictionary

    This function adds shape and style(filled) properties to nodes in a graph.

    Args:
        g: A Graph object that contains nodes and edges.
        recolor_dict: dict of rxn_id/cpd_id[compartment] : hex color code.
    return: a graph object that contains a set of node with defined color.
    """
    cpd_types = ['cpd', 'cpd_biomass_substrate',
                 'cpd_biomass_product', 'cpd_exchange']
    for node in g.nodes:
        node.props['style'] = 'filled'
        if node.props['type'] in cpd_types:
            node.props['fillcolor'] = COMPOUND_COLOR
            node.props['shape'] = 'ellipse'
        elif node.props['type'] == 'rxn':
            node.props['fillcolor'] = REACTION_COLOR
            node.props['shape'] = 'box'
        else:
            if node.props['type'] not in ['bio_rxn', 'Ex_rxn']:
                logger.error('invalid nodes:', type(
                    node.props['entry']), node.props['entry'])
    for r_id in recolor_dict:
        if r_id in g.nodes_original_id_dict:
            for node in g.nodes_original_id_dict[r_id]:
                node.props['fillcolor'] = recolor_dict[r_id]
    return g


def add_node_label(g, cpd_detail, rxn_detail):
    """ Set label of nodes in graph object,

    Set the label of nodes in a graph object based on compound/reaction
    properties provided through the cpd_detail or rxn_detail arguments.

    Args:
        g: A graph object, contain a set of nodes and a dictionary of edges.
        cpd_detail: A list that contains only one
            element, this element is a compound properties name list,
            e.g. detail = [['id', 'name', 'formula']].
        rxn_detail: A list that contains only one
            element, this element is a reaction properties name list,
            e.g. detail = [['id', genes', 'equation']].
    """
    for node in g.nodes:
        if node.props['type'] == 'cpd':
            node.props['label'] = node.props['id']
        elif node.props['type'] == 'rxn':
            node.props['label'] = '\n'.join(
                obj.id for obj in node.props['entry'])

        # update node label based on what users provide in command line
        if cpd_detail is not None:
            if node.props['type'] == 'cpd':
                pre_label = '\n'.join(
                    (node.props['entry'][0].properties.get(value)) for
                    value in cpd_detail[0] if value != 'id' and
                    value in node.props['entry'][0].properties)
                if 'id' in cpd_detail[0]:
                    label = u'{}\n{}'.format(node.props['id'], pre_label)
                else:
                    label = pre_label

                # if all required properties are not in the compound entry,
                # then print compound id
                if label == '':
                    label = node.props['id']

                node.props['label'] = label

        if rxn_detail is not None:
            if node.props['type'] == 'rxn':
                if len(node.props['entry']) == 1:
                    pre_label = u'\n'.join(
                        node.props['entry'][0].properties.get(value)
                        for value in rxn_detail[0] if
                        value != 'equation' and
                        value != 'id' and value in
                        node.props['entry'][0].properties)
                    if 'id' in rxn_detail[0]:
                        label = u'{}\n{}'.format(
                            node.props['entry'][0].properties['id'], pre_label)
                    else:
                        label = pre_label

                    if 'equation' in rxn_detail[0]:
                        label += u'\n{}'.format(
                            node.props['entry'][0].properties.get('equation'))

                    # if all required properties are not in the reaction entry,
                    # then print reaction id
                    if label == '':
                        label = node.props['id']

                    node.props['label'] = label
    return g


def make_cpt_tree(boundaries, extracellular):
    """ This function will create a tree-like dictionary that can be used to
    determine compartment location.

        This function will take a list of compartment boundary tuples
        (eg: [(c, e), (c, p)]) and use this data to construct a tree like
        dictionary of parent compartments to lists of compartments they are
        adjacent to. An extracellular compartment will also be set to use as
        a starting point for the 'outermost' compartment in the model. If
        none is explicitly defined then 'e' will be used by default.
        If an 'e' compartment is not in the model either then a random
        compartment will be used as a starting point.

        args:
        boundaries: a list of tuples of adjacent compartment ids.
        extracellular: the extracellular compartment in the model.
    """
    children = defaultdict(set)
    compartments = set()
    for (j, k) in boundaries:
        compartments.add(j)
        compartments.add(k)
    if extracellular not in compartments:
        etmp = sorted(list(compartments), reverse=True)
        extracellular = etmp[0]
        logger.warning('No extracellular compartment was defined in the '
                       'model.yaml file and no "e" compartment in the model. '
                       'Trying to use {} as the extracellular compartment.'
                       .format(etmp[0]))
    for cpt in compartments:
        for (j, k) in boundaries:
            j = text_type(j)
            k = text_type(k)
            if j == cpt:
                children[cpt].add(k)
            elif k == cpt:
                children[cpt].add(j)
    return children, extracellular


def get_cpt_boundaries(model):
    """This function will determine the compartment boundaries in a model

    This function will take a native model object and determine the
    compartment boundaries either based on the predefined compartments in
    the model.yaml or based on the reaction equations in the model.

    args:
    model: Native model, class 'psamm.datasource.native.NativeModel'.
    """
    if model.extracellular_compartment is not None:
        extracellular = model.extracellular_compartment
    else:
        extracellular = 'e'

    if len(model.compartment_boundaries) != 0:
        boundaries = model.compartment_boundaries
    else:
        boundaries = set()
        for rxn in model.reactions:
            cpd_cpt = set()
            for cpd in rxn.equation.compounds:
                cpd_cpt.add(cpd[0].compartment)
            if len(cpd_cpt) > 1:
                cpd_cpt = list(cpd_cpt)
                boundaries.add(tuple(sorted((cpd_cpt[0], cpd_cpt[1]))))
    return boundaries, extracellular


def add_exchange_rxns(g, rxn_id, reaction, style_flux_dict):
    """ Add exchange reaction nodes and edges to graph object.

    This function is used to add nodes and edges of exchange reactions to
    the graph object. It will return an updated graph object that contains
    nodes representing exchange reactions.

    Args:
        g: A graph object that contains a set of nodes and some edges.
        rxn_id: Exchange reaction id,
        reaction: Exchange reaction object(metabolic model reaction),
            class 'psamm.reaction.Reaction'.
        style_flux_dict: dictionary of reaction ID maps to edge style and
            edge width.
        analysis: "None" type or a string indicates if FBA or FVA file is
            given in command line.
    """
    direction = graph.dir_value(reaction.direction)
    for c, _ in reaction.compounds:
        if text_type(c) in g.nodes_id_dict:
            node_ex = graph.Node({
                'id': text_type(rxn_id),
                'entry': [reaction],
                'shape': 'box',
                'style': 'filled',
                'label': rxn_id,
                'type': 'Ex_rxn',
                'fillcolor': ACTIVE_COLOR,
                'compartment': c.compartment})
            g.add_node(node_ex)

            for c1, _ in reaction.left:
                g.add_edge(graph.Edge(
                    g.get_node(text_type(c1)), node_ex, {
                        'dir': direction,
                        'style': style_flux_dict[rxn_id][0],
                        'penwidth': style_flux_dict[rxn_id][1]
                    }))
            for c2, _ in reaction.right:
                g.add_edge(graph.Edge(
                    node_ex, g.get_node(text_type(c2)), {
                        'dir': direction,
                        'style': style_flux_dict[rxn_id][0],
                        'penwidth': style_flux_dict[rxn_id][1]
                    }))
    return g
