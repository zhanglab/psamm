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

from psamm.expression import boolean
from psamm.lpsolver import generic
from psamm import fluxanalysis
from psamm.lpsolver import glpk, lp
from psamm.graph import make_network_dict, write_network_dict
from psamm.datasource.reaction import Reaction, Compound
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

logger = logging.getLogger(__name__)

REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'

EPSILON = 1e-5

class InteractiveCommand(MetabolicMixin,SolverCommandMixin,
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

        el='C'
        pathway='All'

        pathway_list, rxn_set = get_pathway_list(self._model, "All")
        compounds_list = get_compounds_list(self._model)
        rxns=list(rxn_set)
        rxns_full = []
        for i in self._model.reactions:
            rxns_full.append(i.id)

        # nodes, edges = build_network(nm, rxn_set, network)
        # initialize an empty list. the full is messy
        nodes, edges = [], []

        # Initialize the app
        self._app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        server = self._app.server
        if self._app is not None and hasattr(self, "callbacks"):
            self.callbacks(self._app)

        # Set styles for the app
        styles = {
                'output': {
                        'border': 'thin lightgrey solid',
                        'overflowX': 'scroll'
                },
                'tab': {'height': 'calc(98vh - 115px)'}
        }
        col_swatch = px.colors.qualitative.Dark24
        default_stylesheet = [
                {
                        'selector': 'node',
                        'style': {
                                'background-color': '#BFD7B5',
                                'label': 'data(label)'}},
                {
                        "selector": "edge",
                        "style": {
                                "width": 1,
                                "curve-style": "bezier"}}
        ]

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
                ),],

                brand="Psamm web-based visualization of metabolic models",
                brand_href="#",
                color="dark",
                dark=True,
        )

        # Define the body of the app
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
                                dbc.Col([
                                dbc.Row([

                                        dbc.Col([
                                                        cyto.Cytoscape(id = 'net',
                                                        layout={'name':'cose',
                                                                'animate':True,
                                                                'animationDuration':1000},
                                                        style={'width': '600px', 'height': '600px'},
                                                        elements=nodes+edges,
                                                        stylesheet=[{
                                                                    'selector': 'node',
                                                                    'style': {
                                                                            'background-color': '#BFD7B5',
                                                                            'label': 'data(label)'}},
                                                                    {
                                                                    "selector": "edge",
                                                                    "style": {
                                                                            "width": 1,
                                                                            "curve-style": "bezier",
                                                                            "target-arrow-shape":"triangle"}},
                                                                    {
                                                                    "selector": "[flux > 0]",
                                                                    "style": {
                                                                            "line-color": "blue",
                                                                            "target-arrow-color":"blue",
                                                                            "width": 1,
                                                                            "curve-style": "bezier",
                                                                            "target-arrow-shape":"triangle",
                                                                            "line-dash-pattern": [6,3],
                                                                            "line-dash-offset": 24
                                                                            }
                                                                    },

                                                                    ],
                                                                    minZoom=0.06
                                                                    )],sm=12,md=8,),
                                            dbc.Row([
                                                    dbc.Alert(
                                                            id="node-data",
                                                            children="Click on a node to see its details here",
                                                            color="secondary",
                                                    ),
                                                    dbc.Alert(
                                                            id="edge-data",
                                                            children="Click on an edge to see its details here",
                                                            color="secondary",)
                                            ]),

                                ]),
                                ]),
                                dbc.Col([
                                        html.H5("Pathways:"),
                                        dbc.Row(
                                        [dcc.Dropdown(
                                            id="pathways_dropdown",
                                            options=[
                                            {"label": i,
                                            "value": i,}
                                            for i in list(pathway_list)],
                                        value=pathway_list[0],
                                        multi=True,
                                        placeholder="",
                                        style={"width": "100%"},),]),
                                        html.H5("Element Transfer Networks:"),
                                        dbc.Row([dcc.Dropdown(
                                            id="element_dropdown",
                                            options=[
                                                {"label": i,
                                                "value": i,}
                                                for i in ["C","N","S","P"]],
                                            value=["C","N","S","P"],
                                            multi=False,
                                            style={"width": "100%"},),]),
                                        html.H5("Compound Search:"),
                                        dbc.Row([dcc.Dropdown(
                                            id="compounds_dropdown",
                                            options=[
                                            {"label": i,
                                            "value": i,}
                                        for i in list(compounds_list)],
                                            value=compounds_list,
                                            multi=False,
                                            style={"width": "100%"},),]),
                                        html.H5("Path search"),
                                        html.P("Select two compounds below to see the three shortest paths"),
                                        html.P("Compound 1:"),
                                        dbc.Row([dcc.Dropdown(
                                            id="filter1_dropdown",
                                            options=[{"label": i,
                                                "value": i,}
                                            for i in list(compounds_list)],
                                                value=compounds_list,
                                                multi=False,
                                                style={"width": "100%"},),]),
                                                html.P("Compound 2:"),
                                            dbc.Row([dcc.Dropdown(
                                                id="filter2_dropdown",
                                                options=[{"label": i,
                                                    "value": i,}
                                            for i in list(compounds_list)],
                                                value=compounds_list,
                                                multi=False,
                                                style={"width": "100%"},),]),
                                        html.H5("Flux Analysis:"),
                                            dbc.Row([dcc.Dropdown(
                                                id="fba_dropdown",
                                                options=[{"label": i,
                                                    "value": i,}
                                            for i in list(rxns_full)],
                                                value=rxns,
                                                multi=False,
                                                style={"width": "100%"},),]),
                                        dbc.Alert(id="fba-data",
                                            children="Select a reaction to see the flux here",
                                            color="secondary",),
                                        dbc.Col([
                                            html.Div([
                                                html.Button("Submit", id="btn_sub")
                                                    ]),
                                            html.Div([
                                            html.Button("Download CSV", id="btn_tsv"),
                                            dcc.Download(id="download"),
                                            ]),
                                            html.Div([
                                            html.Button("Download png", id="btn-get-png"),
                                            dcc.Download(id="downloadpng"),
                                            ]),
                                        ]),
                                    ],
                                sm=12,
                                md=4,
                                ),

                                ]),

                        dbc.Row([dcc.Markdown(
                            """
                            \* Data analysis carried out for demonstration of data visualisation purposes only.
                            """
                            )],
                            style={"fontSize": 11, "color": "gray"},),
                ],
                style={"marginTop": 20},
        )

        self._app.layout = html.Div([navbar, body_layout])
        webbrowser.open_new("http://localhost:{}".format("8050"))
        self._app.run_server(debug=True)


    def callbacks(self, _app):

        @_app.callback(
                Output("node-data", "children"),
                [Input("net", "selectedNodeData")]
        )
        def display_nodedata(datalist):
                contents = "Click on a node to see its details here"
                if datalist is not None:
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

                return contents

        @_app.callback(
                Output("edge-data", "children"),
                [Input("net", "selectedEdgeData")]
        )
        def display_nodedata(datalist):
                contents = "Click on an edge to see its details here"
                if datalist is not None:
                        if len(datalist) > 0:
                                data = datalist[-1]
                                contents = []
                                contents.append(html.H5(
                                                "ID: " + data["id"].title()))
                                contents.append(
                                        html.P(
                                                "Name: "
                                                + data["label"]
                                        )
                                )
                                contents.append(
                                        html.P(
                                                "Equation: "
                                                + str(data["equation"])
                                        )
                                )
                                contents.append(
                                        html.P(
                                                "Pathways: "
                                                + data["pathways"]
                                        )
                                )
                                contents.append(
                                        html.P(
                                                "Flux: "
                                                + str(data["flux"])
                                        )
                                )


                return contents

        @_app.callback(
                [Output("net", "elements"),
                Output("fba-data", "children")],
                [
                        Input("btn_sub", "n_clicks"),
                ],
                [
                        State("pathways_dropdown", "value"),
                        State("element_dropdown", "value"),
                        State("compounds_dropdown", "value"),
                        State("fba_dropdown", "value"),
                        State("filter1_dropdown", "value"),
                        State("filter2_dropdown", "value"),
                ],
                prevent_initial_call=True,
        )
        def filter_nodes(n_clicks, pathways_dropdown, element_dropdown, \
                         compounds_dropdown, fba_dropdown, \
                         filter1_dropdown, filter2_dropdown):
            if isinstance(pathways_dropdown, list) \
            or isinstance(element_dropdown, str) \
            or isinstance(compounds_dropdown, str) or isinstance(fba_dropdown, str) \
            or (isinstance(filter1_dropdown, str) and isinstance(filter1_dropdown, str)):
                nm, network = read_model(self._model, self._mm, element_dropdown)
                pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)

                if isinstance(filter1_dropdown, str) and isinstance(filter2_dropdown, str):
                    rxn_list = set()
                    cpd_list = [filter1_dropdown, filter2_dropdown]
                    middle2 = []
                    middle3 = []
                    middle2, rxn_list1 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    middle3, rxn_list2 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    middle4, rxn_list3 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    rxn_list = rxn_list1 + rxn_list2 + rxn_list3
                    if len(rxn_list)==0:
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
                nodes, edges, obj_flux = build_network(self, nm, rxn_list, network, fba_dropdown)
                elements=nodes+edges

                if obj_flux != 0:
                    contents = []
                    contents.append(html.H5("{} flux is {}".format(fba_dropdown, obj_flux)))
                else:
                    contents = []
                    contents.append(hmtl.H5("No flux through the network"))

                return elements, contents


        # @_app.callback(
        #         Output("fba-data", "children"), [Input("fba_dropdown", "value")],
        #         prevent_initial_call=True,
        # )
        # def Run_FBA(fba_dropdown):
        #   if not isinstance(fba_dropdown, list):
        #       objective=str(fba_dropdown)
        #       # First read in the base model
        #       nm, network = read_model(self._model, self._mm, element_dropdown)
        #       # Set up the solver
        #       solver = glpk.Solver()
        #       # Define the flux balance
        #       problem = fluxanalysis.FluxBalanceProblem(nm, solver)
        #       problem.maximize(fba_dropdown)
        #
        #
        #
        #       return(str(problem.get_flux(fba_dropdown)))
        #   else:
        #       return("More than one reaction selected")

        @_app.callback(
            Output("download", "data"),
            [Input("btn_tsv", "n_clicks"),
            Input("pathways_dropdown", "value"),
            Input("element_dropdown", "value"),
            Input("compounds_dropdown", "value"),
            Input("fba_dropdown", "value"),
            Input("filter1_dropdown", "value"),
            Input("filter2_dropdown", "value"),
            ],
            prevent_initial_call=True,
        )
        def func(n_clicks, pathways_dropdown, element_dropdown, compounds_dropdown, fba_dropdown, filter1_dropdown, filter2_dropdown):
            if n_clicks is not None:
                nm, network = read_model(self._model, self._mm, element_dropdown)
                pathway_list, rxn_set = get_pathway_list(nm, pathways_dropdown)

                if isinstance(filter1_dropdown, str) and isinstance(filter2_dropdown, str):
                    rxn_list = set()
                    cpd_list = [filter1_dropdown, filter2_dropdown]
                    middle2 = []
                    middle3 = []
                    middle2, rxn_list1 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    middle3, rxn_list2 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    middle4, rxn_list3 = bfs_compounds(filter1_dropdown, filter2_dropdown, network, rxn_list, rxn_set, middle2, middle3)
                    rxn_list = rxn_list1 + rxn_list2 + rxn_list3
                    if rxn_list is str:
                        raise PreventUpdate
                elif isinstance(compounds_dropdown, str):
                    rxn_list = []
                    for rxn in network[0]:
                        for cpd in network[0][rxn][0]:
                            for i in cpd:
                                if i.name == compounds_dropdown and rxn.id in rxn_set:
                                    rxn_list.append(rxn.id)
                else:
                    rxn_list = rxn_set
                id=[]
                name=[]
                equation=[]
                for rxn in nm.reactions:
                    if rxn.id in rxn_list:
                        id.append(rxn.id)
                        name.append(rxn.name)
                        equation.append(str(rxn.properties["equation"]))
                df = pd.DataFrame({"id": id, "name": name, "equation": equation})
                return(dcc.send_data_frame(df.to_csv, "exported_rxns.csv"))

        @_app.callback(
            Output("net", "generateImage"),
            [Input('btn-get-png', 'n_clicks')
            ],
            prevent_initial_call=True,
        )
        def get_image(n_clicks):
            if n_clicks is not None:
                ftype = "png"
                action = "download"
                return {'type':ftype, 'action':action}

def read_model(model, mm, el):
        if isinstance(el, list) or el == None:
            el = "C"
        excluded_reactions=[]
        for rxn in model.reactions:
            for cpd, v in rxn.equation.compounds:
                if (float(v)).is_integer()==False or float(v) > 10:
                    excluded_reactions.append(rxn.id)
        network = make_network_dict(model, mm, subset=None, method='fpp',
                                    element=el, excluded_reactions=excluded_reactions,
                                    reaction_dict={}, analysis=None)
        return model, network

# Builds a reaction set from a model based on a pathway of interest
def get_pathway_list(nm, pathway):
    pathway_list=set()
    rxn_set=set()
    for path in pathway:
        if path != "All":

            for i in nm.reactions:

                if 'pathways' in i.properties:
                    for j in i.properties['pathways']:
                        pathway_list.add(j)
                    if path in i.properties['pathways']:
                        rxn_set.add(i.id)
                elif  'subsystem' in i.properties:
                    pathway_list.add(i.properties['subsystem'])
                    if path in i.properties['subsystem']:
                        rxn_set.add(i.id)

        else:
            for i in nm.reactions:

                if  'pathways' in i.properties:
                    for j in i.properties['pathways']:
                        pathway_list.add(j)
                elif  'subsystem' in i.properties:
                    pathway_list.add(i.properties['subsystem'])

                rxn_set.add(i.id)
    pathway_list=["All"] + list(pathway_list)
    return pathway_list, rxn_set

def get_compounds_list(nm):
    compounds_list=[]
    for i in nm.compounds:
        compounds_list.append(i.id)
    return compounds_list

def bfs_compounds(start, end, network, rxn_list, rxn_set, middle2, middle3):

    # initialize useful variables
    generic=['h2o', 'h', 'accoa', 'coa', 'nad', 'nadh', 'nadp', 'nadph', 'pi', \
             'co2', 'ppi', 'q8h2', 'q8', 'atp', 'gtp', 'utp', 'ttp']
    depth_start = {(start, "none"):0}
    depth_end={(end, "none"):0}
    i=1

    frontier_start=[(start, "none")]
    frontier_end=[(end, "none")]
    parent_start={}
    parent_start[(start, "none")]="none"
    parent_end={}
    parent_end[(end, "none")]="none"
    frontier_check="False"
    move_start={}
    move_end={}
    final_list={}

    while frontier_check=="False":
        next=[]
        for u in frontier_start:
            adj={}
            for rxn in network[0]:
                left=[]
                for l in rxn._properties['equation'].__dict__['_left']:
                    left.append(l[0].__dict__['_name'])
                right=[]
                for r in rxn._properties['equation'].__dict__['_right']:
                    right.append(r[0].__dict__['_name'])
                for cpd in network[0][rxn][0]:
                    for j in cpd:
                        if j.name==u[0]:
                            for k in cpd:
                                print(left, right)
                                print(j, k)
                                if j.name in left and k.name in right and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Forward" or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
                                    adj[k.name]=rxn.id
                                elif j.name in right and k.name in left and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse" or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
                                    adj[k.name]=rxn.id
            for v, r in adj.items():
                if (v, r) not in depth_start and v not in generic and (v, r) != middle2 and (v, r) != middle3:
                    depth_start[(v, r)]=i
                    parent_start[(v, r)]=u
                    next.append((v, r))
                    move_start[(v, r)]=adj[v]
                if (v, r) in depth_end and v not in generic and (v, r) != middle2 and (v, r) != middle3:
                    frontier_check="True"
                    middle=(v, r)
        frontier_start=next
        next=[]
        for u in frontier_end:
            adj={}
            for rxn in network[0]:
                left=[]
                for l in rxn._properties['equation'].__dict__['_left']:
                    left.append(l[0].__dict__['_name'])
                right=[]
                for r in rxn._properties['equation'].__dict__['_right']:
                    right.append(r[0].__dict__['_name'])
                for cpd in network[0][rxn][0]:
                    for j in cpd:
                        if j.name==u[0]:
                            for k in cpd:
                                if j.name in left and k.name in right and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Forward" or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
                                    adj[k.name]=rxn.id
                                elif j.name in right and k.name in left and (str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse" or str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both"):
                                    adj[k.name]=rxn.id
            for v, r in adj.items():
                if (v, r) not in depth_end and v not in generic and (v, r) != middle2 and (v, r) != middle3:
                    depth_end[(v, r)]=i
                    parent_end[(v, r)]=u
                    next.append((v, r))
                    move_end[(v, r)]=adj[v]
                if (v, r) in depth_start and frontier_check!="True" and v not in generic and (v, r) != middle2  and (v, r) != middle3:
                    frontier_check="True"
                    middle=(v, r)
        frontier_end=next
        i+=1
        if i > 50:
            return([], [])
    # collect the rxns from the start frontier into a final list of rxns
    i=depth_start[middle]

    j=1
    final_list[i]=move_start[middle]
    parent=parent_start[middle]
    while parent!=(start, "none"):
        i-=1
        j+=1
        final_list[i]=move_start[parent]
        parent=parent_start[parent]

    # Checks to make sure the solution isnt a 0 path solution
    if depth_end[middle] > 0:
        # Collects the moves from the end frontier into a final list of moves
        j+=1
        final_list[j]=move_end[middle]
        parent=parent_end[middle]
        while parent!=(end, "none"):
            j+=1
            final_list[j]=move_end[parent]
            parent=parent_end[parent]
    sorted_list=[]
    # Sorts the moves by their index and store them in a final list
    for k in range(1,len(final_list)+1):
        sorted_list.append(final_list[k])
    return(middle, sorted_list)
    #raise NotImplementedError




# Useful function to build the subset network. Returns nodes and edges
# from the network associated with the rxn_set of interest
def build_network(self, nm, rxn_set, network, fba_dropdown):
        name={}
        formula={}
        for i in nm.compounds:
            name[i.id]=i.name
            formula[i.id]=i.formula
        nodes=[]
        edges=[]
        obj = 0
        if not isinstance(fba_dropdown, list):
            objective=str(fba_dropdown)
            mm = nm.create_metabolic_model()

            # Set up the solver
            solver = self._get_solver()
            # Define the flux balance
            problem = fluxanalysis.FluxBalanceProblem(mm, solver)
            problem.maximize(fba_dropdown)
            obj = problem.get_flux(fba_dropdown)
            flux_carrying={}
            for i in nm.reactions:
                flux_carrying[i.id]=problem.get_flux(i.id)
        else:
            flux_carrying={}
            for i in nm.reactions:
                flux_carrying[i.id]=0


        for rxn in network[0]:
            if rxn.id in rxn_set:
                rxn_num=0
                for cpd in network[0][rxn][0]:
                    nodes.append({'data':{'id':str(cpd[0]),
                                  'label': name[str(cpd[0])[0:(str(cpd[0]).find('['))]],
                                  'formula': formula[str(cpd[0])[0:(str(cpd[0]).find('['))]]
                                  }})
                    nodes.append({'data':{'id':str(cpd[1]),
                                  'label': name[str(cpd[1])[0:(str(cpd[1]).find('['))]],
                                  'formula': formula[str(cpd[1])[0:(str(cpd[1]).find('['))]]
                                  }})
                    if  'pathways' in rxn.properties:
                        path = rxn.properties['pathways']
                    elif  'subsystem' in rxn.properties:
                        path = rxn.properties['subsystem']
                    else:
                        path = ['No pathway exists']
                    if rxn.id in edges:
                        if str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Forward":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[0]),
                                'target':str(cpd[1]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                        elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse":
                            flux = 0-float(flux_carrying[rxn.id])
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[1]),
                                'target':str(cpd[0]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                        elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[0]),
                                'target':str(cpd[1]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                            flux = 0-float(flux_carrying[rxn.id])
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[1]),
                                'target':str(cpd[0]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                    else:
                        if str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Forward":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[0]),
                                'target':str(cpd[1]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                        elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Reverse":
                            flux = 0-float(flux_carrying[rxn.id])
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[1]),
                                'target':str(cpd[0]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                        elif str(rxn._properties['equation'].__dict__['_direction']) == "Direction.Both":
                            flux = flux_carrying[rxn.id]
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[0]),
                                'target':str(cpd[1]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1
                            flux = 0-float(flux_carrying[rxn.id])
                            edges.append({'data':{
                                'id':"".join([rxn.id, "_", str(rxn_num)]),
                                'source':str(cpd[1]),
                                'target':str(cpd[0]),
                                'label': "".join([rxn.name, "_", str(rxn_num)]),
                                'equation': str(rxn.properties["equation"]),
                                'pathways':path,
                                'flux':flux
        #                            'equation':rxn.equation
                                }})
                            rxn_num+=1

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
