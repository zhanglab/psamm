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
from pkg_resources import resource_filename
import webbrowser

logger = logging.getLogger(__name__)

# PAGES -------------------------------------------------------------------

#dash.register_page("simulate",  path='/', layout=html.Div('Home Page'))
#dash.register_page("curate", layout=html.Div('Analytics'))

REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'
EPSILON = 1e-5

class HiddenPrints:
    """
    Class that suppresses print statements and avoids cluttering the
    terminal when running commands such as charge and formulacheck
    """
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


class InteractiveCommand(MetabolicMixin, SolverCommandMixin,
                         Command, FilePrefixAppendAction):
    """Launches an interactive interface for the analysis """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', type=convert_to_unicode,
            default=[], action=FilePrefixAppendAction,
            help='Reaction(s) to exclude from metabolite pair prediction')
        super(InteractiveCommand, cls).init_parser(parser)

    def run(self):
        # bootstrap theme CDN
        #bs_theme = "https://cdn.jsdelivr.net/npm/bootswatch@5.2.3/dist/lux/bootstrap.min.css"
        self._app = dash.Dash(
            __name__,  
                external_stylesheets=[dbc.themes.COSMO]
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
        def reaction_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'reaction_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('gene_dialog', 'displayed'),
            Input('gene_help', 'n_clicks')
        )
        def gene_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'gene_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('fba_dialog', 'displayed'),
            Input('fba_help', 'n_clicks')
        )
        def fba_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'fba_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('bfs_dialog', 'displayed'),
            Input('bfs_help', 'n_clicks')
        )
        def bfs_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'bfs_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('compound_dialog', 'displayed'),
            Input('compound_help', 'n_clicks')
        )
        def compound_help(nclicks):
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'compound_help.n_clicks':
                return True
            else:
                return False

        @_app.callback(
            Output('element_dialog', 'displayed'),
            Input('element_help', 'n_clicks')
        )
        def element_help(nclicks):
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
                        self._model._limits[c[0]] = (c[0], Decimal(c[1]),
                                                     Decimal(c[2]))
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
                                "EX_{}[{}]".format(e[0], e[1])) is False:
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
                                          placeholder=float(
                                            self._model._exchange[e][2]),
                                          style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="upper{}".format(str(e.name)),
                                          type="number",
                                          placeholder=float(
                                            self._model._exchange[e][3]),
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
                        with open('{}/{}'.format(path, i['include']), "w") \
                                as f:
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
                    for le in left_out:
                        for c in self._model.compounds:
                            if c.id.upper() == le[0].upper():
                                compound = c
                        compound.charge = int(le[3])
                        compound.formula = le[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(le[0]).in_compartment(
                            le[4])] = int(le[1]) * -1

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
                    unbalanced, count, unchecked, exclude, form_list, \
                            form_dict = formula_check(self)
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
                                r.equation.translated_compounds(
                                    compound_name))))
                    try:
                        contents.append(
                            dbc.Row([
                                dbc.Col([
                                    html.P("Left side total: "
                                           "{}".format(form_dict[r.id][0])),
                                    html.P("Left side off by: "
                                           "{}".format(form_dict[r.id][1])),
                                ]),
                                dbc.Col([
                                    html.P("Right side total: "
                                           "{}".format(form_dict[r.id][2])),
                                    html.P("Right side off by: "
                                           "{}".format(form_dict[r.id][3])),
                                ]),
                            ]),

                        )
                    except KeyError:
                        contents = [
                            "{} was successfully balanced! Select a reaction "
                            "above to curate the charge balance".format(
                                r.id)]
                        left_contents = [
                            "The left side of the reaction will display here"]
                        right_contents = [
                            "The right side of the reaction will display here"]
                        return contents, left_contents, right_contents, [
                            {"label": i, "value": i, } for i in
                            list(form_list)]

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
                                    dcc.Input(id="stoich{}".format(
                                        str(cur_id)),
                                              type="number",
                                              placeholder=str(i[1]),
                                              style={'width': '75%'})],
                                        width=2),
                                dbc.Col([
                                    dcc.Input(id="formula{}".format(
                                        str(cur_id)),
                                              type="text",
                                              placeholder=str(
                                                formula_dict[i[0].name]),
                                              style={'width': '85%'})],
                                        width=3),
                                dbc.Col([
                                    dcc.Input(id="charge{}".format(
                                        str(cur_id)),
                                              type="number",
                                              placeholder=str(
                                                charge_dict[i[0].name]),
                                              style={'width': '75%'})],
                                        width=2),
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
                                          placeholder="0",
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder='',
                                          style={'width': '85%'})], width=3),
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
                                    dcc.Input(id="stoich{}".format(
                                        str(cur_id)),
                                              type="number",
                                              placeholder=str(i[1]),
                                              style={'width': '75%'})],
                                        width=2),
                                dbc.Col([
                                    dcc.Input(id="formula{}".format(str(
                                        cur_id)),
                                              type="text",
                                              placeholder=str(
                                                formula_dict[i[0].name]),
                                              style={'width': '85%'})],
                                        width=3),
                                dbc.Col([
                                    dcc.Input(id="charge{}".format(
                                        str(cur_id)),
                                              type="number",
                                              placeholder=str(
                                                charge_dict[i[0].name]),
                                              style={'width': '75%'})],
                                        width=2),
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
                                          placeholder="0",
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder='',
                                          style={'width': '85%'})], width=3),
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
                        "Select a reaction above to curate the "
                        "formula balance"]
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
                if dash.callback_context.triggered[0]['prop_id'] == \
                        'btn_save_charge.n_clicks':
                    left_out = []
                    eq_dict = {}
                    for i in leftdata:
                        if isinstance(i, dict):
                            left_out_temp = rec_props(
                                i['props']['children'], [])
                            if len(left_out_temp) > 0:
                                if left_out_temp[0] != '':
                                    left_out.append(left_out_temp)
                    for le in left_out:
                        for c in self._model.compounds:
                            if c.id.upper() == le[0].upper():
                                compound = c
                        compound.charge = int(l[3])
                        compound.formula = le[2]
                        self._model.compounds.add_entry(compound)
                        eq_dict[Compound(le[0]).in_compartment(
                            le[4])] = int(le[1]) * -1

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
                unbalanced, count, unchecked, exclude, unbalanced_list, \
                    unbalance_dict = charge_check(self, 1e-6)
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
                        "{} was successfully balanced! Select a reaction "
                        "above to curate the charge balance".format(
                            r.id)]
                    left_contents = [
                        "The left side of the reaction will display here"]
                    right_contents = [
                        "The right side of the reaction will display here"]
                    return contents, left_contents, right_contents, [
                        {"label": i, "value": i, } for
                        i in list(unbalanced_list)]

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
                                          placeholder=str(i[1]),
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder=str(
                                            formula_dict[i[0].name]),
                                          style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(
                                            charge_dict[i[0].name]),
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
                                multi=False)], width=3),
                        dbc.Col([
                            dcc.Input(id="stoich{}".format(str(cur_id)),
                                      type="number",
                                      placeholder="0",
                                      style={'width': '75%'})],
                                width=2),
                        dbc.Col([
                            dcc.Input(id="formula{}".format(str(cur_id)),
                                      type="text",
                                      placeholder='', style={'width': '85%'})],
                                width=3),
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
                                          placeholder=str(i[1]),
                                          style={'width': '75%'})], width=2),
                            dbc.Col([
                                dcc.Input(id="formula{}".format(str(cur_id)),
                                          type="text",
                                          placeholder=str(
                                            formula_dict[i[0].name]),
                                          style={'width': '85%'})], width=3),
                            dbc.Col([
                                dcc.Input(id="charge{}".format(str(cur_id)),
                                          type="number",
                                          placeholder=str(
                                            charge_dict[i[0].name]),
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
                                      placeholder="0",
                                      style={'width': '75%'})], width=2),
                        dbc.Col([
                            dcc.Input(id="formula{}".format(str(cur_id)),
                                      type="text",
                                      placeholder='',
                                      style={'width': '85%'})], width=3),
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
            if dash.callback_context.triggered[0]['prop_id'] == \
                    'btn_update_cur_net.n_clicks':
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
        def display_nodedata(clicks, r, fc):
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

        @_app.callback([Output("download_r", "data"),
                        Output("download_c", "data"),],
                       Input("btn_tsv", "n_clicks"),
                       State("net", "elements"),
                       prevent_initial_call=True,)
        def export_reactions(n_clicks, elements):
            if n_clicks is not None:
                compounds = set()
                reactions = set()
                for e in elements:
                    if 'type' in e['data']:
                        if e['data']['type'] == 'compound':
                            compounds.add(e['data']['orig_id'])
                        elif e['data']['type'] == 'reaction':
                            reactions.add(e['data']['orig_id'])

                rdict = defaultdict(lambda:[])
                props = set()
                for rxn in self._model.reactions:
                    if rxn.id in reactions:
                        for p in rxn.properties:
                            props.add(p)
                for rxn in self._model.reactions:
                    if rxn.id in reactions:
                        for p in props:
                            if p in rxn.properties:
                                rdict[p].append(rxn.properties[p])
                            else:
                                rdict[p].append("NA")

                cdict = defaultdict(lambda:[])
                props = set()
                for cpd in self._model.compounds:
                    if cpd.id in compounds:
                        for p in cpd.properties:
                            props.add(p)
                for cpd in self._model.compounds:
                    if cpd.id in compounds:
                        for p in props:
                            print(p)
                            if p in cpd.properties:
                                cdict[p].append(cpd.properties[p])
                            else:
                                cdict[p].append("NA")
                rdf = pd.DataFrame(rdict)
                print(cdict)
                cdf = pd.DataFrame(cdict)
                return(dcc.send_data_frame(rdf.to_csv, "exported_rxns.csv"),
                       dcc.send_data_frame(cdf.to_csv, "exported_cpds.csv"))

        @_app.callback(Output("net", "generateImage"),
                       [Input('btn-get-png', 'n_clicks')],
                       prevent_initial_call=True,)
        def get_image(n_clicks):
            if n_clicks is not None:
                ftype = "png"
                action = "download"
                return {'type': ftype, 'action': action}

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
    # Define components of the curate tab
    add_tab = dcc.Tab(label='Add', value='add',
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
                      ]),)
    c_n = cyto.Cytoscape(id='net_curate',
                         layout={'name': 'cose',
                                 'animate': True,
                                 'animationDuration': 1000},
                         style={'width': '100%', 'height': '600px'},
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
                                     {"selector": "[flux > 0]",
                                      "style": {"line-color": "blue",
                                                "target-arrow-color":
                                                "blue",
                                                "width": 1,
                                                "curve-style": "bezier",
                                                "target-arrow-shape":
                                                "triangle",
                                                "line-dash-pattern":
                                                [6, 3],
                                                "line-dash-offset": 24}}],
                         minZoom=0.06)
    e_child = html.Div(className='control-tab',
                       children=[
                        dbc.Row([
                            dbc.Col([
                                html.H5("Reaction to edit:"),
                                dbc.Row([
                                    dcc.Dropdown(
                                        id="reaction_dropdown_curate",
                                        options=rxns_full,
                                        value='',
                                        multi=False,
                                        placeholder="",
                                        style={"width": "100%"})]),
                                html.H5("Compound to edit:"),
                                dbc.Row([
                                    dcc.Dropdown(
                                        id="compound_dropdown_curate",
                                        options=[{"label": i,
                                                  "value": i}
                                                 for i in
                                                 list(compounds_list)],
                                        value='',
                                        multi=False,
                                        placeholder="",
                                        style={"width": "100%"}), ]),
                                html.Div([html.Button("Confirm Edit",
                                                      id="btn_update")]),
                                dbc.Alert(id="edit-data",
                                          children="Click on a node",
                                          color="secondary"),
                                html.Div([html.Button("Save",
                                                      id="btn_save")]),
                                dbc.Alert(id="edit-confirm",
                                          children="Click save to confirm "
                                                   "changes to model",
                                          color="secondary")]),
                            dbc.Col([
                                html.H5("Pathways:"),
                                dbc.Row([
                                    dcc.Dropdown(
                                        id="pathways_dropdown_curate",
                                        options=pathway_list,
                                        value=pathway_list[0],
                                        multi=True,
                                        placeholder="",
                                        style={"width": "100%"}, ), ]),
                                html.Div([
                                    html.Button("Submit",
                                                id="btn_update_cur_net")]),
                                c_n,
                                               ]),
                                      ]),
                             ])
    e_t = dcc.Tab(label='Edit',
                  value='edit',
                  children=e_child, )

    chargecheck = dcc.Tab(label='Charge Check', value='chargeCheck',
                          children=html.Div(className='c-tab', children=[
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
                                      children="Select a reaction above "
                                      "to curate the charge balance",
                                      color="secondary",
                                  ),
                                      dbc.Row([
                                          dbc.Col([
                                              dbc.Alert(
                                                  id="chargeLeftCuration",
                                                  children="The left side "
                                                  "of the reaction will "
                                                  "display here",
                                                  color="secondary",
                                              ),
                                              html.Div([
                                                  html.Button(
                                                      "Save Changes",
                                                      id="btn_save_charge")
                                              ]),
                                          ]),
                                          dbc.Col([
                                              dbc.Alert(
                                                  id="chargeRightCuration",
                                                  children="The right "
                                                  "side of the reaction "
                                                  "will display here",
                                                  color="secondary",
                                              ),
                                          ]),
                                      ]),
                                  ]),

                          ]),)

    fc = dcc.Tab(label='Formula Check', value='formulaCheck',
                 children=html.Div(className='control-tab',
                                   children=[dbc.Row(
                                        [html.H5("Unbalanced Reactions:"),
                                            dcc.Dropdown(
                                            id="unbalanced_rxns_formula",
                                            options=[{"label": i,
                                                      "value": i, }
                                                     for i in
                                                     list(form_list)],
                                            value=pathway_list[0],
                                            multi=False,
                                            placeholder="",),
                                            dbc.Alert(
                                            id="formCheckCuration",
                                            children="Select a reaction "
                                            "above to curate the formula "
                                            "balance",
                                            color="secondary",),
                                            dbc.Row([
                                                dbc.Col([
                                                    dbc.Alert(
                                                        id="formLeft"
                                                           "Curation",
                                                        children="The "
                                                                 "left "
                                                                 "side of "
                                                                 "the "
                                                                 "reaction"
                                                                 " will "
                                                                 "display "
                                                                 "here",
                                                        color="secondary"),
                                                    html.Div([
                                                        html.Button(
                                                            "Save Changes",
                                                            id="btn_save"
                                                            "_form")]), ]),
                                                dbc.Col([
                                                    dbc.Alert(
                                                        id="formRight"
                                                           "Curation",
                                                        children="The "
                                                                 "right "
                                                                 "side of "
                                                                 "the "
                                                                 "reaction"
                                                                 " will "
                                                                 "display "
                                                                 "here",
                                                        color="secondary"),
                                                        ]), ]), ]), ]), )


    # Define the main body of the app
    body_layout = dbc.Container(
        [dbc.Row([
            dbc.Row([

                # how do I make the callback only display one of these
                # maybe make them two different

                html.Div(id='modelling', className='modelling', children=[
                    dcc.Tabs(id="parent-tabs", children=[
                        dcc.Tab(id="sim", label='Simulate', style={
                            'width': '30%', 'font-size': '150%', 'height': '50%'}, children=[
                            dbc.Row([
                                dbc.Col([
                                    html.Div(id='model-control-tabs',
                                             className='control-tabs',
                                             style={'width': '69%', 'font-size': '50%', 
                                             'height': '50%'}, 
                                              children=[
                                                 simulate_tabs 
                                             ]), ]),
                                dbc.Col([
                                    s_n,
                                    dbc.Row([
                                        dbc.Alert(
                                            id="node-data",
                                            children="Click on a node to "
                                                     "see its details",
                                            color="secondary",
                                        ),
                                    ]), ]),
                                dbc.Row([dcc.Markdown(
                                    '''
                                    \\* Data analysis carried out for
                                    demonstration of data visualisation
                                    purposes only.
                                    ''')], style={"fontSize": 11,
                                                  "color": "gray"},), ],
                                    style={"marginTop": 20},)]),
                        dcc.Tab(label='Curate', id='cur', style={
                            'width': '100%', 'font-size': '150%', 'height': '50%'},
                                children=[
                                    html.Div(id='curate-control-tabs',
                                             className='curate-tabs',
                                             children=[dcc.Tabs(
                                                id="curate-subtabs",
                                                value="what is",
                                                children=[add_tab,
                                                          e_t,
                                                          chargecheck,
                                                          fc, ]),

                                                       ]),
                                    dbc.Row([
                                        dbc.Alert(
                                            id="save_confirmation",
                                            children="Press Save Model "
                                            "to save the model",
                                            color="secondary",
                                        ),
                                        dbc.Row([
                                            dbc.Col([
                                                html.Button("Save Model",
                                                    id="btn_save_model"),
                                            ]), ]), ],
                                        style={"fontSize": 16,
                                               "color": "gray"},),
                                ]),
                    ]),
                ]), ]),
        ]),
        ], 
        # take up the whole screen
        fluid = True, style={"height": "100vh"} 
        )

    layout = html.Div([navbar, body_layout])
    return layout





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
                                  'orig_id': str(cpd[0])
                                                [0:(str(cpd[0]).
                                                 find('['))],
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
                                  'orig_id': str(cpd[0])
                                                [0:(str(cpd[0]).
                                                 find('['))],
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
