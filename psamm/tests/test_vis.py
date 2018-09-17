#!/usr/bin/env python
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
# Copyright 2018-2018  Ke Zhang <kzhang@my.uri.edu>

from __future__ import unicode_literals

import unittest
from psamm.formula import Formula, Atom, Radical
from psamm.commands import vis
from psamm.reaction import Compound, Reaction, Direction
from psamm.datasource import entry
from collections import defaultdict, OrderedDict
from psamm.database import DictDatabase, ChainedDatabase
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm.datasource import native, context, entry
from psamm.metabolicmodel import MetabolicModel
from psamm import graph
from psamm.datasource.native import NativeModel, ReactionEntry, CompoundEntry
from six import iteritems


class TestMakeFilterDict(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1','equation': parse_reaction('fum_c[c] + h2o_c[c] '
                                                    '<=> mal_L_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'fum_c[c]', 'formula': parse_compound('C4H2O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o_c[c]', 'formula': parse_compound('H2O', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'mal_L_c[c]', 'formula': parse_compound('C4H4O5', 'c')}))


class TestMakeCpairDict(unittest.TestCase):      ##### OK.
    def setUp(self):
        database = DictDatabase()
        database.set_reaction('rxn1', parse_reaction('A[c] + B[c] <=> C[c] + D[c]'))
        database.set_reaction('rxn2', parse_reaction('B[c] <=> (2) D[c]'))
        database.set_reaction('rxn3', parse_reaction('D[c] => D[e]'))
        database.set_reaction('rxn4', parse_reaction('D[e] <=>'))
        # database.set_reaction('rxn5', parse_reaction('A[c] + D[c] <=> C[e] + E'))
        # self.mm.set_reaction('rxn_3', parse_reaction('A => D[e]'))
        # self.mm.set_reaction('rxn_4', parse_reaction('A + 2 B => C'))
        # self.mm.set_reaction('rxn_5', parse_reaction('C => 3 D[e]'))
        self.mm = MetabolicModel.load_model(database)

        self.filter_dict = {'rxn1': [(Compound('A', 'c'), Compound('C', 'c')),
                                     (Compound('B', 'c'), Compound('C', 'c')),
                                     (Compound('B', 'c'), Compound('D', 'c'))],
                            'rxn2': [(Compound('B', 'c'), Compound('D', 'c'))],
                            'rxn3': [(Compound('D', 'c'), Compound('D', 'e'))]}
        self.subset = ['rxn1', 'rxn2', 'rxn3', 'rxn4']
        self.reaction_flux = []   #no fba
        self.method = 'fpp'

    def test_MakeCpairDict_default(self):
        e1, n1 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                  self.reaction_flux, self.method)
        e1_res = defaultdict(lambda: defaultdict(list))
        e1_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'].append('rxn1_1')
        e1_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'].append('rxn1_2')
        # e1_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn1_3','rxn2_1']
        e1_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn1_3', 'rxn2_1']
        e1_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].append('rxn3_1')

        n1_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1','rxn1_3': 'rxn1', 'rxn2_1':'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e1, e1_res)
        self.assertEqual(n1, n1_res)

    def test_MakeCpairDict_nofpp(self):
        e2, n2 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                  self.reaction_flux, 'no-fpp')
        e2_res = defaultdict(lambda: defaultdict(list))
        e2_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'] = ['rxn1']
        e2_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'] = ['rxn1']
        e2_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = \
            ['rxn1', 'rxn2']
        e2_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'] = ['rxn3']

        n2_res = {'rxn1': 'rxn1', 'rxn2': 'rxn2', 'rxn3': 'rxn3'}
        self.assertEqual(e2, e2_res)
        self.assertEqual(n2, n2_res)

    def test_MakeCpairDict_fba(self):
        reaction_flux = {'rxn1': 9.8, 'rxn3': 9.8, 'rxn4':9.8}
        e3, n3 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                     reaction_flux, self.method)
        e3_res = defaultdict(lambda: defaultdict(list))
        e3_res[(Compound('A', 'c'), Compound('C', 'c'))]['forward'].append('rxn1_1')
        e3_res[(Compound('B', 'c'), Compound('C', 'c'))]['forward'].append('rxn1_2')
        e3_res[(Compound('B', 'c'), Compound('D', 'c'))]['forward'] = ['rxn1_3']
        e3_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn2_1']
        e3_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].append('rxn3_1')

        n3_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                  'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e3, e3_res)
        self.assertEqual(n3, n3_res)

    def test_MakeCpairDict_subset_fba(self):
        reaction_flux = {'rxn1': 9.8, 'rxn3': 9.8, 'rxn4':9.8}
        subset = ['rxn1', 'rxn2', 'rxn4']
        e4, n4 = vis.make_cpair_dict(self.mm, self.filter_dict, subset,
                                     reaction_flux, self.method)
        e4_r = defaultdict(lambda: defaultdict(list))
        e4_r[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] = [str('rxn1_1')]
        e4_r[(Compound('B', 'c'), Compound('C', 'c'))]['forward'] = ['rxn1_2']
        e4_r[(Compound('B', 'c'), Compound('D', 'c'))]['forward'] = ['rxn1_3']
        e4_r[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn2_1']

        n4_r = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                'rxn2_1': 'rxn2'}
        print('e4', e4)
        print('en',e4_r)
        self.assertEqual(e4, e4_r)
        self.assertEqual(n4, n4_r)

    def test_MakeCpairDict_filterDict(self):
        filter_dict = {'rxn1': [(Compound('A', 'c'), Compound('C', 'c')),
                                (Compound('B', 'c'), Compound('C', 'c')),
                                (Compound('B', 'c'), Compound('D', 'c')),
                                (Compound('A', 'c'), Compound('D', 'c')),],
                       'rxn2': [(Compound('B', 'c'), Compound('D', 'c'))],
                       'rxn3': [(Compound('D', 'c'), Compound('D', 'e'))]}
        e5, n5 = vis.make_cpair_dict(self.mm, filter_dict, self.subset,
                                  self.reaction_flux, self.method)
        e5_res = defaultdict(lambda: defaultdict(list))
        e5_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'].append('rxn1_1')
        e5_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'].append('rxn1_2')
        e5_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn1_3','rxn2_1']
        e5_res[(Compound('A', 'c'), Compound('D', 'c'))]['both'].append('rxn1_4')
        e5_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].append('rxn3_1')

        n5_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1','rxn1_3': 'rxn1',
                  'rxn1_4': 'rxn1', 'rxn2_1':'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e5, e5_res)
        self.assertEqual(n5, n5_res)


class TestMakeEdgeValues(unittest.TestCase):      ##### OK.
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('FUM', parse_reaction('fum_c[c] + h2o_c[c] <=> mal_L_c[c]'))
        self.database.set_reaction('CYTBD', parse_reaction(
            '(2) h_c[c] + (0.5) o2_c[c] + q8h2_c[c] => h2o_c[c] + '
            '(2) h_e[e] + q8_c[c]'))
        self.database.set_reaction('FRD7', parse_reaction('fum_c[c] + q8h2_c[c] => q8_c[c] + succ_c[c]'))
        self.mm = MetabolicModel.load_model(
            self.database, self.database.reactions)
        self.reaction_flux = {'FUM': 5.06, 'CYTBD': 43.60}
        self.compound_formula = {
            'h2o_c': Formula({Atom('O'): 1, Atom('H'):2}),
            'fum_c': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 2}),
            'h_e': Formula({Atom('O'): 1, Atom('H'): 2}),
            'fum_e': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 2}),
            'o2_c': Formula({Atom('O'): 2}), 'h_c': Formula({Atom('H'): 1}),
            'mal_L_c': Formula({Atom('O'): 5, Atom('C'): 4, Atom('H'): 4}),
            'q8_c': Formula({Atom('O'): 4, Atom('C'): 49, Atom('H'): 74}),
            'q8h2_c': Formula({Atom('O'): 4, Atom('C'): 49, Atom('H'): 76}),
            'succ_c': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 4})
        }
        self.element = 'C'
        self.split_map = False
        self.cpair_dict = defaultdict(lambda: defaultdict(list))
        self.cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'].append('CYTBD_1')
        self.cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'].append('FRD7_1')
        self.cpair_dict[(Compound('fum_c', 'c'), Compound(
            'succ_c', 'c'))]['forward'].append('FRD7_2')
        self.cpair_dict[(Compound('fum_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_1')

        self.new_id_mapping = {'CYTBD_1': 'CYTBD', 'FRD7_1': 'FRD7',
                               'FRD7_2': 'FRD7', 'FUM_1': 'FUM',
                               'FUM_2': 'FUM', 'CYTBD_2': 'CYTBD'}
        self.method = 'fpp'

    def test_edge_values_basic(self):
        e1 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            self.split_map,self.cpair_dict, self.new_id_mapping, self.method)
        e1_res = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
              (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
              (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
              (Compound('fum_c', 'c'), 'FUM'): 5.06
              }
        self.assertEqual(e1, e1_res)

    def test_edge_values_split_map(self):
        e2 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            True, self.cpair_dict, self.new_id_mapping, self.method)
        e2_res = {(Compound('q8_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06
                  }
        self.assertEqual(e2, e2_res)

    def test_edge_values_nofpp(self):
        e3 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, 'no-fpp')
        e3_res = {(Compound('q8_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06
                  }
        self.assertEqual(e3, e3_res)

    def test_edge_values_nofba(self):
        e4 = vis.make_edge_values(
            {}, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, self.method)
        e4_res = {}
        self.assertEqual(e4, e4_res)

    def test_edge_values_elementO(self):
        new_cpair_dict = defaultdict(lambda: defaultdict(list))
        new_cpair_dict[(Compound('o2_c', 'c'), Compound(
            'h2o_c', 'c'))]['forward'].append('CYTBD_1')
        new_cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'] = ['CYTBD_2','FRD7_1']
        new_cpair_dict[(Compound('fum_c', 'c'), Compound(
            'succ_c', 'c'))]['forward'].append('FRD7_2')
        new_cpair_dict[(Compound('fum_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_1')
        new_cpair_dict[(Compound('h2o_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_2')

        e5 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, 'O',
            self.split_map, new_cpair_dict, self.new_id_mapping, self.method)
        e5_res = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06,
                  (Compound('h2o_c', 'c'), 'FUM'): 5.06,
                  (Compound('o2_c', 'c'), 'CYTBD'): 21.80,
                  (Compound('h2o_c', 'c'), 'CYTBD'): 43.60
                  }
        self.assertEqual(e5, e5_res)


class TestAddGraphNodes(unittest.TestCase):
    def setUp(self):
        cpair_dict = defaultdict(lambda: defaultdict(list))
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['both'] = ['rxn1_1','rxn2_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] = ['rxn3_1']
        self.cpair_dict = cpair_dict
        self.new_id_mapping = {'rxn1_1': 'rxn1', 'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.method = 'fpp'
        self.split = False
        self.g = graph.Graph()

    def test_addnodes_default(self):
        g1 = vis.add_graph_nodes(self.g, self.cpair_dict, self.method,
                                 self.new_id_mapping, self.split)
        node_a1 = graph.Node({'id': 'A[c]', 'shape': 'ellipse','style': 'filled',
                             'type': 'cpd', 'original_id': Compound('A', 'c'),
                             'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_c1 = graph.Node({'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
                             'type': 'cpd', 'original_id': Compound('C', 'c'),
                              'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac1_both = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box','style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1', 'rxn2'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        node_ac1_forward = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a1, node_c1, node_ac1_both, node_ac1_forward] for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in [node_a1, node_c1, node_ac1_both, node_ac1_forward]))

    def test_addnodes_split(self):
        g2 = vis.add_graph_nodes(self.g, self.cpair_dict, self.method,
                                 self.new_id_mapping, True)
        node_a2 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_c2 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a2, node_c2, node_rxn1, node_rxn2,
                                  node_rxn3] for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in [
            node_a2, node_c2, node_rxn1, node_rxn2, node_rxn3]))

    def test_addnodes_nofpp(self):
        g3 = vis.add_graph_nodes(self.g, self.cpair_dict, 'no-fpp',
                                 self.new_id_mapping, self.split)
        node_a3 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled','type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_c3 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_rxn1 = graph.Node({
            'id': 'rxn1_1','shape': 'box','style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a3, node_c3, node_rxn1, node_rxn2,
                                  node_rxn3] for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in [
            node_a3, node_c3, node_rxn1, node_rxn2, node_rxn3]))


class TestUpdateNodeColor(unittest.TestCase):
    def setUp(self):
        node_a = graph.Node({'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
                             'type': 'cpd', 'original_id': Compound('A', 'c'),
                             'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_ac = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',  # 'id'?
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        node_ac_comb = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',  # 'id'?
            'original_id': ['rxn1', 'rxn2'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.g = graph.Graph()
        self.node_list = [node_a, node_c, node_ac]
        for node in self.node_list:
            self.g.add_node(node)

        self.g_comb = graph.Graph()
        self.node_list_comb = [node_a, node_c, node_ac_comb]
        for node in self.node_list_comb:
            self.g_comb.add_node(node)

        self.recolor_dict = {'A[c]': '#f4fc55', 'rxn1': '#ef70de'}

    def test1_add_color_emptyrecolordict(self):
        g1 = vis.update_node_color(self.g, {})
        node_a1 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_c1 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a1, node_c1, node_ac1]
                            for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in
                            [node_a1, node_c1, node_ac1]))

    def test2_add_color_recolor(self):
        g2 = vis.update_node_color(self.g, self.recolor_dict)
        node_a2 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#f4fc55'})
        node_c2 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac2 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#ef70de'})

        self.assertTrue(all(i in [node_a2, node_c2, node_ac2] for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in [node_a2, node_c2, node_ac2]))

    def test3_update_color_recolorCombinedRxnNode(self):
        g3 = vis.update_node_color(self.g_comb, self.recolor_dict)
        node_a3 = graph.Node({'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
                             'type': 'cpd', 'original_id': Compound('A', 'c'),
                             'compartment': 'c', 'fillcolor': '#f4fc55'})
        node_c3 = graph.Node({'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
                             'type': 'cpd', 'original_id': Compound('C', 'c'),
                             'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac3 = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1', 'rxn2'], 'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a3, node_c3, node_ac3] for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in [node_a3, node_c3, node_ac3]))


### def add_edges(g, cpairs_dict, method, split):
class TestAddEdges(unittest.TestCase):
    def setUp(self):
        self.node_a = graph.Node({'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
                              'type': 'cpd', 'original_id': Compound('A', 'c'),
                              'compartment': 'c', 'fillcolor':'#ffd8bf'})
        self.node_c = graph.Node({'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
                              'type': 'cpd', 'original_id': Compound('C', 'c'),
                              'compartment': 'c', 'fillcolor':'#ffd8bf'})
        self.node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor':'#c9fccd'})
        self.node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'c', 'fillcolor':'#c9fccd'})
        self.node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor':'#c9fccd'})
        self.node_rxn4 = graph.Node({
            'id': 'rxn4_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn4'], 'compartment': 'c', 'fillcolor':'#c9fccd'})

        self.node_ac_both = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1', 'rxn2'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.node_ac_forw = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor':'#c9fccd'})
        self.node_ac_back = graph.Node({
            'id': 'rxn4_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn4'], 'compartment': 'c', 'fillcolor':'#c9fccd'})

        self.g1 = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_ac_both,
                     self.node_ac_forw, self.node_ac_back]:
            self.g1.add_node(node)          #graph for default setting

        self.g2 = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_rxn1,
                     self.node_rxn2, self.node_rxn3, self.node_rxn4]:
            self.g2.add_node(node)          #graph for split or no-fpp

        cpair_dict = defaultdict(lambda: defaultdict(list))
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['both'] = ['rxn1_1', 'rxn2_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] = ['rxn3_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['back'] = ['rxn4_1']
        self.cpair_dict = cpair_dict
        # self.new_id_mapping = {'rxn1_1': 'rxn1', 'rxn2_1': 'rxn2',
        #                        'rxn3_1': 'rxn3', 'rxn4_1': 'rxn4'}
        self.split = False
        self.method = 'fpp'

    def test1_addedges_default_or_filepath(self):
        g1_default = vis.add_edges(self.g1, self.cpair_dict, self.method, self.split)
        g1_filepath = vis.add_edges(self.g1, self.cpair_dict, 'file_path', self.split)
        edge1 = graph.Edge(self.node_a, self.node_ac_both, {'dir': 'both'})
        edge2 = graph.Edge(self.node_ac_both, self.node_c, {'dir': 'both'})
        edge3 = graph.Edge(self.node_a, self.node_ac_forw, {'dir': 'forward'})
        edge4 = graph.Edge(self.node_ac_forw, self.node_c, {'dir': 'forward'})
        edge5 = graph.Edge(self.node_a, self.node_ac_back, {'dir': 'back'})
        edge6 = graph.Edge(self.node_ac_back, self.node_c, {'dir': 'back'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6]
        self.assertTrue(all(i in edge_list for i in g1_default.edges))
        self.assertTrue(all(i in g1_default.edges for i in edge_list))

        self.assertTrue(all(i in edge_list for i in g1_filepath.edges))
        self.assertTrue(all(i in g1_filepath.edges for i in edge_list))

    def test2_addedges_split_or_nofpp(self):
        g2_split = vis.add_edges(self.g2, self.cpair_dict, self.method, True)

        g2_nofpp = vis.add_edges(self.g2, self.cpair_dict, 'no-fpp', self.split)

        edge1 = graph.Edge(self.node_a, self.node_rxn1, {'dir': 'both'})
        edge2 = graph.Edge(self.node_rxn1, self.node_c, {'dir': 'both'})
        edge3 = graph.Edge(self.node_a, self.node_rxn3, {'dir': 'forward'})
        edge4 = graph.Edge(self.node_rxn3, self.node_c, {'dir': 'forward'})
        edge5 = graph.Edge(self.node_a, self.node_rxn4, {'dir': 'back'})
        edge6 = graph.Edge(self.node_rxn4, self.node_c, {'dir': 'back'})
        edge7 = graph.Edge(self.node_a, self.node_rxn2, {'dir': 'both'})
        edge8 = graph.Edge(self.node_rxn2, self.node_c, {'dir': 'both'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8]

        self.assertTrue(all(i in edge_list for i in g2_split.edges))
        self.assertTrue(all(i in g2_split.edges for i in edge_list))

        self.assertTrue(all(i in edge_list for i in g2_nofpp.edges))
        self.assertTrue(all(i in g2_nofpp.edges for i in edge_list))


# def add_biomass_rxns(g, bio_reaction, biomass_rxn_id):
class TestAddBiomassRnxs(unittest.TestCase):
    def setUp(self):
        self.node_a = graph.Node({'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
                                  'type': 'cpd', 'original_id': Compound('A', 'c'),
                                  'compartment': 'c', 'fillcolor':'#ffd8bf'})
        # self.node_b = graph.Node({'id': 'B[b]', 'shape': 'ellipse', 'style': 'filled',
        #                           'type': 'cpd', 'original_id': Compound('B', 'c'),
        #                           'compartment': 'c'})
        self.node_c = graph.Node({'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
                                  'type': 'cpd', 'original_id': Compound('C', 'c'),
                                  'compartment': 'c', 'fillcolor':'#ffd8bf'})
        self.node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor':'#c9fccd'})
        # self.node_rxn1_2 = graph.Node({
        #     'id': 'rxn1_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
        #     'original_id': ['rxn1'], 'compartment': 'c'})
        self.edge_A_rxn = graph.Edge(self.node_a, self.node_rxn1, {'dir': 'both'})
        self.edge_C_rxn = graph.Edge(self.node_rxn1, self.node_c, {'dir': 'both'})
        self.g = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_rxn1]:
            self.g.add_node(node)
        self.g.add_edge(self.edge_A_rxn)
        self.g.add_edge(self.edge_C_rxn)

        self.biomass_id = 'test_bio'
        self.biomass_rxn = Reaction(
            Direction.Forward, [(Compound('C', 'c'), 1), (Compound('D', 'c'), 2)],
            [(Compound('A', 'c'), 1), (Compound('E', 'c'), 1)])  # C[c] + 2 D[c] => A[c] + E[c]

    def test1_addBiorxn_default(self):
        g1 = vis.add_biomass_rxns(self.g, self.biomass_rxn, self.biomass_id)
        bio_A = graph.Node({
            'id': 'test_bio_2', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor':'#b3fcb8'})
        bio_C = graph.Node({
            'id': 'test_bio_1', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor': '#b3fcb8'})
        edge_bio_A = graph.Edge(bio_A, self.node_a, {'dir': 'forward'})
        edge_bio_C = graph.Edge(self.node_c, bio_C, {'dir': 'forward'})

        # for node in g1.nodes:
        #     print('vis', node.props)
        #
        # for node in [self.node_a, self.node_c, self.node_rxn1,
        #              bio_A, bio_C]:
        #     print('test',node.props)
        # for edge in g1.edges:
        #     print('vis', edge.source, edge.dest, edge.props)
        # for edge in [self.edge_A_rxn, self.edge_C_rxn, edge_bio_A, edge_bio_C]:
        #     print('test',edge.source, edge.dest, edge.props)

        self.assertTrue(all(i in [self.node_a, self.node_c, self.node_rxn1,
                                  bio_A, bio_C] for i in g1.nodes))

        self.assertTrue(all(i in [self.edge_A_rxn, self.edge_C_rxn,
                                  edge_bio_A, edge_bio_C] for i in g1.edges))


# def add_exchange_rxns(g, rxn_id, reaction):
class TestAddExchangeRxns(unittest.TestCase):
    def setUp(self):
        self.a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.a_extracell = graph.Node({
            'id': 'A[a]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'e'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.c_extracell = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'e'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        node_ac = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        node_cc = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',  # rxn2: C[c] <=> C[e]
            'original_id': ['rxn2'], 'compartment': 'e', 'fillcolor': '#c9fccd'})
        node_aa = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',  # rxn3: A[e] => A[c]
            'original_id': ['rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        edge_a_r1 = graph.Edge(self.a, node_ac, {'dir': 'forward'})
        # edge_r1_c = graph.Edge(self.c, node_ac, {'dir': 'forward'})   #why both ok?
        edge_r1_c = graph.Edge(node_ac, self.c, {'dir': 'forward'})
        edge_c_r2 = graph.Edge(self.c, node_cc, {'dir': 'both'})
        edge_r2_c = graph.Edge(node_cc, self.c_extracell, {'dir': 'both'})
        edge_a_r3 = graph.Edge(self.a_extracell, node_aa, {'dir': 'forward'})
        edge_r3_a = graph.Edge(node_aa, self.a, {'dir': 'forward'})

        self.g1 = graph.Graph()
        self.edge_list_c_to_e = [edge_a_r1, edge_r1_c, edge_c_r2, edge_r2_c]
        self.node_list_c_to_e = [self.a, self.c, self.c_extracell, node_ac, node_cc]
        for node in self.node_list_c_to_e:
            self.g1.add_node(node)
        for edge in self.edge_list_c_to_e:
            self.g1.add_edge(edge)
        self.rxn_C = Reaction(Direction.Both, {Compound('C', 'e'): -1})
        self.rxn_A = Reaction(Direction.Both, {Compound('A', 'e'): -1})

        self.g2 = graph.Graph()
        self.edge_list_e_to_c = [edge_a_r1, edge_r1_c, edge_a_r3, edge_r3_a]
        self.node_list_e_to_c = [self.a, self.a_extracell, self.c, node_aa, node_ac]
        for node in self.node_list_e_to_c:
            self.g2.add_node(node)
        for edge in self.edge_list_e_to_c:
            self.g2.add_edge(edge)

    def test1_addExrRxn_cpdC(self):
        g1 = vis.add_exchange_rxns(self.g1, 'test_Ex_C', self.rxn_C)
        node_Ex = graph.Node({
            'id': 'test_Ex_C', 'shape': 'box', 'style': 'filled', 'type': 'Ex_rxn',  # rxn3: A[e] => A[c]
            'original_id': ['test_Ex_C'], 'compartment': 'e', 'fillcolor': '#90f998'})
        edge_Ex = graph.Edge(self.c_extracell, node_Ex, {'dir': 'both'})
        self.node_list_c_to_e.append(node_Ex)
        self.edge_list_c_to_e.append(edge_Ex)
        self.assertTrue(all(i in self.node_list_c_to_e for i in g1.nodes))
        self.assertTrue(all(i in self.edge_list_c_to_e for i in g1.edges))

    def test2_addExrRxn_cpdA(self):
        g2 = vis.add_exchange_rxns(self.g2, 'test_Ex_A', self.rxn_C)

        node_Ex = graph.Node({
            'id': 'test_Ex_A', 'shape': 'box', 'style': 'filled', 'type': 'Ex_rxn',  # rxn3: A[e] => A[c]
            'original_id': ['test_Ex_A'], 'compartment': 'e', 'fillcolor': '#90f998'})
        edge_Ex = graph.Edge(self.a_extracell, node_Ex, {'dir': 'both'})
        self.node_list_e_to_c.append(node_Ex)
        self.edge_list_e_to_c.append(edge_Ex)
        self.assertTrue(all(i in self.node_list_e_to_c for i in g2.nodes))
        self.assertTrue(all(i in self.edge_list_e_to_c for i in g2.edges))


# def add_node_label(g, cpd_detail, rxn_detail, model_compound_entries,
#                    model_reaction_entries, reaction_flux):
class TestAddNodeLabel(unittest.TestCase):
    def setUp(self):
        self.a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.c_extracell = graph.Node({
            'id': 'C[e]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'e'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.node_ac = graph.Node({
            'id': 'rxn1_1,rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1', 'rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.node_cc = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',  # rxn2: C[c] <=> C[e]
            'original_id': ['rxn2'], 'compartment': 'e', 'fillcolor': '#c9fccd'})
        self.node_Ex = graph.Node({
            'id': 'test_Ex_C', 'shape': 'box', 'style': 'filled', 'type': 'Ex_rxn',  # rxn3: A[e] => A[c]
            'original_id': ['test_Ex_C'], 'compartment': 'e', 'fillcolor': '#90f998'})
        self.bio_A = graph.Node({
            'id': 'test_bio_1', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor': '#b3fcb8'})

        edge_a_r13 = graph.Edge(self.a, self.node_ac, {'dir': 'both'})    # in E_coli_core, it is both
        edge_r13_c = graph.Edge(self.node_ac, self.c, {'dir': 'both'})
        edge_c_r2 = graph.Edge(self.c, self.node_cc, {'dir': 'both'})
        edge_r2_c = graph.Edge(self.node_cc, self.c_extracell, {'dir': 'both'})
        edge_Ex = graph.Edge(self.c_extracell, self.node_Ex, {'dir': 'both'})
        edge_bio_A = graph.Edge(self.bio_A, self.a, {'dir': 'forward'})

        self.node_list = [self.a, self.c, self.c_extracell, self.node_ac, self.node_cc, self.node_Ex, self.bio_A]
        self.g = graph.Graph()
        for node in self.node_list:
            self.g.add_node(node)
        for edge in [edge_a_r13, edge_r13_c, edge_c_r2, edge_r2_c, edge_Ex, edge_bio_A]:
            self.g.add_edge(edge)

        self.cpd_detail = None
        self.rxn_detail = None
        # cpd_detail = [['id', 'name', 'charge']]
        # rxn_detail = [['id', 'equation', 'genes']]
        self.reaction_flux = {}

        cpd_A = entry.DictCompoundEntry({
            'id': 'A',
            'name': 'Compound A',
            'formula': 'C6H11O9P',
            'charge': 0})
        cpd_C = entry.DictCompoundEntry({
            'id': 'C',
            'name': 'Compound C',
            'formula': 'C6H11O9P',
            'charge': 0})
        self.cpd_entries = {'A': cpd_A, 'C': cpd_C}

        rxn1 = entry.DictReactionEntry({
            'id': 'rxn1',
            'name': 'Reaction 1',
            'equation': Reaction(Direction.Both, {Compound('A', 'c'): -1, Compound('C', 'c'): 1}),
            'genes': 'gene1 and gene2'
        }),
        rxn2 = entry.DictReactionEntry({
            'id': 'rxn2',
            'name': 'Reaction 2',
            'equation': Reaction(Direction.Both, {
                Compound('C', 'c'): -1, Compound('C', 'e'): 1}),
            'genes': 'gene3'})
        rxn3 = entry.DictReactionEntry({
            'id': 'rxn3',
            'name': 'Reaction 3',
            'equation': Reaction(Direction.Both, {Compound('A', 'c'): -1, Compound('C', 'c'): 1}),
        })
        test_Ex_C = entry.DictReactionEntry({
            'id': 'test_Ex_C',
            'name': 'Exchange Reaction',
            'equation': Reaction(Direction.Both, {Compound('C', 'e'): -1}),
            'genes': 'gene_ex'

        })
        test_bio = entry.DictReactionEntry({
            'id': 'test_bio',
            'name': 'biomass Reaction',
            'equation': Reaction(Direction.Forward, {Compound('A', 'c'): -1, Compound('D', 'c'): -1,
                                                     Compound('E', 'c'): 1,})
        }) # A[c] + D[c] => E[c]

        self.rxn_entries = {'rxn1': rxn1, 'rxn2': rxn2, 'rxn3': rxn3, 'test_Ex_C': test_Ex_C,
                            'test_bio': test_bio}

    def test_detail_default(self):
        g1 = vis.add_node_label(self.g, self.cpd_detail, self.rxn_detail, self.cpd_entries,
                                self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(i in [self.a, self.c, self.c_extracell, self.node_ac,
                                  self.node_cc, self.node_Ex, self.bio_A] for i in g1.nodes))

    def test_detail_fba(self):
        g2 = vis.add_node_label(self.g, self.cpd_detail, self.rxn_detail, self.cpd_entries,
                                self.rxn_entries, {'rxn1': 4.86, 'rxn2': 7.2,
                                                   'rxn3': 5.29, 'test_bio': 0.8})
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3\n10.15'
        self.node_cc.props['label'] = 'rxn2\n7.2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(i in [self.a, self.c, self.c_extracell, self.node_ac,
                                  self.node_cc, self.node_Ex, self.bio_A] for i in g2.nodes))\


    def test_detail_id_name(self): # for both reaction and compound
        g3 = vis.add_node_label(self.g, [['id','name']], [['id','name']], self.cpd_entries,
                                self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]\ncompound A'
        self.c.props['label'] = 'C[c]\ncompound C'
        self.c_extracell.props['label'] = 'C[e]\ncompound C'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2\nReaction 2'
        self.node_Ex.props['label'] = 'test_Ex_C\nExchange Reaction'
        self.bio_A.props['label'] = 'test_bio\nbiomass Reaction'
        self.assertTrue(all(i in [self.a, self.c, self.c_extracell, self.node_ac,
                                  self.node_cc, self.node_Ex, self.bio_A] for i in g3.nodes))

    def test_cpd_detail_id_formula_genes(self):
        g4 = vis.add_node_label(self.g, [['id','formula','genes']], self.rxn_detail, self.cpd_entries,
                                self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]\nC6H11O9P'
        self.c.props['label'] = 'C[c]\nC6H11O9P'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(i in [self.a, self.c, self.c_extracell, self.node_ac,
                                  self.node_cc, self.node_Ex, self.bio_A] for i in g4.nodes))

    def test_rxn_detail_id_formula_genes(self):
        g5 = vis.add_node_label(self.g, self.cpd_detail, [['id','formula','genes']], self.cpd_entries,
                                self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2\gene3'
        self.node_Ex.props['label'] = 'test_Ex_C\ngene_ex'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(i in [self.a, self.c, self.c_extracell, self.node_ac,
                                  self.node_cc, self.node_Ex, self.bio_A] for i in g5.nodes))


# def set_edge_props_withfba(g, edge_values):
class TestEdgePropsWithFBA(unittest.TestCase):
    def setUp(self):
        self.fum = graph.Node({
            'id': 'fum_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('fum_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.mal = graph.Node({
            'id': 'mal_L_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('mal_L_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.succ = graph.Node({
            'id': 'succ_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('succ_c', 'c'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.q8 = graph.Node({
            'id': 'q8_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('q8_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.q8h2 = graph.Node({
            'id': 'q8h2_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('q8h2_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.rxn_FUM = graph.Node({
            'id': 'FUM_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FUM'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.FRD7 = graph.Node({
            'id': 'FRD7_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FRD7'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.CYTBD_FRD7 = graph.Node({
            'id': 'CYTBD_1,FRD7_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['CYTBD', 'FRD7'], 'compartment': 'c', 'fillcolor': '#c9fccd'})  # NO LABEL

        # split graph nodes:
        self.CYTBD = graph.Node({
            'id': 'CYTBD_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['CYTBD'], 'compartment': 'c', 'fillcolor': '#c9fccd'})  # NO LABEL
        self.FRD7_2 = graph.Node({
            'id': 'FRD7_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FRD7'], 'compartment': 'c', 'fillcolor': '#c9fccd'})  # NO LABEL

        self.edge1 = graph.Edge(self.fum, self.rxn_FUM, {'dir': 'both'})  # in E_coli_core, it is both
        self.edge2 = graph.Edge(self.rxn_FUM, self.mal, {'dir': 'both'})
        self.edge3 = graph.Edge(self.fum, self.FRD7, {'dir': 'forward'})
        self.edge4 = graph.Edge(self.FRD7, self.succ, {'dir': 'forward'})
        self.edge5 = graph.Edge(self.q8, self.CYTBD_FRD7, {'dir': 'back'})
        self.edge6 = graph.Edge(self.CYTBD_FRD7, self.q8h2, {'dir': 'back'})

        self.edge7 = graph.Edge(self.q8, self.CYTBD, {'dir': 'back'})
        self.edge8 = graph.Edge(self.CYTBD, self.q8h2, {'dir': 'back'})
        self.edge9 = graph.Edge(self.q8, self.FRD7_2, {'dir': 'back'})
        self.edge10 = graph.Edge(self.FRD7_2, self.q8h2, {'dir': 'back'})

        self.edge_list = [self.edge1, self.edge2, self.edge3, self.edge4, self.edge5, self.edge6]
        self.edge_list_2 = [self.edge1, self.edge2, self.edge3, self.edge4,
                            self.edge7, self.edge8, self.edge9, self.edge10]

        self.g1 = graph.Graph()
        for node in [self.fum, self.mal, self.succ, self.rxn_FUM, self.FRD7,
                     self.CYTBD_FRD7, self.q8, self.q8h2]:
            self.g1.add_node(node)
        for i in self.edge_list:
            self.g1.add_edge(i)

        self.g2 = graph.Graph()
        for node in [self.fum, self.mal, self.succ, self.rxn_FUM, self.FRD7,
                     self.q8, self.q8h2, self.CYTBD, self.FRD7_2]:
            self.g2.add_node(node)
        for i in self.edge_list_2:
            self.g2.add_edge(i)

            # rxn1: 5.0, rxn2: 0, rxn3: 1.7, rxn4: 3.3
        self.edge_values = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
                            (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                            (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
                            (Compound('fum_c', 'c'), 'FUM'): 5.06
                            }

    def test1_nofba(self):
        g1 = vis.set_edge_props_withfba(self.g1, {})
        self.assertTrue(all(i in g1.edges for i in self.edge_list))
        self.assertTrue(all(i in self.edge_list for i in g1.edges))

    def test2_with_fba(self):
        g2 = vis.set_edge_props_withfba(self.g1, self.edge_values)
        edge1 = graph.Edge(self.fum, self.rxn_FUM, {'dir': 'both', 'penwidth': 1.161})  # in E_coli_core, it is bo
        edge2 = graph.Edge(self.rxn_FUM, self.mal, {'dir': 'both', 'penwidth': 1.161})
        edge3 = graph.Edge(self.fum, self.FRD7, {'dir': 'forward', 'style': 'dotted'})
        edge4 = graph.Edge(self.FRD7, self.succ, {'dir': 'forward', 'style': 'dotted'})
        edge5 = graph.Edge(self.q8, self.CYTBD_FRD7, {'dir': 'back', 'penwidth': 10})
        edge6 = graph.Edge(self.CYTBD_FRD7, self.q8h2, {'dir': 'back', 'penwidth': 10})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6]
        self.assertTrue(all(i in edge_list for i in g2.edges))
        self.assertTrue(all(i in g2.edges for i in edge_list))

    def test3_with_fba_splitGraph(self):
        edge_values = {
            (Compound('q8_c', 'c'), 'CYTBD'): 43.60,
            (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
            (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
            (Compound('fum_c', 'c'), 'FUM'): 5.06}
        g3 = vis.set_edge_props_withfba(self.g2, edge_values)
        edge1 = graph.Edge(self.fum, self.rxn_FUM, {'dir': 'both', 'penwidth': 1.161})
        edge2 = graph.Edge(self.rxn_FUM, self.mal, {'dir': 'both', 'penwidth': 1.161})
        edge3 = graph.Edge(self.fum, self.FRD7, {'dir': 'forward', 'style': 'dotted'})
        edge4 = graph.Edge(self.FRD7, self.succ, {'dir': 'forward', 'style': 'dotted'})
        edge5 = graph.Edge(self.q8, self.CYTBD, {'dir': 'back', 'penwidth': 10})
        edge6 = graph.Edge(self.CYTBD, self.q8h2, {'dir': 'back', 'penwidth': 10})
        edge7 = graph.Edge(self.q8, self.FRD7_2, {'dir': 'back', 'style': 'dotted'})
        edge8 = graph.Edge(self.FRD7_2, self.q8h2, {'dir': 'back', 'style': 'dotted'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8]

        self.assertTrue(all(i in edge_list for i in g3.edges))
        self.assertTrue(all(i in g3.edges for i in edge_list))


# def make_cpt_tree(boundaries, extracellular):
class TestMakeCptTree(unittest.TestCase):
    def test1_if_e_inthemodel(self):
        boundaries = [('c', 'e'), ('c', 'p'), ('e', 'p')]
        extra = 'e'
        c1, e1 = vis.make_cpt_tree(boundaries, extra)
        c1_res = defaultdict(set)
        c1_res['c'].add('e')
        c1_res['c'].add('p')
        c1_res['p'].add('e')
        c1_res['p'].add('c')
        c1_res['e'].add('c')
        c1_res['e'].add('p')
        # e1_res = 'e'
        self.assertEqual(c1, c1_res)
        # self.assertEqual(e1, e1_res)

    def test2_if_e_not_inthemodel(self):
        boundaries = [('c', 'p')]
        extra = 'e'
        c2, e2 = vis.make_cpt_tree(boundaries, extra)
        c2_res = defaultdict(set)
        c2_res['c'].add('p')
        c2_res['p'].add('c')
        # e2_res = 'c'
        self.assertEqual(c2, c2_res)
        # self.assertEqual(e2, e2_res)

    def test3_if_e_not_inthemodel_cpmi(self):
        boundaries = [('c', 'p'), ('c', 'mi')]
        extra = 'e'
        c3, e3 = vis.make_cpt_tree(boundaries, extra)
        c3_res = defaultdict(set)
        c3_res['c'].add('p')
        c3_res['c'].add('mi')
        c3_res['mi'].add('c')
        c3_res['p'].add('c')
        # e3_res = 'mi'   failed
        self.assertEqual(c3, c3_res)
        # self.assertEqual(e3, e3_res)

# def make_filter_dict(model, mm, method, element, cpd_formula,
# arg_hide_edges, exclude_rxns):    SUCDi
class TestFilterDict(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({'id': 'rxn1', 'equation': parse_reaction('fum_c[c] + h2o_c[c] <=> mal_L_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'fum_c[c]', 'formula': parse_compound('C4H2O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'h2o_c[c]', 'formula': parse_compound('H2O', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'mal_L_c[c]', 'formula': parse_compound('C4H4O5', 'c')}))
        native_model.reactions.add_entry(ReactionEntry({'id': 'rxn2', 'equation': parse_reaction('q8_c[c] + succ_c[c] => fum_c[c] + q8h2_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'q8_c[c]', 'formula': parse_compound('C49H74O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'q8h2_c[c]', 'formula': parse_compound('C49H76O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'succ_c[c]', 'formula': parse_compound('C4H4O4', 'c')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

        self.cpd_formula = {'fum_c': Formula.parse('C4H2O4'), 'h2o_c': Formula.parse('H2O'),
                            'mal_L_c': Formula.parse('C4H4O5'), 'q8_c': Formula.parse('C49H74O4'),
                            'q8h2_c': Formula.parse('C49H76O4'), 'succ_c': Formula.parse('C4H4O4')}
        self.element = 'C'
        self.method = 'fpp'
        self.exclude_rxns = []
        self.hide_edges = None

    def test1_default_setting(self):
        e1 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e1_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
                  'rxn2': [(Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c'))]}
        print('e1-result', e1)
        self.assertEqual(e1, e1_res)

    def test2_nofpp(self):
        e2 = vis.make_filter_dict(
            self.native, self.mm, 'no-fpp', self.element, self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e2_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
                  'rxn2': [(Compound('q8_c', 'c'), Compound('fum_c', 'c')),
                           (Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e2, e2_res)

    def test3_element_hydrogen(self):
        e3 = vis.make_filter_dict(
            self.native, self.mm, self.method, 'H', self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e3_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c')),
                           (Compound('h2o_c', 'c'), Compound('mal_L_c', 'c'))],
                  'rxn2': [(Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e3, e3_res)

    def test4_exclude_rxn(self):
        e4 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, self.cpd_formula,
            self.hide_edges, ['rxn2'])
        e4_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))]}
        self.assertEqual(e4, e4_res)

    def test5_if_missing_formula_fpp(self):
        cpd_formula = {'fum_c': Formula.parse('C4H2O4'), 'h2o_c': Formula.parse('H2O'),
                       'mal_L_c': Formula.parse('C4H4O5'), 'q8_c': Formula.parse('C49H74O4'),
                       'q8h2_c': Formula.parse('C49H76O4')}
        e5 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e5_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))]}
        self.assertEqual(e5, e5_res)

    def test6_if_missing_formula_nofpp(self):
        cpd_formula = {'fum_c': Formula.parse('C4H2O4'), 'h2o_c': Formula.parse('H2O'),
                       'mal_L_c': Formula.parse('C4H4O5'), 'q8_c': Formula.parse('C49H74O4'),
                       'q8h2_c': Formula.parse('C49H76O4')}
        e6 = vis.make_filter_dict(
            self.native, self.mm, 'no-fpp', self.element, cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e6_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
                  'rxn2': [(Compound('q8_c', 'c'), Compound('fum_c', 'c')),
                           (Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e6, e6_res)

    # def test7_hide_edges(self):
    #     f = open("remove_edges.tsv", "w+")
    #     f.write('q8_c[c]'   'q8h2_c[c]')
    #     f.close()
    #     hide_edges = f
    #
    #     e4 = vis.make_filter_dict(
    #         self.native, self.mm, self.method, self.element, self.cpd_formula,
    #         hide_edges, self.exclude_rxns)
    #     e4_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))]}
    #     self.assertEqual(e4, e4_res)


class TestGetCptBoundaries(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({'id': 'rxn1', 'equation': parse_reaction('A[c] => A[e]')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'A[c]', 'formula': parse_compound('formula_A', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'A[e]', 'formula': parse_compound('formula_A', 'e')}))
        native_model.reactions.add_entry(ReactionEntry({'id': 'rxn2', 'equation': parse_reaction('B[e] <=> B[p]')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'B[e]', 'formula': parse_compound('formula_B', 'e')}))
        native_model.compounds.add_entry(CompoundEntry({'id': 'B[p]', 'formula': parse_compound('formula_B', 'p')}))
        self.native = native_model

    # def test_default_setting(self):
    #     e1_bound,e1_extra = vis.get_cpt_boundaries(self.native)
    #     e1_bound_res = set()
    #     e1_bound_res.add(('c', 'e'))
    #     e1_bound_res.add(('p', 'e'))
    #     e1_extra_res = 'e'
    #     self.assertTrue(all(i in e1_bound for i in e1_bound_res))
    #     self.assertEqual(e1_extra, e1_extra_res)

    # def test2_with_extra_defined(self):
    #     self.native._properties['extracellular'] = 'p'
    #     e2_bound,e2_extra = vis.get_cpt_boundaries(self.native)
    #     e2_bound_res = set()
    #     e2_bound_res.add(('c', 'e'))
    #     e2_bound_res.add(('p', 'e'))
    #     e2_extra_res = 'p'
    #     self.assertTrue(all(i in e2_bound for i in e2_bound_res))
    #     self.assertEqual(e2_extra, e2_extra_res)


if __name__ == '__main__':
    unittest.main()