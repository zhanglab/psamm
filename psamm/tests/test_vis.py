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
# Copyright 2018-2020  Ke Zhang <kzhang@my.uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth@my.uri.edu>

from __future__ import unicode_literals

import unittest
import os
import tempfile

from psamm.formula import Formula, Atom
from psamm.commands import vis
from psamm.reaction import Compound, Reaction, Direction
from collections import defaultdict
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm.datasource import entry
from psamm import graph
from psamm.datasource.native import NativeModel, ReactionEntry, CompoundEntry


# 2019-08-30 new test cases for vis
class TestMakeSubset(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum[c] + h2o[c] <=> mal_L[c]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'fum', 'formula': parse_compound('C4H2O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o', 'formula': parse_compound('H2O', 'c')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'mal_L', 'formula': parse_compound('C4H4O5', 'c')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'q8[c] + succ[c] => fum[c] + q8h2[c]')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'q8', 'formula': parse_compound('C49H74O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'q8h2', 'formula': parse_compound('C49H76O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'succ', 'formula': parse_compound('C4H4O4', 'c')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn3', 'equation': parse_reaction(
                'h2o[c] <=> h2o[e]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o', 'formula': parse_compound('H2O', 'c')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()
        self.sub = None
        self.exclude = []

    def test1_no_subset(self):
        subset = vis.rxnset_for_vis(self.mm, self.sub, self.exclude)
        self.assertEqual(subset, set(['rxn1', 'rxn2', 'rxn3']))

    def test2_subset_reactions(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('rxn1', 'rxn2'))
        sub_file = open(path, 'r')
        subset = vis.rxnset_for_vis(self.mm, sub_file, self.exclude)
        sub_file.close()
        self.assertEqual(subset, set(['rxn1', 'rxn2']))

    def test3_subset_one_compound1(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}'.format('mal_L[c]'))
        sub_file = open(path, 'r')
        subset = vis.rxnset_for_vis(self.mm, sub_file, self.exclude)
        sub_file.close()
        self.assertEqual(subset, set(['rxn1']))

    def test4_subset_one_compound2(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}'.format('fum[c]'))
        sub_file = open(path, 'r')
        subset = vis.rxnset_for_vis(self.mm, sub_file, self.exclude)
        sub_file.close()
        self.assertEqual(subset, set(['rxn1', 'rxn2']))

    def test5_subset_compound3(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}'.format('h2o[c]'))
        sub_file = open(path, 'r')
        subset = vis.rxnset_for_vis(self.mm, sub_file, self.exclude)
        sub_file.close()
        self.assertEqual(subset, set(['rxn1', 'rxn3']))

    def test6_mixed_subset(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('h2o[c]', 'rxn1'))
        sub_file = open(path, 'r')
        self.assertRaises(ValueError, vis.rxnset_for_vis,
                          self.mm, sub_file, self.exclude)

    def test7_not_in_model(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}'.format('h3o'))
        sub_file = open(path, 'r')
        self.assertRaises(ValueError, vis.rxnset_for_vis,
                          self.mm, sub_file, self.exclude)

    def test8_exclude(self):
        exclude = ['rxn3']
        subset = vis.rxnset_for_vis(self.mm, self.sub, exclude)
        self.assertEqual(subset, set(['rxn1', 'rxn2']))

    def test9_exclude_mult(self):
        exclude = ['rxn1', 'rxn3']
        subset = vis.rxnset_for_vis(self.mm, self.sub, exclude)
        self.assertEqual(subset, set(['rxn2']))


class TestAddNodeProps(unittest.TestCase):
    def setUp(self):
        self.cpd_a = CompoundEntry({
            'id': 'A', 'name': 'compound A', 'formula': 'XXXX', 'charge': 0})
        self.cpd_b = CompoundEntry({
            'id': 'B', 'name': 'compound B', 'formula': 'YYYY', 'charge': 0})
        self.cpd_c = CompoundEntry({
            'id': 'C', 'name': 'compound C', 'formula': 'ZZZZ', 'charge': 0})
        self.cpd_d = CompoundEntry({
            'id': 'D', 'name': 'compound D', 'formula': 'MMMM', 'charge': 0})
        self.cpd_e = CompoundEntry({
            'id': 'E', 'name': 'compound E', 'formula': 'NNNN', 'charge': 0})

        # new rxn used to test: rxn1: A => B; rxn2: A => B; rxn3: C <=> B + D
        self.rxn1 = entry.DictReactionEntry({
            'id': 'rxn1', 'name': 'Reaction 1', 'genes': 'gene1 and gene2',
            'equation': Reaction(
                Direction.Forward, {Compound('A', 'c'): -1,
                                    Compound('B', 'c'): 1}),
            })
        self.rxn2 = entry.DictReactionEntry({
            'id': 'rxn2', 'name': 'Reaction 2', 'genes': 'gene3',
            'equation': Reaction(Direction.Forward, {
                Compound('A', 'c'): -1, Compound('B', 'c'): 1})
            })
        self.rxn3 = entry.DictReactionEntry({
            'id': 'rxn3', 'name': 'Reaction 3',
            'equation': Reaction(Direction.Both, {
                Compound('C', 'c'): -1, Compound('B', 'c'): 1,
                Compound('D', 'c'): 1})
        })

        self.node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c'})
        self.node_b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c'})
        self.node_c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c'})
        self.node_d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c'})
        self.node_ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn'})
        self.node_cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3], 'compartment': 'c',
            'type': 'rxn'})

        self.edge1 = graph.Edge(self.node_a, self.node_ab, {'dir': 'forward'})
        self.edge2 = graph.Edge(self.node_ab, self.node_b, {'dir': 'forward'})
        self.edge3 = graph.Edge(self.node_c, self.node_cb, {'dir': 'both'})
        self.edge4 = graph.Edge(self.node_cb, self.node_b, {'dir': 'both'})
        self.edge5 = graph.Edge(self.node_cb, self.node_d, {'dir': 'both'})

        self.g = graph.Graph()
        self.node_list = [self.node_a, self.node_b, self.node_c, self.node_d,
                          self.node_ab, self.node_cb]
        for node in self.node_list:
            self.g.add_node(node)
        self.edge_list = [self.edge1, self.edge2,
                          self.edge3, self.edge4, self.edge5]
        for edge in self.edge_list:
            self.g.add_edge(edge)

        self.recolor_dict = {'A': '#f4fc55', 'rxn1': '#ef70de',
                             'rxn3': '#64b3f4'}

    def test1_DefaultRecolor(self):
        g1 = vis.add_node_props(self.g, {})
        node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd'})
        node_cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd'})
        node_list = [node_a, node_b, node_c, node_d, node_ab, node_cb]
        self.assertTrue(all(i in node_list for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in node_list))

    def test2_Recolor(self):
        g2 = vis.add_node_props(self.g, self.recolor_dict)
        node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#f4fc55'})
        node_b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd'})
        node_cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#64b3f4'})
        node_list = [node_a, node_b, node_c, node_d, node_ab, node_cb]
        self.assertTrue(all(i in node_list for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in node_list))

    def test3_ConflictRecolorDict(self):
        recolor_dict = {'A': '#f4fc55', 'rxn1': '#ef70de', 'rxn2': '#64b3f4'}
        g3 = vis.add_node_props(self.g, recolor_dict)
        node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#f4fc55'})
        node_b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        node_ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd'})
        node_cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd'})
        node_list = [node_a, node_b, node_c, node_d, node_ab, node_cb]
        self.assertTrue(all(i in node_list for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in node_list))


class TestAddNodeLabel(unittest.TestCase):
    def setUp(self):
        self.cpd_a = CompoundEntry({
            'id': 'A', 'name': 'compound A', 'formula': 'XXXX', 'charge': 0})
        self.cpd_b = CompoundEntry({
            'id': 'B', 'name': 'compound B', 'formula': 'YYYY', 'charge': 0})
        self.cpd_c = CompoundEntry({
            'id': 'C', 'name': 'compound C', 'formula': 'ZZZZ', 'charge': 0})
        self.cpd_d = CompoundEntry({
            'id': 'D', 'name': 'compound D', 'formula': 'MMMM', 'charge': 0})
        # self.cpd_e = CompoundEntry({
        #     'id': 'E', 'name': 'compound E', 'formula': 'NNNN', 'charge': 0})

        # new rxn used to test: rxn1: A => B; rxn2: A => B; rxn3: C <=> B + D
        self.rxn1 = entry.DictReactionEntry({
            'id': 'rxn1', 'name': 'Reaction 1', 'genes': 'gene1 and gene2',
            'equation': Reaction(
                Direction.Forward, {Compound('A', 'c'): -1,
                                    Compound('B', 'c'): 1}),
        })
        self.rxn2 = entry.DictReactionEntry({
            'id': 'rxn2', 'name': 'Reaction 2',
            'equation': Reaction(Direction.Forward, {
                Compound('A', 'c'): -1, Compound('B', 'c'): 1})
        })
        self.rxn3 = entry.DictReactionEntry({
            'id': 'rxn3', 'name': 'Reaction 3', 'genes': 'gene3',
            'equation': Reaction(Direction.Both, {
                Compound('C', 'c'): -1, Compound('B', 'c'): 1,
                Compound('D', 'c'): 1})
        })

        self.node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        self.node_b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        self.node_c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        self.node_d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf'})
        self.node_ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd'})
        self.node_cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd'})
        self.edge1 = graph.Edge(self.node_a, self.node_ab, {'dir': 'forward'})
        self.edge2 = graph.Edge(self.node_ab, self.node_b, {'dir': 'forward'})
        self.edge3 = graph.Edge(self.node_c, self.node_cb, {'dir': 'both'})
        self.edge4 = graph.Edge(self.node_cb, self.node_b, {'dir': 'both'})
        self.edge5 = graph.Edge(self.node_cb, self.node_d, {'dir': 'both'})

        self.g = graph.Graph()
        node_list = [self.node_a, self.node_b,
                     self.node_c, self.node_d, self.node_ab, self.node_cb]
        # node_list = [a, b, c, d, ab, cb]
        for node in node_list:
            self.g.add_node(node)
        edge_list = [self.edge1, self.edge2,
                     self.edge3, self.edge4, self.edge5]
        for edge in edge_list:
            self.g.add_edge(edge)

    def test1_detail_default(self):
        """test vis codes when both --cpd-detail and --rxn-detail are None"""
        g1 = vis.add_node_label(self.g, None, None)
        # self.node_a.props['label'] = 'A[c]'
        # self.node_b.props['label'] = 'B[c]'
        # self.node_c.props['label'] = 'C[c]'
        # self.node_d.props['label'] = 'D[c]'
        # self.node_ab.props['label'] = 'rxn1\nrxn2'
        # self.node_cb.props['label'] = 'rxn3\nReaction 3'
        # node_list_1 = [self.node_a, self.node_b, self.node_c, self.node_d,
        #                self.node_ab, self.node_cb]

        a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'A[c]'})
        b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'B[c]'})
        c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'C[c]'})
        d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'D[c]'})
        ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd', 'label': 'rxn1\nrxn2'})
        cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd', 'label': 'rxn3'})
        node_list_1 = [a, b, c, d, ab, cb]
        # for i in g1.nodes:
        #     print('in g1:', i.props)
        #     if i not in node_list_1:
        #         print('in g1 but not in node_list:', i)
        # for j in node_list_1:
        #     print('in node list:', j.props)
        self.assertTrue(all(i in node_list_1 for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in node_list_1))

    def test2_detail_id_name(self):
        """test vis codes when both --cpd-detail
        and --rxn-detail are ['id', 'name']"""
        g2 = vis.add_node_label(self.g, [['id', 'name']], [['id', 'name']])
        # self.node_a.props['label'] = 'A[c]\ncompound A'
        # self.node_b.props['label'] = 'B[c]\ncompound B'
        # self.node_c.props['label'] = 'C[c]\ncompound C'
        # self.node_d.props['label'] = 'D[c]\ncompound D'
        # self.node_ab.props['label'] = 'rxn1\nrxn2'
        # self.node_cb.props['label'] = 'rxn3\nReaction 3'
        # node_list_2 = [self.node_a, self.node_b, self.node_c, self.node_d,
        #                self.node_ab, self.node_cb]

        a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'A[c]\ncompound A'})
        b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'B[c]\ncompound B'})
        c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'C[c]\ncompound C'})
        d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'D[c]\ncompound D'})
        ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd', 'label': 'rxn1\nrxn2'})
        cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd', 'label': 'rxn3\nReaction 3'})
        node_list_2 = [a, b, c, d, ab, cb]
        for i in g2.nodes:
            if i not in node_list_2:
                print('not in node_list :', i, i.props)
        self.assertTrue(all(i in node_list_2 for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in node_list_2))

    def test3_id_formula_genes(self):
        """test vis codes when --cpd-detail and
        --rxn-detail are ['id', 'formula', 'genes']"""
        g3 = vis.add_node_label(self.g, [['id', 'formula', 'genes']],
                                [['id', 'formula', 'genes']])
        # self.node_a.props['label'] = 'A[c]\ncompound A'
        # self.node_b.props['label'] = 'B[c]\ncompound B'
        # self.node_c.props['label'] = 'C[c]\ncompound C'
        # self.node_d.props['label'] = 'D[c]\ncompound D'
        # self.node_ab.props['label'] = 'rxn1\nrxn2'
        # self.node_cb.props['label'] = 'rxn3\nReaction 3'
        # node_list_2 = [self.node_a, self.node_b, self.node_c, self.node_d,
        #                self.node_ab, self.node_cb]

        a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'A[c]\nXXXX'})
        b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'B[c]\nYYYY'})
        c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'C[c]\nZZZZ'})
        d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'D[c]\nMMMM'})
        ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd', 'label': 'rxn1\nrxn2'})
        cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd', 'label': 'rxn3\ngene3'})
        node_list_3 = [a, b, c, d, ab, cb]
        self.assertTrue(all(i in node_list_3 for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in node_list_3))

    def test4_noID(self):
        """test vis codes when 'id' is not in --cpd-detail and --rxn-detail"""
        g4 = vis.add_node_label(self.g, [['name', 'genes']],
                                [['name', 'formula', 'genes']])
        # self.node_a.props['label'] = 'A[c]\ncompound A'
        # self.node_b.props['label'] = 'B[c]\ncompound B'
        # self.node_c.props['label'] = 'C[c]\ncompound C'
        # self.node_d.props['label'] = 'D[c]\ncompound D'
        # self.node_ab.props['label'] = 'rxn1\nrxn2'
        # self.node_cb.props['label'] = 'rxn3\nReaction 3'
        # node_list_2 = [self.node_a, self.node_b, self.node_c, self.node_d,
        #                self.node_ab, self.node_cb]

        a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound A'})
        b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound B'})
        c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound C'})
        d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound D'})
        ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd', 'label': 'rxn1\nrxn2'})
        cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd', 'label': 'Reaction 3\ngene3'})
        node_list_4 = [a, b, c, d, ab, cb]
        self.assertTrue(all(i in node_list_4 for i in g4.nodes))
        self.assertTrue(all(i in g4.nodes for i in node_list_4))

    def test5_(self):
        """test vis codes when obly --cpd-detail is given """
        g5 = vis.add_node_label(self.g, [['name', 'genes']], None)
        # self.node_a.props['label'] = 'A[c]\ncompound A'
        # self.node_b.props['label'] = 'B[c]\ncompound B'
        # self.node_c.props['label'] = 'C[c]\ncompound C'
        # self.node_d.props['label'] = 'D[c]\ncompound D'
        # self.node_ab.props['label'] = 'rxn1\nrxn2'
        # self.node_cb.props['label'] = 'rxn3\nReaction 3'
        # node_list_2 = [self.node_a, self.node_b, self.node_c, self.node_d,
        #                self.node_ab, self.node_cb]

        a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound A'})
        b = graph.Node({
            'id': 'B[c]', 'entry': [self.cpd_b],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound B'})
        c = graph.Node({
            'id': 'C[c]', 'entry': [self.cpd_c],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound C'})
        d = graph.Node({
            'id': 'D[c]', 'entry': [self.cpd_d],
            'type': 'cpd', 'compartment': 'c',
            'style': 'filled', 'shape': 'ellipse',
            'fillcolor': '#ffd8bf', 'label': 'compound D'})
        ab = graph.Node({
            'id': 'rxn1_1, rxn2_1', 'entry': [self.rxn1, self.rxn2],
            'compartment': 'c', 'type': 'rxn', 'style': 'filled',
            'shape': 'box', 'fillcolor': '#c9fccd', 'label': 'rxn1\nrxn2'})
        cb = graph.Node({
            'id': 'rxn3_1', 'entry': [self.rxn3],
            'compartment': 'c', 'type': 'rxn',
            'style': 'filled', 'shape': 'box',
            'fillcolor': '#c9fccd', 'label': 'rxn3'})
        node_list_5 = [a, b, c, d, ab, cb]
        self.assertTrue(all(i in node_list_5 for i in g5.nodes))
        self.assertTrue(all(i in g5.nodes for i in node_list_5))


class TestAddBiomassRnxs(unittest.TestCase):
    def setUp(self):
        cpd_a = CompoundEntry({
            'id': 'A', 'name': 'compound A', 'formula': 'XXXX', 'charge': 0})
        cpd_b = CompoundEntry({
            'id': 'B', 'name': 'compound B', 'formula': 'YYYY', 'charge': 0})
        cpd_c = CompoundEntry({
            'id': 'C', 'name': 'compound C', 'formula': 'ZZZZ', 'charge': 0})
        cpd_d = CompoundEntry({
            'id': 'D', 'name': 'compound D', 'formula': 'MMMM', 'charge': 0})
        rxn1 = entry.DictReactionEntry({
            'id': 'rxn1', 'name': 'Reaction 1', 'genes': 'gene1 and gene2',
            'equation': Reaction(
                Direction.Both, {Compound('A', 'c'): -1,
                                 Compound('C', 'c'): 1}),
        })

        self.node_a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'entry': [cpd_a], 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_b = graph.Node({
            'id': 'B[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'entry': [cpd_b], 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'entry': [cpd_c], 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_d = graph.Node({
            'id': 'D[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'entry': [cpd_d], 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'entry': [rxn1], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.edge_A_rxn = graph.Edge(self.node_a, self.node_rxn1,
                                     {'dir': 'both'})
        self.edge_C_rxn = graph.Edge(self.node_rxn1, self.node_c,
                                     {'dir': 'both'})
        self.g = graph.Graph()
        for node in [self.node_a, self.node_b,
                     self.node_c, self.node_d, self.node_rxn1]:
            self.g.add_node(node)
        self.g.add_edge(self.edge_A_rxn)
        self.g.add_edge(self.edge_C_rxn)

        self.biorxn = ReactionEntry({
            'id': 'test_bio', 'name': 'biomass_equatin', 'equation':
                parse_reaction('A[c] + B[c] => D[c]')})

    def test1_addBiorxn_default(self):
        g1 = vis.add_biomass_rxns(self.g, self.biorxn)
        bio_a = graph.Node({
            'id': 'test_bio_1', 'entry': [self.biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_b = graph.Node({
            'id': 'test_bio_2', 'entry': [self.biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_d = graph.Node({
            'id': 'test_bio_3', 'entry': [self.biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        edge_bio_a = graph.Edge(self.node_a, bio_a, {'dir': 'forward'})
        edge_bio_b = graph.Edge(self.node_b, bio_b, {'dir': 'forward'})
        edge_bio_d = graph.Edge(bio_d, self.node_d, {'dir': 'forward'})

        self.assertTrue(all(i in [self.node_a, self.node_b,
                                  self.node_c, self.node_d,
                                  self.node_rxn1, bio_a, bio_b, bio_d]
                            for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in [
            self.node_a, self.node_b, self.node_c, self.node_d, self.node_rxn1,
            bio_a, bio_b, bio_d]))

        self.assertTrue(all(i in [self.edge_A_rxn, self.edge_C_rxn, edge_bio_a,
                                  edge_bio_b, edge_bio_d] for i in g1.edges))
        self.assertTrue(all(i in g1.edges for i in [
            self.edge_A_rxn, self.edge_C_rxn,
            edge_bio_a, edge_bio_b, edge_bio_d]))

    def test2_TwoCpdsRight(self):
        biorxn = ReactionEntry({
            'id': 'test_bio', 'name': 'biomass_equatin',
            'equation': parse_reaction('A[c] + B[c] => C[c] + D[c]')})
        g2 = vis.add_biomass_rxns(self.g, biorxn)
        bio_a = graph.Node({
            'id': 'test_bio_1', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_b = graph.Node({
            'id': 'test_bio_2', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_c = graph.Node({
            'id': 'test_bio_3', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_d = graph.Node({
            'id': 'test_bio_4', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        edge_bio_a = graph.Edge(self.node_a, bio_a, {'dir': 'forward'})
        edge_bio_b = graph.Edge(self.node_b, bio_b, {'dir': 'forward'})
        edge_bio_c = graph.Edge(bio_c, self.node_c, {'dir': 'forward'})
        edge_bio_d = graph.Edge(bio_d, self.node_d, {'dir': 'forward'})

        final_nodes = [self.node_a, self.node_b, self.node_c, self.node_d,
                       self.node_rxn1, bio_a, bio_b, bio_c, bio_d]
        final_edges = [self.edge_A_rxn, self.edge_C_rxn, edge_bio_a,
                       edge_bio_b, edge_bio_c, edge_bio_d]
        self.assertTrue(all(i in final_nodes for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in final_nodes))
        self.assertTrue(all(i in final_edges for i in g2.edges))
        self.assertTrue(all(i in g2.edges for i in final_edges))

    def test3_EmptyRight(self):
        biorxn = ReactionEntry({
            'id': 'test_bio', 'name': 'biomass_equatin',
            'equation': parse_reaction('A[c] + B[c] => ')})
        g3 = vis.add_biomass_rxns(self.g, biorxn)
        bio_a = graph.Node({
            'id': 'test_bio_1', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        bio_b = graph.Node({
            'id': 'test_bio_2', 'entry': [biorxn], 'shape': 'box',
            'style': 'filled', 'label': 'test_bio', 'type': 'bio_rxn',
            'fillcolor': '#b3fcb8', 'compartment': 'c'})
        edge_bio_a = graph.Edge(self.node_a, bio_a, {'dir': 'forward'})
        edge_bio_b = graph.Edge(self.node_b, bio_b, {'dir': 'forward'})
        final_nodes = [self.node_a, self.node_b, self.node_c, self.node_d,
                       self.node_rxn1, bio_a, bio_b]
        final_edges = [self.edge_A_rxn, self.edge_C_rxn,
                       edge_bio_a, edge_bio_b]
        self.assertTrue(all(i in final_nodes for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in final_nodes))
        self.assertTrue(all(i in final_edges for i in g3.edges))
        self.assertTrue(all(i in g3.edges for i in final_edges))


class TestAddExchange(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum[c] + h2o[c] <=> mal_L[c]')}))
        self.cpd_a = CompoundEntry({
            'id': 'A', 'name': 'compound A', 'formula': 'CO2', 'charge': 0})
        cpd_c = CompoundEntry({
            'id': 'C', 'name': 'compound C', 'formula': 'MNO', 'charge': 0})
        rxn1 = entry.DictReactionEntry({
            'id': 'rxn1', 'name': 'Reaction 1', 'equation': Reaction(
                Direction.Forward, {Compound('A', 'c'): -1,
                                    Compound('C', 'c'): 1})
        })
        rxn2 = entry.DictReactionEntry({
            'id': 'rxn2', 'name': 'Reaction 2', 'equation': Reaction(
                Direction.Forward, {Compound('A', 'e'): -1,
                                    Compound('A', 'c'): 1}),
        })
        rxn3 = entry.DictReactionEntry({
            'id': 'rxn3', 'name': 'Reaction 3', 'equation': Reaction(
                Direction.Forward, {Compound('C', 'c'): 1,
                                    Compound('C', 'e'): 1}),
        })

        node_a = graph.Node({
            'id': 'A[c]', 'entry': [self.cpd_a],
            'compartment': 'c', 'type': 'cpd',
            'style': 'filled', 'shape': 'ellipse', 'fillcolor': '#ffd8bf',
            'label': 'A[c]'})
        node_c = graph.Node({
            'id': 'C[c]', 'entry': [cpd_c], 'compartment': 'c', 'type': 'cpd',
            'style': 'filled', 'shape': 'ellipse', 'fillcolor': '#ffd8bf',
            'label': 'C[c]'})
        self.node_a_extracell = graph.Node({
            'id': 'A[e]', 'entry': [self.cpd_a],
            'compartment': 'e', 'type': 'cpd',
            'style': 'filled', 'shape': 'ellipse', 'fillcolor': '#ffd8bf',
            'label': 'A[e]'})
        self.node_c_extracell = graph.Node({
            'id': 'C[e]', 'entry': [cpd_c], 'compartment': 'e', 'type': 'cpd',
            'style': 'filled', 'shape': 'ellipse', 'fillcolor': '#ffd8bf',
            'label': 'C[e]'})
        node_ac = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'entry': [rxn1], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn1'})
        self.node_aa = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'entry': [rxn2], 'compartment': 'c', 'fillcolor': '#90f998',
            'label': 'rxn2'})
        node_cc = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'entry': [rxn3], 'compartment': 'e', 'fillcolor': '#90f998',
            'label': 'rxn3'})
        edge_a_r1 = graph.Edge(node_a, node_ac, {'dir': 'both'})
        edge_r1_c = graph.Edge(node_ac, node_c, {'dir': 'both'})
        self.edge_a_r2 = graph.Edge(self.node_a_extracell,
                                    self.node_aa, {'dir': 'forward'})
        self.edge_r2_a = graph.Edge(self.node_aa, node_a, {'dir': 'forward'})
        edge_c_r3 = graph.Edge(self.node_c_extracell, node_cc, {'dir': 'back'})
        edge_r3_c = graph.Edge(node_cc, node_c, {'dir': 'back'})

        self.g = graph.Graph()
        self.edge_list = [edge_a_r1, edge_r1_c, self.edge_a_r2, self.edge_r2_a,
                          edge_c_r3, edge_r3_c]
        self.node_list = [node_a, node_c, self.node_a_extracell,
                          self.node_c_extracell, node_ac,
                          self.node_aa, node_cc]
        for node in self.node_list:
            self.g.add_node(node)
        for edge in self.edge_list:
            self.g.add_edge(edge)
        self.EX_A = Reaction(Direction.Both, {Compound('A', 'e'): -5})
        self.EX_C = Reaction(Direction.Both, {Compound('C', 'e'): 5})

    def test1_addExrRxn_cpdA(self):
        g1 = vis.add_exchange_rxns(self.g, 'test_EX_A', self.EX_A,
                                   {'test_EX_A': ('solid', 1)})
        node_ex = graph.Node({
            'id': 'test_EX_A', 'entry': [self.EX_A], 'shape': 'box',
            'style': 'filled', 'label': 'test_EX_A', 'type': 'Ex_rxn',
            'fillcolor': '#90f998', 'compartment': 'e'})
        edge_ex = graph.Edge(self.node_a_extracell, node_ex,
                             {'dir': 'both', 'style': 'solid',
                              'penwidth': 1})
        self.node_list.append(node_ex)
        self.edge_list.append(edge_ex)

        self.assertTrue(all(i in self.node_list for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in self.node_list))

        self.assertTrue(all(i in self.edge_list for i in g1.edges))
        self.assertTrue(all(i in g1.edges for i in self.edge_list))

    def test2_addExrRxn_right(self):
        g2 = vis.add_exchange_rxns(self.g, 'test_EX_C', self.EX_C,
                                   {'test_EX_C': ('solid', 1)})
        node_ex = graph.Node({
            'id': 'test_EX_C', 'entry': [self.EX_C], 'shape': 'box',
            'style': 'filled', 'label': 'test_EX_C', 'type': 'Ex_rxn',
            'fillcolor': '#90f998', 'compartment': 'e'})
        edge_ex = graph.Edge(node_ex, self.node_c_extracell,
                             {'dir': 'both', 'style': 'solid',
                              'penwidth': 1})
        self.node_list.append(node_ex)
        self.edge_list.append(edge_ex)

        self.assertTrue(all(i in self.node_list for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in self.node_list))

        self.assertTrue(all(i in self.edge_list for i in g2.edges))
        self.assertTrue(all(i in g2.edges for i in self.edge_list))

    def test3_addExrCpd(self):
        cpd_formula = {'A': Formula({Atom('O'): 2, Atom('C'): 1})}
        mm_a = Compound('A', 'e')
        new_node_a_extra = graph.Node({'id': 'A[e]', 'entry': [self.cpd_a],
                                       'compartment': 'e', 'type': 'cpd'})
        g_missA = graph.Graph()
        node_list_missA = [n for n in self.node_list if n not in [
            self.node_aa, self.node_a_extracell]]
        edge_list_missA = [e for e in self.edge_list if e not in [
            self.edge_a_r2, self.edge_r2_a]]
        for node in node_list_missA:
            g_missA.add_node(node)
        for edge in edge_list_missA:
            g_missA.add_edge(edge)
        g3 = vis.add_ex_cpd(g_missA, mm_a, self.cpd_a, cpd_formula, 'C')
        node_list_missA.append(new_node_a_extra)

        # for i in g3.nodes:
        #     print(i, i.props)
        # print('fen ge')
        # for i in node_list_missA:
        #     print(i, i.props)

        self.assertTrue(all(i in node_list_missA for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in node_list_missA))

        self.assertTrue(all(i in edge_list_missA for i in g3.edges))
        self.assertTrue(all(i in g3.edges for i in edge_list_missA))

    def test4_addExrCpd_noEle(self):
        cpd_formula = {'A': Formula({Atom('O'): 2, Atom('C'): 1})}
        mm_a = Compound('A', 'e')

        g_missA = graph.Graph()
        node_list_missA = [n for n in self.node_list if n not in [
            self.node_aa, self.node_a_extracell]]
        edge_list_missA = [e for e in self.edge_list if e not in [
            self.edge_a_r2, self.edge_r2_a]]
        for node in node_list_missA:
            g_missA.add_node(node)
        for edge in edge_list_missA:
            g_missA.add_edge(edge)
        g4 = vis.add_ex_cpd(g_missA, mm_a, self.cpd_a, cpd_formula, 'N')

        self.assertTrue(all(i in node_list_missA for i in g4.nodes))
        self.assertTrue(all(i in g4.nodes for i in node_list_missA))

        self.assertTrue(all(i in edge_list_missA for i in g4.edges))
        self.assertTrue(all(i in g4.edges for i in edge_list_missA))

    def test5_addExrCpd_noFormula(self):
        cpd_formula = {}
        mm_a = Compound('A', 'e')

        g_missA = graph.Graph()
        node_list_missA = [n for n in self.node_list if n not in [
            self.node_aa, self.node_a_extracell]]
        edge_list_missA = [e for e in self.edge_list if e not in [
            self.edge_a_r2, self.edge_r2_a]]
        for node in node_list_missA:
            g_missA.add_node(node)
        for edge in edge_list_missA:
            g_missA.add_edge(edge)
        g5 = vis.add_ex_cpd(g_missA, mm_a, self.cpd_a, cpd_formula, 'C')

        self.assertTrue(all(i in node_list_missA for i in g5.nodes))
        self.assertTrue(all(i in g5.nodes for i in node_list_missA))

        self.assertTrue(all(i in edge_list_missA for i in g5.edges))
        self.assertTrue(all(i in g5.edges for i in edge_list_missA))

    def test6_addExrCpd_EleNone(self):
        cpd_formula = {'A': Formula({Atom('O'): 2, Atom('C'): 1})}
        mm_a = Compound('A', 'e')
        new_node_a_extra = graph.Node({'id': 'A[e]', 'entry': [self.cpd_a],
                                       'compartment': 'e', 'type': 'cpd'})
        g_missA = graph.Graph()
        node_list_missA = [n for n in self.node_list if n not in [
            self.node_aa, self.node_a_extracell]]
        edge_list_missA = [e for e in self.edge_list if e not in [
            self.edge_a_r2, self.edge_r2_a]]
        for node in node_list_missA:
            g_missA.add_node(node)
        for edge in edge_list_missA:
            g_missA.add_edge(edge)
        g6 = vis.add_ex_cpd(g_missA, mm_a, self.cpd_a, cpd_formula, None)
        node_list_missA.append(new_node_a_extra)

        self.assertTrue(all(i in node_list_missA for i in g6.nodes))
        self.assertTrue(all(i in g6.nodes for i in node_list_missA))

        self.assertTrue(all(i in edge_list_missA for i in g6.edges))
        self.assertTrue(all(i in g6.edges for i in edge_list_missA))


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
        e1_res = 'e'
        self.assertEqual(c1, c1_res)
        self.assertEqual(e1, e1_res)

    def test2_if_e_not_inthemodel(self):
        boundaries = [('c', 'mi')]
        extra = 'e'
        c2, e2 = vis.make_cpt_tree(boundaries, extra)
        c2_res = defaultdict(set)
        c2_res['c'].add('mi')
        c2_res['mi'].add('c')
        e2_res = 'mi'
        self.assertEqual(c2, c2_res)
        self.assertEqual(e2, e2_res)

    def test3_if_e_not_inthemodel_cpmi(self):
        boundaries = [('c', 'p'), ('c', 'mi')]
        extra = 'e'
        c3, e3 = vis.make_cpt_tree(boundaries, extra)
        c3_res = defaultdict(set)
        c3_res['c'].add('p')
        c3_res['c'].add('mi')
        c3_res['mi'].add('c')
        c3_res['p'].add('c')
        e3_res = 'p'
        self.assertEqual(c3, c3_res)
        self.assertEqual(e3, e3_res)


class TestGetCptBoundaries(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction('A[c] => A[e]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'A[c]', 'formula': parse_compound('formula_A', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'A[e]', 'formula': parse_compound('formula_A', 'e')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction('B[e] <=> B[p]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'B[e]', 'formula': parse_compound('formula_B', 'e')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'B[p]', 'formula': parse_compound('formula_B', 'p')}))
        self.native = native_model

    def test1_default_setting(self):
        e1_bound, e1_extra = vis.get_cpt_boundaries(self.native)
        e1_bound_res = set()
        e1_bound_res.add(('c', 'e'))
        e1_bound_res.add(('e', 'p'))
        e1_extra_res = 'e'
        self.assertTrue(all(i in e1_bound for i in e1_bound_res))
        self.assertEqual(e1_extra, e1_extra_res)

    def test2_with_extra_defined(self):
        self.native._properties['extracellular'] = 'p'
        e2_bound, e2_extra = vis.get_cpt_boundaries(self.native)
        e2_bound_res = set()
        e2_bound_res.add(('c', 'e'))
        e2_bound_res.add(('e', 'p'))
        e2_extra_res = 'p'
        self.assertTrue(all(i in e2_bound for i in e2_bound_res))
        self.assertEqual(e2_extra, e2_extra_res)

    def test3_e_notinmodel(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction('A[c] => A[mi]')}))
        e3_bound, e3_extra = vis.get_cpt_boundaries(native_model)
        e3_bound_res = set()
        e3_bound_res.add(('c', 'mi'))
        e3_extra_res = 'e'
        self.assertEqual(e3_bound, e3_bound_res)
        self.assertEqual(e3_extra, e3_extra_res)

    def test4_cpt_boundaries_defined_in_model(self):
        self.native.compartment_boundaries.add(('c', 'e'))
        self.native.compartment_boundaries.add(('e', 'p'))
        self.native.compartment_boundaries.add(('c', 'mi'))
        e4_bound, e4_extra = vis.get_cpt_boundaries(self.native)
        e4_bound_res = set()
        e4_bound_res.add(('c', 'e'))
        e4_bound_res.add(('e', 'p'))
        e4_bound_res.add(('c', 'mi'))
        e4_extra_res = 'e'
        self.assertTrue(all(i in e4_bound for i in e4_bound_res))
        self.assertEqual(e4_extra, e4_extra_res)


class TestFlux(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum[c] + h2o[c] <=> mal_L[c]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'fum', 'formula': 'C4H2O4'}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o', 'formula': 'H2O'}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'mal_L', 'formula': 'C4H4O5'}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'q8[c] + succ[c] => fum[c] + q8h2[c]')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'q8', 'formula': 'C49H74O4'}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'q8h2', 'formula': 'C49H76O4'}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'succ', 'formula': 'C4H4O4'}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn3', 'equation': parse_reaction(
                'h2o[c] <=> h2o[e]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o', 'formula': 'H2O'}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

    def test1_neg_flux_fba(self):
        reaction_dict = {'rxn1': (-20, 1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fba')
        self.assertEqual(Direction.Reverse,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(5, style_flux_dict['rxn1'][1])

    def test2_pos_flux_fba(self):
        reaction_dict = {'rxn1': (10, 1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fba')
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(5, style_flux_dict['rxn1'][1])

    def test3_0_flux_fba(self):
        reaction_dict = {'rxn1': (0, 1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fba')
        self.assertEqual(Direction.Both,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('dotted', style_flux_dict['rxn1'][0])
        self.assertEqual(1, style_flux_dict['rxn1'][1])

    def test4_neg_flux_fva(self):
        reaction_dict = {'rxn1': (-10, 1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Both,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('dotted', style_flux_dict['rxn2'][0])
        self.assertEqual(5, style_flux_dict['rxn1'][1])

    def test5_pos_flux_neg_fva(self):
        reaction_dict = {'rxn1': (-10, -10)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Reverse,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(5, style_flux_dict['rxn1'][1])

    def test6_pos_flux_pos_fva(self):
        reaction_dict = {'rxn1': (10, 10)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(5, style_flux_dict['rxn1'][1])

    def test7_lower_0_flux_fva(self):
        reaction_dict = {'rxn1': (0, 10)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(1, style_flux_dict['rxn1'][1])

    def test8_upper_0_flux_fva(self):
        reaction_dict = {'rxn1': (-10, 0)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Reverse,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual(1, style_flux_dict['rxn1'][1])

    def test9_both_0_flux_fva(self):
        reaction_dict = {'rxn1': (0, 0)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Both,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual('dotted', style_flux_dict['rxn1'][0])
        self.assertEqual(1, style_flux_dict['rxn1'][1])

    def test10_mult_rxn_fba(self):
        reaction_dict = {'rxn1': (2, 1), 'rxn2': (-4, 1), 'rxn3': (6, 1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fba')
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual(Direction.Reverse,
                         full_pairs_dict[self.native.reactions['rxn2']][1])
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn3']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual('solid', style_flux_dict['rxn2'][0])
        self.assertEqual('solid', style_flux_dict['rxn3'][0])
        self.assertEqual(2.5, style_flux_dict['rxn1'][1])
        self.assertEqual(5, style_flux_dict['rxn2'][1])
        self.assertEqual(7.5, style_flux_dict['rxn3'][1])

    def test11_mult_rxn_fva(self):
        reaction_dict = {'rxn1': (0, 2), 'rxn2': (-4, 1), 'rxn3': (-1, -1)}
        full_pairs_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        self.assertEqual(Direction.Forward,
                         full_pairs_dict[self.native.reactions['rxn1']][1])
        self.assertEqual(Direction.Both,
                         full_pairs_dict[self.native.reactions['rxn2']][1])
        self.assertEqual(Direction.Reverse,
                         full_pairs_dict[self.native.reactions['rxn3']][1])
        self.assertEqual('solid', style_flux_dict['rxn1'][0])
        self.assertEqual('solid', style_flux_dict['rxn2'][0])
        self.assertEqual('solid', style_flux_dict['rxn3'][0])
        self.assertEqual(1, style_flux_dict['rxn1'][1])
        self.assertEqual(8, style_flux_dict['rxn2'][1])
        self.assertEqual(2, style_flux_dict['rxn3'][1])


if __name__ == '__main__':
    unittest.main()
