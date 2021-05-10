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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth1@my.uri.edu>


from __future__ import unicode_literals

from psamm import graph
import os
import tempfile
import unittest
from collections import defaultdict
from psamm.datasource.reaction import parse_reaction, parse_compound, \
    Direction, Compound
from psamm.datasource.native import NativeModel, ReactionEntry, CompoundEntry
from psamm.formula import Formula, Atom, ParseError


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.fum = CompoundEntry({
            'id': 'fum_c', 'formula': 'C4H2O4'})
        self.rxn1 = ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')})
        self.rxn2 = ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')})
        self.g = graph.Graph()
        self.node1 = graph.Node({'id': 'A'})
        self.node2 = graph.Node({'id': 'B', 'color': 'blue'})
        self.node3 = graph.Node({'id': 'C', 'color': 'red'})
        self.node4 = graph.Node({'id': 'D',
                                 'original_id': ['A', 'B'], 'type': 'rxn'})
        self.node5 = graph.Node({'id': 'E',
                                 'original_id': 'cpd_E', 'type': 'cpd'})
        self.node6 = graph.Node({'id': 'Ex_e',
                                 'type': 'Ex_rxn'})
        self.node7 = graph.Node({'id': 'Ex_e',
                                 'type': 'cpd', 'entry': [self.fum]})
        self.node8 = graph.Node({'id': 'rxn1',
                                 'type': 'cpd', 'entry': [self.rxn1]})
        self.node9 = graph.Node({'id': 'rxn1_rxn2',
                                 'type': 'cpd', 'entry':
                                     [self.rxn1, self.rxn2]})
        self.edge1_2 = graph.Edge(self.node1, self.node2)
        self.edge2_3 = graph.Edge(self.node2, self.node3,
                                  props={'id': '2_3'})

    def test_add_node(self):
        self.g.add_node(self.node1)
        self.assertTrue(self.node1 in self.g.nodes)

    def test_node_id_dict(self):
        self.g.add_node(self.node1)
        self.assertEqual(self.node1, self.g._nodes_id[self.node1.props['id']])

    def test_add_multiple_nodes(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.assertTrue(all(i in self.g.nodes
                            for i in [self.node1, self.node2]))

    def test_original_id_EX(self):
        self.g.add_node(self.node6)
        nd = defaultdict(list)
        nd['Ex_e'].append(self.node6)
        self.assertTrue(self.g._nodes_original_id == nd)

    def test_original_id_cpd(self):
        self.g.add_node(self.node7)
        nd = defaultdict(list)
        nd[self.fum.id].append(self.node7)
        self.assertTrue(self.g._nodes_original_id == nd)

    def test_original_id_rxn(self):
        self.g.add_node(self.node8)
        nd = defaultdict(list)
        nd[self.rxn1.id].append(self.node8)
        self.assertTrue(self.g._nodes_original_id == nd)

    def test_original_id_rxn(self):
        self.g.add_node(self.node9)
        nd = defaultdict(list)
        nd['rxn1,rxn2'].append(self.node9)
        self.assertEqual(self.g._nodes_original_id, nd)

    def test_node_count(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.assertEqual(self.g.node_count, 3)

    def test_add_edge(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.g.add_edge(self.edge1_2)
        self.assertTrue(self.edge1_2 in self.g.edges)

    def test_add_edge_with_props(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.assertTrue(self.edge2_3 in self.g.edges)

    def test_add_multiple_edges(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge1_2)
        self.g.add_edge(self.edge2_3)
        self.assertTrue(all(i in self.g.edges
                            for i in [self.edge1_2, self.edge2_3]))

    def test_edge_count(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge1_2)
        self.g.add_edge(self.edge2_3)
        self.assertTrue(self.g.edge_count, 2)

    def test_add_edge_missing_node(self):
        self.g.add_node(self.node1)
        with self.assertRaises(ValueError):
            self.g.add_edge(self.edge1_2)

    def test_edges_for(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge1_2)
        self.g.add_edge(self.edge2_3)
        edges_for_list = []
        for i in self.g.edges_for(self.node1):
            edges_for_list.append(i)
        self.assertTrue(edges_for_list == [(self.node2,
                                            set([self.edge1_2]))])
        self.assertEqual('A', 'A')

    def test_get_node(self):
        self.g.add_node(self.node1)
        test_node = self.g.get_node('A')
        self.assertEqual(test_node, self.node1)

    def test_nodes_id_dict(self):
        self.g.add_node(self.node1)
        self.g.add_node(self.node2)
        self.assertEqual(self.g.nodes_id_dict, {'A': self.node1,
                                                'B': self.node2})

    def test_default_edge_props(self):
        self.g._default_edge_props['style'] = 'dashed'
        d = {'style': 'dashed'}
        self.assertEqual(d, self.g.default_edge_props)

    def test_default_node_props(self):
        self.g._default_node_props['shape'] = 'box'
        d = {'shape': 'box'}
        self.assertEqual(d, self.g.default_node_props)

    def test_write_nodes_table(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_node')
        with open(path, mode='w') as f:
            self.g.write_nodes_tables(f)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file, ['id\tcolor\tlabel\n',
                                     'B\tblue\tB\n', 'C\tred\tC\n'])

    def test_write_nodes_table_with_original_id(self):
        self.g.add_node(self.node2)
        self.node3.props['original_id'] = 'A_1'
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_node')
        with open(path, mode='w') as f:
            self.g.write_nodes_tables(f)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file, ['id\tcolor\tlabel\n',
                                     'B\tblue\tB\n', 'C\tred\tC\n'])

    def test_write_edges_tables(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_edges_tables(f)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file, ['source\ttarget\tid\n', 'B\tC\t2_3\n'])

    def test_write_graphviz(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz(f, 2, 2)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(['digraph {\n', 'size = "2, 2"; ratio = fill;\n',
                          ' "B"[color="blue",id="B"]\n',
                          ' "C"[color="red",id="C"]\n',
                          ' "B" -> "C"[id="2_3"]\n', '}\n'], read_file)

    def test_write_graphviz_graph_props(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g.props['fontsize'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz(f, 2, 2)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(['digraph {\n', 'size = "2, 2"; ratio = fill;\n',
                          ' fontsize="12";\n', ' "B"[color="blue",id="B"]\n',
                          ' "C"[color="red",id="C"]\n',
                          ' "B" -> "C"[id="2_3"]\n', '}\n'], read_file)

    def test_write_graphviz_default_node_props(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g._default_node_props['fontname'] = 'Arial'
        self.g._default_node_props['fontsize'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz(f, 2, 2)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(['digraph {\n', 'size = "2, 2"; ratio = fill;\n',
                          ' node[fontname="Arial",fontsize="12"];\n',
                          ' "B"[color="blue",id="B"]\n',
                          ' "C"[color="red",id="C"]\n',
                          ' "B" -> "C"[id="2_3"]\n', '}\n'], read_file)

    def test_write_graphvizdefault_edge_props(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g._default_edge_props['style'] = 'dashed'
        self.g._default_edge_props['width'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz(f, 2, 2)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(['digraph {\n', 'size = "2, 2"; ratio = fill;\n',
                          ' edge[style="dashed",width="12"];\n',
                          ' "B"[color="blue",id="B"]\n',
                          ' "C"[color="red",id="C"]\n',
                          ' "B" -> "C"[id="2_3"]\n', '}\n'], read_file)

    def test_write_graphviz_default_size(self):
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g._default_node_props['fontname'] = 'Arial'
        self.g._default_node_props['fontsize'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz(f, None, None)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n',
                          ' node[fontname="Arial",fontsize="12"];\n',
                          ' "B"[color="blue",id="B"]\n',
                          ' "C"[color="red",id="C"]\n',
                          ' "B" -> "C"[id="2_3"]\n', '}\n'])

    def test_write_graphviz_compartmentalize(self):
        self.node2.props['compartment'] = 'c'
        self.node3.props['compartment'] = 'e'
        self.edge2_3.props['compartment'] = 'e'
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz_compartmentalized(
                f, {'e': {'c'}, 'c': set()}, 'e', 2, 2)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n',
                          'size="2,2"; ratio = fill;\n',
                          ' subgraph cluster_e {\n',
                          '  style=solid;\n', '  color=black;\n',
                          '  penwidth=4;\n', '  fontsize=35;\n',
                          '  label = "Compartment: e"\n',
                          ' "C"[color="red",compartment="e",id="C"]\n',
                          ' subgraph cluster_c {\n', '  style=dashed;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n', '  label = "Compartment: c"\n',
                          ' "B"[color="blue",compartment="c",id="B"]\n',
                          '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])

    def test_write_graphviz_compartmentalize_default_size(self):
        self.node2.props['compartment'] = 'c'
        self.node3.props['compartment'] = 'e'
        self.edge2_3.props['compartment'] = 'e'
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz_compartmentalized(
                f, {'e': {'c'}, 'c': set()}, 'e', None, None)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n', ' subgraph cluster_e {\n',
                          '  style=solid;\n', '  color=black;\n',
                          '  penwidth=4;\n', '  fontsize=35;\n',
                          '  label = "Compartment: e"\n',
                          ' "C"[color="red",compartment="e",id="C"]\n',
                          ' subgraph cluster_c {\n', '  style=dashed;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n',
                          '  label = "Compartment: c"\n',
                          ' "B"[color="blue",compartment="c",id="B"]\n',
                          '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])

    def test_write_graphviz_compartmentalize_default_node_props(self):
        self.node2.props['compartment'] = 'c'
        self.node3.props['compartment'] = 'e'
        self.edge2_3.props['compartment'] = 'e'
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g._default_node_props['fontname'] = 'Arial'
        self.g._default_node_props['fontsize'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz_compartmentalized(
                f, {'e': {'c'}, 'c': set()}, 'e', None, None)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n',
                          ' node[fontname="Arial",fontsize="12"];\n',
                          ' subgraph cluster_e {\n', '  style=solid;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n', '  label = "Compartment: e"\n',
                          ' "C"[color="red",compartment="e",id="C"]\n',
                          ' subgraph cluster_c {\n', '  style=dashed;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n', '  label = "Compartment: c"\n',
                          ' "B"[color="blue",compartment="c",id="B"]\n',
                          '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])

    def test_write_graphviz_compartmentalize_default_edge_props(self):
        self.node2.props['compartment'] = 'c'
        self.node3.props['compartment'] = 'e'
        self.edge2_3.props['compartment'] = 'e'
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g._default_edge_props['style'] = 'dashed'
        self.g._default_edge_props['width'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz_compartmentalized(
                f, {'e': {'c'}, 'c': set()}, 'e', None, None)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n',
                          ' edge[style="dashed",width="12"];\n',
                          ' subgraph cluster_e {\n', '  style=solid;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n', '  label = "Compartment: e"\n',
                          ' "C"[color="red",compartment="e",id="C"]\n',
                          ' subgraph cluster_c {\n', '  style=dashed;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n', '  label = "Compartment: c"\n',
                          ' "B"[color="blue",compartment="c",id="B"]\n',
                          '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])

    def test_write_graphviz_compartmentalize_default_graph_props(self):
        self.node2.props['compartment'] = 'c'
        self.node3.props['compartment'] = 'e'
        self.edge2_3.props['compartment'] = 'e'
        self.g.add_node(self.node2)
        self.g.add_node(self.node3)
        self.g.add_edge(self.edge2_3)
        self.g.props['fontsize'] = 12
        path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
        with open(path, mode='w') as f:
            self.g.write_graphviz_compartmentalized(
                f, {'e': {'c'}, 'c': set()}, 'e', None, None)
        read_file = []
        f = open(path, 'r')
        for i in f.readlines():
            read_file.append(i)
        f.close()
        self.assertEqual(read_file,
                         ['digraph {\n', ' fontsize="12";\n',
                          ' subgraph cluster_e {\n',
                          '  style=solid;\n', '  color=black;\n',
                          '  penwidth=4;\n', '  fontsize=35;\n',
                          '  label = "Compartment: e"\n',
                          ' "C"[color="red",compartment="e",id="C"]\n',
                          ' subgraph cluster_c {\n', '  style=dashed;\n',
                          '  color=black;\n', '  penwidth=4;\n',
                          '  fontsize=35;\n',
                          '  label = "Compartment: c"\n',
                          ' "B"[color="blue",compartment="c",id="B"]\n',
                          '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])


class TestNodes(unittest.TestCase):
    def setUp(self):
        self.g = graph.Graph()
        self.node1 = graph.Node({'id': 'A'})
        self.node2 = graph.Node({'id': 'B'})
        self.node3 = graph.Node({'id': 'C'})
        self.edge1_2 = graph.Edge(self.node1, self.node2)

    def test_eq_nodes(self):
        self.g.add_node(self.node1)
        self.assertEqual(self.node1 == self.g.get_node('A'),
                         self.node1.props == self.g.get_node('A').props)

    def test_neq_nodes(self):
        self.assertTrue(self.node1 != self.node2)

    def test_neq_node_edge(self):
        self.assertFalse(self.node1 == self.edge1_2)

    def test_node_repr(self):
        rep = self.node1.__repr__()
        self.assertEqual(rep, '<Node id=A>')


class TestOther(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        self.rxn1 = ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')})
        native_model.reactions.add_entry(self.rxn1)

        self.fum = CompoundEntry({
            'id': 'fum_c', 'formula': 'C4H2O4'})
        native_model.compounds.add_entry(self.fum)

        self.h2o = CompoundEntry({
            'id': 'h2o_c', 'formula': 'H2O'})
        native_model.compounds.add_entry(self.h2o)

        self.mal = CompoundEntry(
            {'id': 'mal_L_c', 'formula': 'C4H4O5'})
        native_model.compounds.add_entry(self.mal)

        self.rxn2 = ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'q8_c[c] + succ_c[c] => fum_c[c] + q8h2_c[c]')})
        native_model.reactions.add_entry(self.rxn2)

        self.q8 = CompoundEntry(
            {'id': 'q8_c', 'formula': 'C49H74O4'})
        native_model.compounds.add_entry(self.q8)

        self.q8h2 = CompoundEntry({
            'id': 'q8h2_c', 'formula': 'C49H76O4'})
        native_model.compounds.add_entry(self.q8h2)

        self.succ = CompoundEntry({
            'id': 'succ_c', 'formula': 'C4H4O4'})
        native_model.compounds.add_entry(self.succ)

        self.native = native_model
        self.mm = native_model.create_metabolic_model()

    def test_compound_dict(self):
        cpd_dict = graph.get_compound_dict(self.native)
        test_dict = {'fum_c': Formula({Atom('O'): 4,
                                       Atom('C'): 4, Atom('H'): 2}),
                     'h2o_c': Formula({Atom('O'): 1,
                                       Atom('H'): 2}),
                     'mal_L_c': Formula({Atom('O'): 5,
                                         Atom('C'): 4, Atom('H'): 4}),
                     'q8_c': Formula({Atom('O'): 4,
                                      Atom('C'): 49, Atom('H'): 74}),
                     'q8h2_c': Formula({Atom('O'): 4,
                                        Atom('C'): 49, Atom('H'): 76}),
                     'succ_c': Formula({Atom('O'): 4,
                                        Atom('C'): 4, Atom('H'): 4})}
        self.assertEqual(test_dict, cpd_dict)

    def test_network_dict(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[])

        test_dict = {self.rxn1: ([(Compound(u'fum_c', u'c'),
                                   Compound(u'mal_L_c', u'c')),
                                  (Compound(u'h2o_c', u'c'),
                                   Compound(u'mal_L_c', u'c'))],
                                 Direction.Both),
                     self.rxn2: ([(Compound(u'q8_c', u'c'),
                                   Compound(u'q8h2_c', u'c')),
                                  (Compound(u'succ_c', u'c'),
                                   Compound(u'fum_c', u'c')),
                                  (Compound(u'succ_c', u'c'),
                                   Compound(u'q8h2_c', u'c'))],
                                 Direction.Right)}
        self.assertEqual(net_dict, test_dict)

    def test_network_dict_nofpp(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native, self.mm,
                                    subset=None, method='no-fpp',
                                    element=None, excluded_reactions=[])
        test_dict = {self.rxn1: ([(Compound(u'fum_c', u'c'),
                                   Compound(u'mal_L_c', u'c')),
                                  (Compound(u'h2o_c', u'c'),
                                   Compound(u'mal_L_c', u'c'))],
                                 Direction.Both),
                     self.rxn2: ([(Compound(u'q8_c', u'c'),
                                   Compound(u'fum_c', u'c')),
                                  (Compound(u'q8_c', u'c'),
                                   Compound(u'q8h2_c', u'c')),
                                  (Compound(u'succ_c', u'c'),
                                   Compound(u'fum_c', u'c')),
                                  (Compound(u'succ_c', u'c'),
                                   Compound(u'q8h2_c', u'c')),
                                  ],
                                 Direction.Right)}
        self.assertEqual(net_dict, test_dict)

    def test_reaction_dir_for(self):
        self.assertTrue(graph.dir_value(Direction.Forward) == 'forward')

    def test_reaction_dir_rev(self):
        self.assertTrue(graph.dir_value(Direction.Reverse) == 'back')

    def test_reaction_dir_both(self):
        self.assertTrue(graph.dir_value(Direction.Both) == 'both')


class TestMakeNetworks(unittest.TestCase):
    def setUp(self):
        self.native_model = NativeModel()
        self.atp = CompoundEntry({
            'id': 'atp', 'formula': 'C10H16N5O13P3'})
        self.adp = CompoundEntry({
            'id': 'adp', 'formula': 'C10H15N5O10P2'})
        self.glc = CompoundEntry({
            'id': 'glc', 'formula': 'C6H12O6'})
        self.g6p = CompoundEntry({
            'id': 'g6p', 'formula': 'C6H13O9P'})
        self.f6p = CompoundEntry({
            'id': 'f6p', 'formula': 'C6H13O9P'})
        self.rxn1 = ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'atp + glc => adp + g6p')})
        self.rxn2 = ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'g6p <=> f6p')})
        self.native_model.compounds.add_entry(self.atp)
        self.native_model.compounds.add_entry(self.adp)
        self.native_model.compounds.add_entry(self.g6p)
        self.native_model.compounds.add_entry(self.glc)
        self.native_model.compounds.add_entry(self.f6p)
        self.native_model.reactions.add_entry(self.rxn1)
        self.native_model.reactions.add_entry(self.rxn2)
        self.mm = self.native_model.create_metabolic_model()
        self.node_atp = graph.Node({'id': 'atp', 'entry': Compound('atp')})
        self.node_adp = graph.Node({'id': 'adp', 'entry': Compound('adp')})
        self.node_g6p = graph.Node({'id': 'g6p', 'entry': Compound('g6p')})
        self.node_glc = graph.Node({'id': 'glc', 'entry': Compound('glc')})
        self.node_f6p = graph.Node({'id': 'f6p', 'entry': Compound('f6p')})
        self.edge_1 = graph.Edge(self.node_glc, self.node_g6p,
                                 props={'id': 'g6p_glc_forward',
                                        'dir': 'forward'})
        self.edge_2 = graph.Edge(self.node_g6p, self.node_f6p,
                                 props={'id': 'f6p_g6p_both',
                                        'dir': 'both'})
        self.edge_3 = graph.Edge(self.node_atp, self.node_adp,
                                 props={'id': 'adp_atp_forward',
                                        'dir': 'forward'})
        self.edge_4 = graph.Edge(self.node_atp, self.node_adp,
                                 props={'id': 'adp_atp_forward',
                                        'dir': 'forward'})
        self.edge_5 = graph.Edge(self.node_atp, self.node_g6p,
                                 props={'id': 'atp_g6p_forward',
                                        'dir': 'forward'})
        self.edge_6 = graph.Edge(self.node_glc, self.node_adp,
                                 props={'id': 'adp_glc_forward',
                                        'dir': 'forward'})

        self.model_compound_entries = {}
        for cpd in self.native_model.compounds:
            self.model_compound_entries[cpd.id] = cpd
        self.node_atp_bip = graph.Node(
            {'id': 'atp', 'entry': [self.model_compound_entries['atp']],
             'compartment': None, 'type': 'cpd'})
        self.node_adp_bip = graph.Node(
            {'id': 'adp', 'entry': [self.model_compound_entries['adp']],
             'compartment': None, 'type': 'cpd'})
        self.node_g6p_bip = graph.Node(
            {'id': 'g6p', 'entry': [self.model_compound_entries['g6p']],
             'compartment': None, 'type': 'cpd'})
        self.node_glc_bip = graph.Node(
            {'id': 'glc', 'entry': [self.model_compound_entries['glc']],
             'compartment': None, 'type': 'cpd'})
        self.node_f6p_bip = graph.Node(
            {'id': 'f6p', 'entry': [self.model_compound_entries['f6p']],
             'compartment': None, 'type': 'cpd'})
        self.node_rxn1_1 = graph.Node({'id': 'rxn1_1', 'entry': [self.rxn1],
                                       'compartment': None, 'type': 'rxn'})
        self.node_rxn1_2 = graph.Node({'id': 'rxn1_2', 'entry': [self.rxn1],
                                       'compartment': None, 'type': 'rxn'})
        self.node_rxn2_1 = graph.Node({'id': 'rxn2_1', 'entry': [self.rxn2],
                                       'compartment': None, 'type': 'rxn'})
        self.edge_bip_1 = graph.Edge(self.node_atp_bip, self.node_rxn1_1,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_2 = graph.Edge(self.node_rxn1_1, self.node_adp_bip,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_3 = graph.Edge(self.node_atp_bip, self.node_rxn1_1,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_4 = graph.Edge(self.node_rxn1_1, self.node_g6p_bip,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_5 = graph.Edge(self.node_glc_bip, self.node_rxn1_2,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_6 = graph.Edge(self.node_rxn1_2, self.node_g6p_bip,
                                     props={'dir': 'forward',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_7 = graph.Edge(self.node_g6p_bip, self.node_rxn2_1,
                                     props={'dir': 'both',
                                            'style': 'solid', 'penwidth': 1})
        self.edge_bip_8 = graph.Edge(self.node_rxn2_1, self.node_f6p_bip,
                                     props={'dir': 'both',
                                            'style': 'solid', 'penwidth': 1})
        self.node_list = [
            self.node_atp_bip, self.node_adp_bip, self.node_g6p_bip,
            self.node_glc_bip, self.node_f6p_bip,
            self.node_rxn1_1, self.node_rxn1_2, self.node_rxn2_1]
        self.edge_list = [
            self.edge_bip_1, self.edge_bip_2, self.edge_bip_3,
            self.edge_bip_4, self.edge_bip_5, self.edge_bip_6,
            self.edge_bip_7, self.edge_bip_8]

        self.model_compound_entries = {}
        for cpd in self.native_model.compounds:
            self.model_compound_entries[cpd.id] = cpd

    def test_compound_graph(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        compound_graph = graph.make_compound_graph(net_dict)
        self.assertTrue(all(i in compound_graph.nodes
                            for i in [self.node_atp, self.node_adp,
                                      self.node_g6p, self.node_glc,
                                      self.node_f6p]))
        self.assertTrue(
            all(i in compound_graph.edges
                for i in [self.edge_1, self.edge_2,
                          self.edge_3, self.edge_4, self.edge_5]))
        self.assertTrue(self.edge_6 not in compound_graph.edges)

    def test_compound_graph_filter(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element='C', excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        compound_graph = graph.make_compound_graph(net_dict)
        self.assertTrue(all(i in compound_graph.nodes
                            for i in [self.node_atp, self.node_adp,
                                      self.node_g6p, self.node_glc,
                                      self.node_f6p]))
        self.assertTrue(all(i in compound_graph.edges
                            for i in [self.edge_1, self.edge_2,
                                      self.edge_3, self.edge_4]))
        self.assertTrue(self.edge_5 not in compound_graph.edges)
        self.assertTrue(self.edge_6 not in compound_graph.edges)

    def test_compound_graph_subset(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=['rxn1'], method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        compound_graph = graph.make_compound_graph(net_dict)
        self.assertTrue(
            all(i in compound_graph.nodes
                for i in [self.node_atp, self.node_adp,
                          self.node_g6p, self.node_glc]))
        self.assertTrue(self.node_f6p not in compound_graph.nodes)
        self.assertTrue(all(i in compound_graph.edges
                            for i in [self.edge_1, self.edge_3,
                                      self.edge_4, self.edge_5]))
        self.assertTrue(self.edge_2 not in compound_graph.edges)
        self.assertTrue(self.edge_6 not in compound_graph.edges)

    def test_compound_graph_exclude(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=['rxn2'],
                                    reaction_dict={}, analysis=None)
        compound_graph = graph.make_compound_graph(net_dict)
        self.assertTrue(
            all(i in compound_graph.nodes
                for i in [self.node_atp, self.node_adp,
                          self.node_g6p, self.node_glc]))
        self.assertTrue(self.node_f6p not in compound_graph.nodes)
        self.assertTrue(all(i in compound_graph.edges
                            for i in [self.edge_1, self.edge_3,
                                      self.edge_4, self.edge_5]))
        self.assertTrue(self.edge_2 not in compound_graph.edges)
        self.assertTrue(self.edge_6 not in compound_graph.edges)

    def test_compound_graph_nofpp(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='no-fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        compound_graph = graph.make_compound_graph(net_dict)
        self.assertTrue(all(i in compound_graph.nodes for i in
                            [self.node_f6p, self.node_g6p, self.node_glc,
                             self.node_adp, self.node_atp]))
        self.assertTrue(all(i in compound_graph.edges for i in
                            [self.edge_1, self.edge_2, self.edge_3,
                             self.edge_4, self.edge_5, self.edge_6]))

    def test_make_cpair_dict_combine_0(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_2']
        test_dict[(Compound('g6p'), Compound('f6p'))]['both'] = ['rxn2_1']

        test_new_id = {u'rxn2_1': self.rxn2,
                       u'rxn1_2': self.rxn1, u'rxn1_1': self.rxn1}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_make_cpair_dict_combine_1(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 1, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('g6p'), Compound('f6p'))]['both'] = ['rxn2_1']

        test_new_id = {u'rxn2_1': self.rxn2, u'rxn1_1': self.rxn1}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_make_cpair_dict_combine_2(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 2, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('g6p'), Compound('f6p'))]['both'] = ['rxn2_1']

        test_new_id = {u'rxn2_1': self.rxn2, u'rxn1_1': self.rxn1}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_make_cpair_dict_nofpp(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='no-fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'no-fpp', 0, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1']
        test_dict[(Compound('glc'), Compound('adp'))]['forward'] = ['rxn1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1']
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1']
        test_dict[(Compound('g6p'), Compound('f6p'))]['both'] = ['rxn2']

        test_new_id = {u'rxn2': self.rxn2, u'rxn1': self.rxn1}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_cpair_dict_duplicate_cpair_both(self):
        rxn3 = ReactionEntry({
            'id': 'rxn3', 'equation': parse_reaction(
                'f6p <=> g6p')})
        nm = self.native_model
        nm.reactions.add_entry(rxn3)
        mm = self.native_model.create_metabolic_model()

        net_dict, style_flux_dict = \
            graph.make_network_dict(nm, mm, subset=None,
                                    method='fpp', element=None,
                                    excluded_reactions=[], reaction_dict={},
                                    analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_2']
        test_dict[(Compound('f6p'),
                   Compound('g6p'))]['both'] = ['rxn2_1', 'rxn3_1']

        test_new_id = {u'rxn2_1': self.rxn2, u'rxn1_2': self.rxn1,
                       u'rxn1_1': self.rxn1, 'rxn3_1': rxn3}
        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_cpair_dict_duplicate_cpair_rev(self):
        rxn3 = ReactionEntry({
            'id': 'rxn3', 'equation': parse_reaction(
                'g6p => f6p')})
        nm = self.native_model
        nm.reactions.add_entry(rxn3)
        mm = self.native_model.create_metabolic_model()

        net_dict, style_flux_dict = \
            graph.make_network_dict(nm, mm, subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_2']
        test_dict[(Compound('g6p'), Compound('f6p'))]['both'] = ['rxn2_1']
        test_dict[(Compound('g6p'), Compound('f6p'))]['forward'] = ['rxn3_1']

        test_new_id = {u'rxn1_1': self.rxn1, u'rxn1_2': self.rxn1,
                       u'rxn2_1': self.rxn2, 'rxn3_1': rxn3}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_cpair_dict_duplicate_for(self):
        rxn3 = ReactionEntry({
            'id': 'rxn3', 'equation': parse_reaction(
                'f6p => g6p')})
        nm = self.native_model
        nm.reactions.add_entry(rxn3)
        mm = self.native_model.create_metabolic_model()

        net_dict, style_flux_dict = \
            graph.make_network_dict(nm, mm, subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        test_dict = defaultdict(lambda: defaultdict(list))
        test_dict[(Compound('atp'), Compound('adp'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('atp'), Compound('g6p'))]['forward'] = ['rxn1_1']
        test_dict[(Compound('glc'), Compound('g6p'))]['forward'] = ['rxn1_2']
        test_dict[(Compound('f6p'), Compound('g6p'))]['both'] = ['rxn2_1']
        test_dict[(Compound('f6p'), Compound('g6p'))]['forward'] = ['rxn3_1']

        test_new_id = {u'rxn2_1': self.rxn2, u'rxn1_2': self.rxn1,
                       u'rxn1_1': self.rxn1, 'rxn3_1': rxn3}

        self.assertEqual(cpairs, test_dict)
        self.assertEqual(new_id, test_new_id)

    def test_node_original_id_dict(self):
        g = graph.Graph()
        g.add_node(self.node_atp)
        d = defaultdict(list)
        d['atp'].append(self.node_atp)
        self.assertEqual(g.nodes_original_id_dict, d)

    def test_bipartite_graph(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        self.assertTrue(all(i in bipartite_graph.nodes
                            for i in self.node_list))
        self.assertTrue(all(i in self.node_list
                            for i in bipartite_graph.nodes))
        self.assertTrue(all(i in bipartite_graph.edges
                            for i in self.edge_list))
        self.assertTrue(all(i in self.edge_list
                            for i in bipartite_graph.edges))

    def test_bipartite_graph_filter(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element='C', excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        edge_list = [self.edge_bip_1, self.edge_bip_2, self.edge_bip_5,
                     self.edge_bip_6, self.edge_bip_7, self.edge_bip_8]
        self.assertTrue(all(i in bipartite_graph.nodes
                            for i in self.node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in self.node_list
                            for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_filter2(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element='P', excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        node_list = [i for i in self.node_list
                     if i not in [self.node_glc_bip, self.node_rxn1_2]]
        edge_list = [self.edge_bip_1, self.edge_bip_2, self.edge_bip_3,
                     self.edge_bip_4, self.edge_bip_7, self.edge_bip_8]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_filter3(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        edge_list = [self.edge_bip_1, self.edge_bip_2, self.edge_bip_3,
                     self.edge_bip_4, self.edge_bip_5, self.edge_bip_6,
                     self.edge_bip_7, self.edge_bip_8]
        self.assertTrue(all(i in bipartite_graph.nodes
                            for i in self.node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in self.node_list
                            for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_subset(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=['rxn1'], method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        node_list = [i for i in self.node_list
                     if i not in [self.node_f6p_bip, self.node_rxn2_1]]
        edge_list = [self.edge_bip_1, self.edge_bip_2, self.edge_bip_3,
                     self.edge_bip_4, self.edge_bip_5, self.edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_subset2(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=['rxn2'], method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        node_list = [self.node_g6p_bip, self.node_f6p_bip, self.node_rxn2_1]
        edge_list = [self.edge_bip_7, self.edge_bip_8]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_exclude(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=['rxn2'],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        node_list = [i for i in self.node_list
                     if i not in [self.node_f6p_bip, self.node_rxn2_1]]
        edge_list = [self.edge_bip_1, self.edge_bip_2, self.edge_bip_3,
                     self.edge_bip_4, self.edge_bip_5, self.edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_nofpp(self):
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='no-fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict={}, analysis=None)
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'no-fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'no-fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict)
        node_rxn1 = graph.Node({'id': 'rxn1', 'entry': [self.rxn1],
                                'compartment': None, 'type': 'rxn'})
        node_rxn2 = graph.Node({'id': 'rxn2', 'entry': [self.rxn2],
                                'compartment': None, 'type': 'rxn'})
        edge_bip_1 = graph.Edge(self.node_atp_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 1})
        edge_bip_2 = graph.Edge(self.node_glc_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 1})
        edge_bip_3 = graph.Edge(node_rxn1, self.node_adp_bip,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 1})
        edge_bip_4 = graph.Edge(node_rxn1, self.node_g6p_bip,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 1})
        edge_bip_5 = graph.Edge(self.node_g6p_bip, node_rxn2,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 1})
        edge_bip_6 = graph.Edge(node_rxn2, self.node_f6p_bip,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 1})
        node_list = [self.node_atp_bip, self.node_adp_bip, self.node_g6p_bip,
                     self.node_glc_bip, self.node_f6p_bip,
                     node_rxn1, node_rxn2]
        edge_list = [edge_bip_1, edge_bip_2, edge_bip_3,
                     edge_bip_4, edge_bip_5, edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_fba(self):
        reaction_dict = {'rxn1': (0, 1), 'rxn2': (2, 1)}
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='no-fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict=reaction_dict,
                                    analysis='fba')
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'no-fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'no-fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict, 'fba')
        node_rxn1 = graph.Node({'id': 'rxn1', 'entry': [self.rxn1],
                                'compartment': None, 'type': 'rxn'})
        node_rxn2 = graph.Node({'id': 'rxn2', 'entry': [self.rxn2],
                                'compartment': None, 'type': 'rxn'})
        edge_bip_1 = graph.Edge(self.node_atp_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_2 = graph.Edge(self.node_glc_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_3 = graph.Edge(node_rxn1, self.node_adp_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_4 = graph.Edge(node_rxn1, self.node_g6p_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_5 = graph.Edge(self.node_g6p_bip, node_rxn2,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 5})
        edge_bip_6 = graph.Edge(node_rxn2, self.node_f6p_bip,
                                props={'dir': 'forward',
                                       'style': 'solid', 'penwidth': 5})
        node_list = [self.node_atp_bip, self.node_adp_bip,
                     self.node_g6p_bip, self.node_glc_bip,
                     self.node_f6p_bip, node_rxn1, node_rxn2]
        edge_list = [edge_bip_1, edge_bip_2, edge_bip_3,
                     edge_bip_4, edge_bip_5, edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_fva(self):
        reaction_dict = {'rxn1': (0, 0), 'rxn2': (-1, 1)}
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='no-fpp',
                                    element='C', excluded_reactions=[],
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'no-fpp', 0, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'no-fpp', 0,
                                              self.model_compound_entries,
                                              new_style_flux_dict, 'fva')
        node_rxn1 = graph.Node({'id': 'rxn1', 'entry': [self.rxn1],
                                'compartment': None, 'type': 'rxn'})
        node_rxn2 = graph.Node({'id': 'rxn2', 'entry': [self.rxn2],
                                'compartment': None, 'type': 'rxn'})
        edge_bip_1 = graph.Edge(self.node_atp_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_2 = graph.Edge(self.node_glc_bip, node_rxn1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_3 = graph.Edge(node_rxn1, self.node_adp_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_4 = graph.Edge(node_rxn1, self.node_g6p_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_5 = graph.Edge(self.node_g6p_bip, node_rxn2,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 5})
        edge_bip_6 = graph.Edge(node_rxn2, self.node_f6p_bip,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 5})
        node_list = [self.node_atp_bip, self.node_adp_bip,
                     self.node_g6p_bip, self.node_glc_bip, self.node_f6p_bip,
                     node_rxn1, node_rxn2]
        edge_list = [edge_bip_1, edge_bip_2, edge_bip_3,
                     edge_bip_4, edge_bip_5, edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))

    def test_bipartite_graph_combine_2(self):
        reaction_dict = {'rxn1': (0, 0), 'rxn2': (-1, 1)}
        net_dict, style_flux_dict = \
            graph.make_network_dict(self.native_model, self.mm,
                                    subset=None, method='fpp',
                                    element=None, excluded_reactions=[],
                                    reaction_dict=reaction_dict,
                                    analysis='fva')
        cpairs, new_id, new_style_flux_dict = graph.make_cpair_dict(
            net_dict, 'fpp', 2, style_flux_dict)
        bipartite_graph = \
            graph.make_bipartite_graph_object(cpairs, new_id, 'fpp', 2,
                                              self.model_compound_entries,
                                              new_style_flux_dict, 'fva')
        node_rxn1_1 = graph.Node({'id': 'rxn1_1', 'entry': [self.rxn1],
                                  'compartment': None, 'type': 'rxn'})
        node_rxn2_1 = graph.Node({'id': 'rxn2_1', 'entry': [self.rxn2],
                                  'compartment': None, 'type': 'rxn'})
        edge_bip_1 = graph.Edge(self.node_atp_bip, node_rxn1_1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_2 = graph.Edge(node_rxn1_1, self.node_adp_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_3 = graph.Edge(node_rxn1_1, self.node_g6p_bip,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        edge_bip_4 = graph.Edge(self.node_g6p_bip, node_rxn2_1,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 5})
        edge_bip_5 = graph.Edge(node_rxn2_1, self.node_f6p_bip,
                                props={'dir': 'both',
                                       'style': 'solid', 'penwidth': 5})
        edge_bip_6 = graph.Edge(self.node_glc_bip, node_rxn1_1,
                                props={'dir': 'forward',
                                       'style': 'dotted', 'penwidth': 1})
        node_list = [self.node_atp_bip, self.node_adp_bip,
                     self.node_g6p_bip, self.node_glc_bip, self.node_f6p_bip,
                     node_rxn1_1, node_rxn2_1]
        edge_list = [edge_bip_1, edge_bip_2, edge_bip_3,
                     edge_bip_4, edge_bip_5, edge_bip_6]
        self.assertTrue(all(i in bipartite_graph.nodes for i in node_list))
        self.assertTrue(all(i in bipartite_graph.edges for i in edge_list))
        self.assertTrue(all(i in node_list for i in bipartite_graph.nodes))
        self.assertTrue(all(i in edge_list for i in bipartite_graph.edges))


class TestEdges(unittest.TestCase):
    def setUp(self):
        self.g = graph.Graph()
        self.node1 = graph.Node({'id': 'A'})
        self.node2 = graph.Node({'id': 'B'})
        self.node3 = graph.Node({'id': 'C'})
        self.node4 = graph.Node({'id': 'A'})
        self.edge1_2 = graph.Edge(self.node1, self.node2)
        self.edge2_3 = graph.Edge(self.node2, self.node3)
        self.edge12_21 = graph.Edge(self.node1, self.node2)

    def test_eq_edges(self):
        self.assertTrue(self.edge1_2 == self.edge12_21)

    def test_neq_edges(self):
        self.assertTrue(self.edge1_2 != self.edge2_3)

    def test_neq_edge_node(self):
        self.assertFalse(self.edge1_2 == self.node1)


class TestCompExit(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        self.rxn1 = ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')})
        native_model.reactions.add_entry(self.rxn1)
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

    def test_sys_exit(self):
        with self.assertRaises(SystemExit):
            graph.make_network_dict(self.native, self.mm)

    def test_with_element(self):
        with self.assertRaises(SystemExit):
            graph.make_network_dict(self.native, self.mm,
                                    method='no-fpp', element='C')


if __name__ == '__main__':
    unittest.main()
