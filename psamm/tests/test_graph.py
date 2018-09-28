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
# Copyright 2015-2018  Keith Dufault-Thompson <keitht547@uri.edu>


from __future__ import unicode_literals

from psamm import graph
import os
import tempfile
import unittest
from psamm.datasource.native import NativeModel
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm import formula

from psamm.datasource.native import NativeModel
from psamm.datasource.native import ReactionEntry
from psamm.datasource.native import CompoundEntry
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm.formula import Formula, Atom, Radical


class TestGraph(unittest.TestCase):
	def setUp(self):
		self.g = graph.Graph()
		self.node1 = graph.Node({'id': 'A'})
		self.node2 = graph.Node({'id': 'B', 'color': 'blue'})
		self.node3 = graph.Node({'id': 'C', 'color': 'red'})
		self.node4 = graph.Node({'id': 'D', 'original_id': ['A', 'B'], 'type': 'rxn'})
		self.node5 = graph.Node({'id': 'E', 'original_id': 'cpd_E', 'type': 'cpd'})
		self.edge1_2 = graph.Edge(self.node1, self.node2)
		self.edge2_3 = graph.Edge(self.node2, self.node3, props={'id': '2_3'})

	def test_add_node(self):
		self.g.add_node(self.node1)
		self.assertTrue(self.node1 in self.g.nodes)

	def test_add_node_with_type_cpd(self):
		self.g.add_node(self.node5)
		self.assertEqual(self.g._nodes_original_id['cpd_E'], [self.node5])

	def test_add_node_with_type_rxn(self):
		self.g.add_node(self.node4)
		self.assertEqual(self.g._nodes_original_id['A,B'], [self.node4])

	def test_node_id_dict(self):
		self.g.add_node(self.node1)
		self.assertEqual(self.node1, self.g._nodes_id[self.node1.props['id']])

	def test_add_multiple_nodes(self):
		self.g.add_node(self.node1)
		self.g.add_node(self.node2)
		self.assertTrue(all(i in self.g.nodes for i in [self.node1, self.node2]))

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
		self.assertTrue(all(i in self.g.edges for i in [self.edge1_2, self.edge2_3]))

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
		self.assertTrue(edges_for_list == [(self.node2, set([self.edge1_2]))])
		self.assertEqual('A', 'A')

	def test_get_node(self):
		self.g.add_node(self.node1)
		test_node = self.g.get_node('A')
		self.assertEqual(test_node, self.node1)

	def test_nodes_id_dict(self):
		self.g.add_node(self.node1)
		self.g.add_node(self.node2)
		self.assertEqual(self.g.nodes_id_dict, {'A': self.node1, 'B': self.node2})

	def test_write_cytoscap_nodes(self):
		self.g.add_node(self.node2)
		self.g.add_node(self.node3)
		self.g.add_edge(self.edge2_3)
		path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_node')
		with open(path, mode='w') as f:
			self.g.write_cytoscape_nodes(f)
		read_file = []
		f = open(path, 'r')
		for i in f.readlines():
			read_file.append(i)
		f.close()
		self.assertEqual(read_file, ['id\tcolor\tlabel\n', 'B\tblue\tB\n', 'C\tred\tC\n'])

	def test_write_cytoscape_edges(self):
		self.g.add_node(self.node2)
		self.g.add_node(self.node3)
		self.g.add_edge(self.edge2_3)
		path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
		with open(path, mode='w') as f:
			self.g.write_cytoscape_edges(f)
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
		self.assertEqual(read_file, ['digraph {\n', 'size = "2, 2"; ratio = fill; node[fontname=Arail, fontsize=12]\n', ' "B"[color="blue",id="B"]\n', ' "C"[color="red",id="C"]\n', ' "B" -> "C"[id="2_3"]\n', '}\n'])

	def test_write_graphviz_compartmentalize(self):
		self.node2.props['compartment'] = 'c'
		self.node3.props['compartment'] = 'e'
		self.edge2_3.props['compartment'] = 'e'
		self.g.add_node(self.node2)
		self.g.add_node(self.node3)
		self.g.add_edge(self.edge2_3)
		path = os.path.join(tempfile.mkdtemp(), 'tmp_cyto_edge')
		with open(path, mode='w') as f:
			self.g.write_graphviz_compartmentalized(f, {'e': {'c'}, 'c': set()}, 'e', 2, 2)
		read_file = []
		f = open(path, 'r')
		for i in f.readlines():
			read_file.append(i)
		f.close()
		self.assertEqual(read_file, ['digraph {\n', 'size="2,2"; ratio = fill; node[fontname=Arail, fontsize=12]\n', ' subgraph cluster_e {\n', '  style=solid;\n', '  color=black;\n', '  penwidth=4;\n', '  fontsize=35;\n', '  label = "Compartment: e"\n', ' "C"[color="red",compartment="e",id="C"]\n', ' subgraph cluster_c {\n', '  style=dashed;\n', '  color=black;\n', '  penwidth=4;\n', '  fontsize=35;\n', '  label = "Compartment: c"\n', ' "B"[color="blue",compartment="c",id="B"]\n', '}} "B" -> "C"[compartment="e",id="2_3"]\n', '}\n'])


class TestNodes(unittest.TestCase):
	def setUp(self):
		self.g = graph.Graph()
		self.node1 = graph.Node({'id': 'A'})
		self.node2 = graph.Node({'id': 'B'})
		self.node3 = graph.Node({'id': 'C'})
		self.edge1_2 = graph.Edge(self.node1, self.node2)

	def test_eq_nodes(self):
		self.g.add_node(self.node1)
		self.assertEqual(self.node1 == self.g.get_node('A'), self.node1.props == self.g.get_node('A').props)

	def test_neq_nodes(self):
		self.assertTrue(self.node1 != self.node2)

	def test_neq_node_edge(self):
		self.assertFalse(self.node1 == self.edge1_2)


class TestEdges(unittest.TestCase):
	def setUp(self):
		self.g = graph.Graph()
		self.node1 = graph.Node({'id': 'A'})
		self.node2 = graph.Node({'id': 'B'})
		self.node3 = graph.Node({'id': 'C'})
		self.node4 = graph.Node({'id': 'A'})
		self.edge1_2 = graph.Edge(self.node1, self.node2)
		self.edge2_3 = graph.Edge(self.node2, self.node3)

	def test_eq_edges(self):
		self.assertTrue(self.node1 == self.node4)

	def test_neq_edges(self):
		self.assertTrue(self.edge1_2 != self.edge2_3)

	def test_neq_edge_node(self):
		self.assertFalse(self.edge1_2 == self.node1)


if __name__ == '__main__':
	unittest.main()
