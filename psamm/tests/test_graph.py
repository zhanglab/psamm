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

import unittest
from psamm import graph
from psamm.datasource.native import NativeModel
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm import formula

from psamm.datasource.native import NativeModel
from psamm.datasource.native import ReactionEntry
from psamm.datasource.native import CompoundEntry
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm.formula import Formula, Atom, Radical


class TestFilterDict(unittest.TestCase):
	def setUp(self):
		native_model = NativeModel()
		native_model.reactions.add_entry(ReactionEntry({'id': 'rxn1',
				                                                'equation': parse_reaction('fum_c[c] + h2o_c[c] <=> mal_L_c[c]')}))
		native_model.compounds.add_entry(CompoundEntry({'id': 'fum_c[c]', 'formula': parse_compound('C4H2O4', 'c')}))
		native_model.compounds.add_entry(CompoundEntry({'id': 'h2o_c[c]', 'formula': parse_compound('H2O', 'c')}))
		native_model.compounds.add_entry(CompoundEntry({'id': 'mal_L_c[c]', 'formula': parse_compound('C4H4O5', 'c')}))
		self.native = native_model
		self.mm = native_model.create_metabolic_model()

		self.compounds_dict = {'fum_c' : Formula.parse('C4H2O4'), 'h2o_c' : Formula.parse('H2O'), 'mal_L_c': Formula.parse('C4H4O5')}
		self.element = 'C'
		self.method = 'fpp'

	def test(self):
		print(self.compounds_dict)
		# vis.make_filter_dict(self.native, self.mm, self.method, self.element, self.compounds_dict, None, [])

		self.assertEqual('A', 'A')


class TestMakeModelExample(unittest.TestCase):
	def setUp(self):
		native_model = NativeModel()
		native_model.reactions.add_entry(ReactionEntry({'id': 'rxn1',
		                                                'equation': parse_reaction('A[c] + B[c] => C[c] + D[c]')}))
		self._native = native_model
		self._mm = native_model.create_metabolic_model()

	def test_nm_test(self):
		print(self._native)
		self.assertEqual('A', 'A')


class TestGraph(unittest.TestCase):
	def setUp(self):
		self.g = graph.Graph()
		self.node1 = graph.Node({'id': 'A'})
		self.node2 = graph.Node({'id': 'B'})
		self.node3 = graph.Node({'id': 'C'})
		self.edge1_2 = graph.Edge(self.node1, self.node2)
		self.edge2_3 = graph.Edge(self.node2, self.node3, props={'id': '2_3'})

	def test_add_node(self):
		self.g.add_node(self.node1)
		self.assertTrue(self.node1 in self.g.nodes)

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


if __name__ == '__main__':
	unittest.main()
