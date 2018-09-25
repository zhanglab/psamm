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
# Copyright 2015-2018  Keith Dufault-Thompson <keitht547@uri.edu>

from __future__ import unicode_literals

from itertools import count

from six import iteritems, text_type

from collections import defaultdict


def _graphviz_prop_string(d):
    return ','.join('{}="{}"'.format(k, text_type(v)) for k, v
                    in iteritems(d))


class Entity(object):
    """Base class for graph entities."""
    def __init__(self, props={}):
        self._props = dict(props)

    @property
    def props(self):
        return self._props

    def __repr__(self):
        return '<{} id={}>'.format(
            self.__class__.__name__, self._props.get('id'))


class Graph(Entity):
    """Graph entity representing a collection of nodes and edges."""
    def __init__(self, props={}):
        super(Graph, self).__init__(props)
        self._nodes = set()
        self._edges = {}
        self._node_edges = {}
        self._edge_count = 0
        self._default_node_props = {}
        self._default_edge_props = {}
        self._nodes_id = {}
        self._nodes_original_id = defaultdict(list)

    def add_node(self, node):
        """add node to a Graph entity.
        node: Node entity.
        """
        self._nodes.add(node)
        self._nodes_id[node.props['id']] = node
        if 'type' in node.props:
            self.set_original_id(node)
        else:
            self._nodes_original_id[node.props['id']].append(node)

    def set_original_id(self, node):
        if node.props['type'] == 'cpd':
            original_id_string = text_type(node.props['original_id'])
        else:
            original_id_string = ','.join(node.props['original_id'])
        self._nodes_original_id[original_id_string].append(node)

    def get_node(self, node_id):
        """get Node object.
        Args:
            node_id: text_type(compound object) or text_type(a string that
                contains single or multiple reaction IDs).
        """
        return self._nodes_id[node_id]

    def add_edge(self, edge):
        if edge.source not in self._nodes or edge.dest not in self._nodes:
            raise ValueError('Edge nodes not in graph')

        s_nodes = tuple([edge.source, edge.dest])
        self._edges.setdefault(s_nodes, set()).add(edge)
        self._edge_count += 1

        self._node_edges.setdefault(edge.source, {}).setdefault(
            edge.dest, set()).add(edge)
        self._node_edges.setdefault(edge.dest, {}).setdefault(
            edge.source, set()).add(edge)

    @property
    def nodes(self):
        return iter(self._nodes)

    @property
    def node_count(self):
        return len(self._nodes)

    @property
    def nodes_id_dict(self):
        return self._nodes_id

    @property
    def nodes_original_id_dict(self):
        return self._nodes_original_id

    @property
    def edges(self):
        for pair, edges in iteritems(self._edges):
            for edge in edges:
                yield edge

    @property
    def edge_count(self):
        return self._edge_count

    def edges_for(self, node):
        return iteritems(self._node_edges.get(node, {}))

    def apply_layout(self, layout):
        for node in self.nodes:
            x, y = layout[node]
            node.props['x'] = x
            node.props['y'] = y
            node.props['pos'] = '{},{}!'.format(x, y)

    @property
    def default_node_props(self):
        return self._default_node_props

    @property
    def default_edge_props(self):
        return self._default_edge_props

    def write_graphviz(self, f, width, height):
        """write the nodes and edges information into a dot file.
        Args:
            self: Graph entity, including nodes and edges entities.
            f: An empty file.
            width: Width of final metabolic map.
            height: Height of final metabolic map.
        """
        f.write('digraph {\n')
        if width is not None and height is not None:
            f.write('size = "{}, {}"; ratio = fill; node[fontname=Arail, '
                    'fontsize=12]\n'.format(width, height))
        else:
            f.write('node[fontname=Arail, fontsize=12]\n'.
                    format(width, height))

        if len(self._default_node_props) > 0:
            f.write(' node[{}];\n'.format(
                _graphviz_prop_string(self._default_node_props)))

        if len(self._default_edge_props) > 0:
            f.write(' edge[{}];\n'.format(
                _graphviz_prop_string(self._default_edge_props)))

        for k, v in iteritems(self.props):
            f.write(' {}="{}";\n'.format(k, v))

        next_id = count(0)
        for node in sorted(self.nodes, key=lambda k: k.props['id']):
            if 'id' not in node.props:
                node.props['id'] = 'n{}'.format(next(next_id))

            f.write(' "{}"[{}]\n'.format(
                node.props['id'], _graphviz_prop_string(node.props)))

        for edge in sorted(self.edges, key=lambda k: (k.source.props['id'], k.dest.props['id'], k.props.get('dir'))):
            f.write(' "{}" -> "{}"[{}]\n'.format(
                edge.source.props['id'], edge.dest.props['id'],
                _graphviz_prop_string(edge.props)))

        f.write('}\n')

    def write_graphviz_compartmentalized(self, f, compartment_tree,
                                         extracellular, width, height):
        """Function to write compartmentalized version of dot file
        for graph.
        Args:
            self: Graph entity.
            f: An empty file.
            compartment_tree: a defaultdict of set, each element
                represents a compartment and its adjacent compartments.
            extracellular: the extracellular compartment in the model
            width: Width of final metabolic map.
            height: Height of final metabolic map..
        """
        f.write('digraph {\n')
        if width is not None and height is not None:
            f.write('size="{},{}"; ratio = fill; node[fontname=Arail, '
                    'fontsize=12]\n'.format(width, height))
        else:
            f.write('node[fontname=Arail, fontsize=12]\n'.format(width, height))

        if len(self._default_node_props) > 0:
            f.write(' node[{}];\n'.format(
                _graphviz_prop_string(self._default_node_props)))

        if len(self._default_edge_props) > 0:
            f.write(' edge[{}];\n'.format(
                _graphviz_prop_string(self._default_edge_props)))

        for k, v in iteritems(self.props):
            f.write(' {}="{}";\n'.format(k, v))

        next_id = count(0)
        node_dicts = defaultdict(list)

        for node in self.nodes:
            node_dicts[node.props['compartment']].append(node)

        def edit_labels(string):
            return string.replace('-', '_')

        def write_node_props(f, node_list):
            for node in node_list:
                if 'id' not in node.props:
                    node.props['id'] = 'n{}'.format(next(next_id))
                f.write('  "{}"[{}]\n'.format(
                    node.props['id'], _graphviz_prop_string(node.props)))

        def dfs_recursive(graph, vertex, node_dict, extracellular, f, path=[]):
            path.append(vertex)
            if vertex == extracellular:
                f.write(''.join(
                    [' subgraph cluster_{} '.format(edit_labels(vertex)),
                     '{\n  style=dashed;\n  color=black;\n  penwidth=4;\n  '
                     'fontsize=35;\n', '  label = "Compartment: {}"\n'.format
                     (edit_labels(vertex))]))
                for x in sorted(node_dict[vertex], key=lambda k: k.props['id']):
                    f.write(' "{}"[{}]\n'.format(
                        x.props['id'], _graphviz_prop_string(x.props)))
            elif vertex != extracellular:
                f.write(''.join(
                    [' subgraph cluster_{} '.format(edit_labels(vertex)),
                     '{\n  style=dashed;\n  color=black;\n  penwidth=4;\n  '
                     'fontsize=35;\n', '  label = "Compartment: {}"\n'.format
                     (edit_labels(vertex))]))
                for x in sorted(node_dict[vertex], key=lambda k: k.props['id']):
                    f.write(' "{}"[{}]\n'.format(
                        x.props['id'], _graphviz_prop_string(x.props)))
            for neighbor in graph[vertex]:
                if neighbor not in path:
                    path = dfs_recursive(graph, neighbor, node_dict, path, f)
            f.write('}')
            return path

        dfs_recursive(compartment_tree, extracellular, node_dicts,
                      extracellular, f)

        for edge in sorted(self.edges,  key=lambda k: (k.source.props['id'], k.dest.props['id'], k.props.get('dir'))):
            f.write(' "{}" -> "{}"[{}]\n'.format(
                edge.source.props['id'], edge.dest.props['id'],
                _graphviz_prop_string(edge.props)))

        f.write('}\n')

    def write_cytoscape_nodes(self, f):
        """write a table file (.tsv) that contains nodes information.
        Args:
            self: Graph entity.
            f: An empty file.
        """

        next_id = count(0)
        properties = set()
        for node in self.nodes:
            if 'id' not in node.props:
                node.props['id'] = 'n{}'.format(next(next_id))
            if 'label' not in node.props:
                node.props['label'] = node.props['id']
            properties.update(node.props)

        properties.remove('id')
        properties.remove('label')
        properties.remove('original_id')
        properties = ['id'] + sorted(properties) + ['label']
        f.write('\t'.join(properties) + '\n')
        for node in self.nodes:
            a = '\t'.join(node.props.get(x)
                          for x in properties if x != 'label')
            b = node.props['label'].replace('\n', ',')
            f.write('{}\t{}\n'.format(a, b))

    def write_cytoscape_edges(self, f):
        """ Write a tab separated table that contains edges information,
        including edge source, edge dest, and edge properties.
        Args:
            self: Graph entity.
            f: An empty TSV file.
        """
        properties = set()
        for edge in self.edges:
            properties.update(edge.props)

        properties = sorted(properties)
        header = ['source', 'target'] + properties
        f.write('\t'.join(header) + '\n')
        for edge in self.edges:
            f.write('{}\t{}\t{}\n'.format(
                edge.source.props['id'], edge.dest.props['id'],
                '\t'.join(text_type(edge.props.get(x.encode('utf8')))
                          for x in properties)))


class Node(Entity):
    """Node entity represents a vertex in the graph."""
    __hash__ = Entity.__hash__

    def __init__(self, props={}):
        super(Node, self).__init__(props)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return other._props == self._props
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)


class Edge(Entity):
    """Edge entity represents a connection between nodes."""
    __hash__ = Entity.__hash__

    def __init__(self, source, dest, props={}):
        super(Edge, self).__init__(props)
        self.source = source
        self.dest = dest

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            a = self._props == other._props
            b = self.source == other.source
            c = self.dest == other.dest
            if all(i is True for i in [a, b, c]):
                return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
