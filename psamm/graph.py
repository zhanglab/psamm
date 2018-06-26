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

from itertools import count

from six import iteritems, text_type


def _graphviz_prop_string(d):
     return ','.join('{}="{}"'.format(k, v) for k, v in iteritems(d))


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

    def add_node(self, node):
        self._nodes.add(node)

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

    def write_graphviz(self, f):
        f.write('digraph {\n')

        if len(self._default_node_props) > 0:
            f.write(' node[{}];\n'.format(
                _graphviz_prop_string(self._default_node_props)))

        if len(self._default_edge_props) > 0:
            f.write(' edge[{}];\n'.format(
                _graphviz_prop_string(self._default_edge_props)))

        for k, v in iteritems(self.props):
            f.write(' {}="{}";\n'.format(k, v))

        next_id = count(0)
        for node in self.nodes:
            if 'id' not in node.props:
                node.props['id'] = 'n{}'.format(next(next_id))

            f.write(' "{}"[{}]\n'.format(
                node.props['id'], _graphviz_prop_string(node.props)))

        for edge in self.edges:
            f.write(' "{}" -> "{}"[{}]\n'.format(
                edge.source.props['id'], edge.dest.props['id'],
                _graphviz_prop_string(edge.props)))

        f.write('}\n')

    def write_svg(self, f):
        min_x = max_x = min_y = max_y = None
        for node in self.nodes:
            x = node.props.get('x', 0)
            y = node.props.get('y', 0)
            if min_x is None or x < min_x:
                min_x = x
            if max_x is None or x > max_x:
                max_x = x
            if min_y is None or y < min_y:
                min_y = y
            if max_y is None or y > max_y:
                max_y = y

        width = max_x - min_x + 1
        height = max_y - min_y + 1

        f.write('\n'.join([
            '<?xml version="1.0"?>',
            '<svg xmlns="http://www.w3.org/2000/svg"',
            ' viewBox="{} {} {} {}"'.format(min_x, min_y, width, height),
            ' width="{}cm" height="{}cm">'.format(width, height),
            '']))

        for edge in self.edges:
            x1 = edge.source.props.get('x', 0) + 0.5
            y1 = edge.source.props.get('y', 0) + 0.5
            x2 = edge.dest.props.get('x', 0) + 0.5
            y2 = edge.dest.props.get('y', 0) + 0.5
            f.write('\n'.join([
                '<g>',
                '<title>{}</title>'.format(edge.props.get('id', '')),
                '<path d="M{} {}L{} {}"'.format(x1, y1, x2, y2),
                ' stroke="black" stroke-width="0.1"/>',
                '</g>',
                ''
            ]))

        for node in self.nodes:
            x = node.props.get('x', 0) + 0.25
            y = node.props.get('y', 0) + 0.25
            f.write('\n'.join([
                '<g>',
                '<title>{}</title>'.format(node.props.get('id', '')),
                '<rect x="{}" y="{}"'.format(x, y),
                ' width="0.5" height="0.5" fill="white" stroke="black"',
                ' stroke-width="0.1"/>',
                '</g>',
                '']))

        f.write('</svg>\n')

    def write_cytoscape_nodes(self, f):
        next_id = count(0)
        properties = set()
        for node in self.nodes:
            if 'id' not in node.props:
                node.props['id'] = 'n{}'.format(next(next_id))
            if 'label' not in node.props:
                node.props['label'] = node.props['id']
            properties.update(node.props)

        properties.remove('id')
        properties = ['id'] + sorted(properties)
        f.write('\t'.join(properties) + '\n')
        for node in self.nodes:
            f.write('\t'.join(
                text_type(node.props.get(x)) for x in properties) + '\n')

    def write_cytoscape_edges(self, f):
        properties = set()
        for edge in self.edges:
            properties.update(edge.props)

        properties = sorted(properties)
        header = ['source', 'target'] + properties
        f.write('\t'.join(header) + '\n')
        for edge in self.edges:
            f.write('{}\t{}\t{}\n'.format(
                edge.source.props['edge_id'], edge.dest.props['edge_id'],
                '\t'.join(text_type(edge.props.get(x)) for x in properties)))


class Node(Entity):
    """Node entity represents a vertex in the graph."""
    def __init__(self, props={}):
        super(Node, self).__init__(props)



class Edge(Entity):
    """Edge entity represents a connection between nodes."""
    def __init__(self, source, dest, props={}):
        super(Edge, self).__init__(props)
        self.source = source
        self.dest = dest
