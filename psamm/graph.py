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

from itertools import count
from six import iteritems, text_type
from collections import defaultdict, Counter

from . import findprimarypairs

import logging

from .formula import Formula, Atom, ParseError
from .reaction import Direction, Reaction


logger = logging.getLogger(__name__)


def _graphviz_prop_string(d):
    return ','.join('{}="{}"'.format(k, text_type(v)) for k, v
                    in sorted(iteritems(d)))


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
            original_id_string = ','.join([c.id for c in node.props['entry']])
        else:
            if node.props['type'] == 'Ex_rxn':
                original_id_string = node.props['id']
            else:
                original_id_string = ','.join(
                    [r.id for r in node.props['entry']])
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

    @property
    def default_node_props(self):
        return self._default_node_props

    @property
    def default_edge_props(self):
        return self._default_edge_props

    def write_graphviz(self, f, width, height):
        """ Write the nodes and edges information into a dot file.

        Print a given graph object to a graph file in the dot format.
        This graph can then be converted to an image file using the
        graphviz program.

        Args:
            self: Graph entity, including nodes and edges entities.
            f: An empty file.
            width: Width of final metabolic map.
            height: Height of final metabolic map.
        """
        f.write('digraph {\n')
        if width is not None and height is not None:
            f.write('size = "{}, {}"; ratio = fill;\n'.format(width, height))

        if len(self._default_node_props) > 0:
            f.write(' node[{}];\n'.format(
                    _graphviz_prop_string(self._default_node_props)))

        if len(self._default_edge_props) > 0:
            f.write(' edge[{}];\n'.format(
                    _graphviz_prop_string(self._default_edge_props)))

        for k, v in iteritems(self.props):
            f.write(' {}="{}";\n'.format(k, v))

        for node in sorted(self.nodes, key=lambda k: k.props['id']):
            f.write(' "{}"[{}]\n'.format(node.props['id'],
                                         _graphviz_prop_string(node.props)))

        for edge in sorted(self.edges, key=lambda k: (
                k.source.props['id'], k.dest.props['id'], k.props.get('dir'))):
            f.write(' "{}" -> "{}"[{}]\n'.format(
                    edge.source.props['id'], edge.dest.props['id'],
                    _graphviz_prop_string(edge.props)))

        f.write('}\n')

    def write_graphviz_compartmentalized(self, f, compartment_tree,
                                         extracellular, width, height):
        """Function to write compartmentalized version of dot file
        for graph.

        Print a given graph object to a graph file in the dot format.
        In this graph, reaction nodes will be separated into different
        areas in the visualization based on the defined cellular
        compartments in the GEM.

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
            f.write('size="{},{}"; ratio = fill;\n'.format(width, height))

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

        def dfs_recursive(graph, vertex, node_dict,
                          extracellular, f, path=[]):
            path.append(vertex)
            if vertex == extracellular:
                f.write(''.join(
                        [' subgraph cluster_{} '.format(edit_labels(vertex)),
                         '{\n  style=solid;\n  color=black;\n  penwidth=4;\n  '
                         'fontsize=35;\n',
                         '  label = "Compartment: {}"\n'.format
                         (edit_labels(vertex))]))
                for x in sorted(node_dict[vertex],
                                key=lambda k: k.props['id']):
                    f.write(' "{}"[{}]\n'.format(
                            x.props['id'], _graphviz_prop_string(x.props)))
            elif vertex != extracellular:
                f.write(''.join(
                        [' subgraph cluster_{} '.format(edit_labels(vertex)),
                         '{\n  style=dashed;\n  color=black;\n  '
                         'penwidth=4;\n  fontsize=35;\n',
                         '  label = "Compartment: {}"\n'.format
                         (edit_labels(vertex))]))
                for x in sorted(node_dict[vertex],
                                key=lambda k: k.props['id']):
                    f.write(' "{}"[{}]\n'.format(
                            x.props['id'], _graphviz_prop_string(x.props)))
            for neighbor in graph[vertex]:
                if neighbor not in path:
                    path = dfs_recursive(graph, neighbor, node_dict, path, f)
            f.write('}')
            return path

        dfs_recursive(compartment_tree, extracellular, node_dicts,
                      extracellular, f)

        for edge in sorted(self.edges, key=lambda k: (k.source.props['id'],
                                                      k.dest.props['id'],
                                                      k.props.get('dir'))):
            f.write(' "{}" -> "{}"[{}]\n'.format(
                    edge.source.props['id'], edge.dest.props['id'],
                    _graphviz_prop_string(edge.props)))

        f.write('}\n')

    def write_nodes_tables(self, f):
        """write a table file (.tsv) that contains nodes information.

        Write all node information from a given graph object into a
        tab separated table. This table will include IDs, shapes,
        fill colors, and labels. This table can be used to import
        the graph information to other programs.

        Args:
            self: Graph entity.
            f: An empty file.
        """

        properties = set()
        for node in self.nodes:
            if 'label' not in node.props:
                node.props['label'] = node.props['id']
            properties.update(node.props)
        if 'id' in properties:
            properties.remove('id')
        if 'label' in properties:
            properties.remove('label')
        if 'original_id' in properties:
            properties.remove('original_id')
        if 'entry' in properties:
            properties.remove('entry')
        properties = ['id'] + sorted(properties) + ['label']
        f.write('\t'.join(properties) + '\n')

        for node in sorted(self.nodes, key=lambda k: k.props['id']):
            # a = '\t'.join(text_type(node.props.get(x))
            #               for x in properties if x != 'label')
            a = '\t'.join(text_type(node.props.get(x))
                          for x in properties if x not in ['label', 'entry'])
            b = node.props['label'].replace('\n', ',')
            f.write('{}\t{}\n'.format(a, b))

    def write_edges_tables(self, f):
        """ Write a tab separated table that contains edges information,
        including edge source, edge dest, and edge properties.

        Write all edge information from a given graph object into a
        tab separated table. This table will include IDs, source nodes
        and destination nodes. This table can be used to import
        the graph information to other programs.

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
        for edge in sorted(self.edges, key=lambda k: (edge.source.props['id'],
                                                      edge.dest.props['id'],
                                                      edge.props.get('dir'))):
            f.write('{}\t{}\t{}\n'.format(
                    edge.source.props['id'], edge.dest.props['id'],
                    '\t'.join(
                        text_type(edge.props.get(x)) for x in properties)))


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


def get_compound_dict(model):
    """ Parse through model compounds and return a formula dict.

    This function will parse the compound information
    in a given model and return a dictionary of compound IDs to
    compound formula objects.

    Args:
        model: <class 'psamm.datasource.native.NativeModel'>
    """
    compound_formula = {}
    for compound in model.compounds:
        if compound.formula is not None:
            try:
                f = Formula.parse(compound.formula).flattened()
                if not f.is_variable():
                    compound_formula[compound.id] = f
                else:
                    logger.warning(
                        'Skipping variable formula {}: {}'.format(
                            compound.id, compound.formula))
            except ParseError as e:
                msg = (
                    'Error parsing formula'
                    ' for compound {}:\n{}\n{}'.format(
                        compound.id, e, compound.formula))
                if e.indicator is not None:
                    msg += '\n{}'.format(e.indicator)
                logger.warning(msg)
    return compound_formula


def make_network_dict(nm, mm, subset=None, method='fpp', element=None,
                      excluded_reactions=[], reaction_dict={}, analysis=None):
    """Create a dictionary of reactant/product pairs to reaction directions.

    Returns a dictionary that connects predicted reactant/product pairs
    to their original reaction directions. This dictionary is used when
    generating bipartite graph objects. This can be done either using
    the FindPrimaryPairs method to predict reactant/product pairs or
    by mapping all possible pairs.

    Args:
        nm: <class 'psamm.datasource.native.NativeModel'>, the native model
        mm: <class 'psamm.metabolicmodel.MetabolicModel'>, the metabolic model.
        subset: None or path to a file that contains a list of reactions or
            compound ids. By default, it is None. It defines which reaction
            need to be visualized.
        method: 'fpp' or 'no-fpp', the method used for visualization.
        element: Symbol of chemical atom, such as 'C' ('C' indicates carbon).
        excluded_reactions: a list that contains reactions excluded from
            visualization.
        reaction_dict: dictionary of FBA or FVA results. By default it is an
            empty dictionary.
        analysis: "None" type or a string indicates if FBA or FVA file is
            given in command line.
    """
    compound_formula = get_compound_dict(nm)
    if not compound_formula and (method == 'fpp' or element):
        logger.error(
            'Compound formulas are required for fpp or specific element '
            'visualizations, try --element all to visualize all pathways '
            'without compound formula input.')
        exit(1)

    if subset is not None:
        testing_list_raw = []
        for rxn in mm.reactions:
            if rxn in nm.reactions or mm.is_exchange(rxn):
                if rxn in subset:
                    if rxn not in excluded_reactions:
                        testing_list_raw.append(rxn)
                    else:
                        logger.warning(
                            'Reaction {} is in the subset and exclude file. '
                            'Reaction will be excluded.'.format(rxn))
    else:
        testing_list_raw = [rxn for rxn in mm.reactions
                            if (rxn in nm.reactions or mm.is_exchange(rxn)) and
                            rxn not in excluded_reactions]

    reaction_data = {}
    style_flux_dict = {}

    for rxn in testing_list_raw:
        flux = 0
        style = 'solid' if analysis is None else 'dotted'
        # direction = nm.reactions[rxn].equation.direction
        direction = mm.get_reaction(rxn).direction
        if rxn in reaction_dict:
            style = 'solid'
            reaction = mm.get_reaction(rxn)
            lower = reaction_dict[rxn][0]
            upper = reaction_dict[rxn][1]
            flux = lower * upper

            if flux == 0:
                if analysis == 'fba' or (lower == 0 and upper == 0):
                    style = 'dotted'
                elif lower == 0:
                    direction = Direction.Forward
                else:
                    direction = Direction.Reverse
            elif flux > 0:
                direction = Direction.Forward if (
                    analysis == 'fba' or lower > 0) else Direction.Reverse
            else:
                direction = Direction.Reverse if analysis == 'fba' \
                    else Direction.Both

            r = Reaction(direction, reaction.left, reaction.right)
            r_id = rxn
            mm.remove_reaction(rxn)
            mm.database.set_reaction(r_id, r)
            mm.add_reaction(r_id)
            if rxn in nm.reactions:
                nm.reactions[rxn].equation = r
        if rxn in nm.reactions:
            reaction_data[rxn] = (nm.reactions[rxn], direction)
        style_flux_dict[rxn] = (style, abs(flux))

    flux_list = sorted([f for s, f in style_flux_dict.values()])
    median = 1
    flux_list = list(filter(lambda x: x != 0, flux_list))

    testing_list = [rxn for rxn in testing_list_raw if not mm.is_exchange(rxn)]

    if flux_list:
        mid = len(flux_list) // 2
        median = flux_list[mid] if len(flux_list) % 2 else float(
            flux_list[mid-1] + flux_list[mid]) / 2.0
        if median < 1:
            median = 1

    for rxn, (style, flux) in style_flux_dict.items():
        new_flux = float(5 * flux) / float(median)
        new_flux = max(min(10, new_flux), 1)
        style_flux_dict.update({rxn: (style, new_flux)})

    full_pairs_dict = {}
    if method == 'fpp':
        element_weight = findprimarypairs.element_weight
        testing_list_update = []
        for r in testing_list:
            if all(cpd.name in compound_formula
                   for cpd, _ in nm.reactions[r].equation.compounds):
                testing_list_update.append(r)
            else:
                logger.warning(
                    'Reaction {} is excluded from visualization due to '
                    'missing or undefined compound formula'.format(r))
        reaction_pairs = [(r, nm.reactions[r].equation)
                          for r in testing_list_update]
        fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(
            reaction_pairs, compound_formula, element_weight=element_weight)

        for rxn_id, fpp_pairs in iteritems(fpp_dict):
            compound_pairs = []
            for cpd_pair, transfer in iteritems(fpp_pairs[0]):
                if element is None:
                    compound_pairs.append(cpd_pair)
                else:
                    if any(Atom(element) in k for k in transfer):
                        compound_pairs.append(cpd_pair)
            (rxn_entry, rxn_dir) = reaction_data[rxn_id]
            full_pairs_dict[rxn_entry] = (sorted(compound_pairs), rxn_dir)

    elif method == 'no-fpp':
        for reaction in testing_list:
            compound_pairs = []
            for substrate in nm.reactions[reaction].equation.left:
                for product in nm.reactions[reaction].equation.right:
                    if element is None:
                        compound_pairs.append((substrate[0], product[0]))
                    else:
                        if Atom(element) in \
                                compound_formula[substrate[0].name] and \
                                Atom(element) in \
                                compound_formula[product[0].name]:
                            compound_pairs.append((substrate[0], product[0]))
            full_pairs_dict[nm.reactions[reaction]] = \
                (sorted(compound_pairs),
                 nm.reactions[reaction].equation.direction)

    return full_pairs_dict, style_flux_dict


def write_network_dict(network_dict):
    """Print out network dictionary object to a tab separated table.

    Print information from the dictionary "network_dict". This information
    includes four items: reaction ID, reactant, product, and direction.
    this can table can be used as an input to other graph visualization
    and analysis software.


     Args:
         network_dict: Dictionary object from make_network_dict()
     """
    for key, value in sorted(iteritems(network_dict), key=lambda x: str(x)):
        (cpair_list, dir) = value
        for (c1, c2) in cpair_list:
            print('{}\t{}\t{}\t{}'.format(key.id, c1, c2, dir_value(dir)))


def make_cpair_dict(filter_dict, args_method, args_combine, style_flux_dict,
                    hide_edges=[]):
    """Create a mapping from compound pair to a defaultdict containing
    lists of reactions for the forward, reverse, and both directions.

    Returns a dictionary that connects reactant/product pair to all reactions
    that contain this pair. Those reactions are stored in a dictionary and
    classified by reaction direction. For example:
    {(c1, c2): {'forward': [rxn1], 'back': [rxn3], 'both': [rxn4, rxn5]},
    (c3, c4): {...}, ...}

    Args:
        filter_dict: A dictionary mapping reaction entry to
            compound pairs (inside of the pairs there are cpd
            objects, not cpd IDs)
        args_method: a string, including 'fpp' and 'no-fpp'.
        args_combine: combine level, default = 0, optional: 1 and 2.
        style_flux_dict: a dictionary to set the edge style when fba or
            fva input is given.
        hide_edges: to determine edges between which compound pair need
            to be hidden.
    """

    new_id_mapping = {}
    new_style_flux_dict = {}
    rxn_count = Counter()
    cpair_dict = defaultdict(lambda: defaultdict(list))

    def make_mature_cpair_dict(cpair_dict, hide):
        new_cpair_dict = {}
        cpair_list = []
        for (c1, c2), rxns in sorted(iteritems(cpair_dict)):
            if (c1, c2) not in cpair_list and (text_type(c1),
                                               text_type(c2)) not in hide:
                new_rxns = rxns
                if (c2, c1) in cpair_dict:
                    if len(cpair_dict[(c2, c1)]['forward']) > 0:
                        for r in cpair_dict[(c2, c1)]['forward']:
                            new_rxns['back'].append(r)
                    if len(cpair_dict[(c2, c1)]['back']) > 0:
                        for r in cpair_dict[(c2, c1)]['back']:
                            new_rxns['forward'].append(r)
                    if len(cpair_dict[(c2, c1)]['both']) > 0:
                        for r in cpair_dict[(c2, c1)]['both']:
                            new_rxns['both'].append(r)
                    new_cpair_dict[(c1, c2)] = new_rxns
                    cpair_list.append((c1, c2))
                    cpair_list.append((c2, c1))
                else:
                    new_cpair_dict[(c1, c2)] = new_rxns
                    cpair_list.append((c1, c2))

        rxns_sorted_cpair_dict = defaultdict(lambda: defaultdict(list))
        for (c1, c2), rxns in sorted(iteritems(new_cpair_dict)):
            for direction, rlist in iteritems(rxns):
                rxns_sorted_cpair_dict[(c1, c2)][direction] = sorted(rlist)

        return rxns_sorted_cpair_dict

    if args_method != 'no-fpp':
        if args_combine == 0:
            for rxn, cpairs_dir in iteritems(filter_dict):
                have_visited = set()
                sub_pro = defaultdict(list)
                rxn_mixcpairs = defaultdict(list)
                for (c1, c2) in sorted(cpairs_dir[0]):
                    sub_pro[c1].append(c2)
                for k1, v1 in sorted(iteritems(sub_pro)):
                    if k1 not in have_visited:
                        rxn_count[rxn] += 1
                        have_visited.add(k1)
                        r_id = '{}_{}'.format(rxn.id, rxn_count[rxn])
                        new_id_mapping[r_id] = rxn
                        new_style_flux_dict[r_id] = style_flux_dict[rxn.id]
                        for v in v1:
                            rxn_mixcpairs[r_id].append((k1, v))
                        for k2, v2 in iteritems(sub_pro):
                            if k2 not in have_visited:
                                if k2 != k1:
                                    if v1 == v2:
                                        have_visited.add(k2)
                                        for vtest in v2:
                                            rxn_mixcpairs[r_id].append(
                                                (k2, vtest))
                for rxn_id, cpairs in iteritems(rxn_mixcpairs):
                    for (c1, c2) in cpairs:
                        if cpairs_dir[1] == Direction.Forward:
                            cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                        elif cpairs_dir[1] == Direction.Reverse:
                            cpair_dict[(c1, c2)]['back'].append(rxn_id)
                        else:
                            cpair_dict[(c1, c2)]['both'].append(rxn_id)

        elif args_combine == 1 or args_combine == 2:
            for rxn, cpairs_dir in iteritems(filter_dict):
                cpd_rid = {}
                have_visited = set()
                for (c1, c2) in sorted(cpairs_dir[0]):
                    if c1 not in have_visited:
                        if c2 not in have_visited:
                            rxn_count[rxn] += 1
                            rxn_id = '{}_{}'.format(
                                rxn.id, rxn_count[rxn])
                            new_id_mapping[rxn_id] = rxn
                            new_style_flux_dict[rxn_id] = \
                                style_flux_dict[rxn.id]
                            have_visited.add(c1)
                            have_visited.add(c2)
                            cpd_rid[c1] = rxn_id
                            cpd_rid[c2] = rxn_id
                        else:
                            rxn_id = cpd_rid[c2]
                            have_visited.add(c1)
                            cpd_rid[c1] = rxn_id
                    else:
                        rxn_id = cpd_rid[c1]
                        have_visited.add(c2)
                        cpd_rid[c2] = rxn_id

                    if cpairs_dir[1] == Direction.Forward:
                        cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                    elif cpairs_dir[1] == Direction.Reverse:
                        cpair_dict[(c1, c2)]['back'].append(rxn_id)
                    else:
                        cpair_dict[(c1, c2)]['both'].append(rxn_id)
    else:
        for rxn, cpairs_dir in iteritems(filter_dict):
            for (c1, c2) in cpairs_dir[0]:
                r_id = rxn.id
                new_id_mapping[r_id] = rxn
                new_style_flux_dict[r_id] = style_flux_dict[rxn.id]
                if cpairs_dir[1] == Direction.Forward:
                    cpair_dict[(c1, c2)]['forward'].append(r_id)
                elif cpairs_dir[1] == Direction.Reverse:
                    cpair_dict[(c1, c2)]['back'].append(r_id)
                else:
                    cpair_dict[(c1, c2)]['both'].append(r_id)

    rxns_sorted_cpair_dict = make_mature_cpair_dict(cpair_dict, hide_edges)
    return rxns_sorted_cpair_dict, new_id_mapping, new_style_flux_dict


def make_bipartite_graph_object(cpairs_dict, new_id_mapping, method,
                                args_combine, model_compound_entries,
                                new_style_flux_dict, analysis=None):
    """ Makes a bipartite graph object from a cpair_dict object.

    Start from empty graph() and cpair dict to make a graph object.
    Nodes only have rxn/cpd ID and rxn/cpd entry info in this initial
    graph. The information can be modified to add on other properties
    like color, or names.

    Args:
        cpairs_dict: defaultdict of compound_pair:
            defaultdict of direction: reaction list.
            e.g. {(c1, c2): {'forward":[rx1],
            'both':[rx2}}.
        new_id_mapping: dictionary of rxn_id_suffix: rxn_id.
        method: options=['fpp', 'no-fpp', file_path].
        args_combine: Command line argument, could be 0, 1, 2.
        model_compound_entries: dict of cpd_id:compound_entry.
        new_style_flux_dict: a dictionary to determine the edge style with
            the new reaction IDs.
    return: A graph object that contains basic nodes and edges.
        only ID and rxn/cpd entry are in node properties,
        no features like color, shape.
    """
    g = Graph()
    g._default_node_props['fontname'] = 'Arial'
    g._default_node_props['fontsize'] = 12

    def add_graph_nodes(g, cpairs_dict, method, new_id_mapping, args_combine,
                        model_compound_entries):
        """ Create compound and reaction nodes, adding them to
            empty graph object.

        Args:
            g: An empty Graph object.
            cpairs_dict: defaultdict of compound_pair:
                defaultdict of direction:
                reaction list. e.g. {(c1, c2): {'forward":[rx1],
                    'both':[rx2}}.
            method: Command line argument,
                options=['fpp', 'no-fpp', file_path].
            new_id_mapping: dictionary of rxn_id_suffix: rxn_id.
            args_combine: Command line argument,
                could be 0, 1, 2. By default it is 0.
            model_compound_entries: cpd id map to cpd entry of native model.
        return: A graph object that contains a set of nodes.
        """
        graph_nodes = set()
        for cpair, reactions in sorted(iteritems(cpairs_dict)):
            for c in cpair:
                if c not in graph_nodes:
                    node = Node({
                        'id': text_type(c),
                        'entry': [model_compound_entries[c.name]],
                        'compartment': c.compartment,
                        'type': 'cpd'})
                    g.add_node(node)
                    graph_nodes.add(c)
            for direction, rlist in iteritems(reactions):
                if method != 'no-fpp' and args_combine == 2:
                    real_rxns = [new_id_mapping[r] for r in rlist]
                    rxn_string = text_type(','.join(rlist))
                    if rxn_string not in graph_nodes:
                        rnode = Node({
                            'id': text_type(','.join(rlist)),
                            'entry': real_rxns,
                            'compartment': c.compartment,
                            'type': 'rxn'})
                        g.add_node(rnode)
                        graph_nodes.add(rxn_string)
                else:
                    for sub_rxn in rlist:
                        rxn_ob = new_id_mapping[sub_rxn]
                        if sub_rxn not in graph_nodes:
                            rnode = Node({
                                'id': text_type(sub_rxn),
                                'entry': [rxn_ob],
                                'compartment': c.compartment,
                                'type': 'rxn'})
                            g.add_node(rnode)
                            graph_nodes.add(sub_rxn)
        return g

    def add_edges(g, cpairs_dict, method, args_combine, new_style_flux_dict,
                  analysis):
        """ Add edges to the graph object obtained in last step.

        Args:
            g: A graph object contains a set of nodes.
            cpairs_dict: A defaultdict of compound_pair: defaultdict
                of direction: reaction list. e.g. {(c1, c2):
                {'forward":[rx1], 'both':[rx2}}.
            method: Command line argument, options=
                ['fpp', 'no-fpp', file_path].
            # split: Command line argument, True or False.
                By default split = False.
            args_combine: Command line argument, an
                integer(could be 0, 1 or 2).
            new_style_flux_dict: A dictsionary of reaction is maps to edge
                style and edge width.
            analysis: "None" type or a string indicates if FBA or FVA file is
                given in command line.
        """
        edge_list = []
        for (c1, c2), value in iteritems(cpairs_dict):
            for direction, rlist in iteritems(value):
                new_rlist = ','.join(rlist)
                if args_combine == 0 or args_combine == 1 or \
                        method == 'no-fpp':
                    for sub_rxn in rlist:
                        test1 = c1, sub_rxn
                        if test1 not in edge_list:
                            edge_list.append(test1)
                            g.add_edge(Edge(
                                g.get_node(text_type(c1)),
                                g.get_node(text_type(sub_rxn)),
                                {'dir': direction,
                                 'style': new_style_flux_dict[sub_rxn][0],
                                 'penwidth': new_style_flux_dict[sub_rxn][1]}))

                        test2 = c2, sub_rxn
                        if test2 not in edge_list:
                            edge_list.append(test2)
                            g.add_edge(Edge(
                                g.get_node(text_type(sub_rxn)),
                                g.get_node(text_type(c2)),
                                {'dir': direction,
                                 'style': new_style_flux_dict[sub_rxn][0],
                                 'penwidth': new_style_flux_dict[sub_rxn][1]}))

                else:
                    test1 = c1, new_rlist
                    test2 = c2, new_rlist
                    sub_rxn = list(new_rlist.split(','))

                    style_list = set(new_style_flux_dict[rxn][0] for rxn in
                                     sub_rxn)
                    style = style_list.pop() if len(set(style_list)) == 1 \
                        else 'solid'

                    flux = sum([new_style_flux_dict[rxn][1] for rxn in
                                sub_rxn], 0) if analysis else 0
                    flux = max(min(10, flux), 1)

                    if test1 not in edge_list:
                        edge_list.append(test1)
                        g.add_edge(Edge(
                            g.get_node(text_type(c1)),
                            g.get_node(text_type(new_rlist)),
                            {'dir': direction, 'style': style,
                             'penwidth': flux}))

                    if test2 not in edge_list:
                        edge_list.append(test2)
                        g.add_edge(Edge(
                            g.get_node(text_type(new_rlist)),
                            g.get_node(text_type(c2)),
                            {'dir': direction, 'style': style,
                             'penwidth': flux}))
        return g

    g = add_graph_nodes(g, cpairs_dict, method, new_id_mapping,
                        args_combine, model_compound_entries)
    g = add_edges(g, cpairs_dict, method, args_combine, new_style_flux_dict,
                  analysis)
    return g


def make_compound_graph(network_dictionary):
    g = Graph()
    compound_nodes = []
    edge_list = []
    for reaction, (cpair_list, dir) in iteritems(network_dictionary):
        for (c1, c2) in cpair_list:
            if c1 not in compound_nodes:
                g.add_node(Node({'id': text_type(c1), 'entry': c1}))
                compound_nodes.append(c1)
            if c2 not in compound_nodes:
                g.add_node(Node({'id': text_type(c2), 'entry': c2}))
                compound_nodes.append(c2)
            cpair_sorted = sorted([c1.name, c2.name])
            edge = Edge(g.get_node(text_type(c1)),
                        g.get_node(text_type(c2)),
                        props={'id': '{}_{}_{}'.format(
                            cpair_sorted[0], cpair_sorted[1],
                            dir_value(dir)), 'dir': dir_value(dir)})
            if edge.props['id'] not in edge_list:
                g.add_edge(edge)
                edge_list.append(edge.props['id'])
    return g


def dir_value(direction):
    """ Assign value to different reaction directions"""
    if direction == Direction.Forward:
        return 'forward'
    elif direction == Direction.Reverse:
        return 'back'
    else:
        return 'both'
