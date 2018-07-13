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

import logging
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command, FilePrefixAppendAction, CommandError)
import csv
from ..reaction import Direction
from six import text_type, iteritems, itervalues, iterkeys
from .. import findprimarypairs
from ..formula import Formula, Atom, ParseError
from .. import graph
from collections import Counter
from .tableexport import _encode_value
import argparse
from .. import fluxanalysis
from collections import defaultdict, namedtuple

try:
    from graphviz import render
except ImportError:
    render = None

import subprocess

# from py2cytoscape.data.cyrest_client import CyRestClient

logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#fbb4ae'
ACTIVE_COLOR = '#b3cde3'    # exchange reaction  color
ALT_COLOR = '#f4fc55'       # biomass reaction color
RXN_COMBINED_COLOR = '#fc9a44'
CPD_ONLY_IN_BIO = '#82e593'
CPD_ONLY_IN_EXC = '#5a95f4'


def cpds_properties(cpd, compound, detail): # cpd=Compound object, compound = CompoundEntry object
    """define compound nodes label"""
    compound_set = set()
    compound_set.update(compound.properties)
    if detail is not None:
        cpd_detail = []
        for prop in detail[0]:
            if prop in compound_set:
                cpd_detail.append(str(prop))
        pre_label = '\n'.join(_encode_value(compound.properties[value]) for value in cpd_detail if value != 'id')
        label = '{}\n{}'.format(str(cpd), pre_label)
    else:
        label = str(cpd)
    return label


def rxns_properties(reaction, detail, reaction_flux):
    """define reaction nodes label"""
    reaction_set = set()
    reaction_set.update(reaction.properties)
    if detail is not None:
        rxn_detail = []
        for prop in detail[0]:
            if prop in reaction_set:
                rxn_detail.append(str(prop))
        label = '\n'.join(_encode_value(reaction.properties[value])
                          for value in rxn_detail)
        if len(reaction_flux) > 0:
            if reaction.id in iterkeys(reaction_flux):
                label = '{}\n{}'.format(label, reaction_flux[reaction.id])
    else:
        if len(reaction_flux) > 0:
            if reaction.id in iterkeys(reaction_flux):
                label = '{}\n{}'.format(reaction.id, reaction_flux[reaction.id])
            else:
                label = reaction.id
        else:
            label = reaction.id

    return label


def primary_element(element):
    if element is not None:
        if element in ['c', 'h', 'o', 'n', 'p', 's', 'k', 'b', 'f', 'v', 'y', 'i', 'w']:
            return element.upper()
        else:
            return element


class VisualizationCommand(MetabolicMixin, ObjectiveMixin,
                           SolverCommandMixin, Command, LoopRemovalMixin, FilePrefixAppendAction):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method',type=text_type,
            default='fpp', help='Compound pair prediction method')
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help=('Reaction to exclude (e.g. biomass reactions or'
                  ' macromolecule synthesis)'))
        parser.add_argument(
            '--edge-values', type=argparse.FileType('r'), default=None,
            help='Values for edges, derived from reaction flux')
        parser.add_argument(
            '--fba', action='store_true',
            help='flux balance analysis')
        parser.add_argument(
            '--element', type=text_type, default='C',
            help='primary element flow')
        parser.add_argument(
            '--detail', type=text_type, default=None, action='append', nargs='+',
            help='reaction and compound properties showed on nodes label')
        parser.add_argument(
            '--subset', type=argparse.FileType('r'), default=None,
            help='reactions designated to visualize')
        parser.add_argument(
            '--color', type=argparse.FileType('r'), default=None, nargs='+',
            help='customize color of reaction and compound nodes ')
        parser.add_argument(
            '--Image', type=text_type, default=None, help='generate image file directly')
        parser.add_argument(
            '--exclude-pairs', type=argparse.FileType('r'), default=None,
            help='remove edge of given compound pairs from final graph ')
        parser.add_argument(
            '--split-map', action='store_true',
            help='create dot file for reaction-splitted metabolic network visualization')
        # parser.add_argument(
        #     '--Cytoscape', action='store_true', help='generate graph in Cytoscape')
        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""

        # parse compound id and formula
        compound_formula = {}
        for compound in self._model.compounds:
            if compound.formula is not None:
                try:
                    f = Formula.parse(compound.formula).flattened()
                    if not f.is_variable():
                        compound_formula[compound.id] = f
                    else:
                        logger.warning(
                            'Skipping variable formula {}: {}'.format(
                                compound.id, compound.formula))    #skip generic compounds
                except ParseError as e:
                    msg = (
                        'Error parsing formula'
                        ' for compound {}:\n{}\n{}'.format(
                            compound.id, e, compound.formula))
                    if e.indicator is not None:
                        msg += '\n{}'.format(e.indicator)
                    logger.warning(msg)

        # set edge_values
        reaction_flux = {}
        if self._args.fba is True:
            if self._args.edge_values is None:
                solver = self._get_solver()
                p = fluxanalysis.FluxBalanceProblem(self._mm, solver)
                try:
                    p.maximize(self._get_objective())
                except fluxanalysis.FluxBalanceError as e:
                    self.report_flux_balance_error(e)

                fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

                # Run flux minimization
                flux_var = p.get_flux_var(self._get_objective())
                p.prob.add_linear_constraints(flux_var == p.get_flux(self._get_objective()))
                p.minimize_l1()

                count = 0
                for r_id in self._mm.reactions:
                    flux = p.get_flux(r_id)
                    if abs(flux - fluxes[r_id]) > 1e-6:
                        count += 1
                    if abs(flux) > 1e-6:
                        reaction_flux[r_id] = flux
                logger.info('Minimized reactions: {}'.format(count))

            else:
                logger.warning('Invalid command, the two arguments --flux-analysis and --edge-values '
                               'should not present at the same time')
                quit()

        else:
            if self._args.edge_values is not None:
                for row in csv.reader(self._args.edge_values, delimiter=str(u'\t')):
                    reaction_flux[row[0]] = float(row[1])

        edge_values = None
        if len(reaction_flux) > 0:
            edge_values = {}
            for reaction in self._mm.reactions:
                rx = self._mm.get_reaction(reaction)
                if reaction in reaction_flux:
                    flux = reaction_flux[reaction]
                    if abs(flux) < 1e-9:
                        continue

                    if flux > 0:
                        for compound, value in rx.right:  # value=stoichiometry
                            edge_values[reaction, compound] = (flux * float(value))
                        for compound, value in rx.left:
                            edge_values[compound, reaction] = (flux * float(value))
                    else:
                        for compound, value in rx.left:
                            edge_values[reaction, compound] = (- flux * float(value))
                        for compound, value in rx.right:
                            edge_values[compound, reaction] = (- flux * float(value))

        cpd_object = {}
        for cpd in self._mm.compounds:
            cpd_object[str(cpd)] = cpd

        if self._args.subset is not None:  # a file contains reaction_id in the subset users want to visualize
            subset_reactions = []
            for line in self._args.subset.readlines():
                subset_reactions.append(line.rstrip())
        else:
            subset_reactions = set(self._mm.reactions)

        def iter_reactions():
            for reaction in self._model.reactions:
                if (reaction.id not in self._model.model or
                        reaction.id in self._args.exclude):
                    continue

                if reaction.equation is None:
                    logger.warning(
                        'Reaction {} has no reaction equation, so remove this reaction from final graph'.format(
                            reaction.id))
                    continue

                if any(c.name not in compound_formula
                       for c, _ in reaction.equation.compounds):
                    logger.warning(
                        'Reaction {} has compounds with undefined formula, so remove this reaction from '
                        'final graph'.format(reaction.id))
                    continue

                yield reaction

        exclude_cpairs = []
        if self._args.exclude_pairs is not None:
            for row in csv.reader(self._args.exclude_pairs, delimiter=str('\t')):
                exclude_cpairs.append((cpd_object[row[0]], cpd_object[row[1]]))
                exclude_cpairs.append((cpd_object[row[1]], cpd_object[row[0]]))

        filter_dict = {}
        if self._args.method == 'fpp':
            split_reactions = True

            reaction_pairs = [(r.id, r.equation) for r in iter_reactions()]
            element_weight = findprimarypairs.element_weight
            fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(reaction_pairs, compound_formula,
                                                                           element_weight=element_weight)
            for rxn_id, fpp_pairs in iteritems(fpp_dict):
                compound_pairs = []
                for cpd_pair, transfer in iteritems(fpp_pairs[0]):
                    if cpd_pair not in exclude_cpairs:
                        if self._args.element is None:
                            compound_pairs.append(cpd_pair)
                        else:
                            if any(Atom(primary_element(self._args.element)) in k for k in transfer):
                                compound_pairs.append(cpd_pair)
                filter_dict[rxn_id] = compound_pairs

        elif self._args.method == 'no-fpp':
            split_reactions = False
            for rxn in iter_reactions():
                if rxn.id != self._model.biomass_reaction:
                    rx = self._mm.get_reaction(rxn.id)
                    cpairs = []
                    for c1, _ in rx.left:
                        for c2, _ in rx.right:
                            if (c1, c2) not in exclude_cpairs:
                                if self._args.element is not None:
                                    if Atom(primary_element(self._args.element)) in compound_formula[c1.name]:
                                        if Atom(primary_element(self._args.element)) in compound_formula[c2.name]:
                                            cpairs.append((c1, c2))
                                else:
                                    cpairs.append((c1, c2))
                    filter_dict[rxn.id] = cpairs
        else:
            split_reactions = True
            with open(self._args.method, 'r') as f:
                cpair_list, rxn_list = [], []
                for row in csv.reader(f, delimiter=str(u'\t')):
                    if (cpd_object[row[1]], cpd_object[row[2]]) not in exclude_cpairs:
                        if self._args.element is None:
                            cpair_list.append((cpd_object[row[1]], cpd_object[row[2]]))
                            rxn_list.append(row[0])
                        else:
                            if Atom(primary_element(self._args.element)) in Formula.parse(row[3]).flattened():
                                cpair_list.append((cpd_object[row[1]], cpd_object[row[2]]))
                                rxn_list.append(row[0])

                filter_dict = defaultdict(list)
                for r, cpair in zip(rxn_list, cpair_list):
                    filter_dict[r].append(cpair)

        raw_cpairs_dict = defaultdict(list)     # key=compound pair, value=list of reaction_id
        raw_dict = {k: v for k, v in iteritems(filter_dict) if k in subset_reactions}
        for rxn_id, cpairs in iteritems(raw_dict):
            for pair in cpairs:
                raw_cpairs_dict[pair].append(rxn_id)

        cpairs_dict = {}
        pair_list = []
        for (c1, c2), rxns in iteritems(raw_cpairs_dict):
            if (c1, c2) not in pair_list:
                forward_rxns, back_rxns, bidirectional_rxns = [], [], []
                for r in rxns:
                    reaction = self._mm.get_reaction(r)
                    if reaction.direction == Direction.Forward:
                        forward_rxns.append(r)
                    elif reaction.direction == Direction.Reverse:
                        back_rxns.append(r)
                    else:
                        bidirectional_rxns.append(r)

                if (c2, c1) in iterkeys(raw_cpairs_dict):
                    for r in raw_cpairs_dict[(c2, c1)]:
                        reaction = self._mm.get_reaction(r)
                        if reaction.direction == Direction.Forward:
                            back_rxns.append(r)
                        elif reaction.direction == Direction.Reverse:
                            forward_rxns.append(r)
                        else:
                            bidirectional_rxns.append(r)
                cpair_rxn = namedtuple('cpair_rxn', ['forward', 'back', 'bidirection'])
                cpairs_dict[(c1, c2)] = cpair_rxn._make([forward_rxns, back_rxns, bidirectional_rxns])
            pair_list.append((c1, c2))
            pair_list.append((c2, c1))

        g, g1, g2 = self.create_split_bipartite_graph(self._mm, self._model, filter_dict, cpairs_dict,
                                                      self._args.element, subset_reactions, edge_values,
                                                      compound_formula, reaction_flux, split_graph=split_reactions)
        final_graph = None
        if self._args.method != 'no-fpp':
            if self._args.split_map is True:
                final_graph = g
            else:
                final_graph = g2
        else:
            if self._args.split_map is True:
                logger.warning('--split-map option can\'t be applied on visualization when method is no-fpp,'
                               ' break program')
                quit()
            else:
                final_graph = g

        with open('reactions.dot', 'w') as f:
            final_graph.write_graphviz(f)
        with open('reactions.nodes.tsv', 'w') as f:
            final_graph.write_cytoscape_nodes(f)
        with open('reactions.edges.tsv', 'w') as f:
            final_graph.write_cytoscape_edges(f)

        # raw_network_1 = nx.drawing.nx_agraph.read_dot('./combined-reactions.dot')  # NetworkX graph object
        # raw_network_2 = from_networkx(raw_network_1)    # convert NetworkX graph object to cyjs
        # cy = CyRestClient()
        # network = cy.network.create(name='metabolic network', data=raw_network_2)
        # print(network.get_id())

        if self._args.Image is not None:
            if render is None:
                self.fail(
                    'create image file requires python binding graphviz module'
                    ' ("pip install graphviz")')
            else:
                if len(subset_reactions) > 500:
                    logger.info(
                        'The program is going to create a large graph that contains {} reactions, '
                        'it may take a long time'.format(len(subset_reactions)))
                try:
                    render('dot', self._args.Image, 'reactions.dot')
                except subprocess.CalledProcessError:
                    logger.warning('the graph is too large to create')

    def create_split_bipartite_graph(self, model, nativemodel, predict_results, cpair_dict, element, subset,
                                     edge_values,  cpd_formula, reaction_flux, split_graph=True):
        """create bipartite graph from filter_dict"""
        g = graph.Graph()
        g1 = graph.Graph()
        g2 = graph.Graph()

        cpd_entry = {}
        for cpd in nativemodel.compounds:
            cpd_entry[cpd.id] = cpd

        rxn_entry = {}
        for rxn in nativemodel.reactions:
            rxn_entry[rxn.id] = rxn

        if edge_values is not None and len(edge_values) > 0:
            min_edge_value = min(itervalues(edge_values))
            max_edge_value = max(itervalues(edge_values))
        else:
            min_edge_value = 1
            max_edge_value = 1

        # def final_color(color_arg, node):   # node = str(compound object) or reaction id
        #     color = {}
        #     all_cpds = []
        #     for c in model.compounds:
        #         all_cpds.append(str(c))
        #     if color_arg is not None:
        #         recolor_dict = {}
        #         for f in color_arg:
        #             for row in csv.reader(f, delimiter=str(u'\t')):
        #                 recolor_dict[row[0]] = row[1]  # row[0] =reaction id, row[1] = hex color code, such as #cfe0fc
        #         if node in model.reactions:
        #             if node in recolor_dict:
        #                 color[node] = recolor_dict[node]
        #             else:
        #                 color[node] = REACTION_COLOR
        #         elif node in all_cpds:
        #             if node in recolor_dict:
        #                 color[node] = recolor_dict[node]
        #             else:
        #                 color[node] = COMPOUND_COLOR
        #         else:
        #             logger.warning('Invalid compound or reaction:{}'.format(node))
        #     else:
        #         if node in model.reactions:
        #             color[node] = REACTION_COLOR
        #         elif cpd_object[node] in model.compounds:
        #             color[node] = COMPOUND_COLOR
        #         else:
        #             logger.warning('Invalid compound or reaction:{}'.format(node))
        #
        #     return color[node]

        recolor_dict = {}
        if self._args.color is not None:
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    recolor_dict[row[0]] = row[1]  # row[0] =reaction id or str(cpd object), row[1] = hex color code
        color = {}
        for c in model.compounds:
            if str(c) in recolor_dict:
                color[c] = recolor_dict[str(c)]
            else:
                color[c] = COMPOUND_COLOR
        for r in model.reactions:
            if r in recolor_dict:
                color[r] = recolor_dict[r]
            else:
                color[r] = REACTION_COLOR

        def pen_width(value):
            """calculate edges width"""
            if max_edge_value - min_edge_value == 0:
                return 1
            else:
                alpha = (value - min_edge_value) / (max_edge_value - min_edge_value)
                return 20 * alpha

        def dir_value(direction):
            """assign value to different reaction directions"""
            if direction == Direction.Forward:
                return 'forward'
            elif direction == Direction.Reverse:
                return 'back'
            else:
                return 'both'

        def final_props(reaction, edge1, edge2):
            """set final properties of edges"""
            if edge_values is not None:
                p = {}
                if edge1 in edge_values:
                    p['dir'] = 'forward'
                    p['penwidth'] = pen_width(edge_values[edge1])
                    # print(edge1)
                    # print('forward\t{}\t{}\t{}'.format(reaction, edge1, edge_values[edge1]))
                elif edge2 in edge_values:
                    p['dir'] = 'back'
                    p['penwidth'] = pen_width(edge_values[edge2])
                    # print('reverse\t{}\t{}\t{}'.format(reaction, edge2, edge_values[edge2]))
                else:
                    p['style'] = 'dotted'
                    p['dir'] = dir_value(reaction.direction)
                return p
            else:
                return {'dir': dir_value(reaction.direction)}

        # create standard bipartite graph
        cpds = []  # cpds in predict_results
        rxns = Counter()
        compound_nodes = {}
        edge_list = []
        for rxn_id, cpairs in sorted(iteritems(predict_results)):
            if rxn_id in subset:
                for (c1, c2) in cpairs:
                    if c1 not in cpds:          # c1.name = compound.id, no compartment
                        node = graph.Node({
                            'id': text_type(c1),
                            'edge_id': text_type(c1),
                            'label': cpds_properties(c1, cpd_entry[c1.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': color[c1]})
                        g.add_node(node)
                        g1.add_node(node)
                        compound_nodes[c1] = node
                    if c2 not in cpds:
                        node = graph.Node({
                            'id': text_type(c2),
                            'edge_id': text_type(c2),
                            'label': cpds_properties(c2, cpd_entry[c2.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': color[c2]})
                        g.add_node(node)
                        g1.add_node(node)
                        compound_nodes[c2] = node
                    cpds.append(c1)
                    cpds.append(c2)

                    if split_graph is True:
                        rxns[rxn_id] += 1
                        node_id = '{}_{}'.format(rxn_id, rxns[rxn_id])
                    else:
                        node_id = rxn_id
                    node = graph.Node({
                        'id': node_id,
                        'edge_id': rxn_id,
                        'label': rxns_properties(rxn_entry[rxn_id], self._args.detail, reaction_flux),
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': color[rxn_id]})
                    g.add_node(node)

                    rx = model.get_reaction(rxn_id)

                    edge1 = c1, rxn_id  # forward
                    edge2 = rxn_id, c1
                    if split_graph is True:
                        g.add_edge(graph.Edge(
                            compound_nodes[c1], node, final_props(rx, edge1, edge2)))
                    else:
                        if edge1 and edge2 not in edge_list:
                            edge_list.append(edge1)
                            edge_list.append(edge2)
                            g.add_edge(graph.Edge(
                                compound_nodes[c1], node, final_props(rx, edge1, edge2)))

                    edge1 = rxn_id, c2 # forward
                    edge2 = c2, rxn_id
                    if split_graph is True:
                        g.add_edge(graph.Edge(
                            node, compound_nodes[c2], final_props(rx, edge1, edge2)))
                    else:
                        if edge1 and edge2 not in edge_list:
                            edge_list.append(edge1)
                            edge_list.append(edge2)
                            g.add_edge(graph.Edge(
                                node, compound_nodes[c2], final_props(rx, edge1, edge2)))

                    g1.add_edge(graph.Edge(
                        compound_nodes[c1], compound_nodes[c2], {'dir': dir_value(rx.direction), 'reaction': rxn_id}))

        # create bipartite and reactions-combined graph if --method is nt no-fpp
        cpd_nodes = {}
        cpd_pairs = Counter()
        compound_list = []
        for (c1, c2), rxns in iteritems(cpair_dict):
            if c1 not in compound_list:  # c1.name = compound.id, no compartment
                node = graph.Node({
                    'id': text_type(c1),
                    'label': cpds_properties(c1, cpd_entry[c1.name], self._args.detail),
                    'shape': 'ellipse',
                    'style': 'filled',
                    'fillcolor': color[c1]})
                g2.add_node(node)
                cpd_nodes[c1] = node
            if c2 not in compound_list:  # c1.name = compound.id, no compartment
                node = graph.Node({
                    'id': text_type(c2),
                    'label': cpds_properties(c2, cpd_entry[c2.name], self._args.detail),
                    'shape': 'ellipse',
                    'style': 'filled',
                    'fillcolor': color[c2]})
                g2.add_node(node)
                cpd_nodes[c2] = node
            compound_list.append(c1)
            compound_list.append(c2)

            def final_rxn_color(color_args, rlist):
                if color_args is not None:
                    if len(rlist) == 1:
                        return color[rlist[0]]
                    else:
                        if any(r in recolor_dict for r in rlist):
                            return RXN_COMBINED_COLOR
                        else:
                            return REACTION_COLOR
                else:
                    return REACTION_COLOR

            def condensed_rxn_props(detail, r_list, reaction_flux):
                if len(r_list) == 1:
                    label_comb = rxns_properties(rxn_entry[r_list[0]], detail, reaction_flux)
                else:
                    label_comb = '\n'.join(r for r in r_list)
                return label_comb

            def dir_value_2(r_list):
                if r_list == rxns.forward:
                    return {'dir': 'forward'}
                elif r_list == rxns.back:
                    return {'dir': 'back'}
                else:
                    return {'dir': 'both'}

            def final_props_2(rlist, c, n):
                if edge_values is not None:
                    if len(rlist) == 1:
                        if n == 1:
                            edge1 = c, rlist[0]
                            edge2 = rlist[0], c
                        else:
                            edge1 = rlist[0], c
                            edge2 = c, rlist[0]
                        rx = model.get_reaction(rlist[0])
                        return final_props(rx, edge1, edge2)
                    else:
                        if n == 1:
                            edges1 = [(r, c) for r in rlist]
                            edges2 = [(c, r) for r in rlist]
                        else:
                            edges1 = [(c, r) for r in rlist]
                            edges2 = [(r, c) for r in rlist]
                        if any(i for i in edges1) and any(j for j in edges2) in edge_values:
                            p = {}
                            p['style'] = 'dotted'
                            p['dir'] = dir_value_2(rlist)
                            return p
                        else:
                            return dir_value_2(rlist)
                else:
                    return dir_value_2(rlist)

                # if c1 in reaction.left :
                #     return {'dir': dir_value(reaction.direction)}
                # else:
                #     if reaction.direction == Direction.Forward:
                #         return {'dir': dir_value(Direction.Reverse)}
                #     elif reaction.direction == Direction.Reverse:
                #         return {'dir': dir_value(Direction.Forward)}
                #     else:
                #         return {'dir': dir_value(Direction.Both)}

            for r_list in rxns:
                if len(r_list) > 0:
                    cpd_pairs[(c1, c2)] += 1
                    node = graph.Node({
                        'id': '{}_{}'.format((str(c1), str(c2)), cpd_pairs[(c1, c2)]),
                        'label': condensed_rxn_props(self._args.detail, r_list, reaction_flux),
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': final_rxn_color(self._args.color, r_list)})
                    g2.add_node(node)

                    g2.add_edge(graph.Edge(
                        cpd_nodes[c1], node, final_props_2(r_list,c1, 1)))

                    g2.add_edge(graph.Edge(
                        node, cpd_nodes[c2], final_props_2(r_list,c2, 2)))

        # add exchange reaction nodes
        rxn_set = set()
        for reaction in model.reactions:
            if model.is_exchange(reaction):
                if reaction in subset:
                    raw_exchange_rxn = model.get_reaction(reaction)
                    if element is not None:
                        if any(Atom(primary_element(element)) in cpd_formula[str(c[0].name)] for c in raw_exchange_rxn.compounds):
                            rxn_set.add(reaction)
                    else:
                        rxn_set.add(reaction)
        for r in rxn_set:
            exchange_rxn = model.get_reaction(r)
            label = r
            if len(reaction_flux) > 0:
                if r in iterkeys(reaction_flux):
                    label = '{}\n{}'.format(r, reaction_flux[r])
            node_ex = graph.Node({
                'id': r,
                'edge_id': r,
                'label': label,
                'shape': 'box',
                'style': 'filled',
                'fillcolor': ACTIVE_COLOR})
            g.add_node(node_ex)
            g2.add_node(node_ex)

            for c1, _ in exchange_rxn.left:
                if c1 not in compound_nodes.keys():
                    node_ex_cpd = graph.Node({
                        'id': text_type(c1),
                        'edge_id':text_type(c1),
                        'label': cpds_properties(c1, cpd_entry[c1.name], self._args.detail),
                        'shape': 'ellipse',
                        'style': 'filled',
                        'fillcolor': CPD_ONLY_IN_EXC})
                    g.add_node(node_ex_cpd)
                    compound_nodes[c1] = node_ex_cpd

                    g2.add_node(node_ex_cpd)
                    cpd_nodes[c1] = node_ex_cpd

                edge1 = c1, r
                edge2 = r, c1
                g.add_edge(graph.Edge(
                    compound_nodes[c1], node_ex, final_props(exchange_rxn, edge1, edge2)))
                g2.add_edge(graph.Edge(
                    cpd_nodes[c1], node_ex, final_props(exchange_rxn, edge1, edge2)))

            for c2, _ in exchange_rxn.right:
                if c2 not in compound_nodes.keys():
                    node_ex_cpd = graph.Node({
                        'id': text_type(c2),
                        'edge_id': text_type(c2),
                        'label': cpds_properties(c2, cpd_entry[c2.name], self._args.detail),
                        'shape': 'ellipse',
                        'style': 'filled',
                        'fillcolor': CPD_ONLY_IN_EXC})
                    g.add_node(node_ex_cpd)
                    compound_nodes[c2] = node_ex_cpd

                    g2.add_node(node_ex_cpd)
                    cpd_nodes[c2] = node_ex_cpd

                edge1 = r, c2
                edge2 = c2, r
                g.add_edge(graph.Edge(
                    node_ex, compound_nodes[c2], final_props(exchange_rxn, edge1, edge2)))
                g2.add_edge(graph.Edge(
                    node_ex, cpd_nodes[c2], final_props(exchange_rxn, edge1, edge2)))

        # add biomass reaction nodes
        bio_pair = Counter()
        biomass_cpds = set()
        if nativemodel.biomass_reaction is not None:
            if nativemodel.biomass_reaction in subset:
                biomass_reaction = model.get_reaction(nativemodel.biomass_reaction)
                for c, _ in biomass_reaction.left:
                    if element is not None:
                        if Atom(primary_element(element)) in cpd_formula[str(c.name)]:
                            biomass_cpds.add(c)
                    else:
                        biomass_cpds.add(c)
                for c in biomass_cpds:
                    bio_pair[nativemodel.biomass_reaction] += 1  # bio_pair = Counter({'biomass_rxn_id': 1}), Counter({'biomass_rxn_id': 2})...
                    node_bio = graph.Node({
                        'id': '{}_{}'.format(nativemodel.biomass_reaction, bio_pair[nativemodel.biomass_reaction]),
                        'edge_id': nativemodel.biomass_reaction,
                        'label': nativemodel.biomass_reaction,
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': ALT_COLOR})
                    g.add_node(node_bio)
                    g2.add_node(node_bio)

                    if c not in compound_nodes.keys():
                        node_bio_cpd = graph.Node({
                            'id': text_type(c),
                            'edge_id': text_type(c),
                            'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': CPD_ONLY_IN_BIO})
                        g.add_node(node_bio_cpd)
                        compound_nodes[c] = node_bio_cpd

                        g2.add_node(node_bio_cpd)
                        cpd_nodes[c] = node_bio_cpd

                    edge1 = c, nativemodel.biomass_reaction
                    edge2 = nativemodel.biomass_reaction, c
                    g.add_edge(graph.Edge(
                        compound_nodes[c], node_bio, final_props(biomass_reaction, edge1, edge2)))
                    g2.add_edge(graph.Edge(
                        cpd_nodes[c], node_bio, final_props(biomass_reaction, edge1, edge2)))
        else:
            logger.warning(
                'No biomass reaction in this model')

        return g, g1, g2
