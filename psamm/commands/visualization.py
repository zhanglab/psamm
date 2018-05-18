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
                       MetabolicMixin, Command, FilePrefixAppendAction)
import csv
from ..reaction import Direction
from six import text_type, iteritems, itervalues, string_types, integer_types
from .. import findprimarypairs
from ..formula import Formula, Atom, ParseError
from .. import graph
from collections import Counter
from tableexport import _encode_value
import argparse


logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#fbb4ae'
ACTIVE_COLOR = '#b3cde3'
ALT_COLOR = '#f4fc55'


def cpds_properties(compound, detail):
    """define compound nodes label"""
    if detail == 'basic':
        label = compound.id
    elif detail == 'medium':
        label = '{}\n{}'.format(compound.id, compound.name)
    elif detail == 'all':
        compound_set = set()
        compound_set.update(compound.properties)
        compound_list_sorted = sorted(
            compound_set, key=lambda x: (x != 'id', x != 'name', x))
        label = '\n'.join(_encode_value(compound.properties[value])
                          for value in compound_list_sorted)

    return label


def rxns_properties(reaction, detail):
    """define reaction nodes label"""
    if detail == 'basic':
        label = reaction.id
    elif detail == 'medium':
        label = '{}\n{}'.format(reaction.id, reaction.name)
    elif detail == 'all':
        reaction_set = set()
        reaction_set.update(reaction.properties)
        reaction_list_sorted = sorted(
            reaction_set, key=lambda x: (x != 'id', x != 'equation', x))
        label = '\n'.join(_encode_value(reaction.properties[value])
                          for value in reaction_list_sorted)

    return label


class VisualizationCommand(MetabolicMixin, ObjectiveMixin,
                         SolverCommandMixin, Command, LoopRemovalMixin, FilePrefixAppendAction):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help=('Reaction to exclude (e.g. biomass reactions or'
                  ' macromolecule synthesis)'))
        parser.add_argument(
            '--edge-values', type=text_type, default=None,
            help='Values for edges, derived from reaction flux')
        parser.add_argument(
            '--element', type=text_type, default=None,
            help='primary element flow')
        # parser.add_argument(
        #     '--detail', type=text_type, default=None,
        #     help='reaction and compound information showed on nodes label')
        parser.add_argument(
            '--detail',
            choices= ['basic', 'medium', 'all'],
            default='basic', help='reaction and compound information showed on nodes label')
        parser.add_argument(
            '--subset', type=text_type, default=None,
            help='reactions designated to visualize')
        parser.add_argument(
            '--color', type=argparse.FileType('r'), default=None, nargs='+',
            help='customize color of reaction and compound nodes ')
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
        edge_values = None
        if self._args.edge_values is not None:
            raw_values = {}
            with open(self._args.edge_values, 'r') as f:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    raw_values[row[0]] = float(row[1])

            edge_values = {}
            for reaction in self._mm.reactions:
                rx = self._mm.get_reaction(reaction)
                if reaction in raw_values:
                    flux = raw_values[reaction]
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

        #reaction_pairs = [(r, r.equation) for r in self._mm.reactions if r not in self._args.exclude]
        reaction_pairs = [(r, self._mm.get_reaction(r)) for r in self._mm.reactions if r not in self._args.exclude]
        element_weight = findprimarypairs.element_weight
        fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(reaction_pairs, compound_formula,
                                                                       element_weight=element_weight)

        filter_fpp = {}
        if self._args.element is None:
            filter_fpp = fpp_dict
        else:
            for rxn_id, fpp_pairs in fpp_dict.iteritems():
                pairs_tmp = {}
                for cpair, transfer in fpp_pairs[0].iteritems():
                    if any(Atom(self._args.element) in k for k in transfer):
                        pairs_tmp[cpair] = transfer
                pairs_tmp_2 = (pairs_tmp, {})
                filter_fpp[rxn_id] = pairs_tmp_2

        g = self.create_split_bipartite_graph(self._mm, self._model,
                                              filter_fpp, self._args.element,edge_values, compound_formula)
        with open('reactions.dot', 'w') as f:
            g.write_graphviz(f)
        with open('reactions.nodes.tsv', 'w') as f:
            g.write_cytoscape_nodes(f)
        with open('reactions.edges.tsv', 'w') as f:
            g.write_cytoscape_edges(f)

    def create_split_bipartite_graph(self, model, nativemodel, fpp_results, element, edge_values, cpd_formula):
        """create bipartite graph from fpp_dict"""
        g = graph.Graph()

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

        edge_value_span = max_edge_value - min_edge_value

        color = {}
        # for reaction in model.reactions:
        #     color[reaction] = REACTION_COLOR
        if self._args.color is not None:
            recolor_nodes = []
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    color[row[0]] = row[1]  # row[0] =reaction id, row[1] = hex color code, such as #cfe0fc
                    recolor_nodes.append(row[0])
                    for reaction in model.reactions:
                        if reaction not in recolor_nodes:
                            color[reaction] = REACTION_COLOR
                    for compound in model.compounds:
                        if str(compound.name) not in recolor_nodes:
                            color[compound.name] = COMPOUND_COLOR
                            # print(type(compound.name))
        else:
            for r in model.reactions:
                color[r] = REACTION_COLOR
            for c in model.compounds:
                color[c.name] = COMPOUND_COLOR
        #
        # def cpds_properties(compound, detail):
        #     """define compound nodes label"""
        #     if detail is not None:
        #         label = '\n'.join(_encode_value(compound.properties[value])
        #                           for value in self._args.detail)
        #     else:
        #         label = compound.id
        #
        #     return label
        #
        # def rxns_properties(reaction, detail):
        #     """define reaction nodes label"""
        #     if self._args.edge_values is not None:
        #         if detail is not None:
        #             reaction_set =
        #             label = '\n'.join(_encode_value(reaction.properties[value])
        #                               for value in self._args.detail)
        #         else:
        #             label = '{}\n{}'.format(reaction.id, abs(raw_values[reaction.id]))
        #     else:
        #         if detail is not None:
        #             label = '\n'.join(_encode_value(reaction.properties[value])
        #                               for value in self._args.detail)
        #         else:
        #             label = reaction.id

        if self._args.subset is not None:  # a file contains reaction_id in the subset users want to visualize
            subset_reactions = []
            with open(self._args.subset, 'r') as f:
                for line in f.readlines():
                    subset_reactions.append(line)
        else:
            subset_reactions = set(model.reactions)

        fpp_cpds = []  # cpds in fpp_results
        fpp_rxns = Counter()
        compound_nodes = {}

        for rxn_id, fpp_pairs in sorted(fpp_results.iteritems()):
            if rxn_id in subset_reactions:
                for (c1, c2), transfer in fpp_pairs[0].iteritems():
                    if c1 not in fpp_cpds:          # c1.name = compound.id, no compartment
                        node = graph.Node({
                            'id': text_type(c1),
                            'label': cpds_properties(cpd_entry[c1.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': color[c1.name]})
                        g.add_node(node)
                        compound_nodes[c1] = node
                    if c2 not in fpp_cpds:
                        node = graph.Node({
                            'id': text_type(c2),
                            'label': cpds_properties(cpd_entry[c2.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': color[c2.name]})
                        g.add_node(node)
                        compound_nodes[c2] = node
                    fpp_cpds.append(c1)
                    fpp_cpds.append(c2)

                    fpp_rxns[rxn_id] += 1
                    node = graph.Node({
                        'id': '{}_{}'.format(rxn_id, fpp_rxns[rxn_id]),
                        'label': rxns_properties(rxn_entry[rxn_id], self._args.detail),
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': color[rxn_id]})
                    g.add_node(node)

                    def pen_width(value):
                        """calculate edges width"""
                        if edge_value_span == 0:
                            return 1
                        else:
                            alpha = value / edge_value_span
                            return 19 * alpha + 1

                    def dir_value(direction):
                        """assign value to different reaction directions"""
                        if direction == Direction.Forward:
                            return 'forward'
                        elif direction == Direction.Reverse:
                            return 'back'
                        else:
                            return 'both'

                    rx = model.get_reaction(rxn_id)
                    edge_props = {'dir': dir_value(rx.direction)}

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
                            return dict(edge_props)

                    edge1 = c1, rxn_id  # forward
                    edge2 = rxn_id, c1  # backward
                    g.add_edge(graph.Edge(
                        compound_nodes[c1], node, final_props(rx, edge1, edge2)))

                    edge1 = rxn_id, c2
                    edge2 = c2, rxn_id
                    g.add_edge(graph.Edge(
                        node, compound_nodes[c2], final_props(rx, edge1, edge2)))

        for reaction in model.reactions:
            if model.is_exchange(reaction):
                if reaction in subset_reactions:
                    exchange_rxn = model.get_reaction(reaction)
                    if any(Atom(element) in cpd_formula[str(c[0].name)] for c in exchange_rxn.compounds):
                        node_ex = graph.Node({
                            'id': reaction,
                            'label': reaction,
                            'shape': 'box',
                            'style': 'filled',
                            'fillcolor': ACTIVE_COLOR})
                        g.add_node(node_ex)

                        for c, _ in exchange_rxn.left:
                            if c not in compound_nodes.keys():
                                node_ex = graph.Node({
                                    'id': text_type(c),
                                    'label': cpds_properties(cpd_entry[c.name], self._args.detail),
                                    'shape': 'ellipse',
                                    'style': 'filled',
                                    'fillcolor':'#5a95f4'})
                                g.add_node(node_ex)
                                compound_nodes[c] = node_ex

                        edge1 = reaction, c
                        edge2 = c, reaction
                        g.add_edge(graph.Edge(
                            node_ex, compound_nodes[c], final_props(exchange_rxn, edge1, edge2)))

        if nativemodel.biomass_reaction is not None:
            bio_pair = Counter()
            biomass_reaction = model.get_reaction(nativemodel.biomass_reaction)
            for c, _ in biomass_reaction.left:
                if Atom(element) in cpd_formula[str(c.name)]:
                    bio_pair[nativemodel.biomass_reaction] += 1  # bio_pair = Counter({'biomass_rxn_id': 1}), Counter({'biomass_rxn_id': 2})...
                    node_bio = graph.Node({
                        'id': '{}_{}'.format(nativemodel.biomass_reaction, bio_pair[nativemodel.biomass_reaction]),
                        'label': nativemodel.biomass_reaction,
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': ALT_COLOR})
                    g.add_node(node_bio)

                    if c not in compound_nodes.keys():
                        node_bio_cpd = graph.Node({
                            'id': text_type(c),
                            'label': cpds_properties(cpd_entry[c.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': '#82e593'})
                        g.add_node(node_bio_cpd)
                        compound_nodes[c] = node_bio_cpd

                    edge1 = c, nativemodel.biomass_reaction
                    edge2 = nativemodel.biomass_reaction, c
                    g.add_edge(graph.Edge(
                        compound_nodes[c], node_bio, final_props(biomass_reaction, edge1, edge2)))
        else:
            logger.warning(
                'No biomass reaction in this model')

        return g

