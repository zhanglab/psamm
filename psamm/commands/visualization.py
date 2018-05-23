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
from six import text_type, iteritems, iterkeys, itervalues, string_types, integer_types
from .. import findprimarypairs
from ..formula import Formula, Atom, ParseError
from .. import graph
from collections import Counter
from tableexport import _encode_value
import argparse
from .. import fluxanalysis

logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#fbb4ae'
ACTIVE_COLOR = '#b3cde3'
ALT_COLOR = '#f4fc55'


def cpds_properties(cpd, compound, detail): # cpd=Compound object, compound = CompoundEntry object
    """define compound nodes label"""
    compound_set = set()
    compound_set.update(compound.properties)
    if detail is not None:
        cpd_detail = []
        for prop in detail[0]:
            if prop in compound_set:
                cpd_detail.append(str(prop))
        A = '\n'.join(_encode_value(compound.properties[value])
                          for value in cpd_detail if value != 'id')
        label = '{}\n{}'.format(cpd, A)
    else:
        label = cpd
    return label


def rxns_properties(reaction, detail, flux_analysis, reaction_flux):
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
        if flux_analysis is not None:
            if reaction.id in reaction_flux.iterkeys():
                label = '{}\n{}'.format(label, reaction_flux[reaction.id])
    else:
        if flux_analysis is not None:
            if reaction.id in reaction_flux.iterkeys():
                label = '{}\n{}'.format(reaction.id, reaction_flux[reaction.id])
            else:
                label = reaction.id
        else:
            label = reaction.id

    return label


class VisualizationCommand(MetabolicMixin, ObjectiveMixin,
                         SolverCommandMixin, Command, LoopRemovalMixin, FilePrefixAppendAction):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--prediction-method',
            choices = ['fpp', 'no-prediction'],
            default = 'fpp', help='Compound pair prediction method')
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help=('Reaction to exclude (e.g. biomass reactions or'
                  ' macromolecule synthesis)'))
        parser.add_argument(
            '--edge-values', type=text_type, default=None,
            help='Values for edges, derived from reaction flux')
        parser.add_argument(
            '--flux-analysis', type=text_type, default=None,
            choices = ('None', 'fba'),
            help='flux balance analysis')
        parser.add_argument(
            '--element', type=text_type, default=None,
            help='primary element flow')
        parser.add_argument(
            '--detail', type = text_type, default=None, action='append', nargs='+',
            help='reaction and compound properties showed on nodes label')
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
        reaction_flux = {}
        if self._args.flux_analysis is not None:
            if self._args.edge_values is None:
                reaction = self._get_objective()
                if not self._mm.has_reaction(reaction):
                    self.fail('Specified reaction is not in model: {}'.format(reaction))
                solver = self._get_solver()
                p = fluxanalysis.FluxBalanceProblem(self._mm, solver)
                try:
                    p.maximze(reaction)
                except fluxanalysis.FluxBalanceError as e:
                    self.report_flux_balance_error(e)
                for reaction_id in self._mm.reactions:
                    if abs(p.get_flux(reaction_id)) > 1e-9:
                        reaction_flux[reaction_id] = p.get_flux(reaction_id)
            else:
                logger.warning('Invalid command, the two arguments --flux-analysis and --edge-values should not present at the same time')

        else:
            if self._args.edge_values is not None:
                with open(self._args.edge_values, 'r') as f:
                    for row in csv.reader(f, delimiter=str(u'\t')):
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

        #reaction_pairs = [(r, r.equation) for r in self._model.reactions if r not in self._args.exclude]
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

        if self._args.prediction_method == 'fpp':
            g = self.create_split_bipartite_graph(self._mm, self._model, filter_fpp,
                                                  self._args.element, edge_values, compound_formula, reaction_flux)
            with open('reactions.dot', 'w') as f:
                g.write_graphviz(f)
            with open('reactions.nodes.tsv', 'w') as f:
                g.write_cytoscape_nodes(f)
            with open('reactions.edges.tsv', 'w') as f:
                g.write_cytoscape_edges(f)

        elif self._args.prediction_method == 'no-prediction':
            g1 = self.create_split_bipartite_graph_without_fpp(self._mm, self._model, self._args.element,
                                                               edge_values, compound_formula, reaction_flux)
            with open('reactions.dot', 'w') as f:
                g1.write_graphviz(f)
            with open('reactions.nodes.tsv', 'w') as f:
                g1.write_cytoscape_nodes(f)
            with open('reactions.edges.tsv', 'w') as f:
                g1.write_cytoscape_edges(f)

        else:
            raise CommandError('Invalid prediction method: {}'.format(self._args.prediction_method))

    def create_split_bipartite_graph(self, model, nativemodel, fpp_results, element,
                                     edge_values, cpd_formula, reaction_flux):
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
        else:
            for r in model.reactions:
                color[r] = REACTION_COLOR
            for c in model.compounds:
                color[c.name] = COMPOUND_COLOR

        if self._args.subset is not None:  # a file contains reaction_id in the subset users want to visualize
            subset_reactions = []
            with open(self._args.subset, 'r') as f:
                for line in f.readlines():
                    subset_reactions.append(line.rstrip())
        else:
            subset_reactions = set(model.reactions)

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

        fpp_cpds = []  # cpds in fpp_results
        fpp_rxns = Counter()
        compound_nodes = {}
        for rxn_id, fpp_pairs in sorted(fpp_results.iteritems()):
            if rxn_id in subset_reactions:
                for (c1, c2), transfer in fpp_pairs[0].iteritems():
                    if c1 not in fpp_cpds:          # c1.name = compound.id, no compartment
                        node = graph.Node({
                            'id': text_type(c1),
                            'label': cpds_properties(c1, cpd_entry[c1.name], self._args.detail),
                            'shape': 'ellipse',
                            'style': 'filled',
                            'fillcolor': color[c1.name]})
                        g.add_node(node)
                        compound_nodes[c1] = node
                    if c2 not in fpp_cpds:
                        node = graph.Node({
                            'id': text_type(c2),
                            'label': cpds_properties(c2, cpd_entry[c2.name], self._args.detail),
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
                        'label': rxns_properties(rxn_entry[rxn_id], self._args.detail,
                                                 self._args.flux_analysis, reaction_flux),
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': color[rxn_id]})
                    g.add_node(node)

                    rx = model.get_reaction(rxn_id)

                    edge1 = c1, rxn_id  # forward
                    edge2 = rxn_id, c1  # backward
                    g.add_edge(graph.Edge(
                        compound_nodes[c1], node, final_props(rx, edge1, edge2)))

                    edge1 = rxn_id, c2
                    edge2 = c2, rxn_id
                    g.add_edge(graph.Edge(
                        node, compound_nodes[c2], final_props(rx, edge1, edge2)))

        # add exchange reaction nodes
        rxn_set = set()
        for reaction in model.reactions:
            if model.is_exchange(reaction):
                if reaction in subset_reactions:
                    raw_exchange_rxn = model.get_reaction(reaction)
                    if element is not None:
                        if any(Atom(element) in cpd_formula[str(c[0].name)] for c in raw_exchange_rxn.compounds):
                            rxn_set.add(reaction)
                    else:
                        rxn_set.add(reaction)
        for r in rxn_set:
            exchange_rxn = model.get_reaction(r)
            node_ex = graph.Node({
                'id': r,
                'label': r,
                'shape': 'box',
                'style': 'filled',
                'fillcolor': ACTIVE_COLOR})
            g.add_node(node_ex)

            for c, _ in exchange_rxn.left:
                if c not in compound_nodes.keys():
                    node_ex_cpd = graph.Node({
                        'id': text_type(c),
                        'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
                        'shape': 'ellipse',
                        'style': 'filled',
                        'fillcolor':'#5a95f4'})
                    g.add_node(node_ex_cpd)
                    compound_nodes[c] = node_ex_cpd

            edge1 = r, c
            edge2 = c, r
            g.add_edge(graph.Edge(
                node_ex, compound_nodes[c], final_props(exchange_rxn, edge1, edge2)))

        if nativemodel.biomass_reaction is not None:
            bio_pair = Counter()
            biomass_cpd = set()
            biomass_reaction = model.get_reaction(nativemodel.biomass_reaction)
            if nativemodel.biomass_reaction in subset_reactions:
                for c, _ in biomass_reaction.left:
                    if element is not None:
                        if Atom(element) in cpd_formula[str(c.name)]:
                            bio_pair[nativemodel.biomass_reaction] += 1  # bio_pair = Counter({'biomass_rxn_id': 1})...
                            biomass_cpd.add(c)
                    else:
                        bio_pair[nativemodel.biomass_reaction] += 1
                        biomass_cpd.add(c)

                for c, _ in biomass_reaction.left:
                    bio_pair[nativemodel.biomass_reaction] += 1
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
                            'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
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

    def create_split_bipartite_graph_without_fpp(self, model, nativemodel, element, edge_values, cpd_formula, reaction_flux):
        """create bipartite graph on the model, without any selection"""

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

        subset_reactions = []
        if self._args.subset is not None:  # a file contains reaction_id in the subset users want to visualize
            with open(self._args.subset, 'r') as f:
                for line in f.readlines():
                    subset_reactions.append(line.rstrip())
        else:
            subset_reactions = set(model.reactions)

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

        # edge_props = {'dir': dir_value(rx.direction)}

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

        compound_nodes = {}
        cpd_set = set()
        rxn_set = set()
        for r_id in subset_reactions:
            if r_id != nativemodel.biomass_reaction and not model.is_exchange(r_id):
                rx_raw = model.get_reaction(r_id)
                for c, _ in rx_raw.compounds:
                    if element is not None:
                        if Atom(element) in cpd_formula[c.name]:
                            cpd_set.add(c)
                            rxn_set.add(r_id)
                    else:
                        cpd_set.add(c)
                        rxn_set.add(r_id)

        for rxn_id in rxn_set:
            rx = model.get_reaction(rxn_id)
            for c in cpd_set:
                node = graph.Node({
                    'id': text_type(c.name),
                    'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
                    'shape': 'ellipse',
                    'style': 'filled',
                    'fillcolor': color[c.name]})
                g.add_node(node)
                compound_nodes[c] = node

            node = graph.Node({
                'id': rxn_id,
                'label': rxns_properties(rxn_entry[rxn_id], self._args.detail, self._args.flux_analysis, reaction_flux),
                'shape': 'box',
                'style': 'filled',
                'fillcolor': color[rxn_id]})
            g.add_node(node)

            for c1, _ in rx.left:
                if c1 in cpd_set:
                    edge1 = c1, rxn_id  # forward
                    edge2 = rxn_id, c1  # backward
                    g.add_edge(graph.Edge(
                        compound_nodes[c1], node, final_props(rx, edge1, edge2)))

            for c2, _ in rx.right:
                if c2 in cpd_set:
                    edge1 = rxn_id, c2
                    edge2 = c2, rxn_id
                    g.add_edge(graph.Edge(
                        node, compound_nodes[c2], final_props(rx, edge1, edge2)))

        # add exchange reaction nodes
        EX_reactions = set()
        for r in subset_reactions:
            if model.is_exchange(r):
                ex_rxn = model.get_reaction(r)
                if element is not None:
                    if any(Atom(element) in cpd_formula[str(c[0].name)] for c in ex_rxn.compounds):
                        EX_reactions.add(r)
                else:
                    EX_reactions.add(r)
        for reaction in EX_reactions:
            exchange_rxn = model.get_reaction(reaction)
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
                        'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
                        'shape': 'ellipse',
                        'style': 'filled',
                        'fillcolor': '#5a95f4'})
                    g.add_node(node_ex)
                    compound_nodes[c] = node_ex
            # print(dir_value(exchange_rxn))

            edge1 = reaction, c
            edge2 = c, reaction
            g.add_edge(graph.Edge(
                node_ex, compound_nodes[c], final_props(exchange_rxn, edge1, edge2)))

        # add biomass reaction nodes
        bio_pair = Counter()
        biomass_cpds = set()
        if nativemodel.biomass_reaction is not None:
            if nativemodel.biomass_reaction in subset_reactions:
                biomass_reaction = model.get_reaction(nativemodel.biomass_reaction)
                for c, _ in biomass_reaction.left:
                    if element is not None:
                        if Atom(element) in cpd_formula[str(c.name)]:
                            biomass_cpds.add(c)
                    else:
                        biomass_cpds.add(c)
                for c in biomass_cpds:
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
                            'label': cpds_properties(c, cpd_entry[c.name], self._args.detail),
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