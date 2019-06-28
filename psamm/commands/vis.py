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

import logging
import csv
import argparse
from collections import defaultdict, Counter
from six import text_type, iteritems, itervalues
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command, FilePrefixAppendAction)
from ..reaction import Direction
from ..formula import Atom
from .. import graph
from .. import fluxanalysis

try:
    from graphviz import render
except ImportError:
    render = None

logger = logging.getLogger(__name__)

REACTION_COLOR = '#c9fccd'
COMPOUND_COLOR = '#ffd8bf'
ACTIVE_COLOR = '#90f998'
ALT_COLOR = '#b3fcb8'


class VisualizationCommand(MetabolicMixin, ObjectiveMixin, SolverCommandMixin,
                           Command, LoopRemovalMixin, FilePrefixAppendAction):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method', type=text_type, default='fpp',
            help='Compound pair prediction method,'
                 'choices=[fpp, no-fpp]')
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help='Reaction(s) to exclude from metabolite pair '
                 'prediction and final visualization of the model'
                 '. (e.g. biomass reactions or biosynthesis reactions)')
        parser.add_argument(
            '--fba', action='store_true',
            help='Run and visualize l1min FBA simulation on network image')
        parser.add_argument(
            '--element', type=text_type, default='C',
            help='Element transfer to show on the graph (C, N, O, S).')
        parser.add_argument(
            '--cpd-detail', type=text_type, default=None, action='append',
            nargs='+', help='determine compound properties showed on nodes '
                            'label, e.g. formula, name, charge)')
        parser.add_argument(
            '--rxn-detail', type=text_type, default=None, action='append',
            nargs='+', help='determine reaction properties showed on nodes '
                            'label, e.g. genes, EC, equation)')
        parser.add_argument(
            '--subset', type=argparse.FileType('rU'), default=None,
            help='Subset of reactions to include in graph '
                 '(default=All reactions)')
        parser.add_argument(
            '--color', type=argparse.FileType('rU'), default=None, nargs='+',
            help='Designate non-default colors for nodes in the graph')
        parser.add_argument(
            '--Image', type=text_type, default=None,
            help='Generate image file in designated format '
                 '(e.g. pdf, png, eps)')
        parser.add_argument(
            '--hide-edges', type=argparse.FileType('rU'), default=None,
            help='Remove edges between specific compound pair from network ')
        parser.add_argument(
            '--combine', metavar='Combine Level', type=int, choices=range(3),
            default=0, help='Combined reaction nodes in different two levels.')
        parser.add_argument(
            '--compartment', action='store_true',
            help='Generate visualization of the network with compartments '
                 'shown.')
        parser.add_argument(
            '--output', type=text_type, help='set output name.')
        parser.add_argument(
            '--image-size', type=text_type, default=None,
            help='Set width and height of the graph image. '
                 '(width,height)(inches)')
        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""

        compound_formula = graph.get_compound_dict(self._model)
        print(compound_formula)

        vis_rxns = rxnset_for_vis(self._mm, self._args.subset, self._args.exclude)

        filter_dict = graph.make_network_dict(
            self._model, self._mm, self._args.method, self._args.element, self._args.exclude)

        reaction_flux = {}
        if self._args.fba is True:
            solver = self._get_solver()
            p = fluxanalysis.FluxBalanceProblem(self._mm, solver)
            try:
                p.maximize(self._get_objective())
            except fluxanalysis.FluxBalanceError as e:
                self.report_flux_balance_error(e)

            fluxes = {r: p.get_flux(r) for r in self._mm.reactions}

            flux_var = p.get_flux_var(self._get_objective())
            p.prob.add_linear_constraints(flux_var == p.get_flux(
                self._get_objective()))
            p.minimize_l1()

            count = 0
            for r_id in self._mm.reactions:
                flux = p.get_flux(r_id)
                if abs(flux - fluxes[r_id]) > 1e-6:
                    count += 1
                if abs(flux) > 1e-6:
                    reaction_flux[r_id] = flux

        model_compound_entries, model_reaction_entries = {}, {}
        for c in self._model.compounds:
            model_compound_entries[c.id] = c
        for r in self._model.reactions:
            model_reaction_entries[r.id] = r

        cpair_dict, new_id_mapping = graph.make_cpair_dict(filter_dict, reaction_flux,
            self._args.method, self._args.combine)

        g = graph.make_bipartite_graph_object(
            cpair_dict, new_id_mapping, self._args.method, self._args.combine,
            model_compound_entries)

        recolor_dict = {}
        if self._args.color is not None:
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    recolor_dict[row[0]] = row[1]
        g = add_node_props(g, recolor_dict)

        g = add_node_label(g, self._args.cpd_detail, self._args.rxn_detail)
        if len(reaction_flux) != 0:
            g = update_rxnnode_label(g, reaction_flux)

        bio_cpds_sub = set()
        bio_cpds_pro = set()
        if self._model.biomass_reaction in vis_rxns and \
                self._model.biomass_reaction not in self._args.exclude:
            nm_biomass_rxn = model_reaction_entries[self._model.biomass_reaction]
            g = add_biomass_rxns(g, nm_biomass_rxn, self._model.biomass_reaction)
            for cpd, _ in nm_biomass_rxn.equation.left:
                bio_cpds_sub.add(text_type(cpd))
            for cpd, _ in nm_biomass_rxn.equation.right:
                bio_cpds_pro.add(text_type(cpd))

        exchange_cpds = set()
        for reaction in self._mm.reactions:
            if self._mm.is_exchange(reaction):
                if reaction != self._model.biomass_reaction:
                    exchange_rxn = self._mm.get_reaction(reaction)
                    g = add_exchange_rxns(g, reaction, exchange_rxn)
                    for cpd, _ in exchange_rxn.compounds:
                        exchange_cpds.add(text_type(cpd))

        for node in g.nodes:
            if node.props['id'] in bio_cpds_sub:
                node.props['type'] = 'cpd_biomass_substrate'
            elif node.props['id'] in bio_cpds_pro:
                node.props['type'] = 'cpd_biomass_product'
            elif node.props['id'] in exchange_cpds:
                node.props['type'] = 'cpd_exchange'
            else:
                continue

        edge_values = make_edge_values(
            reaction_flux, self._mm, compound_formula, self._args.element,
            self._args.combine, cpair_dict, new_id_mapping,
            self._args.method)
        if len(edge_values) > 0:
            g = set_edge_props_withfba(g, edge_values)

        if self._args.method == 'no-fpp' and self._args.combine != 0:
            logger.warning('--combine option is not compatible with no-fpp method')

        width = None
        height = None
        if self._args.image_size is not None:
            hw = self._args.image_size.split(',')
            width = hw[0]
            height = hw[1]

        if self._args.output is not None:
            output = self._args.output
        else:
            output = 'reactions'

        if self._args.compartment:
            boundaries, extracellular = get_cpt_boundaries(self._model)
            boundary_tree, extracellular = make_cpt_tree(boundaries,
                                                         extracellular)
            with open('reactions_compartmentalized.dot', 'w') as f:
                g.write_graphviz_compartmentalized(
                    f, boundary_tree, extracellular, width, height)
        else:
            with open('{}.dot'.format(output), 'w') as f:
                g.write_graphviz(f, width, height)
        with open('{}.nodes.tsv'.format(output), 'w') as f:
            g.write_cytoscape_nodes(f)
        with open('{}.edges.tsv'.format(output), 'w') as f:
            g.write_cytoscape_edges(f)

        if self._args.Image is not None:
            if render is None:
                raise ImportError(
                    'Making an image directly requires the '
                    'graphviz python bindings and the graphviz program '
                    'to be installed ("pip install graphviz")')
            else:
                if len(g.nodes_id_dict) > 1000:
                    logger.info(
                        'This graph contains a large number of reactions, '
                        'graphs of this size may take a long time to '
                        'create'.format(len(filter_dict.keys())))
                if self._args.compartment:
                    render('dot', self._args.Image,
                           'reactions_compartmentalized.dot')
                else:
                    # render('dot', self._args.Image, 'reactions.dot')
                    render('dot', self._args.Image, '{}.dot'.format(output))


def rxnset_for_vis(mm, subset_file, exclude):
    """create a collection of reaction IDs that need to be visualized.

    Args:
        mm: Metabolic model, class 'psamm.metabolicmodel.MetabolicModel'.
        subset_file: None or an open file containing a list of reactions
                     or compound ids.
        exclude: rxns to exclude
    """
    all_cpds = set()
    for cpd in mm.compounds:
        all_cpds.add(cpd.name)
    if subset_file is None:
        if len(exclude) == 0:
            final_rxn_set = set(mm.reactions)
        else:
            final_rxn_set = set([rxn for rxn in mm.reactions if rxn not in exclude])
    else:
        final_rxn_set = set()
        cpd_set = set()
        rxn_set = set()
        for l in subset_file.readlines():
            value = l.strip()
            if value in all_cpds:
                cpd_set.add(value)
            elif mm.has_reaction(value):
                rxn_set.add(value)
            else:
                raise ValueError('{} is in subset file but is '
                                 'not a compound or reaction id'.format(value))

        if all(i > 0 for i in [len(cpd_set), len(rxn_set)]):
            raise ValueError('Subset file is a mixture of reactions and '
                             'compounds.')
        else:
            if len(cpd_set) > 0:
                for rx in mm.reactions:
                    rxn = mm.get_reaction(rx)
                    if any(str(c.name) in cpd_set for (c, _) in rxn.compounds):
                        final_rxn_set.add(rx)
            elif len(rxn_set) > 0:
                final_rxn_set = rxn_set

    return final_rxn_set


def make_edge_values(reaction_flux, mm, compound_formula, element, args_combine,
                     cpair_dict, new_id_mapping, method):
    """set edge_values according to reaction fluxes

    Start from reaction equation and reaction flux to create a dictionary
    where the key is an edge(a tuple of (c, r) or (r, c)) and the value is
    a flux between the compound and reaction

    Args:
        reaction_flux: dictionary of reaction id and flux value.
        mm: Metabolic model, class 'psamm.metabolicmodel.MetabolicModel'.
        compound_formula: A dictionary of compound id and
            class 'psamm.formula.Formula'.
        element: Chemical symbol of an element, such as C(carbon),
            N(nitrogen), S(sulfur).
        args_combine: Interger, 0, 1 or 2. Decide whether to combine the
            reaction nodes in different level (2 levels allowed) to get other
            perspectives.
        cpair_dict: A defaultdict of compound pair and another defaultdict of
            list(in the latter defaultdict, key is direction(forward, back or
            both), value is reaction list),{(c1, c2): {'forward':[rxns],...}}.
        new_id_mapping: dictionary of reaction ID with suffix(number) and
            real reaction id, e.g.{r_1:r, r_2: r,...}.
        method: Command line argument, including 3 options:'fpp', 'no-fpp'
            or a file path.
        """

    edge_values = {}
    edge_values_combined = {}
    if len(reaction_flux) > 0:
        for reaction in reaction_flux:
            rx = mm.get_reaction(reaction)
            flux = reaction_flux[reaction]
            if abs(flux) < 1e-9:
                continue

            for cpd, val in rx.compounds:
                c = text_type(cpd)
                if element == 'none':
                    edge_values[c, reaction] = abs(flux * float(val))
                else:
                    if cpd.name in compound_formula:
                        if Atom(element) in compound_formula[cpd.name]:
                            edge_values[c, reaction] = abs(flux * float(val))
                    else:
                        edge_values[c, reaction] = abs(flux * float(val))

        rxn_set = set()
        for (cpd1, cpd2), rxns in iteritems(cpair_dict):
            c1 = text_type(cpd1)
            c2 = text_type(cpd2)
            for direction, rlist in iteritems(rxns):
                rlist_string = ','.join(new_id_mapping[r].id for r in rlist)
                if any(new_id_mapping[r].id in reaction_flux and abs(
                        reaction_flux[new_id_mapping[r].id]) > 1e-9 for
                       r in rlist):
                    x_comb_c1, x_comb_c2 = 0, 0
                    for r in rlist:
                        real_r = new_id_mapping[r].id
                        rxn_set.add(real_r)
                        if real_r in reaction_flux and abs(reaction_flux[real_r]) > 1e-9:
                            x_comb_c1 += edge_values[(c1, real_r)]
                            x_comb_c2 += edge_values[(c2, real_r)]
                    edge_values_combined[(c1, rlist_string)] = x_comb_c1
                    edge_values_combined[(c2, rlist_string)] = x_comb_c2

        for r in reaction_flux:
            if r not in rxn_set:
                for (cpd, rxn) in edge_values:
                    if rxn == r:
                        edge_values_combined[(cpd, r)] = \
                            edge_values[(cpd, rxn)]

    if args_combine == 0 or args_combine == 1 or method == 'no-fpp':
        return edge_values
    else:
        return edge_values_combined


def dir_value(direction):
    """assign value to different reaction directions"""
    if direction == Direction.Forward:
        return 'forward'
    elif direction == Direction.Reverse:
        return 'back'
    else:
        return 'both'


def add_biomass_rxns(g, nm_bio_reaction, biomass_rxn_id):
    """Adds biomass reaction nodes and edges to a graph object.

    Args:
        g: Graph object.
        nm_bio_reaction: Biomass reaction DictReactionEntry.
        biomass_rxn_id: Biomass reaction ID.
    """
    bio_reaction = nm_bio_reaction.equation
    direction = dir_value(bio_reaction.direction)
    reactant_list, product_list = [], []
    for c, _ in bio_reaction.left:
        reactant_list.append(c)
    for c, _ in bio_reaction.right:
        product_list.append(c)
    bio_pair = Counter()
    for c, _ in bio_reaction.compounds:
        if text_type(c) in g.nodes_id_dict:
            bio_pair[biomass_rxn_id] += 1
            node_bio = graph.Node({
                'id': '{}_{}'.format(biomass_rxn_id,
                                     bio_pair[biomass_rxn_id]),
                'entry': [nm_bio_reaction],
                'shape': 'box',
                'style': 'filled',
                'label': biomass_rxn_id,
                'type': 'bio_rxn',
                'fillcolor': ALT_COLOR,
                'original_id': [biomass_rxn_id],
                'compartment': c.compartment})
            g.add_node(node_bio)

            if c in reactant_list:
                g.add_edge(graph.Edge(g.get_node(text_type(c)), node_bio,
                                      {'dir': direction}))
            if c in product_list:
                g.add_edge(graph.Edge(node_bio, g.get_node(text_type(c)),
                                      {'dir': direction}))
    return g


def add_exchange_rxns(g, rxn_id, reaction):
    """Add exchange reaction nodes and edges to graph object.

    Args:
        g: A graph object that contains a set of nodes and some edges.
        rxn_id: Exchange reaction id,
        reaction: Exchange reaction object(metabolic model reaction),
            class 'psamm.reaction.Reaction'.
    """
    for c, _ in reaction.compounds:
        if text_type(c) in g.nodes_id_dict:
            node_ex = graph.Node({
                'id': text_type(rxn_id),
                'entry': [reaction],
                'shape': 'box',
                'style': 'filled',
                'label': rxn_id,
                'type': 'Ex_rxn',
                'fillcolor': ACTIVE_COLOR,
                'original_id': [rxn_id],
                'compartment': c.compartment})
            g.add_node(node_ex)

            direction = dir_value(reaction.direction)
            for c1, _ in reaction.left:
                g.add_edge(graph.Edge(
                    g.get_node(text_type(c1)), node_ex, {'dir': direction}))
            for c2, _ in reaction.right:
                g.add_edge(graph.Edge(
                    node_ex, g.get_node(text_type(c2)), {'dir': direction}))
    return g


def add_node_props(g, recolor_dict):
    """ Update node color in Graph object based on a mapping dictionary

    Args:
        g: A Graph object that contains nodes and edges.
        recolor_dict: dict of rxn_id/cpd_id[compartment] : hex color code.
    return: a graph object that contains a set of node with defined color.
    """
    cpd_types = ['cpd', 'cpd_biomass_substrate', 'cpd_biomass_product', 'cpd_exchange']
    rxn_types = ['rxn', 'bio_rxn', 'Ex_rxn']
    for node in g.nodes:
        node.props['style'] = 'filled'
        if node.props['type'] in cpd_types:
            node.props['fillcolor'] = COMPOUND_COLOR
            node.props['shape'] = 'ellipse'
        elif node.props['type'] == 'rxn':
            node.props['fillcolor'] = REACTION_COLOR
            node.props['shape'] = 'box'
        else:
            if node.props['type'] not in ['bio_rxn', 'Ex_rxn']:
                print('invalid nodes', type(node.props['entry']), node.props['entry'])

    for r_id in recolor_dict:
        if r_id in g.nodes_original_id_dict:
            for node in g.nodes_original_id_dict[r_id]:
                node.props['fillcolor'] = recolor_dict[r_id]
    return g


def add_node_label(g, cpd_detail, rxn_detail):
    """ set label of nodes in graph object,

    Args:
        g: A graph object, contain a set of nodes and a dictionary of edges.
        cpd_detail: Command line argument, a list that contains only one
            element, this element is a compound properties name list,
            e.g. detail = [['id', 'name', 'formula']].
        rxn_detail: Command line argument, a list that contains only one
            element, this element is a reaction properties name list,
            e.g. detail = [['id', genes', 'equation']].
    """
    print(cpd_detail, rxn_detail)
    for node in g.nodes:
        node.props['label'] = '\n'.join(obj.id for obj in node.props['entry'])

        # update node label based on what users provide in command line
        if cpd_detail is not None:
            if node.props['type'] == 'cpd':
                pre_label = '\n'.join(((node.props['entry'][0].properties.get(value).encode(
                'ascii', 'backslashreplace')).decode('ascii')) for value
                                  in cpd_detail[0] if value != 'id')
                if 'id' in cpd_detail[0]:
                    label = '{}\n{}'.format(node.props['entry'][0].properties['id'], pre_label)
                else:
                    label = pre_label
                node.props['label'] = label
        if rxn_detail is not None:
            if node.props['type'] == 'rxn':
                if len(node.props['entry']) == 1:
                    pre_label = '\n'.join((str(node.props['entry'][0].properties.get(value)).encode(
                    'ascii', 'backslashreplace').decode('ascii')) for value
                                          in rxn_detail[0] if value != 'id')
                    if 'id' in rxn_detail[0]:
                        label = '{}\n{}'.format(node.props['entry'][0].properties['id'], pre_label)
                    else:
                        label = pre_label
                    node.props['label'] = label
    return g


def update_rxnnode_label(g, reaction_flux):
    """
    add rxn flux value to the rxn node when --fba is given in command line.
    :param g:
    :param reaction_flux: dict of rxn_id: reaction flux value. Length is not 0.
    :return: an updated graph object, flux values are shown on rxn nodes.
    """
    rxn_types = ['rxn', 'bio_rxn', 'Ex_rxn']
    for node in g.nodes:
        if node.props['type'] in rxn_types:
            if type(node.props['entry']) is list:
                sum_flux = 0
                for r in node.props['entry']:
                    sum_flux += abs(reaction_flux.get(r.id, 0))
            else:
                sum_flux = abs(reaction_flux.get(node.props['entry'].id, 0))
            if sum_flux != 0:
                node.props['label'] = '{}\n{}'.format(
                    node.props['label'], sum_flux)
    return g


def set_edge_props_withfba(g, edge_values):
    """set edge values, including direction, width, style.

    Args:
        g: A graph object contains nodes and edges.
        edge_values: dictionary of (cpd, rxn): flux value.
    return: a complete graph object.
    """

    if len(edge_values) > 0:
        value_list = sorted(edge_values.values())
        ninety_percentile = value_list[int(len(value_list) * 0.9) - 1]
        min_edge_value = min(itervalues(edge_values))
        max_edge_value = ninety_percentile
    else:
        min_edge_value = 1
        max_edge_value = 1

    def pen_width(value):
        """calculate final edges width"""
        if max_edge_value - min_edge_value == 0:
            return 1
        else:
            if value > max_edge_value:
                value = max_edge_value
            alpha = value / max_edge_value
            beta = float("{:.4f}".format(alpha))

            return 10 * beta

    if len(edge_values) > 0:
        for edge in g.edges:
            # print(edge.source.props['id'], edge.dest.props['id'])
            cpd_types = ['cpd', 'cpd_biomass_substrate', 'cpd_biomass_product', 'cpd_exchange']
            if edge.source.props['type'] in cpd_types:
                rxn_string = ','.join([r.id if edge.dest.props['type'] != 'Ex_rxn' else edge.dest.props['id'] for
                                       r in edge.dest.props['entry']])
                edge_test = (edge.source.props['id'], rxn_string)

            else:
                rxn_string = ','.join([r.id if edge.source.props['type'] != 'Ex_rxn' else edge.source.props['id'] for
                                       r in edge.source.props['entry']])
                edge_test = (edge.dest.props['id'], rxn_string)

            if edge_test in edge_values:
                edge.props['penwidth'] = pen_width(edge_values[edge_test])
            else:
                edge.props['style'] = 'dotted'
    return g


def make_cpt_tree(boundaries, extracellular):
    """ This function will create a tree-like dictionary that can be used to
    determine compartment location.

        This function will take a list of compartment boundary tuples
        (eg: [(c, e), (c, p)]) and use this data to construct a tree like
        dictionary of parent compartments to lists of compartments they are
        adjacent to. An extracellular compartment will also be set to use as
        a starting point for the 'outermost' compartment in the model. If
        none is explicitly defined then 'e' will be used by default.
        If an 'e' compartment is not in the model either then a random
        compartment will be used as a starting point.

        args:
        boundaries: a list of tuples of adjacent compartment ids.
        extracellular: the extracellular compartment in the model.
    """
    children = defaultdict(set)
    compartments = set()
    for (j, k) in boundaries:
        compartments.add(j)
        compartments.add(k)
    if extracellular not in compartments:
        etmp = sorted(list(compartments), reverse=True)
        extracellular = etmp[0]
        logger.warning('No extracellular compartment was defined in the '
                       'model.yaml file and no "e" compartment in the model. '
                       'Trying to use {} as the extracellular compartment.'
                       .format(etmp[0]))
    for cpt in compartments:
        for (j, k) in boundaries:
            j = str(j)
            k = str(k)
            if j == cpt:
                children[cpt].add(k)
            elif k == cpt:
                children[cpt].add(j)
    return children, extracellular


def get_cpt_boundaries(model):
    """This function will determine the compartment boundaries in a model

    This function will take a native model object and determine the
    compartment boundaries either based on the predefined compartments in
    the model.yaml or based on the reaction equations in the model.

    args:
    model: Native model, class 'psamm.datasource.native.NativeModel'.
    """
    if model.extracellular_compartment is not None:
        extracellular = model.extracellular_compartment
    else:
        extracellular = 'e'

    if len(model.compartment_boundaries) != 0:
        boundaries = model.compartment_boundaries
    else:
        boundaries = set()
        for rxn in model.reactions:
            cpd_cpt = set()
            for cpd in rxn.equation.compounds:
                cpd_cpt.add(cpd[0].compartment)
            if len(cpd_cpt) > 1:
                cpd_cpt = list(cpd_cpt)
                boundaries.add(tuple(sorted((cpd_cpt[0], cpd_cpt[1]))))
    return boundaries, extracellular
