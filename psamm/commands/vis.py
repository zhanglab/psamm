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
from .tableexport import _encode_value
from .. import findprimarypairs
from ..command import (LoopRemovalMixin, ObjectiveMixin, SolverCommandMixin,
                       MetabolicMixin, Command, FilePrefixAppendAction)
from ..reaction import Direction
from ..formula import Formula, Atom, ParseError
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
                 'choices=[fpp, no-fpp, half-fpp or {file path}]')
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
            '--image-size', type=text_type, default=None,
            help='Set width and height of the graph image. '
                 '(width,height)(inches)')
        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""

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
                                compound.id, compound.formula))
                except ParseError as e:
                    msg = (
                        'Error parsing formula'
                        ' for compound {}:\n{}\n{}'.format(
                            compound.id, e, compound.formula))
                    if e.indicator is not None:
                        msg += '\n{}'.format(e.indicator)
                    logger.warning(msg)

        subset_reactions = make_subset(self._mm, self._args.subset)

        filter_dict = make_filter_dict(
            self._model, self._mm, self._args.method, self._args.element,
            compound_formula, self._args.hide_edges, self._args.exclude)

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

        cpair_dict, new_id_mapping = make_cpair_dict(
            self._mm, filter_dict, subset_reactions, reaction_flux,
            self._args.method, self._args.combine)
        print(len(filter_dict))
        print('length:', len(cpair_dict))

        # cpair_dict_for_table, new_id_mapping_for_table = \
        #     make_cpair_dict(self._mm, filter_dict_for_table,
        #                     subset_reactions, reaction_flux, 'fpp')

        edge_values = make_edge_values(
            reaction_flux, self._mm, compound_formula, self._args.element,
            self._args.combine, cpair_dict, new_id_mapping,
            self._args.method)
        # edge_values_for_table = make_edge_values(
        #     reaction_flux, self._mm, compound_formula, self._args.element,
        #     'True', cpair_dict_for_table, new_id_mapping_for_table,
        #     'fpp')

        model_compound_entries, model_reaction_entries = {}, {}
        for c in self._model.compounds:
            model_compound_entries[c.id] = c
        for r in self._model.reactions:
            model_reaction_entries[r.id] = r

        g = graph.Graph()
        # g_table = graph.Graph()
        g._default_node_props['fontname'] = 'Arial'
        g._default_node_props['fontsize'] = 12
        g = add_graph_nodes(g, cpair_dict, self._args.method,
                            new_id_mapping, self._args.combine)
        # g_table = add_graph_nodes(g_table, cpair_dict_for_table, 'fpp',
        #                           new_id_mapping_for_table, split=True)
        g = add_edges(g, cpair_dict, self._args.method, self._args.combine)
        # g_table = add_edges(g_table, cpair_dict_for_table, 'fpp',
        #                     split=True)

        bio_cpds_sub = set()
        bio_cpds_pro = set()
        if self._model.biomass_reaction in subset_reactions and \
                self._model.biomass_reaction not in self._args.exclude:
            biomass_rxn = self._mm.get_reaction(self._model.biomass_reaction)
            g = add_biomass_rxns(g, biomass_rxn, self._model.biomass_reaction)
            for cpd, _ in biomass_rxn.left:
                bio_cpds_sub.add(text_type(cpd))
            for cpd, _ in biomass_rxn.right:
                bio_cpds_pro.add(text_type(cpd))

        exchange_cpds = set()
        for reaction in self._mm.reactions:
            if self._mm.is_exchange(reaction):
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

        recolor_dict = {}
        if self._args.color is not None:
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    recolor_dict[row[0]] = row[1]

        if len(recolor_dict) > 0:
            g = update_node_color(g, recolor_dict)
            # g_table = update_node_color(g_table, recolor_dict)
        if self._args.cpd_detail is not None or self._args.rxn_detail \
                is not None or self._args.fba is True:
            g = update_node_label(
                g, self._args.cpd_detail, self._args.rxn_detail,
                model_compound_entries, model_reaction_entries, reaction_flux)
            # g_table = update_node_label(
            #     g_table, self._args.cpd_detail, self._args.rxn_detail,
            #     model_compound_entries, model_reaction_entries, reaction_flux)
        g = set_edge_props_withfba(g, edge_values)
        # g_table = set_edge_props_withfba(g_table, edge_values_for_table)

        # if self._args.method == 'no-fpp' and self._args.split_map is True:
        #     logger.warning('--split-map is not compatible with no-fpp option')
        if self._args.method == 'no-fpp' and self._args.combine != 0:
            logger.warning('--combine option is not compatible with no-fpp method')

        width = None
        height = None
        if self._args.image_size is not None:
            hw = self._args.image_size.split(',')
            width = hw[0]
            height = hw[1]

        if self._args.compartment:
            boundaries, extracellular = get_cpt_boundaries(self._model)
            boundary_tree, extracellular = make_cpt_tree(boundaries,
                                                         extracellular)
            with open('reactions_compartmentalized.dot', 'w') as f:
                g.write_graphviz_compartmentalized(
                    f, boundary_tree, extracellular, width, height)
        else:
            with open('reactions.dot', 'w') as f:
                g.write_graphviz(f, width, height)
        with open('reactions.nodes.tsv', 'w') as f:
            g.write_cytoscape_nodes(f)
        with open('reactions.edges.tsv', 'w') as f:
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
                    render('dot', self._args.Image, 'reactions.dot')


def make_subset(mm, subset_file):
    """create a collection of reaction IDs that need to be visualized.

    Args:
        mm: Metabolic model, class 'psamm.metabolicmodel.MetabolicModel'.
        subset_file: None or an open file containing a list of reactions
                     or compound ids.
    """
    if subset_file is None:
        return set(mm.reactions)
    else:
        cpd_set = set()
        rxn_set = set()
        for l in subset_file.readlines():
            id = l.rstrip()
            if mm.has_compound(id):
                cpd_set.add(id)
            elif mm.has_reaction(id):
                rxn_set.add(id)
            else:
                raise ValueError('{} was in subset file but is '
                                 'not a compound or reaction id'.format(id))
        if all(i > 0 for i in [len(cpd_set), len(rxn_set)]):
            raise ValueError('Subset file contains a mix of reactions and '
                             'compounds.')
        if len(cpd_set) > 0:
            for rx in mm.reactions:
                rxn = mm.get_reaction(rx)
                if any(str(c) in cpd_set for (c, _) in rxn.compounds):
                    rxn_set.add(rx)
        return rxn_set


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

            for c, val in rx.compounds:
                if element == 'none':
                    edge_values[c, reaction] = abs(flux * float(val))
                else:
                    if c.name in compound_formula:
                        if Atom(element) in compound_formula[c.name]:
                            edge_values[c, reaction] = abs(flux * float(val))
                    else:
                        edge_values[c, reaction] = abs(flux * float(val))

        rxn_set = set()
        for (c1, c2), rxns in iteritems(cpair_dict):
            for direction, rlist in iteritems(rxns):
                rlist_string = ','.join(new_id_mapping[r] for r in rlist)
                if any(new_id_mapping[r] in reaction_flux and abs(
                        reaction_flux[new_id_mapping[r]]) > 1e-9 for
                       r in rlist):
                    x_comb_c1, x_comb_c2 = 0, 0
                    for r in rlist:
                        real_r = new_id_mapping[r]
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


def make_filter_dict(model, mm, method, element, cpd_formula,
                     arg_hide_edges, exclude_rxns):
    """create a dictionary of reaction id(key) and a list of related
    compound pairs(value) through 4 different methods.

    Args:
        model: Native model, class 'psamm.datasource.native.NativeModel'.
        mm: Metabolic model, class 'psamm.metabolicmodel.MetabolicModel'.
        method: Command line argument, including 'fpp', 'no-fpp' and
            file path.
        element: Chemical symbol of an element, such as C(carbon),
            N(nitrogen), S(sulfur).
        cpd_formula: A dictionary that key is compound id (a string) and
            value is class 'psamm.formula.Formula'.
        arg_hide_edges: Command line argument, by default it is "None",
            if this argument is given in command line, then it will be a file
            that contains two columns(tab-separated) of compounds(in format
            of compound_id[compartment], such as atp[c], akg[c], each row
            represent two reactant/product pairs((c1, c2)and (c2,c1)).
        exclude_rxns: A list of reaction IDs, these reaction will be removed
            when visualize the network.
    """

    cpd_object = {}
    for cpd in mm.compounds:
        cpd_object[str(cpd)] = cpd

    hide_edges = []
    if arg_hide_edges is not None:
        for row in csv.reader(arg_hide_edges, delimiter=str('\t')):
            hide_edges.append((cpd_object[row[0]], cpd_object[row[1]]))
            hide_edges.append((cpd_object[row[1]], cpd_object[row[0]]))

    fpp_rxns, rxns_no_formula = set(), []
    for reaction in model.reactions:
        if (reaction.id in mm.reactions and
                reaction.id not in exclude_rxns):

            if any(c.name not in cpd_formula for c, _ in
                   reaction.equation.compounds):
                rxns_no_formula.append(reaction.id)
                continue

            fpp_rxns.add(reaction)

    if len(rxns_no_formula) > 0:
        logger.warning(
            '{} reactions have compounds with undefined formula. '
            'These reactions contain {}'.format(len(rxns_no_formula),
                                                rxns_no_formula))

    filter_dict = {}
    if method == 'fpp' or method == 'half-fpp':
        if len(fpp_rxns) == 0:
            raise ValueError(
                'All the reactions have compounds with undefined formula or '
                'have no reaction equation, fix them or '
                'add "--method no-fpp --element none" to the command')

        reaction_pairs = [(r.id, r.equation) for r in fpp_rxns
                          if r.id != model.biomass_reaction]
        element_weight = findprimarypairs.element_weight
        fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(
            reaction_pairs, cpd_formula, element_weight=element_weight)

        for rxn_id, fpp_pairs in iteritems(fpp_dict):
            compound_pairs = []
            for cpd_pair, transfer in iteritems(fpp_pairs[0]):
                if cpd_pair not in hide_edges:
                    if element == 'none':
                        compound_pairs.append(cpd_pair)
                    else:
                        if any(Atom(element) in k for k in transfer):
                            compound_pairs.append(cpd_pair)
            if len(compound_pairs) != 0:
                filter_dict[rxn_id] = compound_pairs

    elif method == 'no-fpp':
        for rxn_id in mm.reactions:
            if rxn_id != model.biomass_reaction and rxn_id not in \
                    exclude_rxns:
                rx = mm.get_reaction(rxn_id)
                cpairs = []
                for c1, _ in rx.left:
                    for c2, _ in rx.right:
                        if (c1, c2) not in hide_edges:
                            if element != 'none':
                                if c1.name in cpd_formula and c2.name in \
                                        cpd_formula:
                                    if Atom(element) in cpd_formula[c1.name] \
                                            and Atom(element) in \
                                            cpd_formula[c2.name]:
                                        cpairs.append((c1, c2))
                                else:
                                    cpairs.append((c1, c2))
                            else:
                                cpairs.append((c1, c2))
                if len(cpairs) > 0:
                    filter_dict[rxn_id] = cpairs
    else:
        try:
            with open(method, 'r') as f:
                cpair_list, rxn_list = [], []
                for row in csv.reader(f, delimiter=str(u'\t')):
                    if (cpd_object[row[1]], cpd_object[row[2]]) not in \
                            hide_edges:
                        if row[0] not in exclude_rxns:
                            if element == 'none':
                                cpair_list.append((cpd_object[row[1]],
                                                   cpd_object[row[2]]))
                                rxn_list.append(row[0])
                            else:
                                if Atom(element) in \
                                        Formula.parse(row[3]).flattened():
                                    cpair_list.append((cpd_object[row[1]],
                                                       cpd_object[row[2]]))
                                    rxn_list.append(row[0])

            filter_dict = defaultdict(list)
            for r, cpair in zip(rxn_list, cpair_list):
                filter_dict[r].append(cpair)
        except IOError:
            raise IOError('Invalid file path, no such file or '
                          'directory: {}' .format(method))

    cpairs_ordered_filter_dict = {}
    for r, cpairs in iteritems(filter_dict):
        cpairs_ordered_filter_dict[r] = sorted(cpairs)

    return cpairs_ordered_filter_dict


def make_mature_cpair_dict(cpair_dict):
    new_cpair_dict = {}
    cpair_list = []
    for (c1, c2), rxns in sorted(iteritems(cpair_dict)):
        if (c1, c2) not in cpair_list:
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
    for (c1, c2), rxns in iteritems(new_cpair_dict):
        for direction, rlist in iteritems(rxns):
            rxns_sorted_cpair_dict[(c1, c2)][direction] = sorted(rlist)

    return rxns_sorted_cpair_dict


def make_cpair_dict(mm, filter_dict, subset, reaction_flux, args_method, args_combine):
    """Create a mapping from compound pair to a defaultdict containing
    lists of reactions for the forward, reverse, and both directions.

    Args:
        mm: Class 'psamm.metabolicmodel.MetabolicModel'.
        filter_dict: A dictionary mapping form reaction ID to compound pairs
            that go through this reaction.
        subset: A list of reaction IDs, including all the reactions that
            will be visualized.
        reaction_flux: A dictionary mapping from reaction ID to
            reaction flux value.
        args_method: Command line argument, a string, including 'fpp',
            'no-fpp' and file path.
    """

    new_id_mapping = {}
    rxn_count = Counter()
    cpair_dict = defaultdict(lambda: defaultdict(list))

    if args_method != 'no-fpp':
        if args_combine == 0 or args_combine == 2:
            for rxn, cpairs in iteritems(filter_dict):
                if rxn in subset:
                    rx = mm.get_reaction(rxn)
                    have_visited = set()
                    sub_pro = defaultdict(list)  # substrate map to a list of product
                    rxn_mixcpairs = defaultdict(list)
                    for (c1, c2) in sorted(cpairs):
                        sub_pro[c1].append(c2)
                    for k1, v1 in iteritems(sub_pro):
                        if k1 not in have_visited:
                            rxn_count[rxn] += 1
                            have_visited.add(k1)
                            r_id = str('{}_{}'.format(rxn, rxn_count[rxn]))
                            new_id_mapping[r_id] = rxn
                            for v in v1:
                                rxn_mixcpairs[r_id].append((k1, v))
                            for k2, v2 in iteritems(sub_pro):
                                if k2 not in have_visited:
                                    if k2 != k1:
                                        if v1 == v2:
                                            have_visited.add(k2)
                                            for vtest in v2:
                                                rxn_mixcpairs[r_id].append((k2, vtest))
                    for rxn_id, cpairs in iteritems(rxn_mixcpairs):
                        for (c1, c2) in cpairs:
                            if rx.direction == Direction.Forward:
                                cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                            else:
                                if rxn in reaction_flux:
                                    if reaction_flux[rxn] > 0:
                                        cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                                    else:
                                        cpair_dict[(c1, c2)]['back'].append(rxn_id)
                                else:
                                    cpair_dict[(c1, c2)]['both'].append(rxn_id)

        elif args_combine == 1:
            for rxn, cpairs in iteritems(filter_dict):
                if rxn in subset:
                    rx = mm.get_reaction(rxn)
                    cpd_rid = {}
                    have_visited = set()
                    for (c1, c2) in sorted(cpairs):
                        if c1 not in have_visited:
                            if c2 not in have_visited:
                                rxn_count[rxn] += 1
                                rxn_id = str('{}_{}'.format(rxn, rxn_count[rxn]))
                                new_id_mapping[rxn_id] = rxn
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

                        if rx.direction == Direction.Forward:
                            cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                        else:
                            if rxn in reaction_flux:
                                if reaction_flux[rxn] > 0:
                                    cpair_dict[(c1, c2)]['forward'].append(rxn_id)
                                else:
                                    cpair_dict[(c1, c2)]['back'].append(rxn_id)
                            else:
                                cpair_dict[(c1, c2)]['both'].append(rxn_id)
    else:
        for rxn, cpairs in iteritems(filter_dict):
            if rxn in subset:
                rx = mm.get_reaction(rxn)
                for (c1, c2) in cpairs:
                    r_id = rxn
                    new_id_mapping[r_id] = rxn
                    if rx.direction == Direction.Forward:
                        cpair_dict[(c1, c2)]['forward'].append(r_id)
                    else:
                        if rxn in reaction_flux:
                            if reaction_flux[rxn] > 0:
                                cpair_dict[(c1, c2)]['forward'].append(r_id)
                            else:
                                cpair_dict[(c1, c2)]['back'].append(r_id)
                        else:
                            cpair_dict[(c1, c2)]['both'].append(r_id)

    rxns_sorted_cpair_dict = make_mature_cpair_dict(cpair_dict)

    return rxns_sorted_cpair_dict, new_id_mapping


def add_graph_nodes(g, cpairs_dict, method, new_id_mapping, args_combine):
    """create compound and reaction nodes, adding them to empty graph object.

    Args:
        g: An empty Graph object.
        cpairs_dict: defaultdict of compound_pair: defaultdict of direction:
            reaction list. e.g. {(c1, c2): {'forward":[rx1], 'both':[rx2}}.
        method: Command line argument, options=['fpp', 'no-fpp', file_path].
        new_id_mapping: dictionary of rxn_id_suffix: rxn_id.
        args_combine: Command line argument, could be 0, 1, 2. By default it is 0, .
    return: A graph object that contains a set of nodes.
    """
    graph_nodes = set()
    for cpair, reactions in iteritems(cpairs_dict):
        for c in cpair:
            if c not in graph_nodes:
                node = graph.Node({
                    'id': text_type(c),
                    'shape': 'ellipse',
                    'style': 'filled',
                    'label': text_type(c),
                    'type': 'cpd',
                    'fillcolor': COMPOUND_COLOR,
                    'original_id': c,
                    'compartment': c.compartment})
                g.add_node(node)
                graph_nodes.add(c)
        for direction, rlist in iteritems(reactions):
            if args_combine == 2 and method != 'no-fpp':
                real_rxns = [new_id_mapping[r] for r in rlist]
                rxn_string = text_type(','.join(rlist))
                if rxn_string not in graph_nodes:
                    rnode = graph.Node({
                        'id': text_type(','.join(rlist)),
                        'shape': 'box',
                        'style': 'filled',
                        'label': '\n'.join(real_rxns),
                        'type': 'rxn',
                        'fillcolor': REACTION_COLOR,
                        'original_id': real_rxns,
                        'compartment': c.compartment})
                    g.add_node(rnode)
                    graph_nodes.add(rxn_string)
            else:
                for sub_rxn in rlist:
                    if sub_rxn not in graph_nodes:
                        rnode = graph.Node({
                            'id': text_type(sub_rxn),
                            'shape': 'box',
                            'style': 'filled',
                            'label': new_id_mapping[sub_rxn],
                            'type': 'rxn',
                            'fillcolor': REACTION_COLOR,
                            'original_id': [new_id_mapping[sub_rxn]],
                            'compartment': c.compartment})
                        g.add_node(rnode)
                        graph_nodes.add(sub_rxn)

    return g


def add_edges(g, cpairs_dict, method, args_combine):
    """add edges to the graph object obtained in last step.

    Args:
        g: A graph object contains a set of nodes.
        cpairs_dict: A defaultdict of compound_pair: defaultdict of direction:
            reaction list. e.g. {(c1, c2): {'forward":[rx1], 'both':[rx2}}.
        method: Command line argument, options=['fpp', 'no-fpp', file_path].
        # split: Command line argument, True or False. By default split = False.
        args_combine: Command line argument, an integer(could be 0, 1 or 2).
    """

    edge_list = []
    for (c1, c2), value in iteritems(cpairs_dict):
        for direction, rlist in iteritems(value):
            new_rlist = ','.join(rlist)
            if args_combine == 0 or args_combine == 1 or method == 'no-fpp':
                for sub_rxn in rlist:
                    test1 = c1, sub_rxn
                    if test1 not in edge_list:
                        edge_list.append(test1)
                        g.add_edge(graph.Edge(
                            g.get_node(text_type(c1)),
                            g.get_node(text_type(sub_rxn)), {'dir': direction}))

                    test2 = c2, sub_rxn
                    if test2 not in edge_list:
                        edge_list.append(test2)
                        g.add_edge(graph.Edge(
                            g.get_node(text_type(sub_rxn)),
                            g.get_node(text_type(c2)), {'dir': direction}))
            else:
                test1 = c1, new_rlist
                test2 = new_rlist, c2
                if test1 not in edge_list:
                    edge_list.append(test1)
                    g.add_edge(graph.Edge(
                        g.get_node(text_type(c1)),
                        g.get_node(text_type(new_rlist)), {'dir': direction}))
                if test2 not in edge_list:
                    edge_list.append(test2)
                    g.add_edge(graph.Edge(
                        g.get_node(text_type(new_rlist)),
                        g.get_node(text_type(c2)), {'dir': direction}))
    return g


def dir_value(direction):
    """assign value to different reaction directions"""
    if direction == Direction.Forward:
        return 'forward'
    elif direction == Direction.Reverse:
        return 'back'
    else:
        return 'both'


def add_biomass_rxns(g, bio_reaction, biomass_rxn_id):
    """Adds biomass reaction nodes and edges to a graph object.

    Args:
        g: Graph object.
        bio_reaction: Biomass reaction ReactionEntry.
        biomass_rxn_id: Biomass reaction ID.
    """
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


def update_node_color(g, recolor_dict):
    """ Update node color in Graph object based on a mapping dictionary

    Args:
        g: A Graph object that contains nodes and edges.
        recolor_dict: dict of rxn_id/cpd_id[compartment] : hex color code.
    return: a graph object that contains a set of node with defined color.
    """
    for r_id in recolor_dict:
        if r_id in g.nodes_original_id_dict:
            for node in g.nodes_original_id_dict[r_id]:
                node.props['fillcolor'] = recolor_dict[r_id]
    return g


def update_node_label(g, cpd_detail, rxn_detail, model_compound_entries,
                      model_reaction_entries, reaction_flux):
    """ set label of nodes in graph object,

    Args:
        g: A graph object, contain a set of nodes.
        cpd_detail: Command line argument, a list that contains only one
            element, this element is a compound properties name list,
            e.g. detail = [['id', 'name', 'formula']].
        rxn_detail: Command line argument, a list that contains only one
            element, this element is a reaction properties name list,
            e.g. detail = [['id', genes', 'equation']].
        model_compound_entries: dict of cpd_id:compound_entry.
        model_reaction_entries: dict of rxn_id:reaction_entry.
        reaction_flux: dict of rxn_id: reaction flux value.
    """

    for node in g.nodes:
        if node.props['type'] == 'cpd':
            if cpd_detail is not None:
                cpd_id = node.props['original_id'].name
                props = model_compound_entries[cpd_id].properties
                cpd_detail_list = [i for i in cpd_detail[0] if i in props]
                pre_label = '\n'.join(((props[value].encode(
                    'ascii', 'backslashreplace')).decode('ascii')) for value
                                      in cpd_detail_list if value != 'id')
                if 'id' in cpd_detail[0]:
                    label = '{}\n{}'.format(node.props['id'], pre_label)
                else:
                    label = pre_label
                node.props['label'] = label

        elif node.props['type'] == 'rxn':
            if rxn_detail is not None:
                rxn_id = '\n'.join(node.props['original_id'])
                if len(node.props['original_id']) == 1:
                    props = model_reaction_entries[rxn_id].properties
                    rxn_detail_list = [i for i in rxn_detail[0] if i in props]
                    pre_label = '\n'.join(_encode_value(
                        props[value]) for value in rxn_detail_list
                                          if value != 'id')
                    if 'id' in rxn_detail[0]:
                        label = '{}\n{}'.format(rxn_id, pre_label)
                    else:
                        label = pre_label
                    node.props['label'] = label

            if len(reaction_flux) > 0:
                sum_flux = 0
                for r in node.props['original_id']:
                    if r in reaction_flux:
                        sum_flux += abs(reaction_flux[r])
                if sum_flux != 0:
                    node.props['label'] = '{}\n{}'.format(
                        node.props['label'], sum_flux)

        elif node.props['type'] == 'Ex_rxn':
            if len(reaction_flux) > 0:
                if node.props['original_id'][0] in reaction_flux:
                    node.props['label'] = '{}\n{}'.format(
                        node.props['label'], reaction_flux[
                            node.props['original_id'][0]])
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
            if edge.source.props['type'] == 'cpd':
                rxn_string = ','.join(edge.dest.props['original_id'])
                edge_test = edge.source.props['original_id'], rxn_string
            elif edge.dest.props['type'] == 'cpd':
                rxn_string = ','.join(edge.source.props['original_id'])
                edge_test = edge.dest.props['original_id'], rxn_string

            if edge_values.get(edge_test) is not None:
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
