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

import unittest

logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#fbb4ae'
ACTIVE_COLOR = '#b3cde3'    # exchange reaction  color
ALT_COLOR = '#f4fc55'       # biomass reaction color
RXN_COMBINED_COLOR = '#fc9a44'
CPD_ONLY_IN_BIO = '#82e593'
CPD_ONLY_IN_EXC = '#5a95f4'

from psamm.reaction import Compound


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
            if reaction.id in reaction_flux:
                label = '{}\n{}'.format(reaction.id, reaction_flux[reaction.id])
            else:
                label = reaction.id
        else:
            label = reaction.id

    return label


def primary_element(element):
    """allow both lower and upper case for one-letter element """
    if element != 'none':
        if element in ['c', 'h', 'o', 'n', 'p', 's', 'k', 'b', 'f', 'v', 'y', 'i', 'w']:
            return element.upper()
        else:
            return element


def make_edge_values(reaction_flux, mm, compound_formula, element, split_map, cpair_dict,
                     new_id_mapping, method):
    """set edge_values according to reaction fluxes"""
    edge_values = {}
    if len(reaction_flux) > 0:
        for reaction in mm.reactions:
            rx = mm.get_reaction(reaction)
            if reaction in reaction_flux:
                flux = reaction_flux[reaction]
                if abs(flux) < 1e-9:
                    continue

                if flux > 0:
                    for compound, value in rx.right:  # value=stoichiometry
                        if element == 'none':
                            edge_values[reaction, compound] = (flux * float(value))
                        else:
                            if Atom(element) in compound_formula[compound.name]:
                                edge_values[reaction, compound] = (flux * float(value))
                    for compound, value in rx.left:
                        if element == 'none':
                            edge_values[compound, reaction] = (flux * float(value))
                        else:
                            if Atom(element) in compound_formula[compound.name]:
                                edge_values[compound, reaction] = (flux * float(value))
                else:
                    for compound, value in rx.left:
                        if element == 'none':
                            edge_values[reaction, compound] = (- flux * float(value))
                        else:
                            if Atom(element) in compound_formula[compound.name]:
                                edge_values[reaction, compound] = (- flux * float(value))
                    for compound, value in rx.right:
                        if element == 'none':
                            edge_values[compound, reaction] = (- flux * float(value))
                        else:
                            if Atom(element) in compound_formula[compound.name]:
                                edge_values[compound, reaction] = (- flux * float(value))
        # for (a, b), v in iteritems(edge_values):
        #     print('{}\t{}\t{}'.format(a,b,v))

        if split_map is not True and method != 'no-fpp':
            remove_edges = set()
            for (c1, c2), rxns in iteritems(cpair_dict):
                if len(rxns['forward']) > 1:
                    if all(new_id_mapping[r] not in reaction_flux for r in rxns['forward']):
                        continue
                    else:
                    # for r in rxns.fo
                    # if any(new_id_mapping[r] in reaction_flux for r in rxns.forward):
                        x_forward_c1, x_forward_c2 = 0, 0
                        for r in rxns['forward']:
                            real_r = new_id_mapping[r]
                            if real_r in reaction_flux:
                                x_forward_c1 += edge_values[(c1, real_r)]
                                x_forward_c2 += edge_values[(real_r, c2)]
                        edge_values[(c1, tuple(rxns['forward']))] = x_forward_c1
                        edge_values[(tuple(rxns['forward']), c2)] = x_forward_c2
                if len(rxns['back']) > 1:
                    if all(new_id_mapping[r] not in reaction_flux for r in rxns['back']):
                        continue
                    else:
                        x_back_c1, x_back_c2 = 0, 0
                        for r in rxns['back']:
                            real_r = new_id_mapping[r]
                            if real_r in reaction_flux:
                                x_back_c1 += edge_values[(real_r, c1)]
                                x_back_c2 += edge_values[(c2, real_r)]
                        edge_values[(tuple(rxns['back']), c1)] = x_back_c1
                        edge_values[(c2, tuple(rxns['back']))] = x_back_c2

    # for (a, b), v in iteritems(edge_values):
    #     print('{}\t{}\t{}'.format(a, b, v))

                    # if split_map is not True:
    #     for (c1, c2), rxns in cpair_dict:  # rxns=[ [forward_rxns], [back_rxns], [bidir_rxns] ]
    #         for rlist in rxns:
    #             if len(rlist) > 1:
    #             x_for = 0
    #             for j in forward:
    #                 y = edge_values.get(c1, j)
    #                 if y is not None:
    #                     x_for += y
    #         if len(fowrard) > 1:
    #             x_rev = 0
    #             for j in forward:
    #                 y_2 = edge_values.get(c1, j)
    #                 if y_2 is not None:
    #                     x_rev += y_2
    #         if len(bid) > 1:
    #             x_bid = 0
    #             for j in bid:
    #                 y_3 = edge_values.get(c1, j)
    #                 if y_3 < 0:
    #                     x_rev += y_3
    #
    #                 if y_3 > 0:
    #                     x_for += y_3

    return edge_values


def make_filter_dict(model, mm, method, element, compound_formula, cpd_object, exclude_cpairs, exclude_rxns,
                     subset_reactions):
    """create a dictionary of reaction id(key) and a list of related compound pairs(value)"""
    filter_dict = {}
    if method == 'fpp':
        # def iter_reactions():
        #     """yield reactions that can applied to fpp"""
        fpp_rxns, rxns_no_equation, rxns_no_formula = set(), set(), []
        for reaction in model.reactions:
            if (reaction.id in model.model and
                    reaction.id not in exclude_rxns):
                if reaction.equation is None:
                    rxns_no_equation.add(reaction.id)
                    continue

                if any(c.name not in compound_formula for c, _ in reaction.equation.compounds):
                    rxns_no_formula.append(reaction.id)
                    continue

                fpp_rxns.add(reaction)

        if len(rxns_no_equation) > 0:
            logger.warning(
                '{} reactions have no reaction equation, fix or exclude them.'
                'These reactions contain {}'.format(len(rxns_no_equation), rxns_no_equation))

        if len(rxns_no_formula) > 0:
            logger.warning(
                '{} reactions have compounds with undefined formula'
                # ', fix them or add "--method no-fpp --element none"to '
                'command. These reactions contain {}'.format(len(rxns_no_formula), rxns_no_formula))

        if len(fpp_rxns) == 0:
            logger.warning(
                'All the reactions have compounds with undefined formula or have no reaction equation, fix them or '
                'add "--method no-fpp --element none" to command')
            quit()

        split_reaction = True
        reaction_pairs = [(r.id, r.equation) for r in fpp_rxns]
        element_weight = findprimarypairs.element_weight
        fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(reaction_pairs, compound_formula,
                                                                       element_weight=element_weight)
        # print(fpp_dict)

        for rxn_id, fpp_pairs in iteritems(fpp_dict):
            compound_pairs = []
            for cpd_pair, transfer in iteritems(fpp_pairs[0]):
                if cpd_pair not in exclude_cpairs:
                    if element == 'none':
                        compound_pairs.append(cpd_pair)
                    else:
                        if any(Atom(primary_element(element)) in k for k in transfer):
                            compound_pairs.append(cpd_pair)
            if len(compound_pairs) != 0:
                filter_dict[rxn_id] = compound_pairs

        for r, cpairs in iteritems(filter_dict):
            print(r, cpairs)

    elif method == 'no-fpp':
        split_reaction = False
        for rxn_id in mm.reactions:
            if rxn_id != model.biomass_reaction:
                rx = mm.get_reaction(rxn_id)
                cpairs = []
                for c1, _ in rx.left:
                    for c2, _ in rx.right:
                        if (c1, c2) not in exclude_cpairs:
                            if element != 'none':
                                if Atom(primary_element(element)) in compound_formula[c1.name]:
                                    if Atom(primary_element(element)) in compound_formula[c2.name]:
                                        cpairs.append((c1, c2))
                            else:
                                cpairs.append((c1, c2))
                if len(cpairs) > 0:
                    filter_dict[rxn_id] = cpairs

        for r, cpairs in iteritems(filter_dict):
            print('{}\t{}'.format(r, cpairs))
    else:
        split_reaction = True
        try:
            with open(method, 'r') as f:
                cpair_list, rxn_list = [], []
                for row in csv.reader(f, delimiter=str(u'\t')):
                    if (cpd_object[row[1]], cpd_object[row[2]]) not in exclude_cpairs:
                        if element == 'none':
                            cpair_list.append((cpd_object[row[1]], cpd_object[row[2]]))
                            rxn_list.append(row[0])
                        else:
                            if Atom(primary_element(element)) in Formula.parse(row[3]).flattened():
                                cpair_list.append((cpd_object[row[1]], cpd_object[row[2]]))
                                rxn_list.append(row[0])

                filter_dict = defaultdict(list)
                for r, cpair in zip(rxn_list, cpair_list):
                    filter_dict[r].append(cpair)

            for r, cpairs in iteritems(filter_dict):
                print(r)
                # print('{}\t{}'.format(r, cpairs))

        except:
            if IOError:
                logger.error(' Invalid file path, no such file or directory : {}' .format(method))
            quit()

    return filter_dict, split_reaction


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
            help='Reaction to exclude (e.g. biomass reactions or macromolecule synthesis)')
        parser.add_argument(
            '--fba', action='store_true',
            help='visualize reaction flux')
        parser.add_argument(
            '--element', type=text_type, default='C',
            help='Primary element flow')
        parser.add_argument(
            '--detail', type=text_type, default=None, action='append', nargs='+',
            help='Reaction and compound properties showed on nodes label')
        parser.add_argument(
            '--subset', type=argparse.FileType('r'), default=None,
            help='Reactions designated to visualize')
        parser.add_argument(
            '--color', type=argparse.FileType('r'), default=None, nargs='+',
            help='Customize color of reaction and compound nodes ')
        parser.add_argument(
            '--Image', type=text_type, default=None, help='generate image file directly')
        parser.add_argument(
            '--exclude-pairs', type=argparse.FileType('r'), default=None,
            help='Remove edge of given compound pairs from final graph ')
        parser.add_argument(
            '--split-map', action='store_true',
            help='Create dot file for reaction-splitted metabolic network visualization')
        super(VisualizationCommand, cls).init_parser(parser)

    def run(self):
        """Run visualization command."""

        # Mapping from compound id to formula
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

        # Mapping from string of cpd_id+compartment(eg: pyr_c[c]) to Compound object
        cpd_object = {}
        for cpd in self._mm.compounds:
            cpd_object[str(cpd)] = cpd

        # read exclude_compound_pairs from command-line argument
        exclude_cpairs = []
        if self._args.exclude_pairs is not None:
            for row in csv.reader(self._args.exclude_pairs, delimiter=str('\t')):
                exclude_cpairs.append((cpd_object[row[0]], cpd_object[row[1]]))
                exclude_cpairs.append((cpd_object[row[1]], cpd_object[row[0]]))

        # set of reactions to visualize
        if self._args.subset is not None:
            raw_subset, subset_reactions, mm_cpds = [], set(), []
            for line in self._args.subset.readlines():
                raw_subset.append(line.rstrip())
            for c in self._mm.compounds:
                mm_cpds.append(str(c))
            if set(raw_subset).issubset(set(self._mm.reactions)):
                subset_reactions = raw_subset
            elif set(raw_subset).issubset(set(mm_cpds)):
                for reaction in self._mm.reactions:
                    rx = self._mm.get_reaction(reaction)
                    if any(str(c) in raw_subset for (c, _) in rx.compounds):
                        subset_reactions.add(reaction)
            else:
                logger.warning('Invalid subset file. The file should contain a column of reaction id or a column '
                               'of compound id with compartment, mix of reactions, compounds and other infomation '
                               'in one subset file is not allowed. The function will generate entire metabolic '
                               'network of the model')
                subset_reactions = set(self._mm.reactions)
        else:
            subset_reactions = set(self._mm.reactions)

        # create {rxn_id:[(c1, c2),(c3,c4),...], ...} dictionary, key = rxn id, value = list of compound pairs
        filter_dict, split_reaction = make_filter_dict(self._model, self._mm, self._args.method, self._args.element,
                                                           compound_formula, cpd_object, exclude_cpairs,
                                                           self._args.exclude, subset_reactions)

        # run l1min_fba, get reaction fluxes
        reaction_flux = {}
        if self._args.fba is True:
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

        # edge_values = make_edge_values(reaction_flux, self._mm, compound_formula, self._args.element,
        #                                self._args.split_map, cpairs_dict)


        # create {(c1, c2):[[forward rxns], [back rxns], [bidir rxns]], ...} dictionary, key=cpd_pair, value=rxn list
        raw_cpair_dict = defaultdict(list)     # key=compound pair, value=list of reaction_id
        raw_dict = {k: v for k, v in iteritems(filter_dict) if k in subset_reactions}
        for rxn_id, cpairs in iteritems(raw_dict):
            for pair in cpairs:
                raw_cpair_dict[pair].append(rxn_id)

        # for (c1, c2), v in iteritems(raw_cpair_dict):
        #     print('{}\t{}\t{}'.format(c1, c2, v))

        cpair_dict = {}
        pair_list = []
        rxn_count = Counter()
        new_id_mapping = {} # mapping from reactiion id with suffix (such as 'PGK_2) to reaction id (PGK)
        for (c1, c2), rxns in iteritems(raw_cpair_dict):
            if (c1, c2) not in pair_list:
                forward_rxns, back_rxns, bidirectional_rxns = [], [], []
                for r in rxns:
                    if self._args.method == 'no-fpp':
                        r_id = r
                    else:
                        rxn_count[r] += 1
                        r_id = str('{}_{}'.format(r, rxn_count[r]))
                    new_id_mapping[r_id] = r

                    reaction = self._mm.get_reaction(r)
                    if reaction.direction == Direction.Forward:
                        forward_rxns.append(r_id)
                    elif reaction.direction == Direction.Reverse:
                        back_rxns.append(r_id)
                    else:
                        if self._args.fba is True:
                            a = reaction_flux.get(r)
                            if a is not None:
                                if a > 0:
                                    forward_rxns.append(r_id)
                                else:
                                    back_rxns.append(r_id)
                            else:
                                bidirectional_rxns.append(r_id)
                        else:
                            bidirectional_rxns.append(r_id)
                    # rxn_list.append(r_id)

                if (c2, c1) in raw_cpair_dict:
                    for rx in raw_cpair_dict[(c2, c1)]:
                        if self._args.method == 'no-fpp':
                            rx_id = rx
                        else:
                            rxn_count[rx] += 1
                            # rxn_list.append(r)
                            rx_id = str('{}_{}'.format(rx, rxn_count[rx]))
                        new_id_mapping[rx_id] = rx
                        reaction = self._mm.get_reaction(rx)
                        if reaction.direction == Direction.Forward:
                            back_rxns.append(rx_id)
                        elif reaction.direction == Direction.Reverse:
                            forward_rxns.append(rx_id)
                        else:
                            if self._args.fba is True:
                                a = reaction_flux.get(rx)
                                if a is not None:
                                    if a > 0:
                                        back_rxns.append(rx_id)
                                    else:
                                        forward_rxns.append(rx_id)
                                else:
                                    bidirectional_rxns.append(rx_id)
                            else:
                                bidirectional_rxns.append(rx_id)

                # cpair_rxn = namedtuple('cpair_rxn', ['forward', 'back', 'bidirection'])
                cpair_rxn = {'forward': forward_rxns, 'back': back_rxns, 'bidirectional': bidirectional_rxns}
                cpair_dict[(c1, c2)] = cpair_rxn
            pair_list.append((c1, c2))
            pair_list.append((c2, c1))
        # for k, v in iteritems(new_id_mapping):
        #     print('{}\t{}'.format(k,v))
        # print(len(cpair_dict))
        # for (c1, c2) in cpair_dict:
        #     print(cpair_dict[(c1, c2)]['forward'])

        # for (c1,c2), v in iteritems(cpair_dict):
        #     print(type(v.back))
        # remove_cpairs = []
        # for (c1, c2), rxns in iteritems(cpair_dict):
        #     if (c2, c1) in cpair_dict:
        #         # A = []
        #         # for r in cpair_dict[(c2, c1)].back:
        #         #     A.append(r)
        #         # print(A)
        #
        #         rxns.forward.append(r for r in cpair_dict[(c2, c1)].back)
        #         rxns.back.append(r for r in cpair_dict[(c2, c1)].forward)
        #         rxns.bidirection.append(r for r in cpair_dict[(c2, c1)].bidirection)
        #         remove_cpairs.append((c2, c1))
        # for cpair in remove_cpairs:
        #     del cpair_dict[(cpair)]

        # print(len(filter_dict))
        # for r in filter_dict:
        #     print(r)

        # print(cpair_dict)
        # for (c1,c2), v in iteritems(cpair_dict):
        #     print('{}\t{}\t{}'.format(c1, c2, v))
        # for r, v in iteritems(filter_dict):
        #     print('{}\t{}'.format(r, v))

        # print(len(cpair_dict))
        # for rxn_id, cpairs in sorted(iteritems(raw_dict)):
        #     for (c1, c2) in cpairs:

        edge_values = make_edge_values(reaction_flux, self._mm, compound_formula, self._args.element,
                                       self._args.split_map, cpair_dict, new_id_mapping, self._args.method)
        # for (a, b), v in iteritems(edge_values):
        #     print('{}\t{}\t{}'.format(a,b,v))

        g = create_bipartite_graph(self._mm, self._model, cpair_dict,self._args.element, self._args.split_map,
                                   subset_reactions, edge_values,compound_formula, reaction_flux, self._args.method,
                                   new_id_mapping, self._args.color, self._args.detail)

        # for node in g.nodes:
        #     if node.props['shape'] == 'box':
        #         print(node)
        if self._args.method == 'no-fpp' and self._args.split_map is True:
            logger.warning('--split-map option makes no difference when method is no-fpp')

        with open('reactions.dot', 'w') as f:
            g.write_graphviz(f)
        with open('reactions.nodes.tsv', 'w') as f:
            g.write_cytoscape_nodes(f)
        with open('reactions.edges.tsv', 'w') as f:
            g.write_cytoscape_edges(f)

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


def create_bipartite_graph(mm, model, cpair_dict, element, split_map, subset, edge_values, cpd_formula,
                           reaction_flux, method, new_id_mapping, args_color, args_detail):
    """create bipartite graph of given metabolic network

    Start from a dictionary comprises compound pairs and related reaction ids, Returns a Graph object
    that contains a set of nodes and a dictionary of edges, node and edge properties(such as node color, shape and
    edge direction) are included.

    Args:
    mm: class 'psamm.metabolicmodel.MetabolicModel'.
    model: class 'psamm.datasource.native.NativeModel'.
    predict_results: Dictionary mapping reaction IDs to compound pairs(reactant/product pair that transfers
        specific element,by default the element is carbon(C).
    cpair_dict: Dictionary mapping compound pair to a list of reaction IDs.
    element: a string that represent a specific chemical element, such as C(carbon), S(sulfur), N(nitrogen).
    subset: Set of reactions for visualizing.
    edge_values: Dictionary mapping (reaction ID, compound ID) to values of edge between them.
    cpd_formula: Dictionary mapping compound IDs to
        :class:`psamm.formula.Formula`. Formulas must be flattened.
    reaction_flux: Dictionary mapping reaction ID to reaction flux. Flux is a float number.
    split_graph: An argument used to decide if split node for one reaction. Default is 'True'"""

    g = graph.Graph()

    # Mapping from compound id to DictCompoundEntry object
    cpd_entry = {}
    for cpd in model.compounds:
        cpd_entry[cpd.id] = cpd

    # Mapping from reaction id to DictReactionEntry object
    rxn_entry = {}
    for rxn in model.reactions:
        rxn_entry[rxn.id] = rxn

    # setting node color
    recolor_dict = {}
    if args_color is not None:
        for f in args_color:
            for row in csv.reader(f, delimiter=str(u'\t')):
                recolor_dict[row[0]] = row[1]  # row[0] =reaction id or str(cpd object), row[1] = hex color code
    color = {}
    for c in mm.compounds:
        if str(c) in recolor_dict:
            color[c] = recolor_dict[str(c)]
        else:
            color[c] = COMPOUND_COLOR
    for r in mm.reactions:
        if r in recolor_dict:
            color[r] = recolor_dict[r]
        else:
            color[r] = REACTION_COLOR

    # define reaction node color for rxns-combined graph
    def final_rxn_color(color_args, rlist):
        if color_args is not None:
            if len(rlist) == 1:
                return color[new_id_mapping[rlist[0]]]
            else:
                if any(new_id_mapping[r] in recolor_dict for r in rlist):
                    return RXN_COMBINED_COLOR
                else:
                    return REACTION_COLOR
        else:
            return REACTION_COLOR

    # preparing for scaling width of edges
    if len(edge_values) > 0:
        value_list = sorted(edge_values.values())
        ninty_percentile = value_list[int(len(value_list)*0.9)+1]
        min_edge_value = min(itervalues(edge_values))
        max_edge_value = ninty_percentile
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

            return 10 * alpha

    def dir_value(direction):
        """assign value to different reaction directions"""
        if direction == Direction.Forward:
            return 'forward'
        elif direction == Direction.Reverse:
            return 'back'
        else:
            return 'both'

    # def final_props(reaction, edge1, edge2):
    #     """set final properties of edges"""
    #     if len(edge_values) > 0:
    #         p = {}
    #         if edge1 in edge_values:
    #             p['dir'] = 'forward'
    #             p['penwidth'] = pen_width(edge_values[edge1])
    #         elif edge2 in edge_values:
    #             p['dir'] = 'back'
    #             p['penwidth'] = pen_width(edge_values[edge2])
    #         else:
    #             p['style'] = 'dotted'
    #             p['dir'] = dir_value(reaction.direction)
    #         return p
    #     else:
    #         return {'dir': dir_value(reaction.direction)}

    def condensed_rxn_props(detail, r_list, reaction_flux):
        if len(r_list) == 1:
            r = new_id_mapping[r_list[0]]
            label_comb = rxns_properties(rxn_entry[r], detail, reaction_flux)
        else:
            if len(reaction_flux) > 0:
                sum_flux = 0
                for r in r_list:
                    if new_id_mapping[r] in reaction_flux:
                        sum_flux += abs(reaction_flux[new_id_mapping[r]])
                label_comb = '\n'.join(new_id_mapping[r] for r in r_list)
                if sum_flux != 0:
                    label_comb = '{}\n{}'.format(label_comb, sum_flux)
            else:
                label_comb = '\n'.join(new_id_mapping[r] for r in r_list)
        return label_comb

    # def dir_value_2(r_list):
    #     if r_list == rxns['forward']:
    #         return {'dir': 'forward'}
    #     elif r_list == rxns.back:
    #         return {'dir': 'back'}
    #     else:
    #         return {'dir': 'both'}

    def final_props_2(dir, edge1, edge2):
        if len(edge_values) > 0:
            p = {}
            if edge1 in edge_values:
                p['penwidth'] = pen_width(edge_values[edge1])
            elif edge2 in edge_values:
                p['penwidth'] = pen_width(edge_values[edge2])
            else:
                p['style'] = 'dotted'
            p['dir'] = dir
            return p
        else:
            return {'dir': dir}

    # create compound nodes and add nodes to Graph object
    compound_set = set()
    for r in subset:
        rx = mm.get_reaction(r)
        if r != model.biomass_reaction:
            for c, _ in rx.compounds:
                if element != 'none':
                    if Atom(primary_element(element)) in cpd_formula[str(c.name)]:
                        compound_set.add(c)
                else:
                    compound_set.add(c)
        # if model.biomass_reaction is not None:
        #     if model.biomass_reaction in subset:
        #         for c, _ in rx.compounds:
        #             if element != 'none':
        #                 if Atom(primary_element(element)) in cpd_formula[str(c.name)]:
        #                     compound_set.add(c)
        #             else:
        #                 compound_set.add(c)
        #             if c, _ in mm.get(model.biomass_reaction).right:
        #                 compound_set.remove(c)
        # else:
        #     logger.warning('No biomass reaction in this model')

    compound_nodes = {}
    for cpd in compound_set:    # cpd=compound object, cpd.name=compound id, no compartment
        node = graph.Node({
            'id': text_type(cpd),
            'label': cpds_properties(cpd, cpd_entry[cpd.name], args_detail),
            'shape': 'ellipse',
            'style': 'filled',
            'fillcolor': color[cpd]})
        compound_nodes[cpd] = node
        g.add_node(node)

    # create and add reaction nodes, add edges
    edge_list = []
    for (c1, c2), rxns in iteritems(cpair_dict):
        for dir, r_list in iteritems(rxns):
            if len(r_list) > 0:
                r_node = []
                if split_map is True or method == 'no-fpp':
                    for r in r_list:
                        r_node.append([r])
                else:
                    r_node.append(r_list)
                for i in r_node:
                    # print(i)
                    node = graph.Node({
                        'id': i,
                        'label': condensed_rxn_props(args_detail, i, reaction_flux),
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': final_rxn_color(args_color, i)})
                    g.add_node(node)

                    if len(i) == 1:
                        reac = new_id_mapping[i[0]]  # a single, real rxn id
                        test = i[0]
                    else:
                        reac = tuple(i)  # a list of rxn id
                        test = i

                    edge1 = c1, reac
                    edge2 = reac, c1
                    edge_test_1 = c1, test
                    edge_test_2 = test, c1
                    if edge_test_1 and edge_test_2 not in edge_list:
                        edge_list.append(edge_test_1)
                        edge_list.append(edge_test_2)   # avoid duplicated edges in no-fpp graph
                        g.add_edge(graph.Edge(
                            compound_nodes[c1], node, final_props_2(dir, edge1, edge2)))

                    edge1 = reac, c2
                    edge2 = c2, reac
                    edge_test_1 = test, c2
                    edge_test_2 = c2, test
                    if edge_test_1 and edge_test_2 not in edge_list:
                        edge_list.append(edge_test_1)
                        edge_list.append(edge_test_2)
                        g.add_edge(graph.Edge(
                            node, compound_nodes[c2], final_props_2(dir, edge1, edge2)))
    # for edge in g.edges:
    #     print('{}\t{}\t{}'.format(edge.source, edge.dest, edge.props))

    # add exchange reaction nodes
    rxn_set = set()
    for reaction in subset:
        if mm.is_exchange(reaction):
            raw_exchange_rxn = mm.get_reaction(reaction)
            if element != 'none':
                if any(Atom(primary_element(element)) in cpd_formula[str(c[0].name)]
                       for c in raw_exchange_rxn.compounds):
                        rxn_set.add(reaction)
            else:
                rxn_set.add(reaction)
    for r in rxn_set:
        exchange_rxn = mm.get_reaction(r)
        label = r
        if len(edge_values) > 0:
            if r in iterkeys(edge_values):
                label = '{}\n{}'.format(r, edge_values[r])
        node_ex = graph.Node({
            'id': r,
            # 'edge_id': r,
            'label': label,
            'shape': 'box',
            'style': 'filled',
            'fillcolor': ACTIVE_COLOR})
        g.add_node(node_ex)

        dir = dir_value(exchange_rxn.direction)
        for c1, _ in exchange_rxn.left:
            edge1 = c1, r
            edge2 = r, c1
            g.add_edge(graph.Edge(
                compound_nodes[c1], node_ex, final_props_2(dir, edge1, edge2)))

        for c2, _ in exchange_rxn.right:
            edge1 = r, c2
            edge2 = c2, r
            g.add_edge(graph.Edge(
                node_ex, compound_nodes[c2], final_props_2(dir, edge1, edge2)))

    # add biomass reaction nodes
    bio_pair = Counter()
    if model.biomass_reaction is not None:
        if model.biomass_reaction in subset:
            biomass_reaction = mm.get_reaction(model.biomass_reaction)
            for c, _ in biomass_reaction.left:
                if c in compound_nodes:
                    bio_pair[model.biomass_reaction] += 1
                    # bio_pair = Counter({'biomass_rxn_id': 1}), Counter({'biomass_rxn_id': 2})...
                    node_bio = graph.Node({
                        'id': '{}_{}'.format(model.biomass_reaction, bio_pair[model.biomass_reaction]),
                        # 'edge_id': model.biomass_reaction,
                        'label': model.biomass_reaction,
                        'shape': 'box',
                        'style': 'filled',
                        'fillcolor': ALT_COLOR})
                    g.add_node(node_bio)

                    edge1 = c, model.biomass_reaction
                    edge2 = model.biomass_reaction, c
                    g.add_edge(graph.Edge(
                        compound_nodes[c], node_bio, final_props_2('forward', edge1, edge2)))
    else:
        logger.warning(
            'No biomass reaction in this model')

    return g