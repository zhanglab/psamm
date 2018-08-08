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
import csv
import argparse
from collections import defaultdict, Counter
from six import text_type, iteritems, itervalues, iterkeys

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

import subprocess

logger = logging.getLogger(__name__)

REACTION_COLOR = '#ccebc5'
COMPOUND_COLOR = '#fbb4ae'
ACTIVE_COLOR = '#b3cde3'    # exchange reaction  color
ALT_COLOR = '#f4fc55'       # biomass reaction color
RXN_COMBINED_COLOR = '#fc9a44'


class VisualizationCommand(MetabolicMixin,ObjectiveMixin,SolverCommandMixin,
                           Command, LoopRemovalMixin, FilePrefixAppendAction):
    """Run visualization command on the model."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--method',type=text_type, default='fpp',
            help='Compound pair prediction method,'
                 'choices=[fpp, no-fpp, or file path]')
        parser.add_argument(
            '--exclude', metavar='reaction', type=text_type, default=[],
            action=FilePrefixAppendAction,
            help='Reaction to exclude (e.g. biomass reactions or '
                 'macromolecule synthesis)')
        parser.add_argument(
            '--fba', action='store_true',
            help='visualize reaction flux')
        parser.add_argument(
            '--element', type=text_type, default='C',
            help='Primary element flow, followed by chemical element, '
                 'such as C, N, O,S.')
        parser.add_argument(
            '--detail', type=text_type, default=None, action='append',
            nargs='+', help='Reaction and compound properties showed '
                            'on nodes label(e.g. formula, equation, genes')
        parser.add_argument(
            '--subset', type=argparse.FileType('rU'), default=None,
            help='Reactions designated to visualize')
        parser.add_argument(
            '--color', type=argparse.FileType('rU'), default=None, nargs='+',
            help='Customize color of reaction and compound nodes ')
        parser.add_argument(
            '--Image', type=text_type, default=None,
            help='generate image file directly')
        parser.add_argument(
            '--exclude-cpairs', type=argparse.FileType('rU'), default=None,
            help='Remove edge of given compound pairs from network ')
        parser.add_argument(
            '--split-map', action='store_true',
            help='Create reactions-splitted metabolic network, one node '
                 'only represent one reaction')
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
                                compound.id, compound.formula))
                except ParseError as e:
                    msg = (
                        'Error parsing formula'
                        ' for compound {}:\n{}\n{}'.format(
                            compound.id, e, compound.formula))
                    if e.indicator is not None:
                        msg += '\n{}'.format(e.indicator)
                    logger.warning(msg)

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
                logger.error(
                    'Invalid subset file. The file should contains a column '
                    'of reaction id or a column of compound id with '
                    'compartment, mix of reactions, compounds and other '
                    'infomation in one subset file is not allowed.')
                quit()
        else:
            subset_reactions = set(self._mm.reactions)

        # create {rxn_id:[(c1, c2),(c3,c4),...], ...} dictionary,
        # key = rxn id, value = list of compound pairs
        filter_dict, fpp_rxns = make_filter_dict(
            self._model, self._mm, self._args.method, self._args.element,
            compound_formula, self._args.exclude_cpairs, self._args.exclude)

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
            logger.info('Minimized reactions: {}'.format(count))

        cpair_dict, new_id_mapping = make_cpair_dict(
            self._mm, filter_dict, subset_reactions, reaction_flux,
            self._args.method)

        edge_values = make_edge_values(
            reaction_flux, self._mm, compound_formula, self._args.element,
            self._args.split_map, cpair_dict, new_id_mapping,
            self._args.method)

        # create rxns/cpds recolor dict
        recolor_dict = {}
        if self._args.color is not None:
            for f in self._args.color:
                for row in csv.reader(f, delimiter=str(u'\t')):
                    recolor_dict[row[0]] = row[1]
        print(recolor_dict)

        model_compoundEntries, model_reactionEntries = {}, {}
        for c in self._model.compounds:
            model_compoundEntries[c.id] = c
        for r in self._model.reactions:
            model_reactionEntries[r.id] = r

        # g = create_bipartite_graph(self._mm, self._model, cpair_dict,
        #                            self._args.split_map, edge_values,
        #                            subset_reactions,reaction_flux,
        #                            self._args.method, new_id_mapping,
        #                            self._args.color,self._args.detail)
        g = graph.Graph()
        g = add_graph_nodes(g, cpair_dict, self._args.method,
                            new_id_mapping, split=self._args.split_map)
        g = add_node_color(g, recolor_dict)
        g = add_edges(g, cpair_dict, self._args.method, split = self._args.split_map)

        if self._model.biomass_reaction in subset_reactions:
            biomass_rxn = self._mm.get_reaction(self._model.biomass_reaction)
            g = add_biomass_rxns(g, biomass_rxn, self._model.biomass_reaction)

        for reaction in self._mm.reactions:
            if self._mm.is_exchange(reaction):
                exchange_rxn = self._mm.get_reaction(reaction)
                g = add_Exchange_rxns(g, reaction, exchange_rxn)
        g = add_node_label(g, self._args.detail, model_compoundEntries,
                           model_reactionEntries, reaction_flux)
        g = set_edge_props_withfba(g, edge_values)

        if self._args.method == 'no-fpp' and self._args.split_map is True:
            logger.warning('--split-map is not compatible with no-fpp option')

        with open('reactions.dot', 'w') as f:
            g.write_graphviz(f)
        with open('reactions.nodes.tsv', 'w') as f:
            g.write_cytoscape_nodes(f)
        with open('reactions.edges.tsv', 'w') as f:
            g.write_cytoscape_edges(f)

        if self._args.Image is not None:
            if render is None:
                self.fail(
                    'create image file requires python binding graphviz '
                    'module ("pip install graphviz")')
            else:
                if len(subset_reactions) > 500:
                    logger.info(
                        'This graph contains {} reactions, graphs of this '
                        'size may take a long time to create'.format
                        (len(subset_reactions)))
                render('dot', self._args.Image, 'reactions.dot')
                # except subprocess.CalledProcessError:
                #     logger.warning('the graph is too large to create')


# def cpds_properties(cpd, compound, detail):
#     """define compound nodes label.
#
#     Args:
#     cpd: class 'psamm.reaction.Compound'.
#     compound: class 'psamm.datasource.entry.DictCompoundEntry'.
#     detail: A list that contains only one element, this element is a list of
#         reaction or compound properties name included in the model.
#         e.g. detrail = [['id', 'name', 'formula']].
#     """
#     compound_set = set()
#     compound_set.update(compound.properties)
#     if detail is not None:
#         cpd_detail = []
#         for prop in detail[0]:
#             if prop in compound_set:
#                 cpd_detail.append(prop)
#         pre_label = '\n'.join(_encode_value(compound.properties[value])
#                               for value in cpd_detail if value != 'id')
#         if 'id' in detail[0]:
#             label = '{}\n{}'.format(str(cpd), pre_label)
#         else:
#             if all(prop not in compound_set for prop in detail[0]):
#                 label = str(cpd)
#             else:
#                 label = pre_label
#     else:
#         label = str(cpd)
#     return label
#
#
# def rxns_properties(reaction, detail, reaction_flux):
#     """define reaction nodes label.
#
#     Args:
#         reaction: class 'psamm.datasource.entry.DictReactionEntry'.
#         detail: A list of reaction or compound properties name included
#             in the model. e.g. [id, genes, equation].
#         reaction_flux: Dictionary of reaction ID and flux value.
#     """
#     reaction_set = set()
#     reaction_set.update(reaction.properties)
#     if detail is not None:
#         rxn_detail = []
#         for prop in detail[0]:
#             if prop in reaction_set:
#                 rxn_detail.append(str(prop))
#         if len(rxn_detail) == 0:
#             label = reaction.id
#         else:
#             label = '\n'.join(_encode_value(reaction.properties[value])
#                               for value in rxn_detail)
#         if len(reaction_flux) > 0:
#             if reaction.id in iterkeys(reaction_flux):
#                 label = '{}\n{}'.format(label, reaction_flux[reaction.id])
#     else:
#         if len(reaction_flux) > 0:
#             if reaction.id in reaction_flux:
#                 label ='{}\n{}'.format(reaction.id,reaction_flux[reaction.id])
#             else:
#                 label = reaction.id
#         else:
#             label = reaction.id
#
#     return label


def make_edge_values(reaction_flux, mm, compound_formula, element, split_map,
                     cpair_dict, new_id_mapping, method):
    """set edge_values according to reaction fluxes

    Start from reaction equation and reaction flux to create a dictionary
    that key is an edge(a tuple of (c, r) or (r, c)) and vlaue is flux
    between the compound and reaction

    Args:
        reactions_flux: Dictionary of reaction id and flux value.
        mm: class 'psamm.metabolicmodel.MetabolicModel'.
        compound_formula: Dictionary of compound id and
            class 'psamm.formula.Formula'.
        element: Chemical symbol of an element, such as C(carbon),
            N(nitrogen), S(sulfer).
        slit_map: True or False, determine which kind of gragh to create(
            combine multiple reactions in one node or not).
        cpair_dict: Defauldict of compound pair and another defaultdict of
            list(in the latter defaultdict, key is direction(forward, back or
            both), value is reaction list),{(c1, c2): {'forward':[rxns],...}}.
        new_id_mapping: Dictionary of reaction is with suffix(number) and
            real reaction id, e.g.{r_1:r, r_2: r,...}.
        method: Command line argument, including 3 options:'fpp', 'no-fpp'
            and file path.
        """
    edge_values = {}
    edge_values_combined = {}
    if len(reaction_flux) > 0:
        for reaction in reaction_flux:
            rx = mm.get_reaction(reaction)
            flux = reaction_flux[reaction]
            if abs(flux) < 1e-9:
                continue

            for compound, value in rx.compounds:
                if element == 'none':
                    edge_values[compound, reaction] = abs(flux * float(value))
                else:
                    if compound.name in compound_formula:
                        if Atom(element) in compound_formula[compound.name]:
                            edge_values[compound, reaction] = \
                                abs(flux * float(value))
                    else:
                        edge_values[compound, reaction] = abs(flux * float(value))
        rxn_set = set()
        for (c1, c2), rxns in iteritems(cpair_dict):
            for dir, rlist in iteritems(rxns):
                rlist_string = str(','.join(new_id_mapping[r] for r in rlist))
                if any(new_id_mapping[r] in reaction_flux for r
                       in rlist):
                    x_comb_c1, x_comb_c2 = 0, 0
                    for r in rlist:
                        real_r = new_id_mapping[r]
                        rxn_set.add(real_r)
                        if real_r in reaction_flux:
                            x_comb_c1 += edge_values[(c1, real_r)]
                            x_comb_c2 += edge_values[(c2, real_r)]
                    edge_values_combined[(c1, rlist_string)] = x_comb_c1
                    edge_values_combined[(c2, rlist_string)] = x_comb_c2
        for r in reaction_flux:
            if r not in rxn_set:
                for (_, r) in edge_values:
                    edge_values_combined[(_,r)] = edge_values[(_,r)]

    if split_map or method == 'no-fpp':
        return edge_values
    else:
        return edge_values_combined


def make_filter_dict(model, mm, method, element, cpd_formula,
                     args_exclude_cpairs, exclude_rxns):
    """create a dictionary of reaction id(key) and a list of related
    compound pairs(value) by 3 different methods.

    Args:
        model: class 'psamm.datasource.native.NativeModel'.
        mm: class 'psamm.metabolicmodel.MetabolicModel'.
        method: Command line argument, including 'fpp', 'no-fpp' and
            file path.
        element: Chemical symbol of an element, such as C(carbon),
            N(nitrogen), S(sulfer).
        compound_formula: A dictionary that key is compound id (a string) and
            value is class 'psamm.formula.Formula'.
        args_exclude_cpairs: Command line argument, a file that contains two
            columns(tab-separated) of compounds(in the format of
            compound_id[compartment], such as atp[c], akg[c],each row
            represent two reactant/product pairs((c1, c2)and (c2,c1)).
        exclude_rxns: Command line argument, a file that contains one column
            of reaction IDs, these reaction will be removed for visualization
            when method is fpp."""

    # Mapping from cpd_id+compartment(eg: pyr_c[c]) to Compound object
    cpd_object = {}
    for cpd in mm.compounds:
        cpd_object[str(cpd)] = cpd

    # read exclude_compound_pairs from command-line argument
    exclude_cpairs = []
    if args_exclude_cpairs is not None:
        for row in csv.reader(args_exclude_cpairs, delimiter=str('\t')):
            exclude_cpairs.append((cpd_object[row[0]], cpd_object[row[1]]))
            exclude_cpairs.append((cpd_object[row[1]], cpd_object[row[0]]))

    fpp_rxns, rxns_no_equation, rxns_no_formula = set(), set(), []
    for reaction in model.reactions:
        if (reaction.id in model.model and
                reaction.id not in exclude_rxns):
            if reaction.equation is None:
                rxns_no_equation.add(reaction.id)
                continue

            if any(c.name not in cpd_formula for c, _ in
                   reaction.equation.compounds):
                rxns_no_formula.append(reaction.id)
                continue

            fpp_rxns.add(reaction)

    if len(rxns_no_equation) > 0:
        logger.warning(
            '{} reactions have no reaction equation, fix or exclude them.'
            'These reactions contain {}'.format(len(rxns_no_equation),
                                                rxns_no_equation))

    if len(rxns_no_formula) > 0:
        logger.warning(
            '{} reactions have compounds with undefined formula. '
            'These reactions contain {}'.format(len(rxns_no_formula),
                                                rxns_no_formula))

    filter_dict = {}
    if method == 'fpp':
        if len(fpp_rxns) == 0:
            logger.error(
                'All the reactions have compounds with undefined formula or '
                'have no reaction equation, fix them or '
                'add "--method no-fpp --element none" to the command')
            quit()

        reaction_pairs = [(r.id, r.equation) for r in fpp_rxns]
        element_weight = findprimarypairs.element_weight
        fpp_dict, _ = findprimarypairs.predict_compound_pairs_iterated(
            reaction_pairs, cpd_formula, element_weight=element_weight)

        for rxn_id, fpp_pairs in iteritems(fpp_dict):
            compound_pairs = []
            for cpd_pair, transfer in iteritems(fpp_pairs[0]):
                if cpd_pair not in exclude_cpairs:
                    if element == 'none':
                        compound_pairs.append(cpd_pair)
                    else:
                        if any(Atom(element) in k for k in transfer):
                            compound_pairs.append(cpd_pair)
            if len(compound_pairs) != 0:
                filter_dict[rxn_id] = compound_pairs

    elif method == 'no-fpp':
        for rxn_id in mm.reactions:
            if rxn_id != model.biomass_reaction:
                rx = mm.get_reaction(rxn_id)
                cpairs = []
                for c1, _ in rx.left:
                    for c2, _ in rx.right:
                        if (c1, c2) not in exclude_cpairs:
                            if element != 'none':
                                if c1.name in cpd_formula and c2.name in \
                                        cpd_formula:
                                    if Atom(element) in cpd_formula[c1.name] \
                                            and Atom(element) in \
                                            cpd_formula[c2.name]:
                                        cpairs.append((c1, c2))
                                else:
                                    cpairs.append((c1,c2))
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
                            exclude_cpairs:
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
        except:
            if IOError:
                logger.error('Invalid file path, no such file or directory '
                             ': {}' .format(method))
            quit()

    return filter_dict, fpp_rxns


def make_cpair_dict(mm, filter_dict, subset, reaction_flux, args_method):
    """create a new cpair_dict: mapping from cpair to a
    defaultdict(key is direction, value is rxns list).
    """
    new_id_mapping = {}
    rxn_count = Counter()
    cpair_dict = defaultdict(lambda: defaultdict(list))
    for r, cpairs in iteritems(filter_dict):
        if r in subset:
            rx = mm.get_reaction(r)
            for (c1, c2) in cpairs:
                if args_method == 'no-fpp':
                    r_id = r
                    new_id_mapping[r_id] = r
                else:
                    rxn_count[r] += 1
                    r_id = str('{}_{}'.format(r, rxn_count[r]))
                    new_id_mapping[r_id] = r
                a = tuple(sorted((c1, c2))) == (c1, c2)
                if rx.direction == Direction.Forward:
                    if a:
                        cpair_dict[(c1, c2)]['forward'].append(r_id)
                    else:
                        cpair_dict[(c2, c1)]['back'].append(r_id)
                elif rx.direction == Direction.Reverse:
                    if a:
                        cpair_dict[(c1, c2)]['back'].append(r_id)
                    else:
                        cpair_dict[(c2, c1)]['forward'].append(r_id)
                else:
                    if len(reaction_flux) > 0:
                        if r in reaction_flux:
                            if reaction_flux[r] > 0:
                                if a:
                                    cpair_dict[(c1, c2)]['forward']. \
                                        append(r_id)
                                else:
                                    cpair_dict[(c2, c1)]['back'].append(
                                        r_id)
                            else:
                                if a:
                                    cpair_dict[(c1, c2)]['back'].append(
                                        r_id)
                                else:
                                    cpair_dict[(c2, c1)]['forward']. \
                                        append(r_id)
                        else:
                            cpair_dict[tuple(sorted((c1, c2)))]['both']. \
                                append(r_id)
                    else:
                        cpair_dict[tuple(sorted((c1, c2)))]['both']. \
                            append(r_id)

    # new_cpair_dict = {}
    # for (c1, c2), rxns in iteritems(cpair_dict):
    #     new_rxns = {}
    #     if len(rxns['forward']) == 0:
    #         # new_rxns = defaultdict(list)
    #         if len(rxns['back']) != 0:
    #             new_rxns['forward'] = rxns['back']
    #             new_rxns['back'] = []
    #             new_rxns['both'] = rxns['both']
    #         else:
    #             new_rxns['forward'] = []
    #             new_rxns['back'] = []
    #             new_rxns['both'] = rxns['both']
    #         new_cpair_dict[(c2, c1)]=new_rxns
    #     else:
    #         new_rxns['forward'] = rxns['forward']
    #         new_rxns['back'] = rxns['back']
    #         new_rxns['both'] = rxns['both']
    #         new_cpair_dict[(c1, c2)] = new_rxns   # Note: all use normal dict didn't works

    # new_cpair_dict = defaultdict(lambda: defaultdict(list))
    # for (c1, c2), rxns in iteritems(cpair_dict):
    #     if len(rxns['forward']) == 0:
    #         if len(rxns['back']) != 0:
    #             new_cpair_dict[(c2, c1)]['forward'] = rxns['back']
    #             if len(rxns['both']) != 0:
    #                 new_cpair_dict[(c2, c1)]['both'] = rxns['both']
    #         else:
    #             new_cpair_dict[(c2, c1)]['both'] = rxns['both']
    #     else:
    #         new_cpair_dict[(c1, c2)] = rxns  # Note: all defaultdict, sometimes works, sometimes seems not
    # print(new_cpair_dict)

    new_cpair_dict = {}
    for (c1, c2), rxns in iteritems(cpair_dict):
        if len(rxns['forward']) == 0:
            new_rxns = defaultdict(list)
            if len(rxns['back']) != 0:
                new_rxns['forward'] = rxns['back']
                if len(rxns['both']) != 0:
                    new_rxns['both'] = rxns['both']
            else:
                new_rxns['both'] = rxns['both']
            new_cpair_dict[(c2, c1)]=new_rxns
        else:
            new_cpair_dict[(c1, c2)] = rxns

    # new_cpair_dict = {}
    # for (c1, c2), rxns in iteritems(cpair_dict):
    #     if len(rxns['forward']) == 0:
    #         new_rxns = {}
    #         if len(rxns['back']) != 0:
    #             new_rxns['forward'] = rxns['back']
    #             new_rxns['back'] = []
    #             new_rxns['both'] = rxns['both']
    #         else:
    #             new_rxns['forward'] = []
    #             new_rxns['back'] = []
    #             new_rxns['both'] = rxns['both']
    #         new_cpair_dict[(c2, c1)]=new_rxns
    #     else:
    #         new_cpair_dict[(c1, c2)] = rxns     # Note: part normal dict, part defaultdict, seems works better

    # print(new_cpair_dict)
    # for (c1, c2), rxns in iteritems(new_cpair_dict):
    #     print(str(c1), str(c2), rxns['forward'], rxns['back'], rxns['both'])

    return new_cpair_dict, new_id_mapping

# def create_bipartite_graph(mm, model, cpair_dict, split_map, edge_values,
#                            subset, reaction_flux, method, new_id_mapping,
#                            args_color, args_detail):
#     """create bipartite graph of given metabolic network
#
#     Start from a dictionary comprises compound pairs and related reaction ids,
#     Returns a Graph object that contains a set of nodes and a dictionary of
#     edges, node and edge properties(such as node color, shape and edge
#     direction) are included.
#
#     Args:
#         mm: class 'psamm.metabolicmodel.MetabolicModel'.
#         model: class 'psamm.datasource.native.NativeModel'.
#         cpair_dict: Dictionary mapping compound pair(c1,c2) to another
#             dictionary that mapping from direction to a list of reaction IDs.
#         split_map: Command line argument. True or False.
#         element: Chemical symbol of an element, such as C(carbon),
#             S(sulfur), N(nitrogen).
#         subset: A set of reaction IDs.
#         edge_values: Dictionary mapping an edge (a tuple, (r, c) or (c, r)) to
#             a float number that represent width of edge. r is reaction id, c is
#             compound id, sometimes r is a tuple of multiple reaction ids.
#         reaction_flux: Dictionary mapping from reaction ID to reaction flux.
#             Flux is a float number.
#         method: Command line argument, including 3 choices: 'fpp', 'no-fpp'
#             and a file path.
#         new_id_mapping: A dictionary mapping from reaction id with suffix to
#             real reaction id.
#         args_color: Command line argument, a file comprised by two columns
#             (tab-separated), the first column is reaction id or
#             compound_id[compartment], the second columns is hex color code.
#         args_detail: Command line argument, including text_type of reaction
#             or compound properties which are defined in model, e.g. formula,
#             equation, genes."""
#
#     g = graph.Graph()
#
#     # Mapping from compound id to DictCompoundEntry object
#     cpd_entry = {}
#     for cpd in model.compounds:
#         cpd_entry[cpd.id] = cpd
#
#     # Mapping from reaction id to DictReactionEntry object
#     rxn_entry = {}
#     for rxn in model.reactions:
#         rxn_entry[rxn.id] = rxn
#
#     # setting node color
#     recolor_dict = {}
#     if args_color is not None:
#         for f in args_color:
#             for row in csv.reader(f, delimiter=str(u'\t')):
#                 recolor_dict[row[0]] = row[1]
#     color = {}
#     for c in mm.compounds:
#         if str(c) in recolor_dict:
#             color[c] = recolor_dict[str(c)]
#         else:
#             color[c] = COMPOUND_COLOR
#     for r in mm.reactions:
#         if r in recolor_dict:
#             color[r] = recolor_dict[r]
#         else:
#             color[r] = REACTION_COLOR
#
#     # define reaction node color for rxns-combined graph
#     def final_rxn_color(color_args, rlist):
#         if color_args is not None:
#             if len(rlist) == 1:
#                 return color[new_id_mapping[rlist[0]]]
#             else:
#                 if any(new_id_mapping[r] in recolor_dict for r in rlist):
#                     return RXN_COMBINED_COLOR
#                 else:
#                     return REACTION_COLOR
#         else:
#             return REACTION_COLOR
#
#     # preparing for scaling width of edges
#     if len(edge_values) > 0:
#         value_list = sorted(edge_values.values())
#         ninty_percentile = value_list[int(len(value_list)*0.9)+1]
#         min_edge_value = min(itervalues(edge_values))
#         max_edge_value = ninty_percentile
#     else:
#         min_edge_value = 1
#         max_edge_value = 1
#
#     def pen_width(value):
#         """calculate final edges width"""
#         if max_edge_value - min_edge_value == 0:
#             return 1
#         else:
#             if value > max_edge_value:
#                 value = max_edge_value
#             alpha = value / max_edge_value
#
#             return 10 * alpha
#
#     def dir_value(direction):
#         """assign value to different reaction directions"""
#         if direction == Direction.Forward:
#             return 'forward'
#         elif direction == Direction.Reverse:
#             return 'back'
#         else:
#             return 'both'
#
#     def condensed_rxn_props(detail, r_list, reaction_flux):
#         if len(r_list) == 1:
#             r = new_id_mapping[r_list[0]]
#             label_comb = rxns_properties(rxn_entry[r], detail, reaction_flux)
#         else:
#             if len(reaction_flux) > 0:
#                 sum_flux = 0
#                 for r in r_list:
#                     if new_id_mapping[r] in reaction_flux:
#                         sum_flux += abs(reaction_flux[new_id_mapping[r]])
#                 label_comb = '\n'.join(new_id_mapping[r] for r in r_list)
#                 if sum_flux != 0:
#                     label_comb = '{}\n{}'.format(label_comb, sum_flux)
#             else:
#                 label_comb = '\n'.join(new_id_mapping[r] for r in r_list)
#         return label_comb
#
#     def final_props(dir, edge):
#         if len(edge_values) > 0:
#             p = {}
#             if edge in edge_values:
#                 p['penwidth'] = pen_width(abs(edge_values[edge]))
#             else:
#                 p['style'] = 'dotted'
#             p['dir'] = dir
#             return p
#         else:
#             return {'dir': dir}
#
#     # create compound nodes and add nodes to Graph object
#     compound_set = set()
#     for (c1, c2), rxns in iteritems(cpair_dict):
#         compound_set.add(c1)
#         compound_set.add(c2)
#
#     compound_nodes = {}
#     for cpd in compound_set:  # cpd=cpd object,cpd.name=cpd id,no compartment
#         node = graph.Node({
#             'id': text_type(cpd),
#             'label': cpds_properties(cpd, cpd_entry[cpd.name], args_detail),
#             'shape': 'ellipse',
#             'style': 'filled',
#             'fillcolor': color[cpd]})
#         compound_nodes[cpd] = node
#         g.add_node(node)
#
#     # create and add reaction nodes, add edges
#     edge_list = []
#     for (c1, c2), rxns in iteritems(cpair_dict):
#         for dir, r_list in iteritems(rxns):
#             if len(r_list) > 0:
#                 r_node = []
#                 if split_map is True or method == 'no-fpp':
#                     for r in r_list:
#                         r_node.append([r])
#                 else:
#                     r_node.append(r_list)
#                 for i in r_node:
#                     node = graph.Node({
#                         'id': i,
#                         'label': condensed_rxn_props(args_detail, i,
#                                                      reaction_flux),
#                         'shape': 'box',
#                         'style': 'filled',
#                         'fillcolor': final_rxn_color(args_color, i)})
#                     g.add_node(node)
#
#                     reac = str(','.join(new_id_mapping[r] for r in i))
#                     test = i
#
#                     edge1 = c1, reac
#                     edge_test_1 = c1, test
#                     edge_test_2 = test, c1
#                     if edge_test_1 and edge_test_2 not in edge_list:
#                         edge_list.append(edge_test_1)
#                         edge_list.append(edge_test_2)
#                         g.add_edge(graph.Edge(compound_nodes[c1], node,
#                                               final_props(dir, edge1)))
#                     edge2 = c2, reac
#                     edge_test_1 = test, c2
#                     edge_test_2 = c2, test
#                     if edge_test_1 and edge_test_2 not in edge_list:
#                         edge_list.append(edge_test_1)
#                         edge_list.append(edge_test_2)
#                         g.add_edge(graph.Edge(node, compound_nodes[c2],
#                                           final_props(dir, edge2)))
#     # add exchange reaction nodes
#     rxn_set = set()
#     for reaction in mm.reactions:
#         if mm.is_exchange(reaction):
#             raw_exchange_rxn = mm.get_reaction(reaction)
#             for c, _ in raw_exchange_rxn.compounds:
#                 if c in compound_nodes:
#                     rxn_set.add(reaction)
#     for r in rxn_set:
#         exchange_rxn = mm.get_reaction(r)
#         label = r
#         if len(edge_values) > 0:
#             if r in iterkeys(edge_values):
#                 label = '{}\n{}'.format(r, edge_values[r])
#         node_ex = graph.Node({
#             'id': r,
#             'label': label,
#             'shape': 'box',
#             'style': 'filled',
#             'fillcolor': ACTIVE_COLOR})
#         g.add_node(node_ex)
#
#         dir = dir_value(exchange_rxn.direction)
#         for c1, _ in exchange_rxn.left:
#             edge1 = c1, r
#             g.add_edge(graph.Edge(
#                 compound_nodes[c1], node_ex, final_props(dir, edge1)))
#         for c2, _ in exchange_rxn.right:
#             edge2 = c2, r
#             g.add_edge(graph.Edge(
#                     node_ex, compound_nodes[c2], final_props(dir, edge2)))
#
#     # add biomass reaction nodes
#     bio_pair = Counter()
#     if model.biomass_reaction in subset:
#         biomass_rxn = mm.get_reaction(model.biomass_reaction)
#         dir = dir_value(biomass_rxn.direction)
#         A, B = [], []
#         for c, _ in biomass_rxn.left:
#             A.append(c)
#         for c, _ in biomass_rxn.right:
#             B.append(c)
#         for c, _ in biomass_rxn.compounds:
#             if c in compound_nodes:
#                 bio_pair[model.biomass_reaction] += 1
#                 node_bio = graph.Node({
#                     'id': '{}_{}'.format(model.biomass_reaction,
#                                          bio_pair[model.biomass_reaction]),
#                     'label': model.biomass_reaction,
#                     'shape': 'box',
#                     'style': 'filled',
#                     'fillcolor': ALT_COLOR})
#                 g.add_node(node_bio)
#
#                 edge = c, model.biomass_reaction
#                 if c in A:
#                     g.add_edge(graph.Edge(compound_nodes[c], node_bio,
#                                           final_props(dir, edge)))
#                 if c in B:
#                     g.add_edge(graph.Edge(node_bio, compound_nodes[c],
#                                           final_props(dir, edge)))
#     return g

### dveide create_bipartite_graph() function into several small functions


def add_graph_nodes(g, cpairs_dict, method, new_id_mapping, split): # by default, split = False
    compound_nodes = {}
    reaction_nodes = {}
    for cpair, reactions in iteritems(cpairs_dict):
        for c in cpair:
            node = graph.Node({
                'id': text_type(c),
                'shape': 'ellipse',
                'style': 'filled',
                'type': 'cpd',
                'original_id': c})
            g.add_node(node)
            compound_nodes[c] = node
        for dir, rlist in iteritems(reactions):
            rlist_str = str(','.join(rlist))
            if split or method == 'no-fpp':
                for sub_rxn in rlist:
                    rnode = graph.Node({
                        'id': text_type(sub_rxn),
                        'shape': 'box',
                        'style': 'filled',
                        'type': 'rxn',
                        'original_id': [new_id_mapping[sub_rxn]]})
                    g.add_node(rnode)
                    reaction_nodes[sub_rxn] = rnode
            else:
                real_rxns = [new_id_mapping[r] for r in rlist]
                rnode = graph.Node({
                    'id': text_type(rlist),
                    'shape': 'box',
                    'style': 'filled',
                    'type': 'rxn',
                    'original_id': real_rxns})
                g.add_node(rnode)
                reaction_nodes[rlist_str] = rnode
    return g


def add_node_color(g, recolor_dict):
    for node in g.nodes:
        if node.props['type'] == 'cpd':
            id = str(node.props['original_id'])
            color = recolor_dict.get(id)
        else:
            if len(node.props['original_id']) == 1:
                id = node.props['original_id'][0]
                color = recolor_dict.get(id)
            else:
                if any (i in recolor_dict for i in node.props['original_id']):
                    color = RXN_COMBINED_COLOR
                else:
                    color = REACTION_COLOR
        if color is not None:
            node.props['fillcolor'] = color
        else:
            if node.props['type'] == 'rxn':
                node.props['fillcolor'] = REACTION_COLOR
            elif node.props['type'] == 'cpd':
                node.props['fillcolor'] = COMPOUND_COLOR
    return g


def add_edges(g, cpairs_dict, method, split):
    node_dict = {}
    for node in g.nodes:
        node_dict[node.props['id']] = node
    edge_list = []
    for (c1, c2), value in iteritems(cpairs_dict):
        for dir, rlist in iteritems(value):
            if split or method == 'no-fpp':
                for sub_rxn in rlist:
                    test1 = c1, sub_rxn
                    # test2 = sub_rxn, c1
                    if test1 not in edge_list:
                        edge_list.append(test1)
                        # edge_list.append(test2)
                        g.add_edge(graph.Edge(
                            node_dict[text_type(c1)],
                            node_dict[text_type(sub_rxn)], {'dir': dir}))

                    # test1 = sub_rxn, c2
                    test2 = c2, sub_rxn
                    if test2 not in edge_list:
                        # edge_list.append(test1)
                        edge_list.append(test2)
                        g.add_edge(graph.Edge(
                            node_dict[text_type(sub_rxn)],
                            node_dict[text_type(c2)], {'dir': dir}))
            else:
                g.add_edge(graph.Edge(
                    node_dict[text_type(c1)], node_dict[text_type(rlist)],
                    {'dir': dir}))
                g.add_edge(graph.Edge(node_dict[text_type(rlist)],
                                      node_dict[text_type(c2)], {'dir': dir}))
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
    dir = dir_value(bio_reaction.direction)
    A, B = [], []
    for c, _ in bio_reaction.left:
        A.append(c)
    for c, _ in bio_reaction.right:
        B.append(c)
    cpd_nodes = {}
    for node in g.nodes:
        if node.props['type'] == 'cpd':
            cpd_nodes[node.props['original_id']] = node
    bio_pair = Counter()
    for c, _ in bio_reaction.compounds:
        if c in cpd_nodes:
            bio_pair[biomass_rxn_id] += 1
            node_bio = graph.Node({
                'id': '{}_{}'.format(biomass_rxn_id,
                                     bio_pair[biomass_rxn_id]),
                'label':text_type(biomass_rxn_id),
                'shape': 'box',
                'style': 'filled',
                'fillcolor': ALT_COLOR,
                'original_id': [biomass_rxn_id],
                'type': 'bio_rxn'})
            g.add_node(node_bio)

            if c in A:
                g.add_edge(graph.Edge(cpd_nodes[c], node_bio,
                                      {'dir': dir}))
            if c in B:
                g.add_edge(graph.Edge(node_bio, cpd_nodes[c],
                                      {'dir': dir}))
    return g


def add_Exchange_rxns(g, rxn_id, reaction):
    """add exchange reaction nodes and edges to graph object.

    Args:
        g: a graph object that contains a set of nodes and some edges.
        rxn_id: Exchange reaction id,
        reaction: Exchange reaction object(metabolic model reaction).
    """
    cpd_nodes = {}
    for node in g.nodes:
        if node.props['type'] == 'cpd':
            cpd_nodes[node.props['original_id']] = node

    for c, _ in reaction.compounds:
        if c in cpd_nodes:
            node_ex = graph.Node({
                'id': text_type(rxn_id),
                'shape': 'box',
                'style': 'filled',
                'fillcolor': ACTIVE_COLOR,
                'original_id': [rxn_id],
                'type': 'Ex_rxn'})
            g.add_node(node_ex)

            dir = dir_value(reaction.direction)
            for c1, _ in reaction.left:
                g.add_edge(graph.Edge(
                    cpd_nodes[c1], node_ex, {'dir': dir}))
            for c2, _ in reaction.right:
                g.add_edge(graph.Edge(
                    node_ex, cpd_nodes[c2], {'dir': dir}))
    return g


def add_node_label(graph, detail, model_compoundEntries, model_reactionEntries, reaction_flux):
    """ set label of nodes in graph object,

    Args:
        graph: A graph object, contain a set of nodes.
        detail: Command line option, a list that contains only one element,
            this element is a list of reaction or compound properties name
            defined in the model. e.g. detail = [['id', 'name', 'formula']].
        model_compoundEntries = dict of cpd_id:compound_entries
    """

    for node in graph.nodes:
        if node.props['type'] == 'cpd':
            if detail is not None:
                props =model_compoundEntries[node.props['original_id'].name].properties
                cpd_detail = [i for i in detail[0] if i in props]
                pre_label = '\n'.join(_encode_value(props[value])
                                  for value in cpd_detail if value != 'id')
                if 'id' in detail[0]:
                    label = '{}\n{}'.format(str(node.props['id']), pre_label)
                else:
                    if all(prop not in props for prop in detail[0]):
                        label = str(node.props['id'])
                    else:
                        label = pre_label
            else:
                label = str(node.props['id'])
            node.props['label'] = label

        elif node.props['type'] == 'rxn':
            if len(node.props['original_id']) == 1:
                rxn_id = node.props['original_id'][0]
                if detail is not None:
                    props = model_reactionEntries[rxn_id].properties
                    rxn_detail = list(i for i in detail[0] if i in props)
                    if len(rxn_detail) == 0:
                        label = rxn_id
                    else:
                        label = '\n'.join(_encode_value(props[value])
                                              for value in rxn_detail)
                else:
                    label = rxn_id

                if len(reaction_flux) > 0:
                    if rxn_id in reaction_flux:
                        label = '{}\n{}'.format(label, reaction_flux[rxn_id])
            else:
                if len(reaction_flux) > 0:
                    sum_flux = 0
                    for r in node.props['original_id']:
                        if r in reaction_flux:
                            sum_flux += abs(reaction_flux[r])
                    label = '\n'.join(r for r in node.props['original_id'])
                    if sum_flux != 0:
                        label = '{}\n{}'.format(label, sum_flux)
                else:
                    label = '\n'.join(r for r in node.props['original_id'])
            node.props['label'] = label
        elif node.props['type'] == 'Ex_rxn':
            label = node.props['original_id'][0]
            if len(reaction_flux) > 0:
                if node.props['original_id'][0] in reaction_flux:
                    label = '{}\n{}'.format(
                        label, reaction_flux[node.props['original_id'][0]])
            node.props['label'] = label
    return graph


def set_edge_props_withfba(g, edge_values):
    """set edge values, including direction, width, style."""

    if len(edge_values) > 0:
        value_list = sorted(edge_values.values())
        ninty_percentile = value_list[int(len(value_list) * 0.9) + 1]
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
    if len(edge_values) > 0:
        for edge in g.edges:
            if edge.source.props['type'] == 'cpd':
                rxn_string = ','.join(edge.dest.props['original_id'])
                edge_test = edge.source.props['original_id'], rxn_string
            elif edge.dest.props['type'] == 'cpd':
                rxn_string = ','.join(edge.source.props['original_id'])
                edge_test = edge.dest.props['original_id'], rxn_string

            if edge_test in edge_values:
                edge.props['penwidth'] = pen_width(abs(edge_values[edge_test]))
            else:
                edge.props['style'] = 'dotted'
    return g



