#!/usr/bin/env python
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

import unittest
import os
from psamm.formula import Formula, Atom
from psamm.commands import vis
from psamm.reaction import Compound, Reaction, Direction
from collections import defaultdict
from psamm.datasource.reaction import parse_reaction, parse_compound
from psamm.datasource import entry
from psamm import graph
from psamm.datasource.native import NativeModel, ReactionEntry, CompoundEntry
import tempfile


class TestMakeFilterDict(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'fum_c[c]', 'formula': parse_compound('C4H2O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o_c[c]', 'formula': parse_compound('H2O', 'c')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'mal_L_c[c]', 'formula': parse_compound('C4H4O5', 'c')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'q8_c[c] + succ_c[c] => fum_c[c] + q8h2_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'q8_c[c]', 'formula': parse_compound('C49H74O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'q8h2_c[c]', 'formula': parse_compound('C49H76O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'succ_c[c]', 'formula': parse_compound('C4H4O4', 'c')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

        self.cpd_formula = {'fum_c': Formula.parse('C4H2O4'),
                            'h2o_c': Formula.parse('H2O'),
                            'mal_L_c': Formula.parse('C4H4O5'),
                            'q8_c': Formula.parse('C49H74O4'),
                            'q8h2_c': Formula.parse('C49H76O4'),
                            'succ_c': Formula.parse('C4H4O4')}
        self.element = 'C'
        self.method = 'fpp'
        self.exclude_rxns = []
        self.hide_edges = None

    def test1_default_setting(self):
        e1 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e1_res = {'rxn1': [(Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
                  'rxn2': [(Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c'))]}
        self.assertEqual(e1, e1_res)

    def test2_nofpp(self):
        e2 = vis.make_filter_dict(
            self.native, self.mm, 'no-fpp', self.element, self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e2_res = {'rxn1': [
            (Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
            'rxn2': [
                (Compound('q8_c', 'c'), Compound('fum_c', 'c')),
                (Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e2, e2_res)

    def test3_element_hydrogen(self):
        e3 = vis.make_filter_dict(
            self.native, self.mm, self.method, 'H', self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e3_res = {'rxn1': [
            (Compound('fum_c', 'c'), Compound('mal_L_c', 'c')),
            (Compound('h2o_c', 'c'), Compound('mal_L_c', 'c'))],
            'rxn2': [(Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                     (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                     (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e3, e3_res)

    def test4_exclude_rxn(self):
        e4 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, self.cpd_formula,
            self.hide_edges, ['rxn2'])
        e4_res = {'rxn1': [(Compound('fum_c', 'c'),
                            Compound('mal_L_c', 'c'))]}
        self.assertEqual(e4, e4_res)

    def test5_if_missing_formula_fpp(self):
        cpd_formula = {'fum_c': Formula.parse('C4H2O4'),
                       'h2o_c': Formula.parse('H2O'),
                       'mal_L_c': Formula.parse('C4H4O5'),
                       'q8_c': Formula.parse('C49H74O4'),
                       'q8h2_c': Formula.parse('C49H76O4')}
        e5 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e5_res = {'rxn1': [(Compound('fum_c', 'c'),
                            Compound('mal_L_c', 'c'))]}
        self.assertEqual(e5, e5_res)

    def test6_if_missing_formula_nofpp(self):
        cpd_formula = {'fum_c': Formula.parse('C4H2O4'),
                       'h2o_c': Formula.parse('H2O'),
                       'mal_L_c': Formula.parse('C4H4O5'),
                       'q8_c': Formula.parse('C49H74O4'),
                       'q8h2_c': Formula.parse('C49H76O4')}
        e6 = vis.make_filter_dict(
            self.native, self.mm, 'no-fpp', self.element, cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e6_res = {'rxn1': [
            (Compound('fum_c', 'c'), Compound('mal_L_c', 'c'))],
            'rxn2': [(Compound('q8_c', 'c'), Compound('fum_c', 'c')),
                     (Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                     (Compound('succ_c', 'c'), Compound('fum_c', 'c')),
                     (Compound('succ_c', 'c'), Compound('q8h2_c', 'c'))]}
        self.assertEqual(e6, e6_res)

    def test7_hide_edges(self):
        path = os.path.join(tempfile.mkdtemp(), 'remove_edges')
        with open(path, 'w') as f:
            f.write('{}\t{}\n{}\t{}'.format(
                'q8_c[c]', 'q8h2_c[c]', 'fum_c[c]', 'mal_L_c[c]'))
        e7 = vis.make_filter_dict(
            self.native, self.mm, self.method, self.element, self.cpd_formula,
            open(path), self.exclude_rxns)
        e7_res = {'rxn2': [(Compound('succ_c', 'c'), Compound('fum_c', 'c'))]}
        self.assertEqual(e7, e7_res)

    def test8_file_path(self):
        path = os.path.join(tempfile.mkdtemp(), 'primarypairs_prediction')
        with open(path, 'w') as f:
            row1 = '{}\t{}\t{}\t{}'.format('rxn2', 'q8_c[c]',
                                           'q8h2_c[c]', 'C49H74O4')
            row2 = '{}\t{}\t{}\t{}'.format('rxn2', 'succ_c[c]',
                                           'fum_c[c]', 'C4H2O4')
            row3 = '{}\t{}\t{}\t{}'.format('rxn2', 'succ_c[c]',
                                           'q8h2_c[c]', 'H2')
            f.write('\n'.join([row1, row2, row3]))
        e8 = vis.make_filter_dict(
            self.native, self.mm, path, self.element, self.cpd_formula,
            self.hide_edges, self.exclude_rxns)
        e8_res = {'rxn2': [(Compound('q8_c', 'c'), Compound('q8h2_c', 'c')),
                           (Compound('succ_c', 'c'), Compound('fum_c', 'c'))]}
        self.assertEqual(e8, e8_res)


class TestMakeCpairDict(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'rxn1', 'equation': parse_reaction(
                'A[c] + B[c] <=> C[c] + D[c]')}))
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'rxn2', 'equation': parse_reaction(
                'B[c] <=> (2) D[c]')}))
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'rxn3', 'equation': parse_reaction(
                'D[c] => D[e]')}))
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'rxn4', 'equation': parse_reaction(
                'D[e] <=>')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

        self.filter_dict = {'rxn1': [(Compound('A', 'c'), Compound('C', 'c')),
                                     (Compound('B', 'c'), Compound('C', 'c')),
                                     (Compound('B', 'c'), Compound('D', 'c'))],
                            'rxn2': [(Compound('B', 'c'), Compound('D', 'c'))],
                            'rxn3': [(Compound('D', 'c'), Compound('D', 'e'))]}
        self.subset = ['rxn1', 'rxn2', 'rxn3', 'rxn4']
        self.reaction_flux = []
        self.method = 'fpp'

    def test_MakeCpairDict_default(self):
        e1, n1 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                     self.reaction_flux, self.method)
        e1_res = defaultdict(lambda: defaultdict(list))
        e1_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_1')
        e1_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_2')
        e1_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] \
            = ['rxn1_3', 'rxn2_1']
        e1_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].\
            append('rxn3_1')

        n1_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                  'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e1, e1_res)
        self.assertEqual(n1, n1_res)

    def test_MakeCpairDict_nofpp(self):
        e2, n2 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                     self.reaction_flux, 'no-fpp')
        e2_res = defaultdict(lambda: defaultdict(list))
        e2_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'] = ['rxn1']
        e2_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'] = ['rxn1']
        e2_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = \
            ['rxn1', 'rxn2']
        e2_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'] = ['rxn3']

        n2_res = {'rxn1': 'rxn1', 'rxn2': 'rxn2', 'rxn3': 'rxn3'}
        self.assertEqual(e2, e2_res)
        self.assertEqual(n2, n2_res)

    def test_MakeCpairDict_fba(self):
        reaction_flux = {'rxn1': 9.8, 'rxn3': 9.8, 'rxn4': 9.8}
        e3, n3 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                     reaction_flux, self.method)
        e3_res = defaultdict(lambda: defaultdict(list))
        e3_res[(Compound('A', 'c'), Compound('C', 'c'))]['forward'].\
            append('rxn1_1')
        e3_res[(Compound('B', 'c'), Compound('C', 'c'))]['forward'].\
            append('rxn1_2')
        e3_res[(Compound('B', 'c'), Compound('D', 'c'))]['forward'] \
            = ['rxn1_3']
        e3_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn2_1']
        e3_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].\
            append('rxn3_1')

        n3_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                  'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e3, e3_res)
        self.assertEqual(n3, n3_res)

    def test_MakeCpairDict_subset_fba(self):
        reaction_flux = {'rxn1': 9.8, 'rxn3': 9.8, 'rxn4': 9.8}
        subset = ['rxn1', 'rxn2', 'rxn4']
        e4, n4 = vis.make_cpair_dict(self.mm, self.filter_dict, subset,
                                     reaction_flux, self.method)
        e4_r = defaultdict(lambda: defaultdict(list))
        e4_r[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] = ['rxn1_1']
        e4_r[(Compound('B', 'c'), Compound('C', 'c'))]['forward'] = ['rxn1_2']
        e4_r[(Compound('B', 'c'), Compound('D', 'c'))]['forward'] = ['rxn1_3']
        e4_r[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = ['rxn2_1']

        n4_r = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                'rxn2_1': 'rxn2'}
        self.assertEqual(e4, e4_r)
        self.assertEqual(n4, n4_r)

    def test_MakeCpairDict_filterDict(self):
        filter_dict = {'rxn1': [(Compound('A', 'c'), Compound('C', 'c')),
                                (Compound('B', 'c'), Compound('C', 'c')),
                                (Compound('B', 'c'), Compound('D', 'c')),
                                (Compound('A', 'c'), Compound('D', 'c'))],
                       'rxn2': [(Compound('B', 'c'), Compound('D', 'c'))],
                       'rxn3': [(Compound('D', 'c'), Compound('D', 'e'))]}
        e5, n5 = vis.make_cpair_dict(self.mm, filter_dict, self.subset,
                                     self.reaction_flux, self.method)
        e5_res = defaultdict(lambda: defaultdict(list))
        e5_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_1')
        e5_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_2')
        e5_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] = \
            ['rxn1_3', 'rxn2_1']
        e5_res[(Compound('A', 'c'), Compound('D', 'c'))]['both'].\
            append('rxn1_4')
        e5_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].\
            append('rxn3_1')

        n5_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                  'rxn1_4': 'rxn1', 'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e5, e5_res)
        self.assertEqual(n5, n5_res)

    def test_MakeCpairDict_filepath(self):
        e6, n6 = vis.make_cpair_dict(self.mm, self.filter_dict, self.subset,
                                     self.reaction_flux, 'file_path')
        e6_res = defaultdict(lambda: defaultdict(list))
        e6_res[(Compound('A', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_1')
        e6_res[(Compound('B', 'c'), Compound('C', 'c'))]['both'].\
            append('rxn1_2')
        e6_res[(Compound('B', 'c'), Compound('D', 'c'))]['both'] \
            = ['rxn1_3', 'rxn2_1']
        e6_res[(Compound('D', 'c'), Compound('D', 'e'))]['forward'].\
            append('rxn3_1')

        n6_res = {'rxn1_1': 'rxn1', 'rxn1_2': 'rxn1', 'rxn1_3': 'rxn1',
                  'rxn2_1': 'rxn2', 'rxn3_1': 'rxn3'}
        self.assertEqual(e6, e6_res)
        self.assertEqual(n6, n6_res)


class TestMakeEdgeValues(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'FUM', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')}))
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'CYTBD', 'equation': parse_reaction(
                '(2) h_c[c] + (0.5) o2_c[c] + q8h2_c[c] => '
                'h2o_c[c] + (2) h_e[e] + q8_c[c]')}))
        native_model.reactions.add_entry(
            ReactionEntry({'id': 'FRD7', 'equation': parse_reaction(
                'fum_c[c] + q8h2_c[c] => q8_c[c] + succ_c[c]')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

        self.reaction_flux = {'FUM': 5.06, 'CYTBD': 43.60}
        self.compound_formula = {
            'h2o_c': Formula({Atom('O'): 1, Atom('H'): 2}),
            'fum_c': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 2}),
            'h_e': Formula({Atom('O'): 1, Atom('H'): 2}),
            'fum_e': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 2}),
            'o2_c': Formula({Atom('O'): 2}), 'h_c': Formula({Atom('H'): 1}),
            'mal_L_c': Formula({Atom('O'): 5, Atom('C'): 4, Atom('H'): 4}),
            'q8_c': Formula({Atom('O'): 4, Atom('C'): 49, Atom('H'): 74}),
            'q8h2_c': Formula({Atom('O'): 4, Atom('C'): 49, Atom('H'): 76}),
            'succ_c': Formula({Atom('O'): 4, Atom('C'): 4, Atom('H'): 4})
        }
        self.element = 'C'
        self.split_map = False
        self.cpair_dict = defaultdict(lambda: defaultdict(list))
        self.cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'].append('CYTBD_1')
        self.cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'].append('FRD7_1')
        self.cpair_dict[(Compound('fum_c', 'c'), Compound(
            'succ_c', 'c'))]['forward'].append('FRD7_2')
        self.cpair_dict[(Compound('fum_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_1')

        self.new_id_mapping = {'CYTBD_1': 'CYTBD', 'FRD7_1': 'FRD7',
                               'FRD7_2': 'FRD7', 'FUM_1': 'FUM',
                               'FUM_2': 'FUM', 'CYTBD_2': 'CYTBD'}
        self.method = 'fpp'

    def test_edge_values_defaultSetting_or_filepath(self):
        e1_default = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, self.method)
        e1_filepath = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, 'file_path')
        e1_res = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06}
        self.assertEqual(e1_default, e1_res)
        self.assertEqual(e1_filepath, e1_res)

    def test_edge_values_split_map(self):
        e2 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            True, self.cpair_dict, self.new_id_mapping, self.method)
        e2_res = {(Compound('q8_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06
                  }
        self.assertEqual(e2, e2_res)

    def test_edge_values_nofpp(self):
        e3 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, 'no-fpp')
        e3_res = {(Compound('q8_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06
                  }
        self.assertEqual(e3, e3_res)

    def test_edge_values_nofba(self):
        e4 = vis.make_edge_values(
            {}, self.mm, self.compound_formula, self.element,
            self.split_map, self.cpair_dict, self.new_id_mapping, self.method)
        e4_res = {}
        self.assertEqual(e4, e4_res)

    def test_edge_values_elementO(self):
        new_cpair_dict = defaultdict(lambda: defaultdict(list))
        new_cpair_dict[(Compound('o2_c', 'c'), Compound(
            'h2o_c', 'c'))]['forward'].append('CYTBD_1')
        new_cpair_dict[(Compound('q8h2_c', 'c'), Compound(
            'q8_c', 'c'))]['forward'] = ['CYTBD_2', 'FRD7_1']
        new_cpair_dict[(Compound('fum_c', 'c'), Compound(
            'succ_c', 'c'))]['forward'].append('FRD7_2')
        new_cpair_dict[(Compound('fum_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_1')
        new_cpair_dict[(Compound('h2o_c', 'c'), Compound(
            'mal_L_c', 'c'))]['forward'].append('FUM_2')

        e5 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, 'O',
            self.split_map, new_cpair_dict, self.new_id_mapping, self.method)
        e5_res = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06,
                  (Compound('h2o_c', 'c'), 'FUM'): 5.06,
                  (Compound('o2_c', 'c'), 'CYTBD'): 21.80,
                  (Compound('h2o_c', 'c'), 'CYTBD'): 43.60
                  }
        self.assertEqual(e5, e5_res)

    def test_edge_values_nofpp_split(self):
        e6 = vis.make_edge_values(
            self.reaction_flux, self.mm, self.compound_formula, self.element,
            True, self.cpair_dict, self.new_id_mapping, 'no-fpp')
        e6_res = {(Compound('q8_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                  (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
                  (Compound('fum_c', 'c'), 'FUM'): 5.06
                  }
        self.assertEqual(e6, e6_res)


class TestAddGraphNodes(unittest.TestCase):
    def setUp(self):
        cpair_dict = defaultdict(lambda: defaultdict(list))
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['both'] = \
            ['rxn1_1', 'rxn2_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] \
            = ['rxn3_1']
        self.cpair_dict = cpair_dict
        self.new_id_mapping = {'rxn1_1': 'rxn1', 'rxn2_1': 'rxn2',
                               'rxn3_1': 'rxn3'}
        self.method = 'fpp'
        self.split = False
        self.g = graph.Graph()

    def test_addnodes_default(self):
        g1 = vis.add_graph_nodes(self.g, self.cpair_dict, self.method,
                                 self.new_id_mapping, self.split)
        node_a1 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type': 'cpd',
            'label': 'A[c]', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_c1 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'label': 'C[c]', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac1_both = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled',
            'label': 'rxn1\nrxn2', 'type': 'rxn', 'original_id':
                ['rxn1', 'rxn2'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        node_ac1_forward = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled',
            'label': 'rxn3', 'type': 'rxn', 'original_id': ['rxn3'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a1, node_c1, node_ac1_both,
                                  node_ac1_forward] for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in [
            node_a1, node_c1, node_ac1_both, node_ac1_forward]))

    def test_addnodes_split(self):
        g2 = vis.add_graph_nodes(self.g, self.cpair_dict, self.method,
                                 self.new_id_mapping, True)
        node_a2 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'label': 'A[c]', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_c2 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'label': 'C[c]', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'label': 'rxn1', 'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'label': 'rxn2', 'original_id': ['rxn2'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'label': 'rxn3', 'original_id': ['rxn3'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})

        self.assertTrue(all(i in [node_a2, node_c2, node_rxn1, node_rxn2,
                                  node_rxn3] for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in [
            node_a2, node_c2, node_rxn1, node_rxn2, node_rxn3]))

    def test_addnodes_nofpp(self):
        g3 = vis.add_graph_nodes(self.g, self.cpair_dict, 'no-fpp',
                                 self.new_id_mapping, self.split)
        node_a3 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf', 'label': 'A[c]'})
        node_c3 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf', 'label': 'C[c]'})
        node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn1'})

        node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn2'})

        node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn3'})

        self.assertTrue(all(i in [node_a3, node_c3, node_rxn1, node_rxn2,
                                  node_rxn3] for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in [
            node_a3, node_c3, node_rxn1, node_rxn2, node_rxn3]))


class TestUpdateNodeColor(unittest.TestCase):
    def setUp(self):
        node_a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_d = graph.Node({
            'id': 'D[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('D', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_ac = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        node_ad = graph.Node({
            'id': 'rxn1_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})

        node_ac_comb = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1', 'rxn2'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})

        self.g = graph.Graph()
        self.node_list = [node_a, node_c, node_d, node_ac, node_ad]
        for node in self.node_list:
            self.g.add_node(node)

        self.g_comb = graph.Graph()
        self.node_list_comb = [node_a, node_c, node_d, node_ac_comb, node_ad]
        for node in self.node_list_comb:
            self.g_comb.add_node(node)

        self.recolor_dict = {'A[c]': '#f4fc55', 'rxn1': '#ef70de'}

    def test1_add_color_emptyrecolordict(self):
        g1 = vis.update_node_color(self.g, {})
        self.assertTrue(all(i in self.node_list for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in self.node_list))

    def test2_add_color_recolor_splitGraph(self):
        g2 = vis.update_node_color(self.g, self.recolor_dict)
        node_a2 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#f4fc55'})
        node_c2 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        node_ac2 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#ef70de'})
        node_d2 = graph.Node({
            'id': 'D[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('D', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_ad2 = graph.Node({
            'id': 'rxn1_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#ef70de'})
        node_list = [node_a2, node_c2, node_d2, node_ac2, node_ad2]
        self.assertTrue(all(i in node_list for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in node_list))

    def test3_update_color_recolorCombinedRxnNode(self):
        g3 = vis.update_node_color(self.g_comb, self.recolor_dict)
        node_a3 = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#f4fc55'})
        node_c3 = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_d3 = graph.Node({
            'id': 'D[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('D', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        node_ac3 = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1', 'rxn2'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})
        node_ad3 = graph.Node({
            'id': 'rxn1_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#ef70de'})
        node_list = [node_a3, node_c3, node_ac3, node_d3, node_ad3]
        self.assertTrue(all(i in node_list for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in node_list))

    def test4_add_color_emptyrecolordict_combGraph(self):
        g4 = vis.update_node_color(self.g_comb, {})
        self.assertTrue(all(i in self.node_list_comb for i in g4.nodes))
        self.assertTrue(all(i in g4.nodes for i in self.node_list_comb))


class TestAddEdges(unittest.TestCase):
    def setUp(self):
        self.node_a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.node_rxn2 = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.node_rxn3 = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.node_rxn4 = graph.Node({
            'id': 'rxn4_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn4'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})

        self.node_ac_both = graph.Node({
            'id': 'rxn1_1,rxn2_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['rxn1', 'rxn2'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.node_ac_forw = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.node_ac_back = graph.Node({
            'id': 'rxn4_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn4'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})

        self.g1 = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_ac_both,
                     self.node_ac_forw, self.node_ac_back]:
            self.g1.add_node(node)

        self.g2 = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_rxn1,
                     self.node_rxn2, self.node_rxn3, self.node_rxn4]:
            self.g2.add_node(node)

        cpair_dict = defaultdict(lambda: defaultdict(list))
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['both'] \
            = ['rxn1_1', 'rxn2_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['forward'] \
            = ['rxn3_1']
        cpair_dict[(Compound('A', 'c'), Compound('C', 'c'))]['back'] \
            = ['rxn4_1']
        self.cpair_dict = cpair_dict
        self.split = False
        self.method = 'fpp'

    def test1_addedges_default_or_filepath(self):
        g1_default = vis.add_edges(self.g1, self.cpair_dict, 'fpp', False)
        g1_file = vis.add_edges(self.g1, self.cpair_dict, 'file_path', False)
        edge1 = graph.Edge(self.node_a, self.node_ac_both, {'dir': 'both'})
        edge2 = graph.Edge(self.node_ac_both, self.node_c, {'dir': 'both'})
        edge3 = graph.Edge(self.node_a, self.node_ac_forw, {'dir': 'forward'})
        edge4 = graph.Edge(self.node_ac_forw, self.node_c, {'dir': 'forward'})
        edge5 = graph.Edge(self.node_a, self.node_ac_back, {'dir': 'back'})
        edge6 = graph.Edge(self.node_ac_back, self.node_c, {'dir': 'back'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6]
        self.assertTrue(all(i in edge_list for i in g1_default.edges))
        self.assertTrue(all(i in g1_default.edges for i in edge_list))

        self.assertTrue(all(i in edge_list for i in g1_file.edges))
        self.assertTrue(all(i in g1_file.edges for i in edge_list))

    def test2_addedges_split_or_nofpp(self):
        g2_split = vis.add_edges(self.g2, self.cpair_dict, self.method, True)

        g2_nofpp = vis.add_edges(self.g2, self.cpair_dict, 'no-fpp',
                                 self.split)

        edge1 = graph.Edge(self.node_a, self.node_rxn1, {'dir': 'both'})
        edge2 = graph.Edge(self.node_rxn1, self.node_c, {'dir': 'both'})
        edge3 = graph.Edge(self.node_a, self.node_rxn3, {'dir': 'forward'})
        edge4 = graph.Edge(self.node_rxn3, self.node_c, {'dir': 'forward'})
        edge5 = graph.Edge(self.node_a, self.node_rxn4, {'dir': 'back'})
        edge6 = graph.Edge(self.node_rxn4, self.node_c, {'dir': 'back'})
        edge7 = graph.Edge(self.node_a, self.node_rxn2, {'dir': 'both'})
        edge8 = graph.Edge(self.node_rxn2, self.node_c, {'dir': 'both'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8]

        self.assertTrue(all(i in edge_list for i in g2_split.edges))
        self.assertTrue(all(i in g2_split.edges for i in edge_list))

        self.assertTrue(all(i in edge_list for i in g2_nofpp.edges))
        self.assertTrue(all(i in g2_nofpp.edges for i in edge_list))


class TestAddBiomassRnxs(unittest.TestCase):
    def setUp(self):
        self.node_a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('A', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled', 'type':
                'cpd', 'original_id': Compound('C', 'c'), 'compartment': 'c',
            'fillcolor': '#ffd8bf'})
        self.node_rxn1 = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.edge_A_rxn = graph.Edge(self.node_a, self.node_rxn1,
                                     {'dir': 'both'})
        self.edge_C_rxn = graph.Edge(self.node_rxn1, self.node_c,
                                     {'dir': 'both'})
        self.g = graph.Graph()
        for node in [self.node_a, self.node_c, self.node_rxn1]:
            self.g.add_node(node)
        self.g.add_edge(self.edge_A_rxn)
        self.g.add_edge(self.edge_C_rxn)

        self.biomass_id = 'test_bio'
        self.biomass_rxn = Reaction(
            Direction.Forward, [
                (Compound('C', 'c'), 1), (Compound('D', 'c'), 2)],
            [(Compound('A', 'c'), 1), (Compound('E', 'c'), 1)])

    def test1_addBiorxn_default(self):
        g1 = vis.add_biomass_rxns(self.g, self.biomass_rxn, self.biomass_id)
        bio_a = graph.Node({
            'id': 'test_bio_2', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor': '#b3fcb8'})
        bio_c = graph.Node({
            'id': 'test_bio_1', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor': '#b3fcb8'})
        edge_bio_a = graph.Edge(bio_a, self.node_a, {'dir': 'forward'})
        edge_bio_c = graph.Edge(self.node_c, bio_c, {'dir': 'forward'})

        self.assertTrue(all(i in [self.node_a, self.node_c, self.node_rxn1,
                                  bio_a, bio_c] for i in g1.nodes))

        self.assertTrue(all(i in [self.edge_A_rxn, self.edge_C_rxn,
                                  edge_bio_a, edge_bio_c] for i in g1.edges))


class TestAddExchangeRxns(unittest.TestCase):
    def setUp(self):
        a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'c'), 'label': 'A[c]',
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'c'), 'label': 'C[c]',
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.a_extracell = graph.Node({
            'id': 'A[a]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('A', 'e'), 'label': 'A[e]',
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.c_extracell = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('C', 'e'), 'label': 'C[e]',
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        node_ac = graph.Node({
            'id': 'rxn1_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn1'], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn1'})
        node_cc = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn2'], 'compartment': 'e',
            'fillcolor': '#c9fccd', 'label': 'rxn2'})
        node_aa = graph.Node({
            'id': 'rxn3_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['rxn3'], 'compartment': 'c',
            'fillcolor': '#c9fccd', 'label': 'rxn3'})
        edge_a_r1 = graph.Edge(a, node_ac, {'dir': 'forward'})
        edge_r1_c = graph.Edge(node_ac, c, {'dir': 'forward'})
        edge_c_r2 = graph.Edge(c, node_cc, {'dir': 'both'})
        edge_r2_c = graph.Edge(node_cc, self.c_extracell, {'dir': 'both'})
        edge_a_r3 = graph.Edge(self.a_extracell, node_aa, {'dir': 'forward'})
        edge_r3_a = graph.Edge(node_aa, a, {'dir': 'forward'})

        self.g1 = graph.Graph()
        self.edge_list_c_to_e = [edge_a_r1, edge_r1_c, edge_c_r2, edge_r2_c]
        self.node_list_c_to_e = [a, c, self.c_extracell, node_ac, node_cc]
        for node in self.node_list_c_to_e:
            self.g1.add_node(node)
        for edge in self.edge_list_c_to_e:
            self.g1.add_edge(edge)
        self.rxn_C = Reaction(Direction.Both, {Compound('C', 'e'): -1})
        self.rxn_A = Reaction(Direction.Both, {Compound('A', 'e'): -1})

        self.g2 = graph.Graph()
        self.edge_list_e_to_c = [edge_a_r1, edge_r1_c, edge_a_r3, edge_r3_a]
        self.node_list_e_to_c = [a, self.a_extracell, c, node_aa, node_ac]
        for node in self.node_list_e_to_c:
            self.g2.add_node(node)
        for edge in self.edge_list_e_to_c:
            self.g2.add_edge(edge)

    def test1_addExrRxn_cpdC(self):
        g1 = vis.add_exchange_rxns(self.g1, 'test_Ex_C', self.rxn_C)
        node_ex = graph.Node({
            'id': 'test_Ex_C', 'shape': 'box', 'style': 'filled', 'type':
                'Ex_rxn', 'original_id': ['test_Ex_C'], 'compartment': 'e',
            'fillcolor': '#90f998', 'label': 'test_Ex_C'})
        edge_ex = graph.Edge(self.c_extracell, node_ex, {'dir': 'both'})
        self.node_list_c_to_e.append(node_ex)
        self.edge_list_c_to_e.append(edge_ex)
        self.assertTrue(all(i in self.node_list_c_to_e for i in g1.nodes))
        self.assertTrue(all(i in self.edge_list_c_to_e for i in g1.edges))

    def test2_addExrRxn_cpdA(self):
        g2 = vis.add_exchange_rxns(self.g2, 'test_Ex_A', self.rxn_C)
        node_ex = graph.Node({
            'id': 'test_Ex_A', 'shape': 'box', 'style': 'filled',
            'type': 'Ex_rxn', 'original_id': ['test_Ex_A'], 'compartment':
                'e', 'fillcolor': '#90f998', 'label': 'test_Ex_A'})
        edge_ex = graph.Edge(self.a_extracell, node_ex, {'dir': 'both'})
        self.node_list_e_to_c.append(node_ex)
        self.edge_list_e_to_c.append(edge_ex)
        self.assertTrue(all(i in self.node_list_e_to_c for i in g2.nodes))
        self.assertTrue(all(i in self.edge_list_e_to_c for i in g2.edges))


class TestUpdateNodeLabel(unittest.TestCase):
    def setUp(self):
        self.a = graph.Node({
            'id': 'A[c]', 'shape': 'ellipse', 'style': 'filled',
            'label': 'A[c]', 'type': 'cpd', 'original_id': Compound('A', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.c = graph.Node({
            'id': 'C[c]', 'shape': 'ellipse', 'style': 'filled',
            'label': 'C[c]', 'type': 'cpd', 'original_id': Compound('C', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.c_extracell = graph.Node({
            'id': 'C[e]', 'shape': 'ellipse', 'style': 'filled',
            'label': 'C[e]', 'type': 'cpd', 'original_id': Compound('C', 'e'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.node_ac = graph.Node({
            'id': 'rxn1_1,rxn3_1', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'label': 'rxn1_1\nrxn3_1', 'original_id':
                ['rxn1', 'rxn3'], 'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.node_cc = graph.Node({
            'id': 'rxn2_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'label': 'rxn2_1', 'original_id': ['rxn2'], 'compartment': 'e',
            'fillcolor': '#c9fccd'})
        self.node_Ex = graph.Node({
            'id': 'test_Ex_C', 'shape': 'box', 'style': 'filled',
            'type': 'Ex_rxn', 'label': 'test_Ex_C', 'original_id':
                ['test_Ex_C'], 'compartment': 'e', 'fillcolor': '#90f998'})
        self.bio_A = graph.Node({
            'id': 'test_bio_1', 'label': 'test_bio', 'shape': 'box',
            'style': 'filled', 'type': 'bio_rxn', 'original_id': ['test_bio'],
            'compartment': 'c', 'fillcolor': '#b3fcb8'})

        edge_a_r13 = graph.Edge(self.a, self.node_ac, {'dir': 'both'})
        edge_r13_c = graph.Edge(self.node_ac, self.c, {'dir': 'both'})
        edge_c_r2 = graph.Edge(self.c, self.node_cc, {'dir': 'both'})
        edge_r2_c = graph.Edge(self.node_cc, self.c_extracell,
                               {'dir': 'both'})
        edge_ex = graph.Edge(self.c_extracell, self.node_Ex, {'dir': 'both'})
        edge_bio_a = graph.Edge(self.bio_A, self.a, {'dir': 'forward'})

        self.node_list = [self.a, self.c, self.c_extracell, self.node_ac,
                          self.node_cc, self.node_Ex, self.bio_A]
        self.g = graph.Graph()
        for node in self.node_list:
            self.g.add_node(node)
        for edge in [edge_a_r13, edge_r13_c, edge_c_r2, edge_r2_c, edge_ex,
                     edge_bio_a]:
            self.g.add_edge(edge)

        self.cpd_detail = None
        self.rxn_detail = None
        self.reaction_flux = {}

        cpd_a = entry.DictCompoundEntry({
            'id': 'A',
            'name': 'Compound A',
            'formula': 'C6H11O9P',
            'charge': 0})
        cpd_c = entry.DictCompoundEntry({
            'id': 'C',
            'name': 'Compound C',
            'formula': 'C6H11O9P',
            'charge': 0})
        self.cpd_entries = {'A': cpd_a, 'C': cpd_c}

        rxn1 = entry.DictReactionEntry({
            'id': 'rxn1',
            'name': 'Reaction 1',
            'equation': Reaction(Direction.Both, {Compound('A', 'c'): -1,
                                                  Compound('C', 'c'): 1}),
            'genes': 'gene1 and gene2'
        }),
        rxn2 = entry.DictReactionEntry({
            'id': 'rxn2',
            'name': 'Reaction 2',
            'equation': Reaction(Direction.Both, {
                Compound('C', 'c'): -1, Compound('C', 'e'): 1}),
            'genes': 'gene3'})
        rxn3 = entry.DictReactionEntry({
            'id': 'rxn3',
            'name': 'Reaction 3',
            'equation': Reaction(Direction.Both, {Compound('A', 'c'): -1,
                                                  Compound('C', 'c'): 1}),
        })
        test_ex_c = entry.DictReactionEntry({
            'id': 'test_Ex_C',
            'name': 'Exchange Reaction',
            'equation': Reaction(Direction.Both, {Compound('C', 'e'): -1}),
            'genes': 'gene_ex'

        })
        test_bio = entry.DictReactionEntry({
            'id': 'test_bio',
            'name': 'biomass Reaction',
            'equation': Reaction(Direction.Forward, {
                Compound('A', 'c'): -1, Compound('D', 'c'): -1,
                Compound('E', 'c'): 1})})

        self.rxn_entries = {'rxn1': rxn1, 'rxn2': rxn2, 'rxn3': rxn3,
                            'test_Ex_C': test_ex_c, 'test_bio': test_bio}

    def test_detail_default(self):
        g1 = vis.update_node_label(
            self.g, self.cpd_detail, self.rxn_detail, self.cpd_entries,
            self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(
            i in [self.a, self.c, self.c_extracell, self.node_ac,
                  self.node_cc, self.node_Ex, self.bio_A] for i in g1.nodes))
        self.assertTrue(all(i in g1.nodes for i in [
            self.a, self.c, self.c_extracell, self.node_ac, self.node_cc,
            self.node_Ex, self.bio_A]))

    def test_detail_fba(self):
        g2 = vis.update_node_label(
            self.g, self.cpd_detail, self.rxn_detail, self.cpd_entries,
            self.rxn_entries, {'rxn1': 4.86, 'rxn2': 7.2, 'rxn3': 5.29,
                               'test_bio': 0.8})
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3\n10.15'
        self.node_cc.props['label'] = 'rxn2\n7.2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(
            i in [self.a, self.c, self.c_extracell, self.node_ac,
                  self.node_cc, self.node_Ex, self.bio_A] for i in g2.nodes))
        self.assertTrue(all(i in g2.nodes for i in [
            self.a, self.c, self.c_extracell, self.node_ac, self.node_cc,
            self.node_Ex, self.bio_A]))

    def test_detail_id_name(self):
        g3 = vis.update_node_label(self.g, [['id', 'name']], [['id', 'name']],
                                   self.cpd_entries, self.rxn_entries,
                                   self.reaction_flux)
        self.a.props['label'] = 'A[c]\ncompound A'
        self.c.props['label'] = 'C[c]\ncompound C'
        self.c_extracell.props['label'] = 'C[e]\ncompound C'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2\nReaction 2'
        self.node_Ex.props['label'] = 'test_Ex_C\nExchange Reaction'
        self.bio_A.props['label'] = 'test_bio\nbiomass Reaction'
        self.assertTrue(all(
            i in [self.a, self.c, self.c_extracell, self.node_ac,
                  self.node_cc, self.node_Ex, self.bio_A] for i in g3.nodes))
        self.assertTrue(all(i in g3.nodes for i in [
            self.a, self.c, self.c_extracell, self.node_ac, self.node_cc,
            self.node_Ex, self.bio_A]))

    def test_cpd_detail_id_formula_genes(self):
        g4 = vis.update_node_label(self.g, [['id', 'formula', 'genes']],
                                   self.rxn_detail, self.cpd_entries,
                                   self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]\nC6H11O9P'
        self.c.props['label'] = 'C[c]\nC6H11O9P'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2'
        self.node_Ex.props['label'] = 'test_Ex_C'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(
            i in [self.a, self.c, self.c_extracell, self.node_ac,
                  self.node_cc, self.node_Ex, self.bio_A] for i in g4.nodes))
        self.assertTrue(all(i in g4.nodes for i in [
            self.a, self.c, self.c_extracell, self.node_ac, self.node_cc,
            self.node_Ex, self.bio_A]))

    def test_rxn_detail_id_formula_genes(self):
        g5 = vis.update_node_label(
            self.g, self.cpd_detail, [['id', 'formula', 'genes']],
            self.cpd_entries, self.rxn_entries, self.reaction_flux)
        self.a.props['label'] = 'A[c]'
        self.c.props['label'] = 'C[c]'
        self.c_extracell.props['label'] = 'C[e]'
        self.node_ac.props['label'] = 'rxn1\nrxn3'
        self.node_cc.props['label'] = 'rxn2\gene3'
        self.node_Ex.props['label'] = 'test_Ex_C\ngene_ex'
        self.bio_A.props['label'] = 'test_bio'
        self.assertTrue(all(
            i in [self.a, self.c, self.c_extracell, self.node_ac,
                  self.node_cc, self.node_Ex, self.bio_A] for i in g5.nodes))
        self.assertTrue(all(i in g5.nodes for i in [
            self.a, self.c, self.c_extracell, self.node_ac, self.node_cc,
            self.node_Ex, self.bio_A]))


class TestEdgePropsWithFBA(unittest.TestCase):
    def setUp(self):
        self.fum = graph.Node({
            'id': 'fum_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('fum_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.mal = graph.Node({
            'id': 'mal_L_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('mal_L_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.succ = graph.Node({
            'id': 'succ_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('succ_c', 'c'),
            'compartment': 'e', 'fillcolor': '#ffd8bf'})
        self.q8 = graph.Node({
            'id': 'q8_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('q8_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.q8h2 = graph.Node({
            'id': 'q8h2_c[c]', 'shape': 'ellipse', 'style': 'filled',
            'type': 'cpd', 'original_id': Compound('q8h2_c', 'c'),
            'compartment': 'c', 'fillcolor': '#ffd8bf'})
        self.rxn_FUM = graph.Node({
            'id': 'FUM_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FUM'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.FRD7 = graph.Node({
            'id': 'FRD7_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FRD7'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.CYTBD_FRD7 = graph.Node({
            'id': 'CYTBD_1,FRD7_2', 'shape': 'box', 'style': 'filled',
            'type': 'rxn', 'original_id': ['CYTBD', 'FRD7'],
            'compartment': 'c', 'fillcolor': '#c9fccd'})
        self.CYTBD = graph.Node({
            'id': 'CYTBD_1', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['CYTBD'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})
        self.FRD7_2 = graph.Node({
            'id': 'FRD7_2', 'shape': 'box', 'style': 'filled', 'type': 'rxn',
            'original_id': ['FRD7'], 'compartment': 'c',
            'fillcolor': '#c9fccd'})

        self.edge1 = graph.Edge(self.fum, self.rxn_FUM, {'dir': 'both'})
        self.edge2 = graph.Edge(self.rxn_FUM, self.mal, {'dir': 'both'})
        self.edge3 = graph.Edge(self.fum, self.FRD7, {'dir': 'forward'})
        self.edge4 = graph.Edge(self.FRD7, self.succ, {'dir': 'forward'})
        self.edge5 = graph.Edge(self.q8, self.CYTBD_FRD7, {'dir': 'back'})
        self.edge6 = graph.Edge(self.CYTBD_FRD7, self.q8h2, {'dir': 'back'})

        self.edge7 = graph.Edge(self.q8, self.CYTBD, {'dir': 'back'})
        self.edge8 = graph.Edge(self.CYTBD, self.q8h2, {'dir': 'back'})
        self.edge9 = graph.Edge(self.q8, self.FRD7_2, {'dir': 'back'})
        self.edge10 = graph.Edge(self.FRD7_2, self.q8h2, {'dir': 'back'})

        self.edge_list = [self.edge1, self.edge2, self.edge3, self.edge4,
                          self.edge5, self.edge6]
        self.edge_list_2 = [self.edge1, self.edge2, self.edge3, self.edge4,
                            self.edge7, self.edge8, self.edge9, self.edge10]

        self.g1 = graph.Graph()
        for node in [self.fum, self.mal, self.succ, self.rxn_FUM, self.FRD7,
                     self.CYTBD_FRD7, self.q8, self.q8h2]:
            self.g1.add_node(node)
        for i in self.edge_list:
            self.g1.add_edge(i)

        self.g2 = graph.Graph()
        for node in [self.fum, self.mal, self.succ, self.rxn_FUM, self.FRD7,
                     self.q8, self.q8h2, self.CYTBD, self.FRD7_2]:
            self.g2.add_node(node)
        for i in self.edge_list_2:
            self.g2.add_edge(i)

        self.edge_values = {(Compound('q8_c', 'c'), 'CYTBD,FRD7'): 43.60,
                            (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
                            (Compound('q8h2_c', 'c'), 'CYTBD,FRD7'): 43.60,
                            (Compound('fum_c', 'c'), 'FUM'): 5.06
                            }

    def test1_nofba(self):
        g1 = vis.set_edge_props_withfba(self.g1, {})
        self.assertTrue(all(i in g1.edges for i in self.edge_list))
        self.assertTrue(all(i in self.edge_list for i in g1.edges))

    def test2_with_fba(self):
        g2 = vis.set_edge_props_withfba(self.g1, self.edge_values)
        edge1 = graph.Edge(self.fum, self.rxn_FUM,
                           {'dir': 'both', 'penwidth': 1.161})
        edge2 = graph.Edge(self.rxn_FUM, self.mal,
                           {'dir': 'both', 'penwidth': 1.161})
        edge3 = graph.Edge(self.fum, self.FRD7,
                           {'dir': 'forward', 'style': 'dotted'})
        edge4 = graph.Edge(self.FRD7, self.succ,
                           {'dir': 'forward', 'style': 'dotted'})
        edge5 = graph.Edge(self.q8, self.CYTBD_FRD7,
                           {'dir': 'back', 'penwidth': 10})
        edge6 = graph.Edge(self.CYTBD_FRD7, self.q8h2,
                           {'dir': 'back', 'penwidth': 10})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6]
        self.assertTrue(all(i in edge_list for i in g2.edges))
        self.assertTrue(all(i in g2.edges for i in edge_list))

    def test3_with_fba_splitGraph(self):
        edge_values = {
            (Compound('q8_c', 'c'), 'CYTBD'): 43.60,
            (Compound('mal_L_c', 'c'), 'FUM'): 5.06,
            (Compound('q8h2_c', 'c'), 'CYTBD'): 43.60,
            (Compound('fum_c', 'c'), 'FUM'): 5.06}
        g3 = vis.set_edge_props_withfba(self.g2, edge_values)
        edge1 = graph.Edge(self.fum, self.rxn_FUM,
                           {'dir': 'both', 'penwidth': 1.161})
        edge2 = graph.Edge(self.rxn_FUM, self.mal,
                           {'dir': 'both', 'penwidth': 1.161})
        edge3 = graph.Edge(self.fum, self.FRD7,
                           {'dir': 'forward', 'style': 'dotted'})
        edge4 = graph.Edge(self.FRD7, self.succ,
                           {'dir': 'forward', 'style': 'dotted'})
        edge5 = graph.Edge(self.q8, self.CYTBD,
                           {'dir': 'back', 'penwidth': 10})
        edge6 = graph.Edge(self.CYTBD, self.q8h2,
                           {'dir': 'back', 'penwidth': 10})
        edge7 = graph.Edge(self.q8, self.FRD7_2,
                           {'dir': 'back', 'style': 'dotted'})
        edge8 = graph.Edge(self.FRD7_2, self.q8h2,
                           {'dir': 'back', 'style': 'dotted'})
        edge_list = [edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8]

        self.assertTrue(all(i in edge_list for i in g3.edges))
        self.assertTrue(all(i in g3.edges for i in edge_list))


class TestMakeCptTree(unittest.TestCase):
    def test1_if_e_inthemodel(self):
        boundaries = [('c', 'e'), ('c', 'p'), ('e', 'p')]
        extra = 'e'
        c1, e1 = vis.make_cpt_tree(boundaries, extra)
        c1_res = defaultdict(set)
        c1_res['c'].add('e')
        c1_res['c'].add('p')
        c1_res['p'].add('e')
        c1_res['p'].add('c')
        c1_res['e'].add('c')
        c1_res['e'].add('p')
        e1_res = 'e'
        self.assertEqual(c1, c1_res)
        self.assertEqual(e1, e1_res)

    def test2_if_e_not_inthemodel(self):
        boundaries = [('c', 'mi')]
        extra = 'e'
        c2, e2 = vis.make_cpt_tree(boundaries, extra)
        c2_res = defaultdict(set)
        c2_res['c'].add('mi')
        c2_res['mi'].add('c')
        e2_res = 'mi'
        self.assertEqual(c2, c2_res)
        self.assertEqual(e2, e2_res)

    def test3_if_e_not_inthemodel_cpmi(self):
        boundaries = [('c', 'p'), ('c', 'mi')]
        extra = 'e'
        c3, e3 = vis.make_cpt_tree(boundaries, extra)
        c3_res = defaultdict(set)
        c3_res['c'].add('p')
        c3_res['c'].add('mi')
        c3_res['mi'].add('c')
        c3_res['p'].add('c')
        e3_res = 'p'
        self.assertEqual(c3, c3_res)
        self.assertEqual(e3, e3_res)


class TestGetCptBoundaries(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction('A[c] => A[e]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'A[c]', 'formula': parse_compound('formula_A', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'A[e]', 'formula': parse_compound('formula_A', 'e')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction('B[e] <=> B[p]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'B[e]', 'formula': parse_compound('formula_B', 'e')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'B[p]', 'formula': parse_compound('formula_B', 'p')}))
        self.native = native_model

    def test_default_setting(self):
        e1_bound, e1_extra = vis.get_cpt_boundaries(self.native)
        e1_bound_res = set()
        e1_bound_res.add(('c', 'e'))
        e1_bound_res.add(('e', 'p'))
        e1_extra_res = 'e'
        self.assertTrue(all(i in e1_bound for i in e1_bound_res))
        self.assertEqual(e1_extra, e1_extra_res)

    def test2_with_extra_defined(self):
        self.native._properties['extracellular'] = 'p'
        e2_bound, e2_extra = vis.get_cpt_boundaries(self.native)
        e2_bound_res = set()
        e2_bound_res.add(('c', 'e'))
        e2_bound_res.add(('e', 'p'))
        e2_extra_res = 'p'
        self.assertTrue(all(i in e2_bound for i in e2_bound_res))
        self.assertEqual(e2_extra, e2_extra_res)

    def test3_e_notinmodel(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction('A[c] => A[mi]')}))
        e3_bound, e3_extra = vis.get_cpt_boundaries(native_model)
        e3_bound_res = set()
        e3_bound_res.add(('c', 'mi'))
        e3_extra_res = 'e'
        self.assertEqual(e3_bound, e3_bound_res)
        self.assertEqual(e3_extra, e3_extra_res)


class TestMakeSubset(unittest.TestCase):
    def setUp(self):
        native_model = NativeModel()
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn1', 'equation': parse_reaction(
                'fum_c[c] + h2o_c[c] <=> mal_L_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'fum_c[c]', 'formula': parse_compound('C4H2O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'h2o_c[c]', 'formula': parse_compound('H2O', 'c')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'mal_L_c[c]', 'formula': parse_compound('C4H4O5', 'c')}))
        native_model.reactions.add_entry(ReactionEntry({
            'id': 'rxn2', 'equation': parse_reaction(
                'q8_c[c] + succ_c[c] => fum_c[c] + q8h2_c[c]')}))
        native_model.compounds.add_entry(CompoundEntry(
            {'id': 'q8_c[c]', 'formula': parse_compound('C49H74O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'q8h2_c[c]', 'formula': parse_compound('C49H76O4', 'c')}))
        native_model.compounds.add_entry(CompoundEntry({
            'id': 'succ_c[c]', 'formula': parse_compound('C4H4O4', 'c')}))
        self.native = native_model
        self.mm = native_model.create_metabolic_model()

    def test_subset_reactions(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('rxn1', 'rxn2'))
        sub_file = open(path, 'r')
        subset = vis.make_subset(self.mm, sub_file)
        sub_file.close()
        self.assertEqual(subset, {'rxn1', 'rxn2'})

    def test_subset_one_compound(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}'.format('mal_L_c[c]'))
        sub_file = open(path, 'r')
        subset = vis.make_subset(self.mm, sub_file)
        sub_file.close()
        self.assertEqual(subset, {'rxn1'})

    def test_subset_compounds(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('h2o_c[c]', 'succ_c[c]'))
        sub_file = open(path, 'r')
        subset = vis.make_subset(self.mm, sub_file)
        sub_file.close()
        self.assertEqual(subset, {'rxn1', 'rxn2'})

    def test_mixed_subset(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('h2o_c[c]', 'rxn1'))
        sub_file = open(path, 'r')
        with self.assertRaises(ValueError):
            subset = vis.make_subset(self.mm, sub_file)

    def test_not_in_model(self):
        path = os.path.join(tempfile.mkdtemp(), 'subset')
        with open(path, 'w') as f:
            f.write('{}\n{}'.format('h3o_c[c]', 'rxn44'))
        sub_file = open(path, 'r')
        with self.assertRaises(ValueError):
            subset = vis.make_subset(self.mm, sub_file)

    def test_no_subset(self):
        path = None
        subset = vis.make_subset(self.mm, path)
        self.assertEqual(subset, {'rxn1', 'rxn2'})


if __name__ == '__main__':
    unittest.main()