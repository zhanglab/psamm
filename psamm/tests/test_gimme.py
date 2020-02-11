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
# Copyright 2020  Keith Dufault-Thompson <keitht547@my.uri.edu>


import unittest

from psamm.metabolicmodel import MetabolicModel
from psamm.database import DictDatabase
from psamm.datasource.reaction import parse_reaction
from psamm.expression import boolean
from psamm.lpsolver import generic
from psamm import fluxanalysis
from psamm.commands import gimme


class TestAddReactions(unittest.TestCase):
    def setUp(self):
        self._database = DictDatabase()
        self._database.set_reaction('rxn_1', parse_reaction('A[e] => B[e]'))
        self._database.set_reaction('rxn_2', parse_reaction('B[e] => C[e]'))
        self._database.set_reaction('rxn_3', parse_reaction('B[e] <=> D[e]'))
        self._database.set_reaction('rxn_4', parse_reaction('C[e] <=> E[e]'))
        self._database.set_reaction('rxn_5', parse_reaction('D[e] => E[e]'))
        self._database.set_reaction('rxn_6', parse_reaction('E[e] =>'))
        self._database.set_reaction('ex_A', parse_reaction('A[e] <=>'))

        self._mm = MetabolicModel.load_model(
            self._database, self._database.reactions)
        self._assoc = {
            'rxn_1': str('gene_1'),
            'rxn_2': str('gene_2'),
            'rxn_3': str('gene_5'),
            'rxn_4': str('gene_3 or gene_4'),
            'rxn_5': str('gene_5 and gene_6')
        }
        self._obj_reaction = 'rxn_6'

        try:
            self._solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_reverse_model(self):
        test_associations = {
            'rxn_1': str('gene_1'),
            'rxn_2': str('gene_2'),
            'rxn_3_forward': str('gene_5'),
            'rxn_3_reverse': str('gene_5'),
            'rxn_4_forward': str('gene_3 or gene_4'),
            'rxn_4_reverse': str('gene_3 or gene_4'),
            'rxn_5': str('gene_5 and gene_6')
        }
        mm_irreversible, reversible_gene_assoc, split_rxns = \
            gimme.make_irreversible(self._mm, self._assoc,
            exclude_list=['ex_A', 'rxn_6'])
        self.assertEqual(split_rxns, set([('rxn_4_forward', 'rxn_4_reverse'),
                                          ('rxn_3_forward', 'rxn_3_reverse')]))
        self.assertEqual(reversible_gene_assoc, test_associations)
        self.assertEqual(set([i for i in mm_irreversible.reactions]),
                         {'rxn_1', 'rxn_2', 'rxn_3_forward', 'rxn_3_reverse',
                          'rxn_4_forward', 'rxn_4_reverse', 'rxn_5',
                          'rxn_6', 'ex_A'})

    def test_parse_transcriptome_file(self):
        f = ['gene\texpression', 'gene_1\t15', 'gene_2\t20',
             'gene_3\t15', 'gene_4\t25', 'gene_5\t10', 'gene_6\t25']
        d = gimme.parse_transcriptome_file(f, 20)
        self.assertEqual(d, {'gene_1': 5.0, 'gene_3': 5.0, 'gene_5': 10.0})

    def test_get_rxn_value(self):
        mm_irreversible, reversible_gene_assoc, split_rxns =\
            gimme.make_irreversible(self._mm, self._assoc,
            exclude_list=['ex_A', 'rxn_6'])
        f = ['gene\texpression', 'gene_1\t15', 'gene_2\t20',
             'gene_3\t15', 'gene_4\t25', 'gene_5\t10', 'gene_6\t25']
        d = gimme.parse_transcriptome_file(f, 20)
        root_single = gimme.get_rxn_value(boolean.Expression(
            reversible_gene_assoc['rxn_1'])._root, d)
        self.assertEqual(root_single, 5)

        root_and = gimme.get_rxn_value(boolean.Expression(
            reversible_gene_assoc['rxn_5'])._root, d)
        self.assertEqual(root_and, 10)

        root_none = gimme.get_rxn_value(boolean.Expression(
            reversible_gene_assoc['rxn_2'])._root, d)
        self.assertEqual(root_none, None)

    def test_gimme_model(self):
        f = ['gene\texpression', 'gene_1\t15', 'gene_2\t20',
             'gene_3\t15', 'gene_4\t25', 'gene_5\t10', 'gene_6\t25']
        threshold_dict = gimme.parse_transcriptome_file(f, 20)
        mm_irreversible, reversible_gene_assoc, split_rxns = \
            gimme.make_irreversible(self._mm, self._assoc,
            exclude_list=['ex_A', 'rxn_6'])
        p = fluxanalysis.FluxBalanceProblem(mm_irreversible, generic.Solver())
        final_model, below_threshold_ids, incon_score = \
            gimme.solve_gimme_problem(p, mm_irreversible, 'rxn_6',
            reversible_gene_assoc, split_rxns, threshold_dict, 20)
        self.assertEqual(incon_score, 100)
        self.assertEqual(final_model, set(['rxn_1', 'rxn_2', 'rxn_4']))
        self.assertEqual(below_threshold_ids, set(['rxn_1', 'rxn_3', 'rxn_5']))
