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
# Copyright 2014-2015  Jon Lund Steffensen <jon_steffensen@uri.edu>

import unittest

from psamm.metabolicmodel import MetabolicModel, FlipableModelView
from psamm.database import DictDatabase
from psamm.reaction import Reaction, Compound, Direction
from psamm.datasource.reaction import parse_reaction
from psamm.commands import tmfa
from psamm.lpsolver import generic
from decimal import Decimal
from psamm.lpsolver import lp


from decimal import Decimal
from six import StringIO

class TestTMFA(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('rxn_1', parse_reaction('=> (2) A[c]'))
        self.database.set_reaction('rxn_2', parse_reaction('A[c] <=> B[c]'))
        self.database.set_reaction('rxn_3', parse_reaction('A[c] => D[e]'))
        self.database.set_reaction('rxn_4', parse_reaction('A[c] => C[c]'))
        self.database.set_reaction('rxn_5', parse_reaction('C[c] => D[e]'))
        self.database.set_reaction('rxn_6', parse_reaction('D[e] =>'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)
        self.mm_irreversible, self.split_reversible, new_lump_rxn_dict = tmfa.make_irreversible(self.model, [], {}, False)
        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')

    def test_mm_irrev(self):
        mm_irreversible, split_reversible, new_lump_rxn_dict = tmfa.make_irreversible(self.model, [], {}, False)
        self.assertTrue('rxn_2_forward' in [i for i in mm_irreversible.reactions])
        self.assertTrue('rxn_2_reverse' in [i for i in mm_irreversible.reactions])

    def test_mm_irrev_exclude(self):
        mm_irreversible, split_reversible, new_lump_rxn_dict = tmfa.make_irreversible(self.model, ['rxn_2'], {}, False)
        self.assertTrue('rxn_2_forward' not in [i for i in mm_irreversible.reactions])
        self.assertTrue('rxn_2_reverse' not in [i for i in mm_irreversible.reactions])

    def test_make_tmfa_problem_fluxvar(self):
        prob, v, zi, dgri, xij, cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.assertTrue(v.has_variable('rxn_1'))
        self.assertTrue(v.has_variable('rxn_2_forward'))
        self.assertFalse(v.has_variable('rxn_2'))

    def test_make_tmfa_problem_dgrvar(self):
        prob, v, zi, dgri, xij, cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.assertTrue(dgri.has_variable('rxn_1'))
        self.assertTrue(dgri.has_variable('rxn_2_forward'))
        self.assertFalse(dgri.has_variable('rxn_2'))

    def test_make_tmfa_problem_zivar(self):
        prob, v, zi, dgri, xij, cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.assertTrue(zi.has_variable('rxn_1'))
        self.assertTrue(zi.has_variable('rxn_2_forward'))
        self.assertFalse(zi.has_variable('rxn_2'))

    def test_make_tmfa_problem_xijvar(self):
        prob, v, zi, dgri, xij, cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.assertTrue(xij.has_variable('A[c]'))
        self.assertTrue(xij.has_variable('B[c]'))
        self.assertFalse(xij.has_variable('A[e]'))
        self.assertFalse(xij.has_variable('D[c]'))
        self.assertTrue(xij.has_variable('D[e]'))
        self.assertFalse(xij.has_variable('rxn1'))

    def test_make_tmfa_problem_cplist(self):
        prob, v, zi, dgri, xij, cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.assertEquals(set(cp_list), set(['A[c]', 'B[c]', 'C[c]', 'D[e]']))

    def test_parsedgr(self):
        f = ['rxn_1\t1\t2', 'rxn_2\t-4\t0.5']
        dgr_dict = tmfa.parse_dgr_file(f, self.mm_irreversible)
        self.assertEquals(dgr_dict, {'rxn_1': (Decimal(1), Decimal(2)), 'rxn_2_forward': (Decimal(-4), Decimal(0.5)), 'rxn_2_reverse':(Decimal(4), Decimal(0.5))})

    def test_parsedgrinvalid(self):
        f = ['rxn_1\t1\t2', 'rxn_2\tNA\tNA']
        dgr_dict = tmfa.parse_dgr_file(f, self.mm_irreversible)
        self.assertEquals(dgr_dict, {'rxn_1': (Decimal(1), Decimal(2))})


class TestSolving(unittest.TestCase):
    def setUp(self):
        self.database = DictDatabase()
        self.database.set_reaction('EX_A', parse_reaction('<=> A[e]'))
        self.database.set_reaction('EX_C', parse_reaction('<=> C[e]'))
        self.database.set_reaction('EX_h2o', parse_reaction('<=> h2o[e]'))
        self.database.set_reaction('EX_h', parse_reaction('<=> h[e]'))
        self.database.set_reaction('rxn1', parse_reaction('C[e] + h[e] <=> C[c] + h[c]'))
        self.database.set_reaction('rxn2', parse_reaction('A[e] + h[e] <=> A[c] + h[c]'))
        self.database.set_reaction('rxn3', parse_reaction('A[c] => B[c] + C[c]'))
        self.database.set_reaction('rxn4', parse_reaction('B[c] => D[c]'))
        self.database.set_reaction('rxn5', parse_reaction('D[c] => F[c]'))
        self.database.set_reaction('rxn6', parse_reaction('B[c] => E[c]'))
        self.database.set_reaction('rxn7', parse_reaction('E[c] => F[c]'))
        self.database.set_reaction('rxn8', parse_reaction('h2o[e] <=> h2o[c]'))
        self.database.set_reaction('bio', parse_reaction('F[c] =>'))
        self.model = MetabolicModel.load_model(
            self.database, self.database.reactions)
        try:
            self.solver = generic.Solver()
        except generic.RequirementsError:
            self.skipTest('Unable to find an LP solver for tests')
        self.dgr_dict = {'rxn1_forward': (0,0), 'rxn1_reverse': (0,0), 'rxn2_forward': (0,0), 'rxn2_reverse': (0,0), 'rxn8_forward': (0,0), 'rxn8_reverse': (0,0),
                        'rxn3': (-5, 0), 'rxn4': (-5, 0), 'rxn5': (-5, 0), 'rxn6':(500, 0), 'rxn7': (500,0)}
        self.trans_param = {'rxn1_forward': (Decimal(1), Decimal(1)), 'rxn1_reverse': (-Decimal(1), -Decimal(1)),
                            'rxn2_forward': (Decimal(1), Decimal(1)), 'rxn2_reverse': (-Decimal(1), -Decimal(1)),
                            'rxn8_forward': (Decimal(0), Decimal(0)), 'rxn8_reverse': (-Decimal(0), -Decimal(0))}
        self.mm_irreversible, self.split_reversible, self.new_lump_rxn_dict = tmfa.make_irreversible(self.model, ['EX_A', 'EX_C', 'EX_h', 'EX_h2o', 'bio'], {}, False)
        self.prob, self.v, self.zi, self.dgri, self.xij, self.cp_list = tmfa.make_tmfa_problem(self.mm_irreversible, self.solver)
        self.prob, self.cpd_xij_dict = tmfa.add_conc_constraints(self.xij, self.prob, {}, self.cp_list, ['h2o[c]', 'h2o[e]'], 'h[c]', 'h[e]', '')
        self.prob = tmfa.add_reaction_constraints(self.prob, self.v, self.zi, self.dgri, self.xij, self.mm_irreversible, [], ['EX_A', 'EX_C', 'EX_h', 'EX_h2o', 'bio'], ['EX_A', 'EX_C', 'EX_H', 'EX_h2o', 'bio'], self.dgr_dict,
                                     {}, self.split_reversible, self.trans_param, [i for i in self.mm_irreversible.reactions], None, ['h2o', 'h2o[e]'], 'h[c]',
                                     'h[e]', '', 275, err_est=False, hamilton=False)


    def test_solvebiomass(self):
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.assertEquals(x, 1000)

    def test_rxnfluxes(self):
        flux_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for i in [i for i in self.mm_irreversible.reactions]:
            x= tmfa.get_var_bound(self.prob, self.v(i), lp.ObjectiveSense.Maximize)
            flux_dict[i] = x
        self.assertTrue(flux_dict['bio'] == 1000)
        self.assertTrue(flux_dict['rxn7'] == 0)
        self.assertTrue(flux_dict['rxn1_forward'] == 0)
        self.assertTrue(flux_dict['rxn4'] == 1000)

    def test_dgrrange(self):
        dgr_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for reaction in ['rxn1_forward', 'rxn4', 'rxn6']:
            dgr_dict[reaction] = (tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Maximize))
        self.assertAlmostEqual(dgr_dict['rxn1_forward'][0], (9.999999974752427e-07, 69.28332553115047)[0], places=6)
        self.assertAlmostEqual(dgr_dict['rxn1_forward'][1], (9.999999974752427e-07, 69.28332553115047)[1], places=6)
        self.assertAlmostEqual(dgr_dict['rxn4'][0], (-39.64166326557521, -9.999999903698154e-07)[0], places=6)
        self.assertAlmostEqual(dgr_dict['rxn4'][1], (-39.64166326557521, -9.999999903698154e-07)[1], places=6)
        self.assertAlmostEqual(dgr_dict['rxn6'][0], (465.3583377344248, 499.999999)[0], places=6)
        self.assertAlmostEqual(dgr_dict['rxn6'][1], (465.3583377344248, 499.999999)[1], places=6)


    def test_cpd_range(self):
        cpd_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for cpd in ['A[c]', 'A[e]', 'D[c]', 'E[c]']:
            cpd_dict[cpd] = (tmfa.get_var_bound(self.prob, self.xij(cpd), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.xij(cpd), lp.ObjectiveSense.Maximize))
        self.assertAlmostEqual(cpd_dict['A[c]'][1], (-11.512925464970229, -3.912023005428146)[1], places=6)
        self.assertAlmostEqual(cpd_dict['A[e]'][1], (-11.512925464970229, -3.912023005428146)[1], places=6)
        self.assertAlmostEqual(cpd_dict['D[c]'][1], (-11.512925464970229, -3.912023005428146)[1], places=6)

    def test_constrainedconc(self):
        cpd_dict = {}
        self.prob.add_linear_constraints(self.xij('A[c]') == -10)
        self.prob.add_linear_constraints(self.xij('B[c]') == -6)
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for cpd in ['A[c]', 'A[e]', 'D[c]', 'E[c]', 'B[c]']:
            cpd_dict[cpd] = (tmfa.get_var_bound(self.prob, self.xij(cpd), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.xij(cpd), lp.ObjectiveSense.Maximize))
        dgr_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for reaction in ['rxn1_forward', 'rxn4', 'rxn6']:
            dgr_dict[reaction] = (tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Maximize))
        self.assertAlmostEqual(cpd_dict['D[c]'][1], (-11.51292546497023, -4.90292494313994)[1], places=6)
        self.assertAlmostEqual(cpd_dict['A[c]'][1], (-10.0, -10.0)[1], places=6)
        self.assertAlmostEqual(cpd_dict['B[c]'][1], (-6.0, -6.0)[1], places=6)
        self.assertAlmostEqual(cpd_dict['E[c]'][1], (-11.512925245555174, -6.000000219415031)[1], places=6)
        self.assertAlmostEqual(dgr_dict['rxn4'][1], (-30.125556943039474, -9.999999974752427e-07)[1], places=6)

    def test_reversibledgri(self):
        dgr_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for reaction in ['rxn1_forward', 'rxn1_reverse', 'rxn2_forward', 'rxn2_reverse']:
            dgr_dict[reaction] = (tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.dgri(reaction), lp.ObjectiveSense.Maximize))

        self.assertEqual(dgr_dict['rxn1_forward'], (-dgr_dict['rxn1_reverse'][1], -dgr_dict['rxn1_reverse'][0]))
        self.assertEqual(dgr_dict['rxn2_forward'], (-dgr_dict['rxn2_reverse'][1], -dgr_dict['rxn2_reverse'][0]))

    def test_reversibleflux(self):
        flux_dict = {}
        x = tmfa.get_var_bound(self.prob, self.v('bio'), lp.ObjectiveSense.Maximize)
        self.prob.add_linear_constraints(self.v('bio') == x)
        for reaction in ['rxn1_forward', 'rxn1_reverse', 'rxn2_forward', 'rxn2_reverse']:
            flux_dict[reaction] = (tmfa.get_var_bound(self.prob, self.v(reaction), lp.ObjectiveSense.Minimize), tmfa.get_var_bound(self.prob, self.v(reaction), lp.ObjectiveSense.Maximize))
        self.assertTrue(flux_dict['rxn1_forward'] == (0, 0))
        self.assertTrue(flux_dict['rxn1_reverse'] == (1000, 1000))
        self.assertTrue(flux_dict['rxn2_forward'] == (1000, 1000))
        self.assertTrue(flux_dict['rxn2_reverse'] == (0, 0))
