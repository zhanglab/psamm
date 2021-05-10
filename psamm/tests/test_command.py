# -*- coding: utf-8 -*-
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
# Copyright 2015  Jon Lund Steffensen <jon_steffensen@uri.edu>
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth@my.uri.edu>

from __future__ import unicode_literals

import sys
import os
import argparse
import shutil
import tempfile
from contextlib import contextmanager
import unittest

from six import StringIO, BytesIO, text_type

from psamm.command import (main, Command, MetabolicMixin, SolverCommandMixin,
                           CommandError)
from psamm.lpsolver import generic
from psamm.datasource import native, sbml

from psamm.commands.chargecheck import ChargeBalanceCommand
from psamm.commands.duplicatescheck import DuplicatesCheck
from psamm.commands.excelexport import ExcelExportCommand
from psamm.commands.fastgapfill import FastGapFillCommand
from psamm.commands.fba import FluxBalanceCommand
from psamm.commands.fluxcheck import FluxConsistencyCommand
from psamm.commands.fluxcoupling import FluxCouplingCommand
from psamm.commands.formulacheck import FormulaBalanceCommand
from psamm.commands.fva import FluxVariabilityCommand
from psamm.commands.gapcheck import GapCheckCommand
from psamm.commands.gapfill import GapFillCommand
from psamm.commands.genedelete import GeneDeletionCommand
from psamm.commands.masscheck import MassConsistencyCommand
from psamm.commands.primarypairs import PrimaryPairsCommand
from psamm.commands.randomsparse import RandomSparseNetworkCommand
from psamm.commands.robustness import RobustnessCommand
from psamm.commands.sbmlexport import SBMLExport
from psamm.commands.search import SearchCommand
from psamm.commands.tableexport import ExportTableCommand
from psamm.commands.vis import VisualizationCommand


@contextmanager
def redirected_stdout(target=None):
    stdout = sys.stdout
    if target is None:
        target = StringIO()
    try:
        sys.stdout = target
        yield target
    finally:
        sys.stdout = stdout


class MockCommand(Command):
    """Test command.

    This is a command for testing basic commands.
    """
    init_parser_called = False
    run_called = False
    has_native_model = False

    @classmethod
    def init_parser(cls, parser):
        cls.init_parser_called = True
        parser.add_argument('--test-argument', action='store_true')
        super(MockCommand, cls).init_parser(parser)

    def __init__(self, *args, **kwargs):
        super(MockCommand, self).__init__(*args, **kwargs)

    def run(self):
        self.__class__.run_called = True
        self.__class__.has_native_model = hasattr(self, '_model')


class MockMetabolicCommand(MetabolicMixin, Command):
    """Test metabolic model command."""
    has_metabolic_model = False

    def run(self):
        self.__class__.has_metabolic_model = hasattr(self, '_mm')


class MockSolverCommand(SolverCommandMixin, Command):
    """Test solver command.

    This is a command for testing solver commands.
    """
    def run(self):
        solver = self._get_solver()
        print(solver)


class BaseCommandTest(object):
    """Generic methods used for different test cases.

    This does not inherit from TestCase as it should not be run as a test
    case. Use as a mixin from actual TestCase subclasses. Implement the
    method get_default_model() to have the run_command() methods use a
    default model.
    """
    def assertTableOutputEqual(self, output, table):
        self.assertEqual(
            output, '\n'.join('\t'.join(row) for row in table) + '\n')

    def is_solver_available(self, **kwargs):
        try:
            generic.Solver(**kwargs)
        except generic.RequirementsError:
            return False

        return True

    def skip_test_if_no_solver(self, **kwargs):
        if not self.is_solver_available(**kwargs):
            self.skipTest('No solver available')

    def run_command(self, command_class, args=[], target=None, model=None):
        parser = argparse.ArgumentParser()
        command_class.init_parser(parser)
        parsed_args = parser.parse_args(args)

        if model is None:
            model = self.get_default_model()

        command = command_class(model, parsed_args)

        with redirected_stdout(target=target) as f:
            command.run()

        return f

    def run_solver_command(
            self, command_class, args=[], requirements={}, **kwargs):
        with redirected_stdout() as f:
            if self.is_solver_available(**requirements):
                self.run_command(command_class, args, **kwargs)
            else:
                with self.assertRaises(generic.RequirementsError):
                    self.run_command(command_class, args, **kwargs)

        return f


class TestCommandMain(unittest.TestCase, BaseCommandTest):
    def setUp(self):
        reader = native.ModelReader({
            'name': 'Test model',
            'biomass': 'rxn_1',
            'compartments': [
                {
                    'id': 'e',
                    'name': 'Extracellular',
                    'adjacent_to': ['c']
                },
                {
                    'id': 'c',
                    'name': 'Cytosol',
                    'adjacent_to': 'e'
                }
            ],
            'reactions': [
                {
                    'id': 'rxn_1',
                    'name': 'Reaction 1',
                    'equation': '|A_\u2206[e]| => |B[c]|',
                    'genes': ['gene_1', 'gene_2']
                }, {
                    'id': 'rxn_2_\u03c0',
                    'name': 'Reaction 2',
                    'equation': 'atp[c] + (2) |B[c]| <=> adp[c] + |C[e]|',
                    'genes': 'gene_3 or (gene_4 and gene_5)'
                }, {
                    'id': 'rxn_3',
                    'equation': 'D[c] => E[c]'
                }
            ],
            'compounds': [
                {
                    'id': 'A_\u2206',
                    'name': 'Compound A',
                    'charge': 0
                },
                {
                    'id': 'B',
                    'name': 'Compound B',
                    'formula': 'H2O',
                    'charge': -1
                },
                {
                    'id': 'C',
                    'charge': -1,
                    'formula': 'O2'
                },
                {
                    'id': 'atp',
                    'name': '\u2192 ATP ',
                    'charge': -4,
                    'formula': 'C10H12N5O13P3'
                },
                {
                    'id': 'adp',
                    'name': ' ADP',
                    'charge': -3,
                    'formula': 'C10H12N5O10P2'
                }
            ],
            'exchange': [
                {
                    'compartment': 'e',
                    'compounds': [
                        {'id': 'A_\u2206'},
                        {'id': 'C'}
                    ]
                }
            ],
            'limits': [
                {
                    'reaction': 'rxn_2_\u03c0',
                    'upper': 100
                }
            ]
        })
        self._model = reader.create_model()

        reader = native.ModelReader({
            'biomass': 'rxn_1',
            'reactions': [
                {
                    'id': 'rxn_1',
                    'equation': 'A[c] => B[c]'
                }
            ],
            'limits': [
                {
                    'reaction': 'rxn_1',
                    'fixed': 10
                }
            ]
        })
        self._infeasible_model = reader.create_model()

    def get_default_model(self):
        return self._model

    def test_invoke_version(self):
        with redirected_stdout():
            with self.assertRaises(SystemExit) as context:
                main(args=['--version'])
            self.assertEqual(context.exception.code, 0)

    def test_run_chargecheck(self):
        self.run_command(ChargeBalanceCommand)

    def test_run_duplicatescheck(self):
        self.run_command(DuplicatesCheck)

    def test_run_duplicatescheck_compare_stoichiometry(self):
        self.run_command(DuplicatesCheck, ['--compare-stoichiometry'])

    def test_run_duplicatescheck_compare_direction(self):
        self.run_command(DuplicatesCheck, ['--compare-direction'])

    def test_run_excelexport(self):
        dest = tempfile.mkdtemp()
        try:
            dest_path = os.path.join(dest, 'model.xlsx')
            self.run_command(ExcelExportCommand, [dest_path])
            self.assertTrue(os.path.isfile(dest_path))
        finally:
            shutil.rmtree(dest)

    def test_run_fastgapfill(self):
        # Skip test if solver is GLPK because of issue #61.
        try:
            solver = generic.Solver()
            if solver.properties['name'] == 'glpk':
                self.skipTest('Test has known issue with GLPK')
        except generic.RequirementsError:
            pass

        self.run_solver_command(FastGapFillCommand)

    def test_run_fba(self):
        self.run_solver_command(FluxBalanceCommand)

    def test_run_fba_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxBalanceCommand, model=self._infeasible_model)

    def test_run_fba_show_all_reactions(self):
        self.run_solver_command(FluxBalanceCommand, ['--all-reactions'])

    def test_run_fba_with_tfba(self):
        self.run_solver_command(
            FluxBalanceCommand, ['--loop-removal', 'tfba'], {'integer': True})

    def test_run_fba_with_tfba_infeasible(self):
        self.skip_test_if_no_solver(integer=True)
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxBalanceCommand, ['--loop-removal=tfba'], {'integer': True},
                model=self._infeasible_model)

    def test_run_fba_with_l1min(self):
        self.run_solver_command(
            FluxBalanceCommand, ['--loop-removal', 'l1min'])

    def test_run_fba_with_l1min_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxBalanceCommand, ['--loop-removal=l1min'],
                model=self._infeasible_model)

    def test_run_fluxcheck(self):
        self.run_solver_command(FluxConsistencyCommand)

    def test_run_fluxcheck_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxConsistencyCommand, model=self._infeasible_model)

    def test_run_fluxcheck_with_unrestricted(self):
        self.run_solver_command(FluxConsistencyCommand, ['--unrestricted'])

    def test_run_fluxcheck_with_fastcore(self):
        self.run_solver_command(FluxConsistencyCommand, ['--fastcore'])

    def test_run_fluxcheck_with_fastcore_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxConsistencyCommand, ['--fastcore'],
                model=self._infeasible_model)

    def test_run_fluxcheck_with_tfba(self):
        self.run_solver_command(
            FluxConsistencyCommand, ['--loop-removal=tfba'], {'integer': True})

    def test_run_fluxcheck_with_reduce_lp(self):
        self.run_solver_command(FluxConsistencyCommand, ['--reduce-lp'])

    def test_run_fluxcheck_with_reduce_lp_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxConsistencyCommand, ['--reduce-lp'],
                model=self._infeasible_model)

    def test_run_fluxcheck_with_both_tfba_and_fastcore(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(CommandError):
            self.run_solver_command(
                FluxConsistencyCommand, ['--loop-removal=tfba', '--fastcore'],
                {'integer': True})

    def test_run_fluxcoupling(self):
        self.run_solver_command(FluxCouplingCommand)

    def test_run_fluxcoupling_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxCouplingCommand, model=self._infeasible_model)

    def test_run_formulacheck(self):
        self.run_command(FormulaBalanceCommand)

    def test_run_fva(self):
        self.run_solver_command(FluxVariabilityCommand)

    def test_run_fva_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                FluxVariabilityCommand, model=self._infeasible_model)

    def test_run_fva_with_tfba(self):
        self.run_solver_command(
            FluxVariabilityCommand, ['--loop-removal=tfba'], {'integer': True})

    def test_run_gapcheck_prodcheck(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=prodcheck'])

    def test_run_gapcheck_prodcheck_without_implicit_sinks(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=prodcheck', '--no-implicit-sinks'])

    def test_run_gapcheck_prodcheck_without_extracellular(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=prodcheck', '--exclude-extracellular'])

    def test_run_gapcheck_prodcheck_with_unrestricted_exchange(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=prodcheck', '--unrestricted-exchange'])

    def test_run_gapcheck_sinkcheck(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=sinkcheck'])

    def test_run_gapcheck_sinkcheck_without_implicit_sinks(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=sinkcheck', '--no-implicit-sinks'])

    def test_run_gapcheck_gapfind(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=gapfind'], {'integer': True})

    def test_run_gapcheck_gapfind_without_implicit_sinks(self):
        self.run_solver_command(
            GapCheckCommand, ['--method=gapfind', '--no-implicit-sinks'],
            {'integer': True})

    def test_run_gapfill(self):
        self.run_solver_command(
            GapFillCommand, requirements={'integer': True})

    def test_run_gapfill_with_blocked(self):
        self.run_solver_command(
            GapFillCommand, ['--compound', 'D[c]'], {'integer': True})

    def test_run_genedelete(self):
        self.run_solver_command(GeneDeletionCommand, ['--gene', 'gene_1'])

    def test_run_genedelete_with_fba(self):
        self.run_solver_command(
            GeneDeletionCommand, ['--gene=gene_1', '--method=fba'])

    def test_run_genedelete_with_lin_moma(self):
        self.run_solver_command(
            GeneDeletionCommand, ['--gene=gene_1', '--method=lin_moma'])

    def test_run_genedelete_with_lin_moma2(self):
        self.run_solver_command(
            GeneDeletionCommand, ['--gene=gene_1', '--method=lin_moma2'])

    def test_run_genedelete_with_moma(self):
        self.run_solver_command(
            GeneDeletionCommand, ['--gene=gene_1', '--method=moma'],
            {'quadratic': True})

    def test_run_genedelete_with_moma2(self):
        self.run_solver_command(
            GeneDeletionCommand, ['--gene=gene_1', '--method=moma2'],
            {'quadratic': True})

    def test_run_genedelete_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                GeneDeletionCommand, ['--gene=gene_1'],
                model=self._infeasible_model)

    def test_run_masscheck_compounds(self):
        self.run_solver_command(MassConsistencyCommand, ['--type', 'compound'])

    def test_run_masscheck_reactions(self):
        self.run_solver_command(MassConsistencyCommand, ['--type', 'reaction'])

    def test_run_masscheck_reaction_with_checked(self):
        self.run_solver_command(
            MassConsistencyCommand, ['--type=reaction', '--checked=rxn_3'])

    def test_run_primarypairs_with_fpp(self):
        self.run_command(PrimaryPairsCommand, ['--method', 'fpp'])

    def test_run_primarypairs_with_fpp_and_report_element(self):
        self.run_command(
            PrimaryPairsCommand, ['--method', 'fpp', '--report-element', 'C'])

    def test_run_primarypairs_with_fpp_and_report_all_transfers(self):
        self.run_command(
            PrimaryPairsCommand, ['--method', 'fpp', '--report-all-transfers'])

    def test_run_primarypairs_with_fpp_and_weight(self):
        self.run_command(
            PrimaryPairsCommand, [
                '--method', 'fpp', '--weights', 'C=1,H=0,R=0.5,*=0.4'])

    def test_run_primarypairs_with_mapmaker(self):
        self.run_solver_command(
            PrimaryPairsCommand, ['--method', 'mapmaker'], {'integer': True})

    def test_run_randomsparse_reactions(self):
        self.run_solver_command(
            RandomSparseNetworkCommand, ['--type=reactions', '50%'])

    def test_run_randomsparse_genes(self):
        self.run_solver_command(
            RandomSparseNetworkCommand, ['--type=genes', '50%'])

    def test_run_randomsparse_exchange(self):
        self.run_solver_command(
            RandomSparseNetworkCommand, ['--type=exchange', '50%'])

    def test_run_randomsparse_reactions_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                RandomSparseNetworkCommand, ['--type=reactions', '50%'],
                model=self._infeasible_model)

    def test_run_robustness(self):
        self.run_solver_command(RobustnessCommand, ['rxn_2_\u03c0'])

    def test_run_robustness_all_reactions(self):
        self.run_solver_command(
            RobustnessCommand, ['--all-reaction-fluxes', 'rxn_2_\u03c0'])

    def test_run_robustness_with_tfba(self):
        self.run_solver_command(
            RobustnessCommand, ['--loop-removal=tfba', 'rxn_2_\u03c0'],
            {'integer': True})

    def test_run_robustness_with_l1min(self):
        self.run_solver_command(
            RobustnessCommand, ['--loop-removal=l1min', 'rxn_2_\u03c0'])

    def test_run_robustness_with_fva(self):
        self.run_solver_command(
            RobustnessCommand, ['--fva', 'rxn_2_\u03c0']
        )

    def test_run_robustness_with_invalid_objective(self):
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                RobustnessCommand, ['--objective=rxn_4', 'rxn_2_\u03c0'],
                model=self._model)

    def test_run_robustness_with_neg_steps(self):
        with self.assertRaises(CommandError):
            self.run_solver_command(
                RobustnessCommand, ['--steps=0', 'rxn_2_\u03c0'])

    def test_run_robustness_with_loop_removal_fva(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                RobustnessCommand, ['--loop-removal=tfba', '--fva',
                                    'rxn_2_\u03c0'])

    def test_run_robustness_min_greater(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(CommandError):
            self.run_solver_command(
                RobustnessCommand, ['--minimum=50', '--maximum=20',
                                    'rxn_2_\u03c0'])

    def test_run_robustness_with_infeasible(self):
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                RobustnessCommand, ['rxn_2_\u03c0'],
                model=self._infeasible_model)

    def test_run_sbmlexport(self):
        dest = tempfile.mkdtemp()
        try:
            dest_path = os.path.join(dest, 'model.xml')
            self.run_command(SBMLExport, [dest_path])
            self.assertTrue(os.path.isfile(dest_path))
        finally:
            shutil.rmtree(dest)

    def test_run_search_compound(self):
        self.run_command(SearchCommand, ['compound', '--id', 'A_\u2206'])

    def test_run_search_compound_name(self):
        self.run_command(SearchCommand, ['compound', '--name', 'Compound A'])

    def test_run_search_reaction(self):
        self.run_command(SearchCommand, ['reaction', '--id', 'rxn_1'])

    def test_run_search_reaction_by_compound(self):
        self.run_command(SearchCommand, ['reaction', '--compound=A'])

    def test_run_tableexport_reactions(self):
        f = self.run_command(ExportTableCommand, ['reactions'])

        self.assertTableOutputEqual(f.getvalue(), [
            ['id', 'equation', 'genes', 'name', 'in_model'],
            ['rxn_1', 'A_\u2206[e] => B[c]', '["gene_1", "gene_2"]',
                'Reaction 1', 'true'],
            ['rxn_2_\u03c0', 'atp[c] + (2) B[c] <=> adp[c] + C[e]',
                'gene_3 or (gene_4 and gene_5)', 'Reaction 2', 'true'],
            ['rxn_3', 'D[c] => E[c]', '', '', 'true'],
        ])

    def test_run_tableexport_compounds(self):
        f = self.run_command(ExportTableCommand, ['compounds'])
        self.assertTableOutputEqual(f.getvalue(), [
            ['id', 'name', 'charge', 'formula', 'in_model'],
            ['A_\u2206', 'Compound A', '0', '', 'true'],
            ['B', 'Compound B', '-1', 'H2O', 'true'],
            ['C', '', '-1', 'O2', 'true'],
            ['atp', '\u2192 ATP ', '-4', 'C10H12N5O13P3', 'true'],
            ['adp', ' ADP', '-3', 'C10H12N5O10P2', 'true']
        ])

    def test_run_tableexport_exchange(self):
        f = self.run_command(ExportTableCommand, ['exchange'])

        self.assertTableOutputEqual(f.getvalue(), [
            ['Compound ID', 'Reaction ID', 'Lower Limit', 'Upper Limit'],
            ['A_\u2206[e]', '', '-1000', '1000'],
            ['C[e]', '', '-1000', '1000']
        ])

    def test_run_tableexport_limits(self):
        f = self.run_command(ExportTableCommand, ['limits'])

        self.assertTableOutputEqual(f.getvalue(), [
            ['Reaction ID', 'Lower Limits', 'Upper Limits'],
            ['rxn_2_\u03c0', '', '100']
        ])

    def test_run_tableexport_metadata(self):
        f = self.run_command(ExportTableCommand, ['metadata'])

        self.assertTableOutputEqual(f.getvalue(), [
            ['Model Name', 'Test model'],
            ['Biomass Reaction', 'rxn_1'],
            ['Default Flux Limits', '1000']
        ])

    def test_run_vis(self):
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_element_all(self):
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                "--element", "all"])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_hide_edges(self):
        path = os.path.join(tempfile.mkdtemp(), 'hide_edges.csv')
        with open(path, 'w') as f:
            f.write('{}\t{}'.format('D[c]', 'E[c]'))

        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                "--hide-edges", path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_recolor(self):
        path = os.path.join(tempfile.mkdtemp(), 'color.csv')
        with open(path, 'w') as f:
            f.write('{}\t{}'.format('B', '#f4fc55'))
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                "--color", path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_output(self):
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'test')
        dest_path_dot = os.path.join(dest, 'test.dot')
        dest_path_nodes = os.path.join(dest, 'test.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'test.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_compartment(self):
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                '--compartment'])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_fba(self):
        path = os.path.join(tempfile.mkdtemp(), 'fba.tsv')
        with open(path, 'w') as f:
            f.write('{}\t{}'.format('rxn_1', -0.000001))
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                '--fba', path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_fba_invalid_flux(self):
        path = os.path.join(tempfile.mkdtemp(), 'fba.tsv')
        with open(path, 'w') as f:
            f.write('{}\t{}'.format('rxn_1', 'a'))
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                '--fba', path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_fva(self):
        path = os.path.join(tempfile.mkdtemp(), 'fva.tsv')
        with open(path, 'w') as f:
            f.write('{}\t{}\t{}'.format('rxn_1', -0.000001, 0.000001))
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                '--fva', path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_run_vis_fva_invalid_flux(self):
        path = os.path.join(tempfile.mkdtemp(), 'fva.tsv')
        with open(path, 'w') as f:
            f.write('{}\t{}\t{}'.format('rxn_1', 'a', -0.000001))
        dest = tempfile.mkdtemp()
        dest_path = os.path.join(dest, 'reactions')
        dest_path_dot = os.path.join(dest, 'reactions.dot')
        dest_path_nodes = os.path.join(dest, 'reactions.nodes.tsv')
        dest_path_edges = os.path.join(dest, 'reactions.edges.tsv')
        self.run_command(VisualizationCommand, ['--output', dest_path,
                                                '--fva', path])
        self.assertTrue(os.path.isfile(dest_path_dot))
        self.assertTrue(os.path.isfile(dest_path_nodes))
        self.assertTrue(os.path.isfile(dest_path_edges))

    def test_command_main(self):
        self.run_command(MockCommand)
        self.assertTrue(MockCommand.run_called)
        self.assertTrue(MockCommand.init_parser_called)
        self.assertTrue(MockCommand.has_native_model)

    def test_solver_command_main(self):
        with self.assertRaises(generic.RequirementsError):
            self.run_command(
                MockSolverCommand, ['--solver', 'name=not-an-actual-solver'])

    def test_metabolic_command_main(self):
        self.run_command(MockMetabolicCommand)
        self.assertTrue(MockMetabolicCommand.has_metabolic_model)

    def test_command_fail(self):
        mock_command = MockCommand(self._model, None)
        with self.assertRaises(SystemExit):
            mock_command.fail('Run time error')

    def test_command_argument_error(self):
        mock_command = MockCommand(self._model, None)
        with self.assertRaises(CommandError):
            mock_command.argument_error(msg='Argument error')


class TestSBMLCommandMain(unittest.TestCase, BaseCommandTest):
    def setUp(self):
        doc = BytesIO('''<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"
      xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version1"
      xmlns:html="http://www.w3.org/1999/xhtml"
      level="3" version="1"
      fbc:required="false">
 <model id="test_model" name="Test model">
  <listOfCompartments>
   <compartment id="C_c" name="cell" constant="true"/>
   <compartment id="C_b" name="boundary" constant="true"/>
  </listOfCompartments>
  <listOfSpecies>
   <species id="M_Glucose" name="Glucose" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false" fbc:charge="0" fbc:chemicalFormula="C6H12O6"/>
   <species id="M_Glucose_6_P" name="Glucose-6-P" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false" fbc:charge="-2" fbc:chemicalFormula="C6H11O9P"/>
   <species id="M_H2O" name="H2O" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false" fbc:charge="0" fbc:chemicalFormula="H2O"/>
   <species id="M_Phosphate" name="Phosphate" compartment="C_c" constant="false" boundaryCondition="false" hasOnlySubstanceUnits="false" fbc:charge="-2" fbc:chemicalFormula="HO4P"/>
   <species id="M_Biomass" name="Biomass" compartment="C_b" constant="false" boundaryCondition="true" hasOnlySubstanceUnits="false"/>
  </listOfSpecies>
  <listOfReactions>
   <reaction id="R_G6Pase" reversible="true" fast="false">
    <listOfReactants>
     <speciesReference species="M_Glucose" stoichiometry="2" constant="true"/>
     <speciesReference species="M_Phosphate" stoichiometry="2" constant="true"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_H2O" stoichiometry="2" constant="true"/>
     <speciesReference species="M_Glucose_6_P" stoichiometry="2" constant="true"/>
    </listOfProducts>
   </reaction>
   <reaction id="R_Biomass" reversible="false" fast="false">
    <listOfReactants>
     <speciesReference species="M_Glucose_6_P" stoichiometry="0.56" constant="true"/>
    </listOfReactants>
    <listOfProducts>
     <speciesReference species="M_Biomass" stoichiometry="1" constant="true"/>
    </listOfProducts>
   </reaction>
  </listOfReactions>
  <fbc:listOfObjectives fbc:activeObjective="obj1">
   <fbc:objective fbc:id="obj1" fbc:name="Objective 1" fbc:type="maximize">
    <fbc:listOfFluxObjectives>
     <fbc:fluxObjective fbc:reaction="R_Biomass" fbc:coefficient="1"/>
    </fbc:listOfFluxObjectives>
   </fbc:objective>
  </fbc:listOfObjectives>
  <fbc:listOfFluxBounds>
   <fbc:fluxBound fbc:reaction="R_G6Pase" fbc:operation="greaterEqual" fbc:value="-10"/>
   <fbc:fluxBound fbc:reaction="R_G6Pase" fbc:operation="lessEqual" fbc:value="1000"/>
   <fbc:fluxBound fbc:reaction="R_Biomass" fbc:operation="greaterEqual" fbc:value="0"/>
   <fbc:fluxBound fbc:reaction="R_Biomass" fbc:operation="lessEqual" fbc:value="1000"/>
  </fbc:listOfFluxBounds>
 </model>
</sbml>'''.encode('utf-8'))

        reader = sbml.SBMLReader(doc)
        self._model = reader.create_model()

    def get_default_model(self):
        return self._model

    def test_run_fba(self):
        self.run_solver_command(FluxBalanceCommand)


if __name__ == '__main__':
    unittest.main()
