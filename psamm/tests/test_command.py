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

from __future__ import unicode_literals

import sys
import os
import argparse
import codecs
import shutil
import tempfile
from contextlib import contextmanager
import unittest

from six import StringIO, BytesIO

from psamm.command import (main, Command, MetabolicMixin, SolverCommandMixin,
                           CommandError)
from psamm.lpsolver import generic
from psamm.datasource.native import NativeModel

from psamm.commands.chargecheck import ChargeBalanceCommand
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
from psamm.commands.randomsparse import RandomSparseNetworkCommand
from psamm.commands.robustness import RobustnessCommand
from psamm.commands.sbmlexport import SBMLExport
from psamm.commands.search import SearchCommand
from psamm.commands.tableexport import ExportTableCommand


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


class TestCommandMain(unittest.TestCase):
    def setUp(self):
        self._model = NativeModel({
            'name': 'Test model',
            'biomass': 'rxn_1',
            'reactions': [
                {
                    'id': 'rxn_1',
                    'name': 'Reaction 1',
                    'equation': '|A_\u2206[e]| => |B[c]|',
                    'genes': ['gene_1', 'gene_2']
                }, {
                    'id': 'rxn_2_\u03c0',
                    'name': 'Reaction 2',
                    'equation': '|B[c]| <=> |C[e]|',
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
                }
            ],
            'media': [
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

        self._infeasible_model = NativeModel({
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
            model = self._model

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

    def test_invoke_version(self):
        with redirected_stdout():
            with self.assertRaises(SystemExit) as context:
                main(args=['--version'])
            self.assertEqual(context.exception.code, 0)

    def test_run_chargecheck(self):
        self.run_command(ChargeBalanceCommand)

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

    def test_run_gapcheck(self):
        self.run_solver_command(
            GapCheckCommand, requirements={'integer': True})

    def test_run_gapfill(self):
        self.run_solver_command(
            GapFillCommand, requirements={'integer': True})

    def test_run_gapfill_with_blocked(self):
        self.run_solver_command(
            GapFillCommand, ['--compound', 'D[c]'], {'integer': True})

    def test_run_genedelete(self):
        self.run_solver_command(GeneDeletionCommand, ['--gene', 'gene_1'])

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

    def test_run_robustness_with_infeasible(self):
        self.skip_test_if_no_solver()
        with self.assertRaises(SystemExit):
            self.run_solver_command(
                RobustnessCommand, ['rxn_2_\u03c0'],
                model=self._infeasible_model)

    def test_run_sbmlexport(self):
        self.run_command(SBMLExport, target=BytesIO())

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
            ['rxn_1', '|A_\u2206[e]| => |B[c]|', '["gene_1", "gene_2"]',
                'Reaction 1', 'true'],
            ['rxn_2_\u03c0', '|B[c]| <=> |C[e]|',
                'gene_3 or (gene_4 and gene_5)', 'Reaction 2', 'true'],
            ['rxn_3', '|D[c]| => |E[c]|', '', '', 'true'],
        ])

    def test_run_tableexport_compounds(self):
        f = self.run_command(ExportTableCommand, ['compounds'])

        self.assertTableOutputEqual(f.getvalue(), [
            ['id', 'name', 'charge', 'formula', 'in_model'],
            ['A_\u2206', 'Compound A', '0', '', 'true'],
            ['B', 'Compound B', '-1', 'H2O', 'true'],
            ['C', '', '-1', 'O2', 'true']
        ])

    def test_run_tableexport_medium(self):
        f = self.run_command(ExportTableCommand, ['medium'])

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
