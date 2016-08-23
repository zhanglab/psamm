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
        self._model_dir = tempfile.mkdtemp()
        path = os.path.join(self._model_dir, 'model.yaml')
        with codecs.open(path, 'w', 'utf-8') as f:
            f.write('''---
                name: Test model
                biomass: rxn_1
                reactions:
                  - id: rxn_1
                    equation: '|A_\u2206[e]| => |B[c]|'
                    genes:
                      - gene_1
                      - gene_2
                  - id: rxn_2_\u03c0
                    equation: '|B[c]| => |C[e]|'
                    genes: 'gene_3 or (gene_4 and gene_5)'
                compounds:
                  - id: A_\u2206
                  - id: B
                  - id: C
                media:
                  - compartment: e
                    compounds:
                      - id: A_\u2206
                      - id: C
                limits:
                  - reaction: rxn_2_\u03c0
                    upper: 100
                ''')
        self._model = NativeModel.load_model_from_path(self._model_dir)

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def is_solver_available(self, **kwargs):
        try:
            generic.Solver(**kwargs)
        except generic.RequirementsError:
            return False

        return True

    def run_solver_command(self, args, requirements):
        with redirected_stdout() as f:
            if self.is_solver_available(**requirements):
                main(args=args)
            else:
                with self.assertRaises(generic.RequirementsError):
                    main(args=args)

        return f

    def test_invoke_version(self):
        with redirected_stdout():
            with self.assertRaises(SystemExit) as context:
                main(args=['--version'])
            self.assertEqual(context.exception.code, 0)

    def test_run_chargecheck(self):
        with redirected_stdout():
            main(args=['--model', self._model_dir, 'chargecheck'])

    def test_run_excelexport(self):
        dest_path = os.path.join(self._model_dir, 'model.xlsx')
        with redirected_stdout():
            main(args=['--model', self._model_dir, 'excelexport', dest_path])

        self.assertTrue(os.path.isfile(dest_path))

    def test_run_fastgapfill(self):
        # Skip test if solver is GLPK because of issue #61.
        try:
            solver = generic.Solver()
            if solver.properties['name'] == 'glpk':
                self.skipTest('Test has known issue with GLPK')
        except generic.RequirementsError:
            pass

        self.run_solver_command([
            '--model', self._model_dir, 'fastgapfill'], {})

    def test_run_fba(self):
        self.run_solver_command(['--model', self._model_dir, 'fba'], {})

    def test_run_fba_with_tfba(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fba', '--loop-removal', 'tfba'],
            {'integer': True})

    def test_run_fba_with_l1min(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fba', '--loop-removal', 'l1min'], {})

    def test_run_fluxcheck(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fluxcheck'], {})

    def test_run_fluxcheck_with_fastcore(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fluxcheck', '--fastcore'], {})

    def test_run_fluxcheck_with_tfba(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fluxcheck', '--tfba'],
            {'integer': True})

    def test_run_fluxcheck_with_reduce_lp(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fluxcheck', '--reduce-lp'], {})

    def test_run_fluxcheck_with_both_tfba_and_fastcore(self):
        with self.assertRaises(SystemExit):
            self.run_solver_command([
                '--model', self._model_dir, 'fluxcheck',
                '--tfba', '--fastcore'], {})

    def test_run_fluxcoupling(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fluxcoupling'], {})

    def test_run_formulacheck(self):
        with redirected_stdout():
            main(args=['--model', self._model_dir, 'formulacheck'])

    def test_run_fva(self):
        self.run_solver_command([
            '--model', self._model_dir, 'fva'], {})

    def test_run_gapfill(self):
        self.run_solver_command([
            '--model', self._model_dir, 'gapfill'], {'integer': True})

    def test_run_genedelete(self):
        self.run_solver_command([
            '--model', self._model_dir, 'genedelete', '--gene', 'gene_1'], {})

    def test_run_masscheck_compounds(self):
        self.run_solver_command([
            '--model', self._model_dir, 'masscheck', '--type', 'compound'], {})

    def test_run_masscheck_reactions(self):
        self.run_solver_command([
            '--model', self._model_dir, 'masscheck', '--type', 'reaction'], {})

    def test_run_randomsparse_reactions(self):
        self.run_solver_command([
            '--model', self._model_dir, 'randomsparse', '--type=reactions',
            '50%'], {})

    def test_run_randomsparse_genes(self):
        self.run_solver_command([
            '--model', self._model_dir, 'randomsparse', '--type=genes', '50%'],
            {})

    def test_run_randomsparse_exchange(self):
        self.run_solver_command([
            '--model', self._model_dir, 'randomsparse', '--type=exchange',
            '50%'], {})

    def test_run_robustness(self):
        self.run_solver_command([
            '--model', self._model_dir, 'robustness', 'rxn_2_\u03c0'], {})

    def test_run_sbmlexport(self):
        with redirected_stdout(BytesIO()):
            main(args=['--model', self._model_dir, 'sbmlexport'])

    def test_run_search_compound(self):
        with redirected_stdout():
            main(args=[
                '--model', self._model_dir, 'search', 'compound', '--id', 'A'])

    def test_run_search_reaction(self):
        with redirected_stdout():
            main(args=[
                '--model', self._model_dir, 'search', 'reaction',
                '--id', 'rxn_1'])

    def test_run_tableexport_reactions(self):
        with redirected_stdout() as f:
            main(args=[
                '--model', self._model_dir, 'tableexport', 'reactions'])

        self.assertEqual(f.getvalue(), '\n'.join([
            'id\tequation\tgenes\tin_model',
            "rxn_1\t|A_\u2206[e]| => |B[c]|\t['gene_1', 'gene_2']\tTrue",
            'rxn_2_\u03c0\t|B[c]| => |C[e]|\tgene_3 or'
            ' (gene_4 and gene_5)\tTrue',
            ''
        ]))

    def test_run_tableexport_compounds(self):
        with redirected_stdout() as f:
            main(args=[
                '--model', self._model_dir, 'tableexport', 'compounds'])

        self.assertEqual(f.getvalue(), '\n'.join([
            'id\tin_model', 'A_\u2206\tTrue', 'B\tTrue', 'C\tTrue', ''
        ]))

    def test_run_tableexport_medium(self):
        with redirected_stdout() as f:
            main(args=[
                '--model', self._model_dir, 'tableexport', 'medium'])

        self.assertEqual(f.getvalue(), '\n'.join([
            'Compound ID\tReaction ID\tLower Limit\tUpper Limit',
            'A_\u2206[e]\tNone\t-1000\t1000',
            'C[e]\tNone\t-1000\t1000',
            ''
        ]))

    def test_run_tableexport_limits(self):
        with redirected_stdout() as f:
            main(args=[
                '--model', self._model_dir, 'tableexport', 'limits'])

        self.assertEqual(f.getvalue(), '\n'.join([
            'Reaction ID\tLower Limits\tUpper Limits',
            'rxn_2_\u03c0\tNone\t100',
            ''
        ]))

    def test_run_tableexport_metadata(self):
        with redirected_stdout() as f:
            main(args=[
                '--model', self._model_dir, 'tableexport', 'metadata'])

        self.assertEqual(f.getvalue(), '\n'.join([
            'Model Name: Test model',
            'Biomass Reaction: rxn_1',
            'Default Flux Limits: 1000',
            ''
        ]))

    def test_command_main(self):
        main(MockCommand, args=[
            '--model', self._model_dir, '--test-argument'])
        self.assertTrue(MockCommand.run_called)
        self.assertTrue(MockCommand.init_parser_called)
        self.assertTrue(MockCommand.has_native_model)

    def test_solver_command_main(self):
        with self.assertRaises(generic.RequirementsError):
            main(MockSolverCommand, args=[
                '--model', self._model_dir,
                '--solver', 'name=not-an-actual-solver'])

    def test_metabolic_command_main(self):
        main(MockMetabolicCommand, args=[
            '--model', self._model_dir])
        self.assertTrue(MockMetabolicCommand.has_metabolic_model)

    def test_command_fail(self):
        mock_command = MockCommand(self._model, None)
        with self.assertRaises(SystemExit):
            mock_command.fail('Run time error')

    def test_command_argument_error(self):
        mock_command = MockCommand(self._model, None)
        with self.assertRaises(CommandError):
            mock_command.argument_error(msg='Argument error')
