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

import sys
import os
import shutil
import tempfile
import unittest

from six import StringIO

from psamm.command import main, Command, SolverCommandMixin
from psamm.lpsolver.generic import RequirementsError


class MockCommand(Command):
    """Test command.

    This is a command for testing basic commands.
    """
    init_parser_called = False
    run_called = False
    has_native_model = False
    has_metabolic_model = False

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
        with open(os.path.join(self._model_dir, 'model.yaml'), 'w') as f:
            f.write('''---
                name: Test model
                reactions:
                  - id: rxn_1
                    equation: '|A| => |B|'
                compounds:
                  - id: A
                  - id: B
                ''')

    def tearDown(self):
        shutil.rmtree(self._model_dir)

    def test_invoke_version(self):
        _stdout = sys.stdout
        try:
            sys.stdout = StringIO()
            with self.assertRaises(SystemExit) as context:
                main(args=['--version'])
            self.assertEqual(context.exception.code, 0)
        finally:
            sys.stdout = _stdout

    def test_command_main(self):
        main(MockCommand, args=[
            '--model', self._model_dir, '--test-argument'])
        self.assertTrue(MockCommand.run_called)
        self.assertTrue(MockCommand.init_parser_called)
        self.assertTrue(MockCommand.has_native_model)
        self.assertTrue(MockCommand.has_metabolic_model)

    def test_solver_command_main(self):
        with self.assertRaises(RequirementsError):
            main(MockSolverCommand, args=[
                '--model', self._model_dir,
                '--solver', 'name=not-an-actual-solver'])
