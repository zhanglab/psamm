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

"""Command line interface.

Each command in the command line interface is implemented as a subclass of
:class:`Command`. Commands are also referenced from ``setup.py`` using the
entry point mechanism which allows the commands to be automatically
discovered.

The :func:`.main` function is the entry point of command line interface.
"""

import os
import argparse
import logging
import abc

import pkg_resources
from six import add_metaclass, iteritems

from . import __version__ as package_version
from .datasource.native import NativeModel
from .metabolicmodel import MetabolicModel
from .database import DictDatabase
from .lpsolver import generic
from . import util

logger = logging.getLogger(__name__)


class CommandError(Exception):
    """Error from running a command.

    This should be raised from a ``Command.run()`` if any arguments are
    misspecified. When the command is run and the ``CommandError`` is raised,
    the caller will exit with an error code and print appropriate usage
    information.
    """


@add_metaclass(abc.ABCMeta)
class Command(object):
    """Represents a command in the interface, operating on a model.

    The constructor will be given the NativeModel and the command line
    namespace. The subclass must implement :meth:`run` to handle command
    execution. The doc string will be used as documentation for the command
    in the command line interface.

    In addition, :meth:`init_parser` can be implemented as a classmethod which
    will allow the command to initialize an instance of
    :class:`argparse.ArgumentParser` as desired. The resulting argument
    namespace will be passed to the constructor.
    """

    def __init__(self, model, args):
        self._model = model
        self._args = args

        name = self._model.get_name()
        if name is None:
            name = str(self._model.context)
        logger.info('Model: {}'.format(name))

        version = util.git_try_describe(self._model.context.basepath)
        if version is not None:
            logger.info('Model Git version: {}'.format(version))

        # Create metabolic model
        database = DictDatabase()
        for reaction in model.parse_reactions():
            if reaction.equation is not None:
                database.set_reaction(reaction.id, reaction.equation)

        # Warn about undefined compounds
        compounds = set()
        for compound in model.parse_compounds():
            compounds.add(compound.id)

        undefined_compounds = set()
        for reaction in database.reactions:
            for compound, _ in database.get_reaction_values(reaction):
                if compound.name not in compounds:
                    undefined_compounds.add(compound.name)

        for compound in sorted(undefined_compounds):
            logger.warning(
                'The compound {} was not defined in the list'
                ' of compounds'.format(compound))

        self._mm = MetabolicModel.load_model(
            database, model.parse_model(), model.parse_medium(),
            model.parse_limits(), v_max=model.get_default_flux_limit())

    @classmethod
    def init_parser(cls, parser):
        """Initialize command line parser (:class:`argparse.ArgumentParser`)"""

    @abc.abstractmethod
    def run(self):
        """Execute command"""


class SolverCommandMixin(object):
    """Mixin for commands that use an LP solver.

    This adds a ``--solver`` parameter to the command that the user can use to
    select a specific solver. It also adds the method :meth:`_get_solver` which
    will return a solver with the specified default requirements. The user
    requirements will override the default requirements.
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--solver', action='append', type=str,
            help='Specify solver requirements (e.g. "rational=yes")')
        super(SolverCommandMixin, cls).init_parser(parser)

    def __init__(self, *args, **kwargs):
        super(SolverCommandMixin, self).__init__(*args, **kwargs)
        self._solver_args = {}
        if self._args.solver is not None:
            for s in self._args.solver:
                key, value = generic.parse_solver_setting(s)
                self._solver_args[key] = value

    def _get_solver(self, **kwargs):
        """Return a new :class:`psamm.lpsolver.lp.Solver` instance"""
        solver_args = dict(kwargs)
        solver_args.update(self._solver_args)
        return generic.Solver(**solver_args)


def main(command_class=None, args=None):
    """Run the command line interface with the given :class:`Command`.

    If no command class is specified the user will be able to select a specific
    command through the first command line argument. If the ``args`` are
    provided, these should be a list of strings that will be used instead of
    ``sys.argv[1]``. This is mostly useful for testing.
    """

    # Set up logging for the command line interface
    if 'PSAMM_DEBUG' in os.environ:
        level = getattr(logging, os.environ['PSAMM_DEBUG'].upper(), None)
        if level is not None:
            logging.basicConfig(level=level)
    else:
        logging.basicConfig(
            level=logging.INFO, format='%(levelname)s: %(message)s')

    title = 'Metabolic modeling tools'
    if command_class is not None:
        title, _, _ = command_class.__doc__.partition('\n\n')

    parser = argparse.ArgumentParser(description=title)
    parser.add_argument('--model', metavar='file', default='.',
                        help='Model definition')
    parser.add_argument(
        '-V', '--version', action='version',
        version='%(prog)s ' + package_version)

    if command_class is not None:
        # Command explicitly given, only allow that command
        command_class.init_parser(parser)
        parser.set_defaults(command=command_class)
    else:
        # Discover all available commands
        commands = {}
        for entry in pkg_resources.iter_entry_points('psamm.commands'):
            canonical = entry.name.lower()
            if canonical not in commands:
                command_class = entry.load()
                commands[canonical] = command_class
            else:
                logger.warning('Command {} was found more than once!'.format(
                    canonical))

        # Create parsers for subcommands
        subparsers = parser.add_subparsers(title='Commands', metavar='command')
        for name, command_class in sorted(iteritems(commands)):
            title, _, _ = command_class.__doc__.partition('\n\n')
            subparser = subparsers.add_parser(
                name, help=title.rstrip('.'),
                description=command_class.__doc__)
            subparser.set_defaults(command=command_class)
            command_class.init_parser(subparser)

    parsed_args = parser.parse_args(args)

    # Load model definition
    model = NativeModel(parsed_args.model)

    # Instantiate command with model and run
    command = parsed_args.command(model, parsed_args)
    try:
        command.run()
    except CommandError as e:
        parser.error(str(e))
