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

from __future__ import division

import os
import argparse
import logging
import abc
from itertools import islice
import multiprocessing as mp

import pkg_resources
from six import add_metaclass, iteritems, itervalues

from . import __version__ as package_version
from .datasource.native import NativeModel
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

    @classmethod
    def init_parser(cls, parser):
        """Initialize command line parser (:class:`argparse.ArgumentParser`)"""

    @abc.abstractmethod
    def run(self):
        """Execute command"""


class MetabolicMixin(object):
    """Mixin for commands that use a metabolic model representation."""

    def __init__(self, *args, **kwargs):
        super(MetabolicMixin, self).__init__(*args, **kwargs)

        self._mm = self._model.create_metabolic_model()


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


class ParallelTaskMixin(object):
    """Mixin for commands that run parallel computation tasks."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--parallel', help='Set number of parallel processes (0=auto)',
            type=int, default=0)
        super(ParallelTaskMixin, cls).init_parser(parser)

    def _create_executor(self, handler, args, cpus_per_worker=1):
        """Return a new :class:`.Executor` instance."""
        if self._args.parallel > 0:
            workers = self._args.parallel
        else:
            try:
                workers = mp.cpu_count() // cpus_per_worker
            except NotImplementedError:
                workers = 1

        if workers != 1:
            logger.info('Using {} parallel worker processes...'.format(
                workers))
            executor = ProcessPoolExecutor(
                processes=workers, handler_init=handler, handler_args=args)
        else:
            logger.info('Using single worker...')
            executor = SequentialExecutor(
                handler_init=handler, handler_args=args)

        return executor


class FilePrefixAppendAction(argparse.Action):
    """Action that appends one argument or multiple from a file.

    If the argument starts with a character in ``fromfile_prefix_chars``
    the remaining part of the argument is taken to be a file path. The file
    is read and every line is appended. Otherwise, the argument is simply
    appended.
    """
    def __init__(self, option_strings, dest, nargs=None,
                 fromfile_prefix_chars='@', **kwargs):
        if nargs is not None:
            raise ValueError('nargs not allowed')
        self.fromfile_prefix_chars = fromfile_prefix_chars
        super(FilePrefixAppendAction, self).__init__(
            option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        arguments = getattr(namespace, self.dest)
        if arguments is None:
            arguments = []
            setattr(namespace, self.dest, arguments)

        if len(values) > 0 and values[0] in self.fromfile_prefix_chars:
            filepath = values[1:]
            try:
                with open(filepath, 'r') as f:
                    for line in f:
                        arguments.append(line.strip())
            except IOError:
                parser.error('Unable to read arguments from file: {}'.format(
                    filepath))
        else:
            arguments.append(values)


class _ErrorMarker(object):
    """Signals error in the child process."""


class ExecutorError(Exception):
    """Error running tasks on executor."""


class _ExecutorProcess(mp.Process):
    def __init__(self, task_queue, result_queue, handler_init,
                 handler_args=()):
        super(_ExecutorProcess, self).__init__()
        self._task_queue = task_queue
        self._result_queue = result_queue
        self._handler_init = handler_init
        self._handler_args = handler_args

    def run(self):
        try:
            handler = self._handler_init(*self._handler_args)
            for tasks in iter(self._task_queue.get, None):
                results = [
                    (task, handler.handle_task(*task)) for task in tasks]
                self._result_queue.put(results)
        except BaseException:
            self._result_queue.put(_ErrorMarker())
            raise


class Executor(object):
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def close(self):
        pass


class ProcessPoolExecutor(Executor):
    def __init__(self, handler_init, handler_args=(), processes=None):
        if processes is None:
            try:
                processes = mp.cpu_count()
            except NotImplementedError:
                processes = 1

        self._process_count = processes
        self._processes = []
        self._task_queue = mp.Queue()
        self._result_queue = mp.Queue()

        for _ in range(self._process_count):
            p = _ExecutorProcess(
                self._task_queue, self._result_queue, handler_init,
                handler_args)
            p.start()
            self._processes.append(p)

    def apply(self, task):
        self._task_queue.put([task])
        results = self._result_queue.get()
        if isinstance(results, _ErrorMarker):
            raise ExecutorError('Child process failed')

        return next(itervalues(results))

    def imap_unordered(self, iterable, chunksize=1):
        def iter_chunks():
            while True:
                chunk = list(islice(iterable, chunksize))
                if len(chunk) == 0:
                    break
                yield chunk

        it = iter_chunks()
        workers = 0
        for i in range(self._process_count):
            tasks = next(it, None)
            if tasks is None:
                break

            self._task_queue.put(tasks)
            workers += 1

        while workers > 0:
            results = self._result_queue.get()
            if isinstance(results, _ErrorMarker):
                raise ExecutorError('Child process failed')

            tasks = next(it, None)
            if tasks is None:
                workers -= 1

            self._task_queue.put(tasks)

            for task, result in results:
                yield task, result

    def close(self):
        for i in range(self._process_count):
            self._task_queue.put(None)

    def join(self):
        for p in self._processes:
            p.join()


class SequentialExecutor(Executor):
    def __init__(self, handler_init, handler_args=()):
        self._handler = handler_init(*handler_args)

    def apply(self, task):
        return self._handler.handle_task(*task)

    def imap_unordered(self, iterable, chunksize=1):
        for task in iterable:
            yield task, self._handler.handle_task(*task)

    def join(self):
        pass


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
        logging.basicConfig(level=logging.INFO)
        base_logger = logging.getLogger('psamm')
        if len(base_logger.handlers) == 0:
            handler = logging.StreamHandler()
            handler.setFormatter(
                logging.Formatter(u'%(levelname)s: %(message)s'))
            base_logger.addHandler(handler)
            base_logger.propagate = False

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
