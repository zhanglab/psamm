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

from __future__ import unicode_literals

from ..command import Command, MetabolicMixin


class ConsoleCommand(MetabolicMixin, Command):
    """Start an interactive Python console with the model loaded."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            '--type', choices=('python', 'ipython', 'ipython-kernel'),
            default='python', help='type of console to open')

    def open_python(self, message, namespace):
        """Open interactive python console"""

        # Importing readline will in some cases print weird escape
        # characters to stdout. To avoid this we only import readline
        # and related packages at this point when we are certain
        # they are needed.
        from code import InteractiveConsole
        import readline
        import rlcompleter

        readline.set_completer(rlcompleter.Completer(namespace).complete)
        readline.parse_and_bind('tab: complete')
        console = InteractiveConsole(namespace)
        console.interact(message)

    def open_ipython(self, message, namespace):
        from IPython.terminal.embed import InteractiveShellEmbed
        console = InteractiveShellEmbed(user_ns=namespace, banner2=message)
        console()

    def open_ipython_kernel(self, message, namespace):
        from IPython import embed_kernel
        embed_kernel(local_ns=namespace)

    def run(self):
        message = ('Native model has been loaded into: "model"\n' +
                   'Metabolic model has been loaded into: "mm"')
        namespace = {'model': self._model, 'mm': self._mm}
        console_type = self._args.type

        if console_type == 'python':
            self.open_python(message, namespace)
        elif console_type == 'ipython':
            self.open_ipython(message, namespace)
        elif console_type == 'ipython-kernel':
            self.open_ipython_kernel(message, namespace)
