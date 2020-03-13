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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>

import sys
import argparse

from ..datasource import sbml
from ..command import MetabolicMixin, Command


class SBMLExport(MetabolicMixin, Command):
    """Export model as SBML level 3 file."""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            'file', type=argparse.FileType('w'), nargs='?', default=sys.stdout,
            help='File path for writing the SBML file'
                 ' (use "-" for standard output)')
        parser.add_argument(
            '--pretty', action='store_true',
            help='Format output for readability')
        super(SBMLExport, cls).init_parser(parser)

    def run(self):
        writer = sbml.SBMLWriter()
        writer.write_model(self._args.file, self._model, self._args.pretty)
