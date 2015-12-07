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
# Copyright 2015  Keith Dufault-Thompson <keitht547@my.uri.edu>

from __future__ import unicode_literals

import logging

from ..command import Command

logger = logging.getLogger(__name__)


class ExportTableCommand(Command):
    """Export parts of the model as tab separated tables.

    This command can export various parts of the model as a TSV table.
    The output is written to standard output. The type argument is used
    to determine which part of the model to export:

    - reactions: Export reactions and reaction metadata
    - compounds: Export compounds and compound metadata
    - medium: Export the list of exchange compounds/reactions
    - limits: Export list of internal flux limits
    - metadata: Export general model metadata
    """

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            'export', metavar='export_type',
            choices=['reactions', 'compounds', 'medium', 'limits', 'metadata'],
            help='Type of model data to export')

    def run(self):
        export_type = self._args.export
        if export_type == 'reactions':
            self._reaction_export()
        elif export_type == 'compounds':
            self._compound_export()
        elif export_type == 'medium':
            self._media_export()
        elif export_type == 'limits':
            self._limits_export()
        elif export_type == 'metadata':
            self._metadata_export()

    def _reaction_export(self):
        property_set = set()
        for j in self._model.parse_reactions():
            property_set.update(j.properties)

        property_list_sorted = sorted(
            property_set, key=lambda x: (x != 'id', x != 'equation', x))

        print('\t'.join([str(x) for x in property_list_sorted]))
        for i in self._model.parse_reactions():
            print('\t'.join(str(i.properties.get(property))
                            for property in property_list_sorted))

    def _compound_export(self):
        compound_set = set()
        for compound in self._model.parse_compounds():
            compound_set.update(compound.properties)

        compound_list_sorted = sorted(
            compound_set, key=lambda x: (x != 'id', x != 'name', x))

        print('\t'.join([str(x) for x in compound_list_sorted]))
        for compound in self._model.parse_compounds():
            print('\t'.join(str(compound.properties.get(property))
                            for property in compound_list_sorted))

    def _media_export(self):
        print('{}\t{}\t{}\t{}'.format('Compound ID', 'Reaction ID',
                                      'Lower Limit', 'Upper Limit'))
        default_flux = self._model.get_default_flux_limit()
        if default_flux is None:
            default_flux = 1000

        for compound, reaction, lower, upper in self._model.parse_medium():
            if lower is None:
                lower = -1 * default_flux

            if upper is None:
                upper = default_flux

            print('{}\t{}\t{}\t{}'.format(compound, reaction, lower, upper))

    def _limits_export(self):
        print('{}\t{}\t{}'.format(
            'Reaction ID', 'Lower Limits', 'Upper Limits'))
        for reaction, lower, upper in self._model.parse_limits():
            print('{}\t{}\t{}'.format(reaction, lower, upper))

    def _metadata_export(self):
        print('Model Name: {}'.format(self._model.get_name()))
        print('Biomass Reaction: {}'.format(
            self._model.get_biomass_reaction()))
        print('Default Flux Limits: {}'.format(
            self._model.get_default_flux_limit()))
