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

import re
import json
import logging

from six import text_type, string_types, integer_types

from ..command import Command
from .. import util

logger = logging.getLogger(__name__)


_JSON_TYPES = string_types + integer_types + (
    dict, list, tuple, float, bool, type(None))
_QUOTE_REGEXP = re.compile(r'[\t\n]')


def _encode_value(value):
    if value is None:
        return ''
    if not isinstance(value, _JSON_TYPES):
        value = text_type(value)
    if (isinstance(value, string_types) and not _QUOTE_REGEXP.search(value) and
            value != ''):
        return value
    return json.dumps(value)


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
        for reaction in self._model.parse_reactions():
            property_set.update(reaction.properties)

        property_list_sorted = sorted(
            property_set, key=lambda x: (x != 'id', x != 'equation', x))

        print('\t'.join(
            [text_type(x) for x in property_list_sorted] + ['in_model']))

        model_reactions = set(self._model.parse_model())
        for reaction in self._model.parse_reactions():
            line_content = [reaction.properties.get(property)
                            for property in property_list_sorted]
            in_model = (not self._model.has_model_definition() or
                        reaction.id in model_reactions)
            line_content.append(in_model)
            print('\t'.join(_encode_value(value) for value in line_content))

    def _compound_export(self):
        compound_set = set()
        for compound in self._model.parse_compounds():
            compound_set.update(compound.properties)

        compound_list_sorted = sorted(
            compound_set, key=lambda x: (x != 'id', x != 'name', x))

        print('\t'.join(
            [text_type(x) for x in compound_list_sorted] + ['in_model']))

        metabolic_model = self._model.create_metabolic_model()
        model_compounds = set(x.name for x in metabolic_model.compounds)
        for compound in self._model.parse_compounds():
            line_content = [compound.properties.get(property)
                            for property in compound_list_sorted]
            in_model = (not self._model.has_model_definition() or
                        compound.id in model_compounds)
            line_content.append(in_model)
            print('\t'.join(_encode_value(value) for value in line_content))

    def _media_export(self):
        print('{}\t{}\t{}\t{}'.format('Compound ID', 'Reaction ID',
                                      'Lower Limit', 'Upper Limit'))
        default_flux = self._model.default_flux_limit

        for compound, reaction, lower, upper in self._model.parse_medium():
            if lower is None:
                lower = -1 * default_flux

            if upper is None:
                upper = default_flux

            print('\t'.join(_encode_value(value) for value in (
                compound, reaction, lower, upper)))

    def _limits_export(self):
        print('{}\t{}\t{}'.format(
            'Reaction ID', 'Lower Limits', 'Upper Limits'))
        for reaction, lower, upper in self._model.parse_limits():
            print('\t'.join(_encode_value(value) for value in (
                reaction, lower, upper)))

    def _metadata_export(self):
        git_version = None
        if self._model.context is not None:
            git_version = util.git_try_describe(self._model.context.basepath)

        print('Model Name\t{}'.format(
            _encode_value(self._model.name)))
        print('Biomass Reaction\t{}'.format(
            _encode_value(self._model.biomass_reaction)))
        print('Default Flux Limits\t{}'.format(
            _encode_value(self._model.default_flux_limit)))

        if git_version is not None:
            print('Git version\t{}'.format(_encode_value(git_version)))
