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

from ..command import Command, CommandError

try:
    import xlsxwriter
except ImportError:
    xlsxwriter = None

logger = logging.getLogger(__name__)


class ExcelExportCommand(Command):
    """Export the metabolic model as an Excel workbook"""

    @classmethod
    def init_parser(cls, parser):
        parser.add_argument(
            'file', type=str, help='File path for writing the Excel workbook')

    def run(self):
        model = self._model
        if xlsxwriter is None:
            raise CommandError(
                'Excel export requires the XlsxWriter python module')
        workbook = xlsxwriter.Workbook(self._args.file)
        reaction_sheet = workbook.add_worksheet(name='Reactions')

        property_set = set()
        for j in model.parse_reactions():
            property_set.update(j.properties)
        property_list = list(property_set)
        property_list_sorted = sorted(property_list,
                                      key=lambda x: (x != 'id',
                                                     x != 'equation', x))
        for z, i in enumerate(property_list_sorted):
            reaction_sheet.write_string(0, z, str(i))
        for x, i in enumerate(model.parse_reactions()):
            for y, j in enumerate(property_list_sorted):
                reaction_sheet.write_string(x+1, y, str(i.properties.get(j)))

        compound_sheet = workbook.add_worksheet(name='Compounds')

        compound_set = set()
        for compound in model.parse_compounds():
            compound_set.update(compound.properties)

        compound_list_sorted = sorted(compound_set,
                                      key=lambda x: (x != 'id',
                                                     x != 'name', x))

        for z, i in enumerate(compound_list_sorted):
            compound_sheet.write_string(0, z, str(i))
        for x, i in enumerate(model.parse_compounds()):
            for y, j in enumerate(compound_list_sorted):
                compound_sheet.write_string(x+1, y, str(i.properties.get(j)))

        media_sheet = workbook.add_worksheet(name='Medium')

        media_sheet.write_string(0, 0, 'Compound ID')
        media_sheet.write_string(0, 1, 'Reaction ID')
        media_sheet.write_string(0, 2, 'Lower Limit')
        media_sheet.write_string(0, 3, 'Upper Limit')

        default_flux = model.get_default_flux_limit()
        if default_flux is None:
            default_flux = 1000

        for x, (compound, reaction, lower, upper) in enumerate(
                model.parse_medium()):
            if lower is None:
                lower = -1 * default_flux

            if upper is None:
                upper = default_flux

            media_sheet.write(x+1, 0, str(compound))
            media_sheet.write(x+1, 1, str(reaction))
            media_sheet.write(x+1, 2, str(lower))
            media_sheet.write(x+1, 3, str(upper))

        limits_sheet = workbook.add_worksheet(name='Limits')

        limits_sheet.write_string(0, 0, 'Reaction ID')
        limits_sheet.write_string(0, 1, 'Lower Limit')
        limits_sheet.write_string(0, 2, str('Upper Limit'))

        for x, i in enumerate(model.parse_limits()):
            limits_sheet.write(x, 0, (i[0]))
            limits_sheet.write(x, 1, (i[1]))
            limits_sheet.write(x, 2, (i[2]))

        model_info = workbook.add_worksheet(name='Model Info')

        model_info.write(0, 0, ('Model Name: {}'.
                                format(model.get_name())))
        model_info.write(1, 0, ('Biomass Reaction: {}'.
                                format(model.get_biomass_reaction())))
        model_info.write(2, 0, ('Default Flux Limits: {}'.
                                format(model.get_default_flux_limit())))

        workbook.close()
