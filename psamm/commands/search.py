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

from __future__ import print_function

import re

from ..command import Command
from ..reaction import Compound


class SearchCommand(Command):
    """Search for reactions and compounds in the model."""

    @classmethod
    def init_parser(cls, parser):
        """Initialize argument parser"""
        subparsers = parser.add_subparsers(title='Search domain')

        # Compound subcommand
        parser_compound = subparsers.add_parser(
            'compound', help='Search in compounds')
        parser_compound.set_defaults(which='compound')
        parser_compound.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            default=[], action='append', help='Compound ID')
        parser_compound.add_argument(
            '--name', '-n', dest='name', metavar='name', type=str,
            default=[], action='append', help='Name of compound')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser(
            'reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument(
            '--id', '-i', dest='id', metavar='id', type=str,
            default=[], action='append', help='Reaction ID')
        parser_reaction.add_argument(
            '--compound', '-c', dest='compound', metavar='compound', type=str,
            default=[], action='append',
            help='Comma-separated list of compound IDs')

    def run(self):
        """Run search command."""

        def parse_compound(s):
            """Parse a compound specification with optional compartment"""
            m = re.match(r'([^\[]+)\[(\w+)\]', s)
            if m is not None:
                return Compound(m.group(1), compartment=m.group(2))
            return Compound(s)

        def filter_search_term(s):
            return re.sub(r'[^a-z0-9]+', '', s.lower())

        which_command = self._args.which
        if which_command == 'compound':
            selected_compounds = set()

            for compound in self._model.parse_compounds():
                if len(self._args.id) > 0:
                    if any(c == compound.id for c in self._args.id):
                        selected_compounds.add(compound)
                        continue

                if len(self._args.name) > 0:
                    names = set()
                    if 'name' in compound.properties:
                        names.add(compound.properties['name'])
                    names.update(compound.properties.get('names', []))
                    names = set(filter_search_term(n) for n in names)
                    if any(filter_search_term(n) in names
                            for n in self._args.name):
                        selected_compounds.add(compound)
                        continue

            # Show results
            for compound in selected_compounds:
                print('ID: {}'.format(compound.id))
                if 'name' in compound.properties:
                    print('Name: {}'.format(compound.properties['name']))
                if 'names' in compound.properties:
                    print('Additional names: {}'.format(
                        ', '.join(compound.properties['names'])))
                if 'formula' in compound.properties:
                    print('Formula: {}'.format(compound.properties['formula']))
                if compound.filemark is not None:
                    print('Parsed from: {}'.format(compound.filemark))
                print()
        elif which_command == 'reaction':
            selected_reactions = set()

            # Prepare translation table from compound id to name
            compound_name = {}
            for compound in self._model.parse_compounds():
                if 'name' in compound.properties:
                    compound_name[compound.id] = compound.properties['name']

            # Prepare sets of matched compounds
            search_compounds = []
            for compound_list in self._args.compound:
                compound_set = set()
                for compound_spec in compound_list.split(','):
                    compound_set.add(parse_compound(compound_spec.strip()))
                search_compounds.append(compound_set)

            for reaction in self._model.parse_reactions():
                if len(self._args.id) > 0:
                    if any(r == reaction.id for r in self._args.id):
                        selected_reactions.add(reaction)
                        continue

                if len(search_compounds) > 0:
                    compounds = set(c for c, _ in reaction.equation.compounds)
                    if any(c.issubset(compounds) for c in search_compounds):
                        selected_reactions.add(reaction)
                        continue

            # Show results
            for reaction in selected_reactions:
                print('ID: {}'.format(reaction.id))
                print('Reaction (IDs): {}'.format(
                    reaction.equation))
                translated_equation = reaction.equation.translated_compounds(
                    lambda x: compound_name.get(x, x))
                if reaction.equation != translated_equation:
                    print('Reaction (names): {}'.format(translated_equation))
                if reaction.filemark is not None:
                    print('Parsed from: {}'.format(reaction.filemark))
                print()
