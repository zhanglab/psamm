# coding=utf-8
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
# Copyright 2015-2020  Keith Dufault-Thompson <keitht547@my.uri.edu>
# Copyright 2020-2020  Elysha Sameth <esameth1@my.uri.edu>

from __future__ import print_function, unicode_literals
import re

from six import text_type

from ..command import Command, FilePrefixAppendAction, convert_to_unicode
from ..datasource.reaction import parse_compound


def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())


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
            '--id', '-i', dest='id', metavar='id',
            action=FilePrefixAppendAction, type=convert_to_unicode, default=[],
            help='Compound ID')
        parser_compound.add_argument(
            '--name', '-n', dest='name', metavar='name',
            action=FilePrefixAppendAction, type=convert_to_unicode, default=[],
            help='Name of compound')
        parser_compound.add_argument(
            '--key', dest='key', metavar='key',
            type=convert_to_unicode, default=None,
            help='String to search for within compound '
                 'properties. (case insensitive). Only one --key allowed')
        parser_compound.add_argument(
            '--exact', action='store_true',
            help='Match full compound property')
        parser_compound.add_argument(
            '--inmodel', action='store_true',
            help='Only search compounds in the model reactions')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser(
            'reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument(
            '--id', '-i', dest='id', metavar='id',
            action=FilePrefixAppendAction, type=convert_to_unicode, default=[],
            help='Reaction ID')
        parser_reaction.add_argument(
            '--compound', '-c', dest='compound', metavar='compound',
            action=FilePrefixAppendAction, type=convert_to_unicode, default=[],
            help='Comma-separated list of compound IDs')
        parser_reaction.add_argument(
            '--key', dest='key', metavar='key',
            type=convert_to_unicode, default=None,
            help='String to search for within reaction properties. '
                 '(case insensitive). Only one --key allowed')
        parser_reaction.add_argument(
            '--exact', action='store_true',
            help='Match full reaction property')
        parser_reaction.add_argument(
            '--inmodel', action='store_true',
            help='Only search reactions in the model')

    def run(self):
        """Run search command."""
        self._mm = self._model.create_metabolic_model()
        which_command = self._args.which
        if which_command == 'compound':
            self._search_compound()
        elif which_command == 'reaction':
            self._search_reaction()

    def _search_compound(self):
        selected_compounds = set()
        for compound in self._model.compounds:
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

            # find compounds that contains any of given properties
            if self._args.key is not None:
                # prepare s list of all compound properties
                compound_prop_list = []
                for cpd_property in compound.properties.values():
                    if isinstance(cpd_property, list):
                        for i in cpd_property:
                            compound_prop_list.append(convert_to_unicode(
                                text_type(i)).lower())
                    else:
                        compound_prop_list.append(convert_to_unicode(
                            text_type(cpd_property)).lower())

                # find compound entry based on given property argument
                if self._args.exact:
                    if self._args.key.lower() in compound_prop_list:
                        selected_compounds.add(compound)
                else:
                    if self._args.key.lower() in \
                            ('|'.join(compound_prop_list)):
                        selected_compounds.add(compound)

        # Show results
        if self._args.inmodel:
            model_compounds = set(x.name for x in self._mm.compounds)
            final_compounds = [i for i in selected_compounds
                               if i.id in model_compounds]
        else:
            final_compounds = selected_compounds

        for compound in final_compounds:
            props = set(compound.properties) - {'id'}
            print('id: {}'.format(compound.id))
            for prop in sorted(props):
                print('{}: {}'.format(prop, compound.properties[prop]))
            if compound.filemark is not None:
                print('Defined in {}'.format(compound.filemark))
            print()

    def _search_reaction(self):
        selected_reactions = set()

        # Prepare translation table from compound id to name
        compound_name = {}
        for compound in self._model.compounds:
            if 'name' in compound.properties:
                compound_name[compound.id] = compound.properties['name']

        # Prepare sets of matched compounds
        search_compounds = []
        for compound_list in self._args.compound:
            compound_set = set()
            for compound_spec in compound_list.split(','):
                compound_set.add(parse_compound(compound_spec.strip()))
            search_compounds.append(compound_set)

        for reaction in self._model.reactions:
            if len(self._args.id) > 0:
                if any(r == reaction.id for r in self._args.id):
                    selected_reactions.add(reaction)
                    continue

            if len(search_compounds) > 0:
                compounds = set(c for c, _ in reaction.equation.compounds)
                compounds.update(c.in_compartment(None) for c, _ in
                                 reaction.equation.compounds)
                if any(c.issubset(compounds) for c in search_compounds):
                    selected_reactions.add(reaction)
                    continue

            if self._args.key is not None:
                # prepare a list of all reaction properties
                raw_reaction_prop_list = [
                    reaction.properties[key] for key in reaction.properties]
                reaction_prop_list = []
                for rxn_property in raw_reaction_prop_list:
                    if isinstance(rxn_property, list):
                        for i in rxn_property:
                            reaction_prop_list.append(convert_to_unicode(
                                text_type(i)).lower())

                    else:
                        reaction_prop_list.append(convert_to_unicode(
                            text_type(rxn_property)).lower())

                # find reaction based on given property argument
                if self._args.exact:
                    if self._args.key.lower() in reaction_prop_list:
                        selected_reactions.add(reaction)
                        continue
                else:
                    if self._args.key.lower() in '|'.join(reaction_prop_list):
                        selected_reactions.add(reaction)
        if self._args.inmodel:
            final_reactions = [i for i in selected_reactions
                               if i.id in self._mm.reactions]
        else:
            final_reactions = selected_reactions

        # Show results
        for reaction in final_reactions:
            props = set(reaction.properties) - {'id', 'equation'}
            print('id: {}'.format(reaction.id))
            print('equation: {}'.format(
                reaction.equation))
            translated_equation = reaction.equation.translated_compounds(
                lambda x: compound_name.get(x, x))
            if reaction.equation != translated_equation:
                print('equation (compound names): {}'.format(
                    translated_equation))
            for prop in sorted(props):
                print('{}: {}'.format(prop, reaction.properties[prop]))
            if reaction.filemark is not None:
                print('Defined in {}'.format(reaction.filemark))
            print()
