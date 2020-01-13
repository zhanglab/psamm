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

from __future__ import print_function, unicode_literals

import re

from six import text_type

from ..command import Command, FilePrefixAppendAction
from ..datasource.reaction import parse_compound

from psamm.expression import boolean


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
            action=FilePrefixAppendAction, type=text_type, default=[],
            help='Compound ID')
        parser_compound.add_argument(
            '--name', '-n', dest='name', metavar='name',
            action=FilePrefixAppendAction, type=text_type, default=[],
            help='Name of compound')
        parser_compound.add_argument(
            '--props', '-c', dest='props', metavar='props',
            nargs='+', type=str, default=None,
            help='Space-separated list of compound properties, such as '
                 'compound formula, compound charge, molecular weight...')
        parser_compound.add_argument(
            '--match-type', '-m', dest='match_type', metavar='match_type',
            type=str, choices=['exact', 'vague'], default='vague',
            help='chose the map type when using --props to find compound. '
                 'exact means completely match, vague means partially match')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser(
            'reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument(
            '--id', '-i', dest='id', metavar='id',
            action=FilePrefixAppendAction, type=str, default=[],
            help='Reaction ID')
        parser_reaction.add_argument(
            '--compound', '-c', dest='compound', metavar='compound',
            action=FilePrefixAppendAction, type=str, default=[],
            help='Comma-separated list of compound IDs')
        parser_reaction.add_argument(
            '--props', '-p', dest='props', metavar='props',
            nargs='+', type=str, default=None,
            help='Space-separated list of reaction properties, such as '
                 'reaction name, ec number, pathway or subsystem name or '
                 'gene association')
        parser_reaction.add_argument(
            '--match-type', '-m', dest='match_type', metavar='match_type',
            type=str, choices=['exact', 'vague'], default='exact',
            help='chose the map type when using --props to find reaction. '
                 'exact means completely match, vague means partially match')

    def run(self):
        """Run search command."""

        which_command = self._args.which
        if which_command == 'compound':
            self._search_compound()
        elif which_command == 'reaction':
            self._search_reaction()

        for reaction in self._model.reactions:
            print(reaction)
            if reaction.genes is not None:
                e = boolean.Expression(reaction.genes)
                # e_1 = e.substitute(lambda v: target_genes_l.get(v.symbol, v))
                print(e)
                print(reaction.genes)

    def _search_compound(self):
        selected_compounds = set()
        print(self._args.props)
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
            if self._args.props is not None:

                # prepare s list of all compound properties
                compound_prop_list = []
                for cpd_property in compound.properties.values():
                    if isinstance(cpd_property, list):
                        for i in cpd_property:
                            compound_prop_list.append(str(i).lower())
                    else:
                        compound_prop_list.append(str(cpd_property).lower())

                # find compound entry based on given property argument
                props = set()
                if self._args.match_type == 'exact':
                    for property in self._args.props:
                        props.add(property.lower())
                    if any(prop in props for prop in compound_prop_list):
                        selected_compounds.add(compound)
                        continue
                elif self._args.match_type == 'vague':
                    for property_list in self._args.props:
                        for property in property_list.split(','):
                            props.add(property.lower())
                    if any(prop in ('|'.join(compound_prop_list)) for
                           prop in props):
                        selected_compounds.add(compound)

        # Show results
        for compound in selected_compounds:
            props = set(compound.properties) - {'id'}
            print('id: {}'.format(compound.id))
            for prop in sorted(props):
                print('{}: {}'.format(prop, compound.properties[prop]))
            if compound.filemark is not None:
                print('Defined in {}'.format(compound.filemark))
            print()

    def _search_reaction(self):
        selected_reactions = set()
        print(self._args.props)
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

            if self._args.props is not None:
                props = set()

                # prepare s list of all reaction properties
                # except reaction equation
                raw_reaction_prop_list = [
                    reaction.properties[key] for key in reaction.properties
                    if key != 'equation']
                reaction_prop_list = []
                for rxn_property in raw_reaction_prop_list:
                    if isinstance(rxn_property, list):
                        for i in rxn_property:
                            reaction_prop_list.append(str(i).lower())
                    else:
                        reaction_prop_list.append(str(rxn_property).lower())

                # find reaction based on given property argument
                if self._args.match_type == 'exact':
                    for property in self._args.props:
                        props.add(property.lower())
                    if any(prop in props for prop in reaction_prop_list):
                        selected_reactions.add(reaction)
                        continue
                elif self._args.match_type == 'vague':
                    for property in self._args.props:
                        props.add(property.lower())
                    print(props)
                    if any(prop in '|'.join(reaction_prop_list) for
                           prop in props):
                        selected_reactions.add(reaction)

        # Show results
        for reaction in selected_reactions:
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
