#!/usr/bin/env python

'''Search the database of reactions or compounds'''

import re

from metnet import command
from metnet import modelseed
from metnet.reaction import Compound
from metnet.formula import Formula

def parse_compound(s):
    '''Parse a compound specification with optional compartment'''
    m = re.match(r'([^\[]+)\[(\w+)\]', s)
    if m is not None:
        return Compound(m.group(1), compartment=m.group(2))
    return Compound(s)

def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())

class SearchCommand(command.Command):
    '''Defines the search command'''

    def __init__(self):
        super(SearchCommand, self).__init__('Search the database of reactions or compounds')

    def init_parser(self, parser):
        '''Initialize argument parser'''
        subparsers = parser.add_subparsers(title='Search domain')

        # Compound subcommand
        parser_compound = subparsers.add_parser('compound', help='Search in compounds')
        parser_compound.set_defaults(which='compound')
        parser_compound.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                        default=[], action='append', help='Compound ID')
        parser_compound.add_argument('--name', '-n', dest='name', metavar='name', type=str,
                                        default=[], action='append', help='Name of compound')

        # Reaction subcommand
        parser_reaction = subparsers.add_parser('reaction', help='Search in reactions')
        parser_reaction.set_defaults(which='reaction')
        parser_reaction.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                        default=[], action='append', help='Reaction ID')
        parser_reaction.add_argument('--compound', '-c', dest='compound', metavar='compound', type=str,
                                        default=[], action='append', help='Comma-separated list of compound IDs')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run search command'''

        # Load compound information
        compound_name = {}
        compound_synonyms = {}
        compound_formula = {}
        compound_for_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id
                compound_synonyms[compound.id] = compound.names

                for n in compound.names:
                    n = filter_search_term(n)
                    compound_for_name[n] = compound.id

                if compound.formula is not None and not '.' in compound.formula:
                    compound_formula[compound.id] = Formula.parse(compound.formula)

        which_command = kwargs['which']
        if which_command == 'compound':
            compound_ids = []
            for name in kwargs['name']:
                search_name = filter_search_term(name)
                if search_name in compound_for_name:
                    compound_ids.append(compound_for_name[search_name])
            for compound_id in kwargs['id']:
                if compound_id in compound_name:
                    compound_ids.append(compound_id)

            # Show results
            for compound_id in compound_ids:
                print 'ID: {}'.format(compound_id)
                if compound_id in compound_name:
                    print 'Name: {}'.format(compound_name[compound_id])
                if compound_id in compound_synonyms:
                    print 'Synonyms: {}'.format(', '.join(compound_synonyms[compound_id]))
                if compound_id in compound_formula:
                    print 'Formula: {}'.format(compound_formula[compound_id])
                print
        elif which_command == 'reaction':
            reaction_ids = []
            for reaction_id in kwargs['id']:
                if reaction_id in database.reactions:
                    reaction_ids.append(reaction_id)
            for compound_list in kwargs['compound']:
                matched_reactions = None
                for compound_spec in compound_list.split(','):
                    compound = parse_compound(compound_spec.strip())
                    if matched_reactions is None:
                        matched_reactions = set(database.get_compound_reactions(compound))
                    else:
                        matched_reactions &= database.compound_reactions.get(compound, set())
                if matched_reactions is not None:
                    reaction_ids.extend(matched_reactions)

            # Show results
            for reaction_id in reaction_ids:
                print 'ID: {}'.format(reaction_id)

                reaction = database.get_reaction(reaction_id)
                print 'Reaction (IDs): {}'.format(database.get_reaction(reaction_id))
                translated_reaction = reaction.translated_compounds(lambda x: compound_name.get(x, x))
                if reaction != translated_reaction:
                    print 'Reaction (names): {}'.format(translated_reaction)
                print

if __name__ == '__main__':
    command.main(SearchCommand())
