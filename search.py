#!/usr/bin/env python

'''Search the database of reactions or compounds'''

import argparse
import csv
import re

from metnet.reaction import Compound
from metnet.metabolicmodel import MetabolicDatabase
from metnet.formula import Formula

def parse_compound(s):
    '''Parse a compound specification with optional compartment'''
    m = re.match(r'([^\[]+)\[(\w+)\]', s)
    if m is not None:
        return Compound(m.group(1), compartment=m.group(2))
    return Compound(s)

def filter_search_term(s):
    return re.sub(r'[^a-z0-9]+', '', s.lower())

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Search the database of reactions or compounds')
    parser.add_argument('--database', metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to use as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Compound information table')
    parser.add_argument('--model', metavar='reactionlist',
                        type=argparse.FileType('r'), help='Model definition')
    subparsers = parser.add_subparsers(title='Search domain', description='Where to search')

    # Compound subcommand
    parser_compound = subparsers.add_parser('compound')
    parser_compound.set_defaults(which='compound')
    parser_compound.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                    default=[], action='append', help='Compound ID')
    parser_compound.add_argument('--name', '-n', dest='name', metavar='name', type=str,
                                    default=[], action='append', help='Name of compound')

    # Reaction subcommand
    parser_reaction = subparsers.add_parser('reaction')
    parser_reaction.set_defaults(which='reaction')
    parser_reaction.add_argument('--id', '-i', dest='id', metavar='id', type=str,
                                    default=[], action='append', help='Reaction ID')
    parser_reaction.add_argument('--compound', '-c', dest='compound', metavar='compound', type=str,
                                    default=[], action='append', help='Comma-separated list of compound IDs')
    args = parser.parse_args()

    # Load reaction database
    database = MetabolicDatabase.load_from_files(*args.database)

    # Load compound information
    compounds = {}
    compound_synonyms = {}
    compound_formula = {}
    compound_for_name = {}
    for compound_table in args.compounds:
        compound_table.readline() # Skip header
        for row in csv.reader(compound_table, delimiter='\t'):
            compound_id, names, formula = row[:3]

            synonyms = names.split(',<br>')
            name = synonyms.pop()

            compounds[compound_id] = name
            compound_synonyms[compound_id] = synonyms

            for n in synonyms + [name]:
                n = filter_search_term(n)
                compound_for_name[n] = compound_id

            if formula != '' and '*' not in formula:
                compound_formula[compound_id] = Formula.parse(formula)

    if args.which == 'compound':
        compound_ids = []
        for name in args.name:
            search_name = filter_search_term(name)
            if search_name in compound_for_name:
                compound_ids.append(compound_for_name[search_name])
        for compound_id in args.id:
            if compound_id in compounds:
                compound_ids.append(compound_id)

        # Show results
        for compound_id in compound_ids:
            print 'ID: {}'.format(compound_id)
            if compound_id in compounds:
                print 'Name: {}'.format(compounds[compound_id])
            if compound_id in compound_synonyms:
                print 'Synonyms: {}'.format(', '.join(compound_synonyms[compound_id]))
            if compound_id in compound_formula:
                print 'Formula: {}'.format(compound_formula[compound_id])
            print
    elif args.which == 'reaction':
        reaction_ids = []
        for reaction_id in args.id:
            if reaction_id in database.reactions:
                reaction_ids.append(reaction_id)
        for compound_list in args.compound:
            matched_reactions = None
            for compound_spec in compound_list.split(','):
                compound = parse_compound(compound_spec.strip())
                if matched_reactions is None:
                    matched_reactions = set(database.compound_reactions.get(compound, set()))
                else:
                    matched_reactions &= database.compound_reactions.get(compound, set())
            if matched_reactions is not None:
                reaction_ids.extend(matched_reactions)

        # Show results
        for reaction_id in reaction_ids:
            print 'ID: {}'.format(reaction_id)

            reaction = database.get_reaction(reaction_id)
            print 'Reaction (IDs): {}'.format(database.get_reaction(reaction_id))
            translated_reaction = reaction.translated_compounds(lambda x: compounds.get(x, x))
            if reaction != translated_reaction:
                print 'Reaction (names): {}'.format(translated_reaction)
            print
