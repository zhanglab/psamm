#!/usr/bin/env python

import argparse
import csv

from metnet.metabolicmodel import DictDatabase
from metnet.gapfill import gapfind, gapfill

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run GapFind and GapFill on a metabolic model')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to use as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    args = parser.parse_args()

    database = DictDatabase.load_from_files(*args.database)
    model = database.load_model_from_file(args.reactionlist)

    # Load compound information
    compounds = {}
    for compound_table in args.compounds:
        compound_table.readline() # Skip header
        for row in csv.reader(compound_table, delimiter='\t'):
            cpdid, names = row[:2]
            synonyms = names.split(',<br>')
            name = synonyms.pop()
            compounds[cpdid] = name

    model_compartments = { None, 'e' }

    # Run GapFind on model
    print 'Searching for blocked compounds...'
    blocked = set(compound for compound in gapfind(model) if compound.compartment is not 'e')
    if len(blocked) > 0:
        print 'Blocked:'
        for compound in blocked:
            print compound.translate(lambda x: compounds.get(x, x))

    if len(blocked) > 0:
        # Add exchange and transport reactions to database
        model_complete = model.copy()
        print 'Adding database, exchange and transport reactions...'
        model_complete.add_all_database_reactions(model_compartments)
        model_complete.add_all_exchange_reactions()
        model_complete.add_all_transport_reactions()

        print 'Searching for reactions to fill gaps...'
        added_reactions, reversed_reactions = gapfill(model_complete, model.reactions, blocked)

        for rxnid in added_reactions:
            rx = database.get_reaction(rxnid).translated_compounds(lambda x: compounds.get(x, x))
            print '{}\t{}\t{}'.format(rxnid, 'Add', rx)

        for rxnid in reversed_reactions:
            rx = database.get_reaction(rxnid).translated_compounds(lambda x: compounds.get(x, x))
            print '{}\t{}\t{}'.format(rxnid, 'Reverse', rx)
    else:
        print 'No blocked compounds found'
