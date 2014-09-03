#!/usr/bin/env python

import argparse
import csv
from itertools import chain

from metnet.metabolicmodel import MetabolicDatabase
from metnet.reaction import Reaction, Compound
from metnet.fastcore import Fastcore
from metnet.fluxanalysis import flux_balance
from metnet import lpsolver

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run FastGapFill on a metabolic model')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to use as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('reaction', help='Reaction to maximize')
    args = parser.parse_args()

    database = MetabolicDatabase.load_from_files(*args.database)
    model = database.load_model_from_file(args.reactionlist)

    # Create fastcore object
    solver = lpsolver.CplexSolver(None)
    fastcore = Fastcore(solver)

    # Load compound information
    compounds = {}
    for compound_table in args.compounds:
        compound_table.readline() # Skip header
        for row in csv.reader(compound_table, delimiter='\t'):
            cpdid, names = row[:2]
            synonyms = names.split(',<br>')
            name = synonyms.pop()
            compounds[cpdid] = name

    epsilon = 1e-5
    model_compartments = { None, 'e' }

    # Add exchange and transport reactions to database
    model_complete = model.copy()
    print 'Adding database, exchange and transport reactions...'
    model_complete.add_all_database_reactions(model_compartments)
    model_complete.add_all_exchange_reactions()
    model_complete.add_all_transport_reactions()

    # Run Fastcore and print the induced reaction set
    print 'Calculating Fastcore induced set on model...'
    core = set(model.reaction_set)

    induced = fastcore.fastcore(model_complete, core, epsilon)
    print 'Result: |A| = {}, A = {}'.format(len(induced), induced)
    added_reactions = induced - core
    print 'Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions)

    # Load bounds on exchange reactions
    #model.load_exchange_limits()

    print 'Flux balance on induced model maximizing {}...'.format(args.reaction)
    model_induced = model.copy()
    for rxnid in induced:
        model_induced.add_reaction(rxnid)
    for rxnid, flux in sorted(flux_balance(model_induced, args.reaction)):
        reaction_class = 'Dbase'
        if rxnid in model.reaction_set:
            reaction_class = 'Model'
        reaction = database.get_reaction(rxnid).translated_compounds(lambda x: compounds.get(x, x))
        print '{}\t{}\t{}\t{}'.format(rxnid, reaction_class, flux, reaction)

    print 'Calculating Fastcc consistent subset of induced model...'
    consistent_core = fastcore.fastcc_consistent_subset(model_induced, epsilon)
    print 'Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core)
    removed_reactions = model_induced.reaction_set - consistent_core
    print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)
