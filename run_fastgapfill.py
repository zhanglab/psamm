#!/usr/bin/env python

import argparse
import csv
from itertools import chain

from metnet.metabolicmodel import DictDatabase
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
    parser.add_argument('--limits', metavar='limitsfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional limits on flux of reactions')
    parser.add_argument('--penalty', metavar='penaltyfile', type=argparse.FileType('r'),
                        help='List of penalty scores for database reactions')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('reaction', help='Reaction to maximize')
    args = parser.parse_args()

    database = DictDatabase.load_from_files(*args.database)
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
    db_added = model_complete.add_all_database_reactions(model_compartments)
    ex_added = model_complete.add_all_exchange_reactions()
    tp_added = model_complete.add_all_transport_reactions()

    # Add penalty weights on reactions
    weights = {}
    #weights.update((rxnid, 25) for rxnid in db_added)
    weights.update((rxnid, 50) for rxnid in tp_added)
    weights.update((rxnid, 250) for rxnid in ex_added)

    if args.penalty:
        for line in args.penalty:
            line, _, comment = line.partition('#')
            line = line.strip()
            if line == '':
                continue
            rxnid, weight = line.split(None, 1)
            weights[rxnid] = float(weight)

    # Run Fastcore and print the induced reaction set
    print 'Calculating Fastcore induced set on model...'
    core = set(model.reactions)

    induced = fastcore.fastcore(model_complete, core, epsilon, weights=weights)
    print 'Result: |A| = {}, A = {}'.format(len(induced), induced)
    added_reactions = induced - core
    print 'Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions)

    # Load bounds on exchange reactions
    for limits_table in args.limits:
        model.load_reaction_limits(limits_table)

    print 'Flux balance on induced model maximizing {}...'.format(args.reaction)
    model_induced = model.copy()
    for rxnid in induced:
        model_induced.add_reaction(rxnid)
    for rxnid, flux in sorted(flux_balance(model_induced, args.reaction)):
        reaction_class = 'Dbase'
        weight = weights.get(rxnid, 1)
        if model.has_reaction(rxnid):
            reaction_class = 'Model'
            weight = 0
        reaction = database.get_reaction(rxnid).translated_compounds(lambda x: compounds.get(x, x))
        print '{}\t{}\t{}\t{}\t{}'.format(rxnid, reaction_class, weight, flux, reaction)

    print 'Calculating Fastcc consistent subset of induced model...'
    consistent_core = fastcore.fastcc_consistent_subset(model_induced, epsilon)
    print 'Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core)
    removed_reactions = set(model_induced.reactions) - consistent_core
    print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)
