#!/usr/bin/env python

import argparse
import csv
from itertools import chain

from metabolicmodel import MetabolicDatabase
from reaction import Reaction, Compound
import fastcore as fc
import fluxanalysis

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run FastGapFill on a metabolic model')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to usa as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    args = parser.parse_args()

    database = MetabolicDatabase.load_from_files(*args.database)
    model = database.load_model_from_file(args.reactionlist)

    # Create fastcore object
    fastcore = fc.Fastcore()

    # Load compound information
    compounds = {}
    for compound_table in args.compounds:
        compound_table.readline() # Skip header
        for row in csv.reader(compound_table, delimiter='\t'):
            cpdid, names = row[:2]
            synonyms = names.split(',<br>')
            name = synonyms.pop()
            compounds[cpdid] = name

    # Run Fastcc and print the inconsistent set
    print 'Calculating Fastcc consistent subset...'
    consistent_core = fastcore.fastcc_consistent_subset(model, 0.001)
    print 'Result: |A| = {}, |A| = {}'.format(len(consistent_core), consistent_core)

    all_reactions = {}
    for rxnid in database.reactions:
        rx = database.get_reaction(rxnid)
        all_reactions[rx] = rxnid

    model_compartments = set([None, 'e'])

    # Add exchange and transport reactions to database
    print 'Adding exchange and transport reactions...'
    for cpdid, comp in database.compounds:
        if comp in model_compartments:
            rxnid_ex = 'rxnex_'+cpdid
            reaction_ex = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [])
            if reaction_ex not in all_reactions:
                print '{}\t{}'.format(rxnid_ex, reaction_ex)
                database.set_reaction(rxnid_ex, reaction_ex)
                all_reactions[reaction_ex] = rxnid_ex
            else:
                print '{}\tAlready in database! ({})'.format(rxnid_ex, all_reactions[reaction_ex])

            rxnid_tp = 'rxntp_'+cpdid
            reaction_tp = Reaction('<=>', [(Compound(cpdid), 1, 'e')], [(Compound(cpdid), 1, None)])
            if reaction_tp not in all_reactions:
                print '{}\t{}'.format(rxnid_tp, reaction_tp)
                database.set_reaction(rxnid_tp, reaction_tp)
                all_reactions[reaction_tp] = rxnid_tp
            else:
                print '{}\tAlready in database! ({})'.format(rxnid_tp, all_reactions[reaction_tp])

    # Add reaction from database to model. Do not add reactions from compartments
    # that are not in the model.
    model_complete = model.copy()
    for rxnid in database.reactions:
        reaction = database.get_reaction(rxnid)
        if any(comp not in model_compartments for compound, value, comp in chain(reaction.left, reaction.right)):
            continue
        model_complete.add_reaction(rxnid)

    #print 'Calculating Fastcc consistent subset of database...'
    #database_consistent = fastcore.fastcc_consistent_subset(model_complete, 0.001)
    #print 'Result: |A| = {}, A = {}'.format(len(database_consistent), database_consistent)
    #removed_reactions = model_complete.reaction_set - database_consistent
    #print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)

    # Run Fastcore and print the induced reaction set
    print 'Calculating Fastcore induced set on model...'
    core = model.reaction_set

    induced = fastcore.fastcore(model_complete, core, 0.001)
    print 'Result: |A| = {}, A = {}'.format(len(induced), induced)
    added_reactions = induced - core
    print 'Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions)

    # Load bounds on exchange reactions
    #model.load_exchange_limits()

    print 'Flux balance on original model maximizing growth...'
    for rxnid, flux in sorted(fluxanalysis.flux_balance(model, 'Growth')):
        print '{}\t{}'.format(rxnid, flux)

    print 'Flux balance on induced model maximizing growth...'
    model_induced = model.copy()
    for rxnid in induced:
        model_induced.add_reaction(rxnid)
    for rxnid, flux in sorted(fluxanalysis.flux_balance(model_induced, 'Growth')):
        reaction_class = 'Dbase'
        if rxnid in consistent_core:
            reaction_class = 'Core'
        elif rxnid in model.reaction_set:
            reaction_class = 'Model'
        reaction = database.get_reaction(rxnid).translated_compounds(lambda x: compounds.get(x, x))
        print '{}\t{}\t{}\t{}'.format(rxnid, reaction_class, flux, reaction)

    print 'Calculating Fastcc consistent subset of induced model...'
    consistent_core = fastcore.fastcc_consistent_subset(model_induced, 0.001)
    print 'Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core)
    removed_reactions = model_induced.reaction_set - consistent_core
    print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)
