#!/usr/bin/env python

import sys
import argparse
import csv

from metnet.metabolicmodel import MetabolicDatabase
from metnet.fluxanalysis import flux_balance, flux_minimization

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run flux analysis on a metabolic model')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to usa as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('reaction', help='Reaction to maximize')
    args = parser.parse_args()

    database = MetabolicDatabase.load_from_files(*args.database)
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

    reaction = args.reaction
    if reaction not in model.reaction_set:
        sys.stderr.write('Specified reaction is not in model: {}\n'.format(reaction))
        sys.exit(-1)

    model.load_exchange_limits()

    epsilon = 1e-5

    # Run FBA on model
    fba_fluxes = dict(flux_balance(model, reaction))
    optimum = fba_fluxes[reaction]

    # Run flux minimization
    fmin_fluxes = dict(flux_minimization(model, { reaction: optimum }))
    count = 0
    for rxnid, flux in sorted(fmin_fluxes.iteritems()):
        if fba_fluxes[rxnid] - epsilon > flux:
            count += 1
        rx = database.get_reaction(rxnid)
        print '{}\t{}\t{}\t{}'.format(rxnid, fba_fluxes[rxnid], flux, rx.translated_compounds(lambda x: compounds.get(x, x)))
    print 'Maximum flux: {}'.format(optimum)
    print 'Minimized reactions: {}'.format(count)