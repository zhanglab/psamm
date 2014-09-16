#!/usr/bin/env python

'''Run robustness analysis on metabolic model

Given a reaction to maximize and a reaction to vary,
the robustness analysis will run FBA while fixing the
reaction to vary at each iteration. The reaction will
be fixed at the specified number of steps between the
minimum and maximum flux value specified in the model.'''

import sys
import argparse
import csv

from metnet.metabolicmodel import MetabolicDatabase
from metnet.fluxanalysis import flux_balance

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run robustness analysis on a metabolic model')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to usa as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('reaction', help='Reaction to maximize')
    parser.add_argument('varying', help='Reaction to vary')
    parser.add_argument('--steps', metavar='N', help='Number of flux value steps for varying reaction',
                        type=int, default=10)
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

    varying_reaction = args.varying
    if varying_reaction not in model.reaction_set:
        sys.stderr.write('Specified reaction is not in model: {}\n'.format(varying_reaction))
        sys.exit(-1)

    steps = args.steps
    if steps <= 0:
        sys.stderr.write('Invalid number of steps: {}\n'.format(steps))
        sys.exit(-1)

    # Run FBA on model at different fixed flux values
    flux_min = model.limits[varying_reaction].lower
    flux_max = model.limits[varying_reaction].upper

    model.load_exchange_limits()

    for i in xrange(steps):
        fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
        test_model = model.copy()
        test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

        try:
            fba_fluxes = dict(flux_balance(test_model, reaction))
            for reaction, flux in fba_fluxes.iteritems():
                rx = database.get_reaction(reaction)
                print '{}\t{}\t{}\t{}'.format(reaction, fixed_flux, flux, rx.translated_compounds(lambda x: compounds.get(x, x)))
        except Exception:
            pass
