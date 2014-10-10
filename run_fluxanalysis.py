#!/usr/bin/env python

import sys
import argparse

from metnet.metabolicmodel import DictDatabase
from metnet.fluxanalysis import flux_balance, flux_minimization
from metnet import modelseed

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run flux analysis on a metabolic model')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to usa as database')
    parser.add_argument('--compounds', metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('--limits', metavar='limitsfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional limits on flux of reactions')
    parser.add_argument('reactionlist', type=argparse.FileType('r'), help='Model definition')
    parser.add_argument('reaction', help='Reaction to maximize')
    args = parser.parse_args()

    database = DictDatabase.load_from_files(*args.database)
    model = database.load_model_from_file(args.reactionlist)

    # Load compound information
    compounds = {}
    for compound_table in args.compounds:
        for compound in modelseed.parse_compound_file(compound_table):
            compounds[compound.id] = compound.name if compound.name is not None else compound.id

    reaction = args.reaction
    if not model.has_reaction(reaction):
        sys.stderr.write('Specified reaction is not in model: {}\n'.format(reaction))
        sys.exit(-1)

    for limits_table in args.limits:
        model.load_reaction_limits(limits_table)

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