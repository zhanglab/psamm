#!/usr/bin/env python

'''Run robustness analysis on metabolic model

Given a reaction to maximize and a reaction to vary,
the robustness analysis will run FBA while fixing the
reaction to vary at each iteration. The reaction will
be fixed at the specified number of steps between the
minimum and maximum flux value specified in the model.'''

import sys
import argparse

from metnet.metabolicmodel import DictDatabase
from metnet.fluxanalysis import flux_balance, flux_minimization
from metnet import modelseed

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run robustness analysis on a metabolic model')
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
    parser.add_argument('--steps', metavar='N', help='Number of flux value steps for varying reaction',
                        type=int, default=10)
    parser.add_argument('--minimum', metavar='V', help='Minumum flux value of varying reacton',
                        type=float)
    parser.add_argument('--maximum', metavar='V', help='Maximum flux value of varying reacton',
                        type=float)
    parser.add_argument('varying', help='Reaction to vary')
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

    varying_reaction = args.varying
    if not model.has_reaction(varying_reaction):
        sys.stderr.write('Specified reaction is not in model: {}\n'.format(varying_reaction))
        sys.exit(-1)

    steps = args.steps
    if steps <= 0:
        sys.stderr.write('Invalid number of steps: {}\n'.format(steps))
        sys.exit(-1)

    # Run FBA on model at different fixed flux values
    flux_min = model.limits[varying_reaction].lower if args.minimum is None else args.minimum
    flux_max = model.limits[varying_reaction].upper if args.maximum is None else args.maximum

    if flux_min > flux_max:
        sys.stderr.write('Invalid flux range: {}, {}'.format(flux_min, flux_max))
        sys.exit(-1)

    # Load reaction limits
    for limits_table in args.limits:
        model.load_reaction_limits(limits_table)

    for i in xrange(steps):
        fixed_flux = flux_min + i*(flux_max - flux_min)/float(steps-1)
        test_model = model.copy()
        test_model.limits[varying_reaction].bounds = fixed_flux, fixed_flux

        try:
            fba_fluxes = dict(flux_balance(test_model, reaction))
            optimum = fba_fluxes[reaction]
            fmin_fluxes = dict(flux_minimization(test_model, { reaction: optimum }))
            for other_reaction, flux in fmin_fluxes.iteritems():
                print '{}\t{}\t{}'.format(other_reaction, fixed_flux, flux)
        except Exception:
            pass
