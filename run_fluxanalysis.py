#!/usr/bin/env python

import sys
import argparse

from metnet import command
from metnet.metabolicmodel import DictDatabase
from metnet.fluxanalysis import flux_balance, flux_minimization
from metnet import modelseed

class FluxAnalysisCommand(command.Command):
    def __init__(self):
        super(FluxAnalysisCommand, self).__init__('Run flux analysis on a metabolic model')

    def init_parser(self, parser):
        parser.add_argument('--limits', metavar='limitsfile', action='append',
                            type=argparse.FileType('r'), default=[],
                            help='Optional limits on flux of reactions')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run flux analysis command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        reaction = kwargs['reaction']
        if not model.has_reaction(reaction):
            sys.stderr.write('Specified reaction is not in model: {}\n'.format(reaction))
            sys.exit(-1)

        for limits_table in kwargs['limits']:
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
            print '{}\t{}\t{}\t{}'.format(rxnid, fba_fluxes[rxnid], flux,
                                            rx.translated_compounds(lambda x: compound_name.get(x, x)))
        print 'Maximum flux: {}'.format(optimum)
        print 'Minimized reactions: {}'.format(count)

if __name__ == '__main__':
    command.main(FluxAnalysisCommand())
