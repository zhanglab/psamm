#!/usr/bin/env python

import sys
import argparse

from metnet import command
from metnet.reaction import Reaction, Compound
from metnet.fastcore import Fastcore
from metnet.fluxanalysis import flux_balance
from metnet import lpsolver, modelseed

class FastGapFillCommand(command.Command):
    def __init__(self):
        super(FastGapFillCommand, self).__init__('Run FastGapFill on a metabolic model')

    def init_parser(self, parser):
        parser.add_argument('--limits', metavar='file', action='append',
                            type=argparse.FileType('r'), default=[],
                            help='Optional limits on flux of reactions')
        parser.add_argument('--penalty', metavar='file', type=argparse.FileType('r'),
                            help='List of penalty scores for database reactions')
        parser.add_argument('reaction', help='Reaction to maximize')

    def __call__(self, database, model, compounds=[], **kwargs):
        '''Run FastGapFill command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Create fastcore object
        solver = lpsolver.CplexSolver(None)
        fastcore = Fastcore(solver)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

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

        if kwargs.get('penalty') is not None:
            for line in kwargs['penalty']:
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
        for limits_table in kwargs['limits']:
            model.load_reaction_limits(limits_table)

        maximized_reaction = kwargs['reaction']
        print 'Flux balance on induced model maximizing {}...'.format(maximized_reaction)
        model_induced = model.copy()
        for rxnid in induced:
            model_induced.add_reaction(rxnid)
        for rxnid, flux in sorted(flux_balance(model_induced, maximized_reaction)):
            reaction_class = 'Dbase'
            weight = weights.get(rxnid, 1)
            if model.has_reaction(rxnid):
                reaction_class = 'Model'
                weight = 0
            reaction = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
            print '{}\t{}\t{}\t{}\t{}'.format(rxnid, reaction_class, weight, flux, reaction)

        print 'Calculating Fastcc consistent subset of induced model...'
        consistent_core = fastcore.fastcc_consistent_subset(model_induced, epsilon)
        print 'Result: |A| = {}, A = {}'.format(len(consistent_core), consistent_core)
        removed_reactions = set(model_induced.reactions) - consistent_core
        print 'Removed: |R| = {}, R = {}'.format(len(removed_reactions), removed_reactions)


if __name__ == '__main__':
    command.main(FastGapFillCommand())
