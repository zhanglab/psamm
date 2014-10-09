#!/usr/bin/env python

from metnet import command
from metnet.gapfill import gapfind, gapfill
from metnet import modelseed

class GapFillCommand(command.Command):
    def __init__(self):
        super(GapFillCommand, self).__init__('Run GapFind and GapFill on a metabolic model')

    def __call__(self, database, model, compounds=[]):
        '''Run GapFill command'''

        if model is None:
            sys.stderr.write('Please specify a model. Use --help for more information.\n')
            sys.exit(-1)

        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        model_compartments = { None, 'e' }

        # Run GapFind on model
        print 'Searching for blocked compounds...'
        blocked = set(compound for compound in gapfind(model) if compound.compartment is not 'e')
        if len(blocked) > 0:
            print 'Blocked:'
            for compound in blocked:
                print compound.translate(lambda x: compound_name.get(x, x))

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
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Add', rx)

            for rxnid in reversed_reactions:
                rx = model_complete.get_reaction(rxnid).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(rxnid, 'Reverse', rx)
        else:
            print 'No blocked compounds found'


if __name__ == '__main__':
    command.main(GapFillCommand())
