#!/usr/bin/env python

from metnet import command
from metnet.massconsistency import MassConsistencyCheck
from metnet import modelseed

class MassConsistencyCommand(command.Command):
    def __init__(self):
        super(MassConsistencyCommand, self).__init__('Run mass consistency check on a database')

    def init_parser(self, parser):
        parser.add_argument('--exclude', metavar='reaction', action='append',
                            type=str, default=[], help='Exclude reaction from mass consistency')

    def __call__(self, database, model, compounds=[], **kwargs):
        # Load compound information
        compound_name = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                compound_name[compound.id] = compound.name if compound.name is not None else compound.id

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for reaction_id in database.reactions:
            rx = database.get_reaction(reaction_id)
            if len(rx.left) == 0 or len(rx.right) == 0:
                exchange.add(reaction_id)

        # Other reactions to exclude from consistency check
        exclude = set(kwargs['exclude'])

        # Create set of compounds allowed to have mass zero
        zeromass = set()
        zeromass.add('cpd11632') # Photon
        zeromass.add('cpd12713') # Electron

        mass_consistency = MassConsistencyCheck()
        known_inconsistent = exclude | exchange

        print 'Mass consistency on database...'
        epsilon = 1e-5
        compound_iter = mass_consistency.check_compound_consistency(database, known_inconsistent, zeromass)

        print 'Compound consistency...'
        good = 0
        total = 0
        for compound, mass in sorted(compound_iter, key=lambda x: (x[1], x[0]), reverse=True):
            if mass >= 1-epsilon or compound.name in zeromass:
                good += 1
            total += 1
            print '{}: {}'.format(compound.translate(lambda x: compound_name.get(x, x)), mass)
        print 'Consistent compounds: {}/{}'.format(good, total)

        print 'Is consistent? {}'.format(mass_consistency.is_consistent(database, known_inconsistent, zeromass))

        print 'Reaction consistency...'
        reaction_iter, compound_iter = mass_consistency.check_reaction_consistency(database, known_inconsistent, zeromass)
        for reaction_id, residual in sorted(reaction_iter, key=lambda x: abs(x[1]), reverse=True):
            if abs(residual) >= epsilon:
                reaction = database.get_reaction(reaction_id).translated_compounds(lambda x: compound_name.get(x, x))
                print '{}\t{}\t{}'.format(reaction_id, residual, reaction)


if __name__ == '__main__':
    command.main(MassConsistencyCommand())
