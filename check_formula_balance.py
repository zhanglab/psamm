#!/usr/bin/env python

'''Check whether reactions in a given database or model are balanced

Balanced reactions are those reactions where the number of atoms
is consistent on the left and right side of the reaction equation.
Reactions that are not balanced will be printed out.'''

import operator

from metnet import command
from metnet import modelseed
from metnet.formula import Formula, Radical

class FormulaBalanceCommand(command.Command):
    def __init__(self):
        super(FormulaBalanceCommand, self).__init__('Check formula balance on a model or database')

    def __call__(self, database, model, compounds):
        '''Run formula balance command'''

        # Mapping from compound id to formula
        compound_formula = {}
        for compound_table in compounds:
            for compound in modelseed.parse_compound_file(compound_table):
                # Create pseudo-radical group for compounds with
                # missing formula, so they don't match up. Only
                # cpd11632 (Photon) is allowed to have an empty formula.
                if compound.formula is None or '.' in compound.formula:
                    if compound.id != 'cpd11632':
                        f = Formula({Radical('R'+compound.id): 1})
                    else:
                        f = Formula()
                else:
                    try:
                        f = Formula.parse(compound.formula).flattened()
                    except ValueError as e:
                        print 'Error parsing {}: {}'.format(compound.formula, e)
                compound_formula[compound.id] = f

        # Create a set of known mass-inconsistent reactions
        exchange = set()
        for rxnid in database.reactions:
            rx = database.get_reaction(rxnid)
            if len(rx.left) == 0 or len(rx.right) == 0:
                exchange.add(rxnid)

        def multiply_formula(compound_list):
            for compound, count in compound_list:
                yield count * compound_formula.get(compound.name, Formula())

        for reaction in database.reactions:
            if reaction not in exchange:
                rx = database.get_reaction(reaction)
                left_form = reduce(operator.or_, multiply_formula(rx.left), Formula())
                right_form = reduce(operator.or_, multiply_formula(rx.right), Formula())

                if right_form != left_form:
                    right_missing, left_missing = Formula.balance(right_form, left_form)
                    print '{}\t{}\t{}\t{}\t{}'.format(reaction, left_form, right_form, left_missing, right_missing)

if __name__ == '__main__':
    command.main(FormulaBalanceCommand())
