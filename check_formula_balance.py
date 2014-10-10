#!/usr/bin/env python

'''Check whether reactions in a given database or model are balanced

Balanced reactions are those reactions where the number of atoms
is consistent on the left and right side of the reaction equation.
Reactions that are not balanced will be printed out.'''

import argparse
import re
import operator

from metnet.metabolicmodel import DictDatabase
from metnet.formula import Formula, Radical
from metnet.reaction import ModelSEED
from metnet import modelseed

# Main program
if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Check formula balance on a model or database')
    parser.add_argument('--database', required=True, metavar='reactionfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Reaction definition list to use as database')
    parser.add_argument('--compounds', required=True, metavar='compoundfile', action='append',
                        type=argparse.FileType('r'), default=[],
                        help='Optional compound information table')
    parser.add_argument('reactionlist', nargs='?', type=argparse.FileType('r'),
                        help='Model definition')
    args = parser.parse_args()

    database = DictDatabase.load_from_files(*args.database)

    # Load model from file if given, otherwise run on full database
    if args.reactionlist:
        model = database.load_model_from_file(args.reactionlist)
    else:
        model = database.get_model(database.reactions)

    # Mapping from compound id to formula
    compound_formula = {}
    for compound_table in args.compounds:
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
    for rxnid in model.reactions:
        rx = database.get_reaction(rxnid)
        if len(rx.left) == 0 or len(rx.right) == 0:
            exchange.add(rxnid)

    def multiply_formula(compound_list):
        for compound, count in compound_list:
            yield count * compound_formula.get(compound.name, Formula())

    for reaction in model.reactions:
        if reaction not in exchange:
            rx = database.get_reaction(reaction)
            left_form = reduce(operator.or_, multiply_formula(rx.left), Formula())
            right_form = reduce(operator.or_, multiply_formula(rx.right), Formula())

            if right_form != left_form:
                right_missing, left_missing = Formula.balance(right_form, left_form)
                print '{}\t{}\t{}\t{}\t{}'.format(reaction, left_form, right_form, left_missing, right_missing)
