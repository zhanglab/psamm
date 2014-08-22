#!/usr/bin/env python

'''Check whether reactions in a given database are balanced

Balanced reactions are those reactions where the number of atoms
is consistent on the left and right side of the reaction equation.
Reactions that are not balanced will be printed out.'''

import csv
import argparse

from metnet.formula import Formula, Radical
from metnet.reaction import ModelSEED

# Main program
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert reaction table to GapFind input format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('cpdfile', type=argparse.FileType('r'), help='Compound table file')

    args = parser.parse_args()

    rxn_file = args.rxnfile
    cpd_file = args.cpdfile

    # Mapping from compound id to formula
    compound_formula = {}

    readerc = csv.reader(cpd_file,delimiter='\t')
    for rowc in readerc:
        SEED_cid, cpd_names, formula, Mass,KEGG_maps, KEGG_cid = rowc

        # Create pseudo-radical group for compounds with
        # missing formula, so they don't match up. Only
        # cpd11632 (Photon) is allowed to have an empty formula.
        if (formula.strip() == '' and SEED_cid != 'cpd11632') or '*' in formula:
            f = Formula({Radical('R'+SEED_cid): 1})
        else:
            f = Formula.parse(formula)

        compound_formula[SEED_cid] = f

    readerr = csv.reader(rxn_file,delimiter='\t')
    for rowr in readerr:
        rxn_id,Equation_cpdid = rowr[:2]

        if Equation_cpdid.strip() == '':
            continue

        rx = ModelSEED.parse(Equation_cpdid).normalized()

        def multiply_formula(compound_list):
            for compound, count, comp in compound_list:
                if compound.name in compound_formula:
                    yield count * compound_formula[compound.name]

        left_form = sum(multiply_formula(rx.left), Formula())
        right_form = sum(multiply_formula(rx.right), Formula())

        if right_form != left_form:
            right_missing, left_missing = Formula.balance(right_form, left_form)
            print '{}\t{}\t{}\t{}\t{}'.format(rxn_id, left_form, right_form, left_missing, right_missing)

    rxn_file.close()
    cpd_file.close()