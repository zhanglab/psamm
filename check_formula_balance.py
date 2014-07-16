#!/usr/bin/env python

import csv
from formula import Formula
import reaction
import argparse
import pprint

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

        compound_formula[SEED_cid] = Formula.parse(formula)

    readerr = csv.reader(rxn_file,delimiter='\t')
    for rowr in readerr:
        rxn_id,Equation_cpdid = rowr[:2]

        if Equation_cpdid.strip() == '':
            continue

        rx = reaction.ModelSEED.parse(Equation_cpdid).normalized()

        left_form = sum((count * compound_formula[cpdid] for cpdid, count, comp in rx.left), Formula())
        right_form = sum((count * compound_formula[cpdid] for cpdid, count, comp in rx.right), Formula())

        if right_form != left_form:
                print '{}\t{}\t{}'.format(rxn_id, pprint.pformat(left_form),
                                            pprint.pformat(right_form))

    rxn_file.close()
    cpd_file.close()