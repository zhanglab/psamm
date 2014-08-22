#!/usr/bin/env python

import csv
import argparse

from metnet.reaction import Reaction, SimpleSuden

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert Suden reaction table to standard format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')

    args = parser.parse_args()

    rxn_table = args.rxnfile

    rxn_table.readline() # Skip header
    for row in csv.reader(rxn_table, delimiter='\t'):
        SEED_rid, rxn_name, equation_cpd, equation_name, Subsystem, Suden_protien_ID, gene_association, protein_id, ec_num, reference, Reversibility = row[:11]

        rx = SudenSimple.parse(equation_cpd).normalized()
        if Reversibility == 'N':
            rx = Reaction('=>', rx.left, rx.right)

        if not SEED_rid.islower() and SEED_rid[:3].lower() == 'rxn':
            SEED_rid = 'simplesuden_' + SEED_rid.lower()
        else:
            SEED_rid = SEED_rid.lower()

        print '{}\t{}'.format(SEED_rid, rx)

    rxn_table.close()
