#!/usr/bin/env python

import csv
import argparse
import reaction

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Convert Suden reaction table to standard format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    
    args = parser.parse_args()

    rxn_table = args.rxnfile
    
    rxn_table.readline() # Skip header
    for row in csv.reader(rxn_table, delimiter='\t'):
        SEED_rid, rxn_name, equation_cpd, equation_name, Subsystem, Suden_protien_ID, gene_association, protein_id, ec_reference, Reversibility = row[:10]
        
        rx = reaction.SudenSimple.parse(equation_cpd).normalized()
        if Reversibility == 'N':
            rx.direction = '=>'
        SEED_rid = SEED_rid.lower()
        
        print '{}\t{}'.format(SEED_rid, rx)
    
    rxn_table.close()
    