#!/usr/bin/env python
'''This program compares a Changed Reaction List to and Original Reaction List. It lists what compounds have
been added and removed from the Original List to make the Changed List.'''

import csv
import argparse

# Main program
if __name__ == '__main__':
    # sys.argv[1] contains the file name to open
    
    parser = argparse.ArgumentParser(description='Compares Gene IDs')
    parser.add_argument('ssrxnlist', type=argparse.FileType('r'), help='Changed List')
    parser.add_argument('csrxnlist', type=argparse.FileType('r'), help='Orignial List')

    args = parser.parse_args()

    ss_rxn = args.ssrxnlist
    cs_rxn = args.csrxnlist

    original_set = set()   # Empty Set
    changed_set = set()
    
    # Original List
    cs_rxn.readline() # Skip header
    for rowcs in csv.reader(cs_rxn,delimiter='\t'):
        SEED_rid = rowcs[0]
        
        SEED_rid = SEED_rid.lower()
        original_set.add(SEED_rid)
     
    # Changed List    
    ss_rxn.readline() # Skip header
    for rowss in csv.reader(ss_rxn,delimiter='\t'):
        rxn_id = rowss[0]

        rxn_id = rxn_id.lower()
        changed_set.add(rxn_id)
        
        if rxn_id not in original_set:
            print '{}\tAdded to Original List'.format(rxn_id)

    for SEED_rid in sorted(original_set):
        if SEED_rid.lower() not in changed_set:
            print '{}\tRemoved from Original List'.format(SEED_rid)

        #if rxn_id.startswith('rxn'):
        #    continue
        #else:
        #    for gene_id in suden_protein_id.split(','):
        #        if gene_id in gene_map:
        #            print '{} possibly the same as {}'.format(rxn_id, gene_map[gene_id])

    ss_rxn.close()
    cs_rxn.close()