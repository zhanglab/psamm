#!/usr/bin/env python
'''This program compares SimpleSuden to the ModelSuden. It lists what compounds have
been added and removed from ModelSuden and the changes in direction made to the equations
It also lists which of the 'changed' compounds have been only changed or changed and added.'''

import csv
import argparse

from metnet.reaction import Reaction, ModelSEED, SimpleSuden

# Main program
if __name__ == '__main__':
    # sys.argv[1] contains the file name to open

    parser = argparse.ArgumentParser(description='Compares Gene IDs')
    parser.add_argument('ssrxnlist', type=argparse.FileType('r'), help='SimpleSuden Reaction List')
    parser.add_argument('csrxnlist', type=argparse.FileType('r'), help='Raw ModelSEED Suden')

    args = parser.parse_args()

    ss_rxn = args.ssrxnlist
    cs_rxn = args.csrxnlist

    #gene_map = {}   # Empty Dictionary
    ModelSEED_reactions = {}   # Empty Dictionary
    SimpleSuden_reactions = {}

    cs_rxn.readline() # Skip header
    for rowcs in csv.reader(cs_rxn,delimiter='\t'):
        SEED_rid, rxn_name, EC, equation_cpdname, equation_cpdid, KEGG_rid, KEGG_maps, gene_ids = rowcs[:8]

        ModelSEED_reactions[SEED_rid] = ModelSEED.parse(equation_cpdid).normalized()

        #for gene_id in gene_ids.split(','):
        #    gene_map[gene_id] = SEED_rid

    ss_rxn.readline() # Skip header
    for rowss in csv.reader(ss_rxn,delimiter='\t'):
        rxn_id, annotation, equation_cpd, equation_name, subsystem, suden_protein_id, gene_attributes, protein_id, ec_num, reference, reversibility = rowss[:11]

        rx = SudenSimple.parse(equation_cpd).normalized()
        if reversibility == 'N':
            rx = Reaction('=>', rx.left, rx.right)

        SimpleSuden_reactions[rxn_id] = rx

        #if rxn_id.lower() not in reactions:
        if rxn_id not in ModelSEED_reactions:

            rxn_id_8 = rxn_id[:8]

            if rxn_id_8.lower() in ModelSEED_reactions:
                print '{}\tChanged from ModelSuden'.format(rxn_id)
            else:
                #if rxn_id.startswith('rxn'):
                if rxn_id_8 != rxn_id_8.lower():
                    rxn_id_lower = rxn_id_8.lower()
                    if rxn_id_lower.startswith('r'):
                        print '{}\tChanged and added to SimpleSuden'.format(rxn_id)
                    else:
                        print '{}\tAdded to SimpleSuden'.format(rxn_id)
                else:
                    print '{}\tAdded to SimpleSuden'.format(rxn_id)

        else:
            #other_rx = reactions[rxn_id.lower()]
            other_rx = ModelSEED_reactions[rxn_id]
            if rx != other_rx:
                rx_copy = rx.copy()
                rx_copy.direction = other_rx.direction
                if rx_copy == other_rx:
                    print '{}\tEquations have different directions\t{}to{}'.format(rxn_id,other_rx.direction,rx.direction)
                else:
                    print '{}\tEquations are different: {} vs. {}'.format(rxn_id,rx, other_rx)

    for SEED_rid in sorted(ModelSEED_reactions):
        if SEED_rid not in SimpleSuden_reactions:
            print '{}\tRemoved from ModelSuden'.format(SEED_rid)

        #if rxn_id.startswith('rxn'):
        #    continue
        #else:
        #    for gene_id in suden_protein_id.split(','):
        #        if gene_id in gene_map:
        #            print '{} possibly the same as {}'.format(rxn_id, gene_map[gene_id])

    ss_rxn.close()
    cs_rxn.close()