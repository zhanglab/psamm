#!/usr/bin/env python

import csv
import formula
import reaction
import argparse
import pprint

# Main program
if __name__ == '__main__':
    # sys.argv[1] contains the file name to open
    
    parser = argparse.ArgumentParser(description='Convert reaction table to GapFind input format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('cpdfile', type=argparse.FileType('r'), help='Compound table file')
    
    args = parser.parse_args()

    rxn_file = args.rxnfile
    cpd_file = args.cpdfile
    
    compound_map = {}   #empty dictionary
    
    readerc = csv.reader(cpd_file,delimiter='\t')
    for rowc in readerc:
        SEED_cid, cpd_names, Formula, Mass,KEGG_maps, KEGG_cid = rowc  

        compound_map[SEED_cid] = Formula
        #setting the output of SEED_cid equal to Formula   
    
    readerr = csv.reader(rxn_file,delimiter='\t')
    for rowr in readerr:
        rxn_id,Equation_cpdid = rowr[:2]
        
        if Equation_cpdid.strip() == '':
            continue

        rx = reaction.ModelSEED.parse(Equation_cpdid).normalized()
        
        left_form = {} #empty dictionary
        for cpdid, count, comp in rx.left:
            form = compound_map[cpdid]
            f = formula.parse(form)
            left_form = formula.sum(left_form, formula.multiply(f, count))
        #print left_form
            
        right_form = {} #empty dictionary
        for cpdid, count, comp in rx.right:
            form = compound_map[cpdid]
            f = formula.parse(form)
            right_form = formula.sum(right_form, formula.multiply(f, count))   
        #print right_form    
            
        if right_form != left_form:
                print '{}\t{}\t{}'.format(SEED_rid, pprint.pformat(left_form),
                                            pprint.pformat(right_form))        
            
    rxn_file.close()
    cpd_file.close()  