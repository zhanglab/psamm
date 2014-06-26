#!/usr/bin/env python

'''Parse GapFill output and print added reactions

The GapFill file should be produced with parameter "ps=0".
'''

import csv
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse GapFill output')
    parser.add_argument('gapfill', type=argparse.FileType('r'), help='GapFill output file (ps=0)')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('compoundfile', type=argparse.FileType('r'), help='Compound Map file')
    
    args = parser.parse_args()

    reaction_map = {}
    compound_map = {}

    rxnfile = args.rxnfile
    rxnfile.readline() # Skip header
    for row in csv.reader(rxnfile, delimiter='\t'):
        seed_rid, name, equation = row[:3]
        reaction_map[seed_rid] = name, equation
    rxnfile.close()

    compoundfile = args.compoundfile
    compoundfile.readline() # Skip header
    for row in csv.reader(compoundfile, dialect='excel'):
        SEED_cid, Formula, Mass, KEGG_cid, cpd_names = row
        main_name = cpd_names.split(',<br>')[0]
        compound_map[SEED_cid] = Formula, main_name
    compoundfile.close()

    # Parse GapFind output
    count = 0
    count_not_produced = 0
    count_in_kegg = 0
    flag = False
    gffile = args.gapfill
    for line in gffile:
        if line.startswith('---- VAR yd'):
            flag = True
            for i in range(3):
                next(gffile)
        elif flag:
            if line == '\n':
                break

            fields = line.split()
            
            rxn_id = fields[0]
            enabled = fields[2] != '.'
            if enabled:
                if rxn_id in reaction_map:
                    name, equation = reaction_map[rxn_id]
                    print '{}\t{}\t{}'.format(rxn_id, name, equation)
                elif rxn_id.startswith('rxntp_'):
                    cpdid = rxn_id[6:14]
                    formula, cpdname = compound_map[cpdid]
                    print '{}\t{}\t{}'.format(rxn_id, 'Transport of ' + cpdname, '|{}[e]| => |{}|'.format(cpdname, cpdname))
                elif rxn_id.startswith('rxnex_'):
                    cpdid = rxn_id[6:14]
                    formula, cpdname = compound_map[cpdid]
                    print '{}\t{}\t{}'.format(rxn_id, 'Exchange of ' + cpdname, '|{}[e]| <=>'.format(cpdname, cpdname))
            
    gffile.close()

