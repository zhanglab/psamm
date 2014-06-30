#!/usr/bin/env python

'''Convert ModelSEED reaction table to standard format'''

import csv
import reaction
import argparse
import re

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert ModelSEED reaction table to standard format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('cpdfile', type=argparse.FileType('r'), help='Compound table file')
    
    args = parser.parse_args()

    rxn_file = args.rxnfile
    cpd_file = args.cpdfile
    
    compound_map = {}

    cpd_file.readline() # Skip header
    for row in csv.reader(cpd_file, delimiter='\t'):
        seed_cid, cpdnames, formula, mass, kegg_maps, kegg_cid = row[:6]
        for cpdname in cpdnames.split(',<br>'):
            compound_map[cpdname] = seed_cid
    cpd_file.close()

    rxn_file.readline() # Skip header
    for row in csv.reader(rxn_file, delimiter='\t'):
        seed_rid, rxn_name, equation_cpdname, roles, subsystems, kegg_maps, enzyme, kegg_rid = row[:8]

        if equation_cpdname.strip() == '':
            continue

        def translate(name):
            m = re.match(r'cdp(\d+)', name) # [sic]
            if m is not None:
                return 'cpd' + m.group(1)
            return compound_map[name]

        rx = reaction.ModelSEED.parse(equation_cpdname).normalized().translated_compounds(translate)

        print '{}\t{}'.format(seed_rid,rx)

    rxn_file.close()