#!/usr/bin/env python

'''Enter a list of reaction ID's and the program gives you the equation_name and equation_cpd'''

import csv
import argparse
import re

from metnet.reaction import ModelSEED

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prints Equation_cpd')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('cpdfile', type=argparse.FileType('r'), help='Compound table file')
    parser.add_argument('rxnlist', type=argparse.FileType('r'), help='List of Rxns')

    args = parser.parse_args()

    rxn_file = args.rxnfile
    cpd_file = args.cpdfile
    rxn_list = args.rxnlist

    compound_map = {}
    reaction_map = set()

    for row in csv.reader(rxn_list, delimiter='\t'):
        rxn_id = row[0]

        reaction_map.add(rxn_id)

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

        rx = ModelSEED.parse(equation_cpdname).normalized().translated_compounds(translate)

        if seed_rid in reaction_map:
            print '{}\t{}\t{}'.format(seed_rid,equation_cpdname,rx)


    rxn_file.close()