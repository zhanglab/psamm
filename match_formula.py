#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Match blocked metabolites with available isomers')
    parser.add_argument('cpdfile', type=argparse.FileType('r'), help='Compound table file')
    args = parser.parse_args()

    # Read formulas of compound list
    compound_map = {}
    for row in csv.reader(args.cpdfile, dialect='excel'):
        seed_cid, formula, mass, kegg_cid, cpd_names = row[:5]
        main_name = cpd_names.split(',<br>')[0].replace('&#39;', '\'')
        compound_map[seed_cid] = main_name, formula

    # Read list of cytosolic compounds
    compound_c = set()
    with open('cytosol_metabolites.txt', 'r') as f:
        for line in f:
            cpd = line.strip()
            compound_c.add(cpd)

    # Read list of blocked compounds
    blocked = set()
    with open('blocked.txt', 'r') as f:
        for line in f:
            cpd = line.strip()
            blocked.add(cpd)

    # Record formula of available compounds
    available_formula = defaultdict(set)
    for cpd in compound_c - blocked:
        cpd_name, formula = compound_map[cpd]
        available_formula[formula].add(cpd)

    # Find isomers of blocked compounds
    count = 0
    for cpd in sorted(blocked):
        cpd_name, formula = compound_map[cpd]
        other = available_formula[formula]
        if len(other) > 0:
            count += 1
            print '{} ({}) [{}] possibly available as isomer(s):'.format(cpd_name, cpd, formula)
            for other_id in other:
                print ' - {} ({})'.format(compound_map[other_id][0], other_id)

    print '{} blocked compounds can be resolved'.format(count)