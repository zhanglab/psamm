#!/usr/bin/env python

'''Parse GapFind output and print no-production metabolites

The GapFind file should be produced with parameter "ps=0".
'''

import csv
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse GapFind output and generate no-production list')
    parser.add_argument('gapfind', type=argparse.FileType('r'), help='GapFind output file (ps=0)')
    parser.add_argument('--cpdfile', type=argparse.FileType('r'), help='Compound table file')
    args = parser.parse_args()

    blocked_file = open('blocked.txt', 'w')

    model_cpds = set()
    with open('model_cpds.txt', 'r') as f:
        for line in f:
            cpdid = line.strip()
            model_cpds.add(cpdid)

    compound_e = set()
    with open('extracellular_metabolites.txt', 'r') as f:
        for line in f:
            cpdid = line.strip()
            compound_e.add(cpdid)

    # Load compounds
    compound_map = {}
    if args.cpdfile is not None:
        cpdfile = args.cpdfile
        cpdfile.readline() # Skip header
        for row in csv.reader(cpdfile, dialect='excel'):
            seed_cid, formula, mass, kegg_cid, cpd_name = row[:5]
            compound_map[seed_cid] = formula, kegg_cid, cpd_name
        cpdfile.close()

    # Parse GapFind output
    count = 0
    count_not_produced = 0
    count_in_kegg = 0
    flag = False
    gffile = args.gapfind
    for line in gffile:
        if line.startswith('---- VAR xp'):
            flag = True
            for i in range(3):
                next(gffile)
        elif flag:
            if line == '\n':
                break

            fields = line.split()
            cpdid = fields[0]
            produced = fields[2] != '.'
            comp = None

            if cpdid not in compound_e and cpdid in model_cpds:
                count += 1
                if not produced:
                    blocked_file.write('{}\n'.format(cpdid))
                    count_not_produced += 1

                    if cpdid in compound_map:
                        Formula, KEGG_cid, CPD_name = compound_map[cpdid]
                        if KEGG_cid != 'None':
                            count_in_kegg += 1
                        print '{}\t{}\t{}\t{}'.format(cpdid, KEGG_cid, Formula, CPD_name)
                    else:
                        print cpdid
    gffile.close()
    blocked_file.close()

    print 'Not produced: {}/{}'.format(count_not_produced, count)
    print 'In KEGG: {}/{}'.format(count_in_kegg, count_not_produced)
