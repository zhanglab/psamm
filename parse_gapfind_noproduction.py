#!/usr/bin/env python

'''Parse GapFind output and print no-production metabolites

The GapFind file should be produced with parameter "ps=0".
'''

import sys
import csv
import re

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'USAGE: {} GAPFIND CPDLIST'.format(sys.argv[0])
        sys.exit(-1)

    f = open(sys.argv[1], 'r')
    c = open(sys.argv[2], 'r')
    blocked_file = open('blocked.txt', 'w')

    compound_map = {}

    c.readline() # Skip header
    readerc = csv.reader(c, dialect='excel')
    for rowc in readerc:
        SEED_cid, Formula, Mass, KEGG_cid, CPD_name = rowc[:5]
        compound_map[SEED_cid] = Formula, KEGG_cid, CPD_name

    count = 0
    count_not_produced = 0
    count_in_kegg = 0
    flag = False
    for line in f:
        if line.startswith('---- VAR xp'):
            flag = True
            for i in range(3):
                next(f)
        elif flag:
            if line == '\n':
                break

            fields = line.split()
            cpdid = fields[0]
            produced = fields[2] != '.'
            comp = None

            m = re.match(r'(.+?)_(.+)$', cpdid)
            if m is not None:
                cpdid = m.group(1)
                comp = m.group(2)

            if comp is None and cpdid in compound_map:
                count += 1
                if not produced:
                    Formula, KEGG_cid, CPD_name = compound_map[cpdid]
                    count_not_produced += 1
                    if KEGG_cid != 'None':
                        count_in_kegg += 1
                    print '{}\t{}\t{}\t{}'.format(cpdid, KEGG_cid, Formula, CPD_name)
                    blocked_file.write('{}\n'.format(cpdid))

    print 'Not produced: {}/{}'.format(count_not_produced, count)
    print 'In KEGG: {}/{}'.format(count_in_kegg, count_not_produced)

    blocked_file.close()