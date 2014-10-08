#!/usr/bin/env python

'''Convert SBML model file to standard format'''

import csv
import argparse
import re

from metnet import sbml

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert SBML model file to standard format')
    parser.add_argument('sbmlfile', type=argparse.FileType('r'), help='SBML file')
    args = parser.parse_args()

    database = sbml.SBMLDatabase(args.sbmlfile)

    for reaction in database.reactions:
        rx = database.get_reaction(reaction).normalized()
        print '{}\t{}'.format(reaction, rx)