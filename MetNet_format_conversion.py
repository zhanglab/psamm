#!/usr/bin/env python

'''Convert MetNet reaction table into standard tab-separated reaction table'''

import argparse
import csv
import reaction

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert MetNet reaction table to standardized format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    args = parser.parse_args()

    rxnfile = args.rxnfile
    rxnfile.readline()
    for row in csv.reader(rxnfile, dialect='excel'):
        Reaction_abbreviation, Reaction_name, Equation, EC_number, Gene, Protein, H_Pylori_ortholog, Pathway, Source, Notes = row
        rxnid = 'rxn_' + Reaction_abbreviation.replace('(', '__').replace(')', '__')

        rx = reaction.MetNet.parse(Equation).normalized()

        print '{}\t{}'.format(rxnid, rx)

    rxnfile.close()