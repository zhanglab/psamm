#!/usr/bin/env python

'''This program lists the FLux values of the Fluxbalance.lst.
'''

import csv
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This program lists the FLux values of the Fluxbalance.lst')
    parser.add_argument('Fluxbalance', type=argparse.FileType('r'), help='Fluxbalance file (ps=0)')
    args = parser.parse_args()

    fbfile = args.Fluxbalance
    
    flag = False 

    for line in fbfile:
        if line.startswith('---- VAR v'):
            flag = True
            for i in range(3):
                next(fbfile)
        elif flag:
            if line == '\n':
                break

            fields = line.split()
            rxnid = fields[0]
            enabled = fields[2] != '.'
            if enabled:
                print '{}\t{}'.format(rxnid, float(fields[2]))