#!/usr/bin/env python

'''List the flux values of the FluxBalance GAMS output'''

import csv
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This program lists the FLux values of the Fluxbalance.lst')
    parser.add_argument('fluxbalance', type=argparse.FileType('r'), help='Fluxbalance file (ps=0)')
    parser.add_argument('rxnlist', nargs='?', type=argparse.FileType('r'), help='Model reactions')
    args = parser.parse_args()

    fbfile = args.fluxbalance

    # Load list of model reactions
    model_reactions = None
    if args.rxnlist:
        model_reactions = set()
        for line in args.rxnlist:
            rxnid = line.strip()
            model_reactions.add(rxnid)

    # Parse GAMS output file
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
            flux = float(fields[2]) if enabled else 0
            if model_reactions is None or rxnid in model_reactions:
                print '{}\t{}'.format(rxnid, flux)