#!/usr/bin/env python

'''This program lists the Flux values of the Fluxbalance.lst that 
have a value that is close to either 1000 or -10000.
'''

import csv
import re
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This program lists the FLux values of the Fluxbalance.lst')
    parser.add_argument('Fluxbalance', type=argparse.FileType('r'), help='Fluxbalance file (ps=0)')
    args = parser.parse_args()

    fbfile = args.Fluxbalance
    
    value_1000 = open('1000_fluxvalue.txt', 'w')
    
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
            
            if flux == 1000.000:
                value_1000.write('{}\t{}\n'.format(rxnid, flux))
            elif flux == -1000.000:
                value_1000.write('{}\t{}\n'.format(rxnid, flux))
            elif flux > 800:
                value_1000.write('{}\t{}\n'.format(rxnid, flux))
            elif flux < -800:
                 value_1000.write('{}\t{}\n'.format(rxnid, flux))

    value_1000.close()