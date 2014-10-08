#!/usr/bin/env python

'''Convert SBML model file to standard format'''

import sys
import argparse

from metnet import sbml

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert SBML model file to standard format')
    parser.add_argument('sbmlfile', type=argparse.FileType('r'), help='SBML file')
    parser.add_argument('--reactions', dest='reactions', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Reaction table output file')
    parser.add_argument('--compounds', dest='compounds', nargs='?', type=argparse.FileType('w'),
                        help='Compound table output file')
    args = parser.parse_args()

    database = sbml.SBMLDatabase(args.sbmlfile)

    # Write reaction table
    for reaction in database.reactions:
        rx = database.get_reaction(reaction).normalized()
        args.reactions.write('{}\t{}\n'.format(reaction, rx))

    # Write compound table
    if args.compounds:
        args.compounds.write('id\tname\n')
        for compound in database.compounds:
            name = database.get_compound_name(compound)
            name = name if name is not None else compound.name
            args.compounds.write('{}\t{}\n'.format(compound.name, name))