#!/usr/bin/env python

'''This program lists (in separate files) the reaction names, the names of the
reactions that are reversible, the cpdid of the compounds used in all of the
reactions, writes the compounds and stoichiometric values in matrix format'''

import argparse
import csv
import reaction
import re

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert reaction table to GapFind input format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
    parser.add_argument('--exchange', action='store_true', help='Add exchange reactions to database list')
    parser.add_argument('--transport', action='store_true', help='Add transport reactions to database list')
    parser.add_argument('--database', metavar='rxnfile', action='append', type=argparse.FileType('r'),
                        help='Supply reaction table to add to database list')
    args = parser.parse_args()

    rxn_file = args.rxnfile

    # Opens files to write in
    w = open('rxnnames.txt', 'w')
    database_list = open('databaserxn.txt', 'w')
    model_list = open('modelrxn.txt', 'w')
    rr = open('rev.txt', 'w')
    cl = open('metabolites.txt', 'w')
    m = open('mat.txt', 'w')
    cm = open('cytosol_metabolites.txt', 'w')
    exm = open('extracellular_metabolites.txt', 'w')
    rnp = open('root_no_production.txt', 'w')
    model_cpds = open('model_cpds.txt', 'w')

    compound = set()
    compound_e = set()
    compound_c = set()
    compound_produced = set()
    compound_model = set()
    model_reactions = set()

    for line in rxn_file:
        rxn_id = line.strip()

        # Lists all the reaction names
        model_list.write('{}\n'.format(rxn_id))
        model_reactions.add(rxn_id)

    # Load database reaction tables
    for db_rxn_file in args.database:
        for row in csv.reader(db_rxn_file, delimiter='\t'):
            rxn_id, equation = row[:2]
            rx = reaction.ModelSEED.parse(equation).normalized()

            # Lists all the reaction names
            w.write('{}\n'.format(rxn_id))
            database_list.write('{}\n'.format(rxn_id))

            # Lists the reverse reactions
            if rx.direction == '<=>' or rx.direction.strip() == '' or rx.direction == '?':
                rr.write('{}\n'.format(rxn_id))

            # Add compound names to the set
            for cpdid, value, comp in rx.left + rx.right:
                id = cpdid if comp is None else cpdid + '_' + comp
                compound.add(id)

                if rxn_id in model_reactions:
                    compound_model.add(id)

                # Lists the compounds in two seprate files by compartartment
                if comp is None:
                    compound_c.add(id)
                elif comp == 'e':
                    compound_e.add(id)

            # Lists the matrix
            for cpdid, value, comp in rx.left:
                id = cpdid if comp is None else cpdid + '_' + comp
                m.write('{}.{}\t{}\n'.format(id, rxn_id, -value))

            for cpdid, value, comp in rx.right:
                id = cpdid if comp is None else cpdid + '_' + comp
                m.write('{}.{}\t{}\n'.format(id, rxn_id, value))

        db_rxn_file.close()

    # Optionally create transport reactions in database
    if args.transport:
        for cpdid in sorted(compound_c):
            if cpdid in compound_model:
                rxnid = 'rxntp_' + cpdid
                w.write('{}\n'.format(rxnid))
                database_list.write('{}\n'.format(rxnid))

                # Write to matrix
                compound.add(cpdid + '_e')
                compound_e.add(cpdid + '_e')
                m.write('{}_e.{}\t{}\n'.format(cpdid, rxnid, -1))
                m.write('{}.{}\t{}\n'.format(cpdid, rxnid, 1))

    # Optionally create exchange reactions in database
    if args.exchange:
        for cpdid in sorted(compound_e):
            rxnid = 'rxnex_' + cpdid
            w.write('{}\n'.format(rxnid))
            database_list.write('{}\n'.format(rxnid))

            # Write to matrix
            m.write('{}.{}\t{}\n'.format(cpdid, rxnid, -1))
            rr.write('{}\n'.format(rxnid))

    # Lists all the compound names in the set
    for cpdid in sorted(compound):
        cl.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_c):
        cm.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_e):
        exm.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_model):
        model_cpds.write('{}\n'.format(cpdid))

    compound_not_produced = compound_c - compound_produced
    for cpdid in sorted(compound_not_produced):
        rnp.write('{}\n'.format(cpdid))

    rxn_file.close()

    database_list.close()
    model_list.close()

    rr.close()
    w.close()
    cl.close()
    m.close()
    cm.close()
    exm.close()
    rnp.close()

    model_cpds.close()