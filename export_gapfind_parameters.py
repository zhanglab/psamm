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
    parser.add_argument('--database', nargs=2, metavar=('rxnfile', 'cpdfile'),
                        action='append', type=argparse.FileType('r'),
                        help='Supply reaction table to add to database list')
    args = parser.parse_args()

    rxn_table = args.rxnfile

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
    reaction_model = set()

    rxn_table.readline() # Skip header
    for row in csv.reader(rxn_table, dialect='excel'):
        SEED_rid, RXN_name, EC, Equation_cpdname, Equation_cpdid, KEGG_rid, KEGG_maps, Gene_ids = row[:8]
        rx = reaction.ModelSEED.parse(Equation_cpdid).normalized()

        # Lists all the reaction names
        w.write('{}\n'.format(SEED_rid))
        model_list.write('{}\n'.format(SEED_rid))
        reaction_model.add(SEED_rid)

        # Lists the reverse reactions
        if rx.direction == '<=>':
            rr.write('{}\n'.format(SEED_rid))
            for cpdid, value, comp in rx.left + rx.right:
                id = cpdid if comp is None else cpdid + '_' + comp
                compound_produced.add(id)
        else:
            for cpdid, value, comp in rx.right:
                id = cpdid if comp is None else cpdid + '_' + comp
                compound_produced.add(id)

        # Add compound names to the set
        for cpdid, value, comp in rx.left + rx.right:
            id = cpdid if comp is None else cpdid + '_' + comp
            compound.add(id)

            # Lists the compounds in two seprate files by compartartment
            if comp is None:
                compound_c.add(id)
            elif comp == 'e':
                compound_e.add(id)

        # Lists the matrix
        for cpdid, value, comp in rx.left:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, SEED_rid, -value))

        for cpdid, value, comp in rx.right:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, SEED_rid, value))

    # Write out list of compounds in the model
    for cpdid in sorted(compound):
        model_cpds.write('{}\n'.format(cpdid))

    # Optionally create transport reactions in database
    if args.transport:
        for cpdid in sorted(compound_c):
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

    # Load database reaction tables
    if args.database:
        for db_rxn_file, db_cpd_file in args.database:
            compound_map = {}

            db_cpd_file.readline() # Skip header
            for row in csv.reader(db_cpd_file, delimiter='\t'):
                seed_cid, cpdnames, formula, mass, kegg_maps, kegg_cid = row[:6]
                for cpdname in cpdnames.split(',<br>'):
                    compound_map[cpdname] = seed_cid
            db_cpd_file.close()

            db_rxn_file.readline() # Skip header
            for row in csv.reader(db_rxn_file, delimiter='\t'):
                seed_rid, rxn_name, equation_cpdname, roles, subsystems, kegg_maps, enzyme, kegg_rid = row[:8]

                if seed_rid in reaction_model:
                    continue
                if equation_cpdname.strip() == '':
                    continue

                def translate(name):
                    m = re.match(r'cdp(\d+)', name) # [sic]
                    if m is not None:
                        return 'cpd' + m.group(1)
                    return compound_map[name]

                rx = reaction.ModelSEED.parse(equation_cpdname).normalized().translated_compounds(translate)

                # Lists all the reaction names
                w.write('{}\n'.format(seed_rid))
                database_list.write('{}\n'.format(seed_rid))

                # Lists the reverse reactions
                if rx.direction == '<=>' or rx.direction.strip() == '':
                    rr.write('{}\n'.format(seed_rid))

                # Add compound names to the set
                for cpdid, value, comp in rx.left + rx.right:
                    id = cpdid if comp is None else cpdid + '_' + comp
                    compound.add(id)

                    # Lists the compounds in two seprate files by compartartment
                    if comp is None:
                        compound_c.add(id)
                    elif comp == 'e':
                        compound_e.add(id)

                # Lists the matrix
                for cpdid, value, comp in rx.left:
                    id = cpdid if comp is None else cpdid + '_' + comp
                    m.write('{}.{}\t{}\n'.format(id, seed_rid, -value))

                for cpdid, value, comp in rx.right:
                    id = cpdid if comp is None else cpdid + '_' + comp
                    m.write('{}.{}\t{}\n'.format(id, seed_rid, value))

            db_rxn_file.close()

    # Lists all the compound names in the set
    for cpdid in sorted(compound):
        cl.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_c):
        cm.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_e):
        exm.write('{}\n'.format(cpdid))

    compound_not_produced = compound_c - compound_produced
    for cpdid in sorted(compound_not_produced):
        rnp.write('{}\n'.format(cpdid))

    rxn_table.close()

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