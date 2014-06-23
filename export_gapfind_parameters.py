#!/usr/bin/env python

'''This program lists (in separate files) the reaction names, the names of the 
reactions that are reversible, the cpdid of the compounds used in all of the 
reactions, writes the compounds and stoichiometric values in matrix format'''

import argparse
import csv
import reaction

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert reaction table to GapFind input format')
    parser.add_argument('rxnfile', type=argparse.FileType('r'), help='Reaction table file')
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

    compound = set()
    compound_e = set()
    compound_c = set()
    compound_produced = set() 
    
    rxn_table.readline() # Skip header
    for row in csv.reader(rxn_table, dialect='excel'):
        SEED_rid, RXN_name, EC, Equation_cpdname, Equation_cpdid, KEGG_rid, KEGG_maps, Gene_ids = row[:8]

        direction, left, right = reaction.normalize(reaction.parse(Equation_cpdid))

        # Lists all the reaction names
        w.write('{}\n'.format(SEED_rid))
        model_list.write('{}\n'.format(SEED_rid))

        # Lists the reverse reactions
        if direction == '<=>':
            rr.write('{}\n'.format(SEED_rid))
            for cpdid, value, comp in left + right:
                id = cpdid if comp is None else cpdid + '_' + comp
                compound_produced.add(id)
        else:
            for cpdid, value, comp in right:
                id = cpdid if comp is None else cpdid + '_' + comp
                compound_produced.add(id)
                
        # Add compound names to the set
        for cpdid, value, comp in left + right:
            id = cpdid if comp is None else cpdid + '_' + comp
            compound.add(id)

            # Lists the compounds in two seprate files by compartartment
            if comp is None:
                compound_c.add(id)
            elif comp == 'e':
                compound_e.add(id)

        # Lists the matrix
        for cpdid, value, comp in left:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, SEED_rid, -value))

        for cpdid, value, comp in right:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, SEED_rid, value))

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
