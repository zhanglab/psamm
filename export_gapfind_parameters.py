
'''This program lists (in separate files) the reaction names, the names of the 
reactions that are reversible, the cpdid of the compounds used in all of the 
reactions, writes the compounds and stoichiometric values in matrix format'''

import sys
import csv
import reaction

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'USAGE: {} RXNFILE'.format(sys.argv[0])
        sys.exit(-1)

    r = open(sys.argv[1], 'r')

    # Opens files to write in 
    w = open('rxnnames.txt', 'w')
    rr = open('rev.txt', 'w')
    cl = open('metabolites.txt', 'w')
    m = open('mat.txt', 'w')
    cm = open('cytosol_metabolites.txt', 'w')
    exm = open('extracellular_metabolites.txt', 'w')

    compound = set()

    r.readline()
    readerr = csv.reader(r, dialect='excel')
    for rowr in readerr:
        SEED_rid, RXN_name, EC, Equation_cpdname, Equation_cpdid, KEGG_rid, KEGG_maps, Gene_ids, Visited, Pathway, Additional_paths, Inconsistencies, Gene_id_and_rn_used, Actual_ec, Reaction_technicalities = rowr

        direction, left, right = reaction.parse(Equation_cpdid)

        # Lists all the reaction names
        w.write('{}\n'.format(SEED_rid))

        # Lists the reverse reactions
        if direction == '<=>':
            rr.write('{}\n'.format(SEED_rid))


        # Add compound names to the set
        for cpdid, value, comp in left + right:
            id = cpdid if comp is None else cpdid + '_' + comp
            compound.add(id)

            # Lists the compounds in two seprate files by compartartment
            if comp is None:
                cm.write('{}\n'.format(id))
            elif comp == 'e':
                exm.write('{}\n'.format(id))

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

    r.close()
    rr.close()
    w.close()
    cl.close()
    m.close()
    cm.close()
    exm.close()
