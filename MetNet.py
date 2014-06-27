#!/usr/bin/env python

import sys
import csv
from decimal import Decimal
import reaction

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'USAGE: {} RXNFILE CPDFILE'.format(sys.argv[0])
        sys.exit(-1)
        
    # c = compound table file
    c = open(sys.argv[2],'r')
    c.readline()
    
    ## reader = compound table reader
    reader = csv.reader(c, dialect='excel')
    for row in reader:
        Metabolite_abbreviation, Name, Formula, Charge = row
        
    # Close compound table    
    c.close()

    r = open(sys.argv[1],'r')

    w = open('rxnnames.txt', 'w')
    model_list = open('modelrxn.txt', 'w')
    rr = open('rev.txt', 'w')
    cl = open('metabolites.txt', 'w')
    m = open('mat.txt', 'w')
    cm = open('cytosol_metabolites.txt', 'w')
    exm = open('extracellular_metabolites.txt', 'w')
    model_cpds = open('model_cpds.txt', 'w')

    compound = set()
    compound_e = set()
    compound_c = set()
    
    r.readline()
    reader1 = csv.reader(r, dialect='excel')
    for row1 in reader1:
        Reaction_abbreviation, Reaction_name, Equation, EC_number, Gene, Protein, H_Pylori_ortholog, Pathway, Source, Notes= row1
        rxnid = 'rxn_' + Reaction_abbreviation.replace('(', '__').replace(')', '__')

        rx = reaction.MetNet.parse(Equation).normalized()

        # Lists all the reaction names
        w.write('{}\n'.format(rxnid))
        model_list.write('{}\n'.format(rxnid))

        # Lists the reverse reactions
        if rx.direction == '<=>':
            rr.write('{}\n'.format(rxnid))


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
            m.write('{}.{}\t{}\n'.format(id, rxnid, -value))

        for cpdid, value, comp in rx.right:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, rxnid, value))

    # Lists all the compound names in the set
    for cpdid in sorted(compound):
        cl.write('{}\n'.format(cpdid))
        model_cpds.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_c):
        cm.write('{}\n'.format(cpdid))

    for cpdid in sorted(compound_e):
        exm.write('{}\n'.format(cpdid))
    
    r.close()
    rr.close()
    w.close()
    cl.close()
    m.close()
    cm.close()
    exm.close() 

   
