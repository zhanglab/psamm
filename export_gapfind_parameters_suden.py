#!/usr/bin/env python

'''This program lists (in separate files) the reaction names, the names of the 
reactions that are reversible, the cpdid of the compounds used in all of the 
reactions, writes the compounds and stoichiometric values in matrix format'''

import sys
import csv
#import reaction
from decimal import Decimal

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
    compound_e = set()
    compound_c = set()
    compound_produced = set() 
    
    r.readline()
    readerr = csv.reader(r, delimiter='\t')
    for rowr in readerr:
        rxn_id, rxn_name, equation_cpd, equation_name, Subsystem, Suden_protien_ID, gene_association, protein_id, ec_reference, Reversibility, Ref_directionality, category, SEED_reversibility, Reversibility_without_seed, Pfam_domains, Pfam_Consensus_Score, DIFF_list = rowr[:17]
        
        
        left = []
        right = []

        def parse_compound_list(s):
            cpds = []
            for cpd in s.split('+'):
                if cpd == '':
                    continue
                count, spec = cpd.strip().split(' ')
                spec_split = spec.split('[')
                if len(spec_split) == 2:
                    # compartment
                    cpdid = spec_split[0]
                    comp = spec_split[1].rstrip(']')
                else:
                    cpdid = spec
                    comp = None
                d = Decimal(count)
                if d % 1 == 0:
                    d = int(d)
                cpds.append((cpdid, d, comp))
            return cpds

        cpd_left, cpd_right = equation_cpd.split('<=>')
        left = parse_compound_list(cpd_left)
        right = parse_compound_list(cpd_right)

        direction = '<=>'
        if Reversibility == 'N':
            direction = '=>'

        #direction, left, right = reaction.normalize(reaction.parse(Equation_cpdid))
    
        # Lists all the reaction names
        w.write('{}\n'.format(rxn_id))

        # Lists the reverse reactions
        if direction == '<=>':
            rr.write('{}\n'.format(rxn_id))
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
            m.write('{}.{}\t{}\n'.format(id, rxn_id, -value))
        
        for cpdid, value, comp in right:
            id = cpdid if comp is None else cpdid + '_' + comp
            m.write('{}.{}\t{}\n'.format(id, rxn_id, value))

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
 
    
 
 
    r.close()
    rr.close()
    w.close()
    cl.close()
    m.close()
    cm.close()
    exm.close()
