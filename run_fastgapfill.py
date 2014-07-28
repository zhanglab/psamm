#!/usr/bin/env python

from metabolicmodel import MetabolicDatabase
import fastcore
import fluxanalysis

if __name__ == '__main__':
    database = MetabolicDatabase.load_from_file()
    model = database.load_model_from_file()

    # Run Fastcc and print the inconsistent set
    print 'Calculating Fastcc consistent subset...'
    consistent_core = fastcore.fastcc(model, 0.001)
    print 'Result: |A| = {}, |A| = {}'.format(len(consistent_core), consistent_core)

    # Run Fastcore and print the induced reaction set
    model_complete = model.copy()
    for rxnid in database.reactions:
        model_complete.add_reaction(rxnid)
    print 'Calculating Fastcore induced set with consistent core...'
    core = consistent_core | { 'Biomass' }
    induced = fastcore.fastcore(model_complete, core, 0.001)
    print 'Result: |A| = {}, A = {}'.format(len(induced), induced)
    added_reactions = induced - core
    print 'Extended: |E| = {}, E = {}'.format(len(added_reactions), added_reactions)

    # Load bounds on exchange reactions
    model.load_exchange_limits()

    print 'Flux balance on original model maximizing Biomass...'
    for rxnid, flux in sorted(fluxanalysis.flux_balance(model, 'Biomass')):
        print '{}\t{}'.format(rxnid, flux)

    print 'Flux balance on induced model maximizing Biomass...'
    for rxnid in induced:
        model.add_reaction(rxnid)
    for rxnid, flux in sorted(fluxanalysis.flux_balance(model, 'Biomass')):
        print '{}\t{}'.format(rxnid, flux)