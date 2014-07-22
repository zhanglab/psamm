#!/usr/bin/env python

import argparse
from metabolicmodel import MetabolicDatabase
import fluxanalysis
import random

FLUX_THRESHOLD = 0.75

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find random minimal set of reactions that keep biomass above a threshold')
    parser.add_argument('reaction', help='Name of the biomass reaction')
    args = parser.parse_args()

    biomass_reaction = args.reaction

    database = MetabolicDatabase.load_from_file()
    model = database.load_model_from_file()
    model.load_exchange_limits()

    # Obtain normal biomass flux
    fluxes = dict(fluxanalysis.flux_balance(model, biomass_reaction))
    threshold = fluxes[biomass_reaction] * FLUX_THRESHOLD
    
    # number of experiments
    model_test = model.copy()
    essential = { biomass_reaction }
    deleted = set()
    test_set = model_test.reaction_set - essential

    while len(test_set) > 0:
        testing_rxn = random.sample(test_set, 1)[0]
        test_set.remove(testing_rxn)
        saved_bounds = model_test.limits[testing_rxn].bounds
        model_test.limits[testing_rxn].bounds = 0, 0

        try:
            fluxes = dict(fluxanalysis.flux_balance(model, biomass_reaction))
        except:
            for rxnid, bounds in model.limits.iteritems():
                print rxnid, bounds
            raise

        if fluxes[biomass_reaction] < threshold:
            model_test.limits[testing_rxn].bounds = saved_bounds
            essential.add(testing_rxn)
        else:
            deleted.add(testing_rxn)

    for rxnid in model_test.reaction_set:
        value = 0 if rxnid in deleted else 1
        print '{}\t{}'.format(rxnid, value)