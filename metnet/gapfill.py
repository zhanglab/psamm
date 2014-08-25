
'''Identify blocked metabolites and possible reconstructions

This implements a variant of the algorithms described by Kumar,
Vinay Satish, Madhukar S. Dasika, and Costas D. Maranas.
"Optimization based automated curation of metabolic reconstructions."'''

import cplex

from . import lpsolver

def cpdid_str(compound):
    cpdid, comp = compound
    if comp is None:
        return cpdid
    return cpdid+'_'+comp

def gapfind(model, epsilon=1e-5, v_max=1000, solver=lpsolver.CplexSolver()):
    pass

def gapfill(model, core, blocked, epsilon=1e-5, v_max=1000, solver=lpsolver.CplexSolver()):
    pass
