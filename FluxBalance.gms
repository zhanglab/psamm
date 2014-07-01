*** FluxBalance
* This program will run Flux Balance Analysis on a given model, trying to
* optimize the output of a certain compound or pseudocompound (e.g biomass)
* under the constraints given by the reactions in the model.
*
* Based on FBA_CoreTextbookModel.gms by
* Joshua J. Hamilton, Joohnoon Kim, and Jennifer L. Reed
* Department of Chemical and Biological Engineering
* University of Wisconsin-Madison, Madison, Wisconsin, USA
* http://reedlab.che.wisc.edu/
*
* Computational Tools for Analyzing Microbial Metabolism Workshop
* American Society of Microbiology 2013 Annual Meeting
* Denver, Colorado, USA
* May 18, 2013
* This workshop was supported in part by an NSF CAREER Award to J. Reed
* (NSF 1053712) and an NSF Graduate Research Fellowship to J. Hamilton
* (DGE-0718123).


*** Read sets from external files
Set i metabolites in S (m) /
$include "metabolites.txt"
/;

Set j reactions in S (n) /
$include "rxnnames.txt"
/;

Set mod(j) reactions of the model organism /
$include "modelrxn.txt"
/;

Set rev(j) reversible reactions /
$include "rev.txt"
/;

Set exch(j) exchange reactions /
$include "exchangerxn.txt"
/;

*** Write solver options to file
option lp = cplex;
$onecho > cplex.opt
eprhs 1E-9
epopt 1E-9
epint 0
numericalemphasis 1
$offecho

Parameters
        UpperLimits(j)          maximum value flux can take
        LowerLimits(j)          minimum value flux can take
        c(j)                    indicates the reaction which is to be maximized
        Vmax    /1000/
;

Parameter S(i,j) stoichiometric coefficients (S matrix) /
$include "mat.txt"
/;

*** Set lower and upper bounds on fluxes accordingly
* By default, all fluxes are allowed to vary between Vmax and -Vmax
* The S matrix contains the set of reversible directions. Reactions outside this
* set are required to have a lower bound of 0
* We also set the lower bound of all exchange reactions to zero and specify
* the media composition below.
UpperLimits(j) = Vmax;
LowerLimits(j) = 0;
LowerLimits(j)$rev(j) = -Vmax;
LowerLimits(j)$exch(j) = 0;

Variables
        v(j)                    flux value through reaction
        z                       objective
;

Equations
        massbalance(i)			mass balance equations for each metabolite
        obj						calculates the dot product of the c vector and the flux vector
;

massbalance(i).. sum( j,S(i,j)*v(j) ) =e= 0;
obj.. z =e= sum( j,c(j)*v(j) );

*** Define the reaction to maximize
c('Biomass') = 1;


Model FBA /
        obj
        massbalance
/;

*** Define the media
* Place limits on the exchange fluxes based on the minimal media.
* By default the upper limits of the exchange fluxes are all set to Vmax,
* indicating that the extracellular metabolites can flow away from the cell.
* A negative flux through the exchange reactions implies that the metabolite
* are being made available to the cell.

*** Carbon sources
*LowerLimits('EX_ac_e') = -5;
*LowerLimits('EX_akg_e') = -5;
*LowerLimits('EX_etoh_e') = -5;
*LowerLimits('EX_for_e') = -5;
*LowerLimits('EX_fum_e') = -5;
LowerLimits('EX_glc_e') = -5;
*LowerLimits('EX_lacD_e') = -5;
*LowerLimits('EX_pyr_e') = -5;
*LowerLimits('EX_succ_e') = -5;

*** Electron Acceptor
LowerLimits('EX_o2_e') = -Vmax;

*** Essential nutrients
LowerLimits('EX_co2_e') = -Vmax;
LowerLimits('EX_h2o_e') = -Vmax;
LowerLimits('EX_h_e') = -Vmax;
LowerLimits('EX_pi_e') = -Vmax;

*** Use the LowerLimits and UpperLimits parameters to establish global bounds
* on the fluxes. Reactions that are not in the model are fixed to 0.
v.lo(j) = LowerLimits(j);
v.up(j) = UpperLimits(j);
v.fx(j)$( not mod(j) ) = 0;

*** Solve model
FBA.optfile = 1;
solve FBA using lp maximizing z;
