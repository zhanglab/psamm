*** LoopFluxReduction
* This program will run a progam on a given model that is similar to
* FBA, but it will keep one reaction fixed at a certain flux. At the
* same time, the L1-norm of all fluxes is minimized. This is useful
* for detecting reactions that are not necessary but still have flux
* under normal FBA.

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

Parameter exchlimit(j) lower bound of exchange reactions (max uptake) /
$include "exchangelimit.txt";
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

Positive variables
	t(j)		        bound on fluxes
;

Equations
        massbalance(i)		mass balance equations for each metabolite
	fixedflux	        fix flux for certain reaction
	boundedflux_lower(j)	bound fluxes in the t interval
	boundedflux_upper(j)	bound fluxes in the t interval
        obj			L1 norm of fluxes
;

massbalance(i).. sum( j,S(i,j)*v(j) ) =e= 0;
fixedflux.. v('Biomass') =g= 0.4897;
boundedflux_lower(j).. v(j) =g= -t(j);
boundedflux_upper(j).. v(j) =l= t(j);
obj.. z =e= sum( j,t(j) );


Model FBA /
        massbalance
	fixedflux
	boundedflux_lower
	boundedflux_upper
        obj
/;

*** Define the media
* Place limits on the exchange fluxes based on the minimal media.
* By default the upper limits of the exchange fluxes are all set to Vmax,
* indicating that the extracellular metabolites can flow away from the cell.
* A negative flux through the exchange reactions implies that the metabolite
* are being made available to the cell.
LowerLimits(j)$exch(j) = exchlimit(j);

*** Use the LowerLimits and UpperLimits parameters to establish global bounds
* on the fluxes. Reactions that are not in the model are fixed to 0.
v.lo(j) = LowerLimits(j);
v.up(j) = UpperLimits(j);
v.fx(j)$( not mod(j) ) = 0;

*** Solve model
FBA.optfile = 1;
solve FBA using lp minimize z;
