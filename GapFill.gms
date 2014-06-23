*** GapFill
* This program determines the minimal number of modifications required to restore
* connectivity of each blocked metabolite in the model. A full database of possible
* reactions must be provided. The binary variable ym will indicate a model reaction that
* needs to be reversed and yd will indicate a database reaction that needs to be added.
*
* This implements the GapFill model described by Kumar, Vinay Satish, Madhukar S. Dasika,
* and Costas D. Maranas. "Optimization based automated curation of metabolic reconstructions."
* BMC bioinformatics 8.1 (2007): 212.


*** Allow empty sets to be defined
$onempty

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

Set database(j)	reactions from the external database /
$include "databaserxn.txt"
/;

Set rev(j) reversible reactions /
$include "rev.txt"
/;

Set cytosol(i) cytosolic metabolites /
$include "cytosol_metabolites.txt"
/;

Set extracellular(i) extracellular metabolites /
$include "extracellular_metabolites.txt"
/;

Set blocked(i) metabolites to be resolved /
$include "blocked.txt"
/;

Parameters
        UpperLimits(j)          maximum value flux can take
        LowerLimits(j)          minimum value flux can take
        Vmax    /100/
        delta   /1/
        M       /1000/
;

Parameter S(i,j) stoichiometric coefficients (S matrix) /
$include "mat.txt"
/;

*** Set lower and upper bounds on fluxes accordingly
* Database reactions are assumed to be always reversible
UpperLimits(j)$( mod(j) ) = Vmax;
LowerLimits(j)$( mod(j) ) = -Vmax;
LowerLimits(j)$( mod(j) and not rev(j) ) = 0;
UpperLimits(j)$( database(j) ) = Vmax;
LowerLimits(j)$( database(j) ) = -Vmax;

Variables
        v(j)			flux value through reaction
	z			objective
;

Binary variables
	yd(j)			whether database reaction j is active
        ym(j)                   whether model reaction j is reversed
	w(i,j)                  whether reaction j that produces metabolite i is active
;

Equations
	massbalance(i)          mass balance for cytosolic metabolites
	massbalance_comp(i)     mass balance for metabolites in internal compartments
	modelfluxrev_min(j)     minimum flux for irreversible model reactions
	modelfluxirrev_min(j)   minimum flux for reversible model reactions
	modelflux_max(j)        maximum flux for model reactions
	databaseflux_min(j)     minimum flux for database reactions
	databaseflux_max(j)     maximum flux for database reactions
	prodcons_min(i,j)       minimum production binary constraint
	prodcons_max(i,j)       maximum production binary constraint
	binarycons(i)           metabolite production binary constraint
	obj                     objective function
;

*** The mass balance constraints
* There is an implicit metabolite sink for cytosolic metabolites because of
* the diluting effect from cell division, export from diffusion, etc., which means
* that the production of these metabolites only has to be positive. For metabolites
* in internal compartments the uptake and secretion has to even out.
massbalance(i)$( cytosol(i) and not extracellular(i) )..  sum(j$(S(i,j)),S(i,j)*v(j)) =g= 0;
massbalance_comp(i)$( not cytosol(i) and not extracellular(i) )..  sum(j$(S(i,j)),S(i,j)*v(j)) =e= 0;

*** Constraints on flux of model reactions
* Reversing a one-directional reaction requires the binary variable y to be 
modelfluxrev_min(j)$( mod(j) and rev(j) ).. v(j) =g= LowerLimits(j);
modelfluxirrev_min(j)$( mod(j) and not rev(j) ).. v(j) =g= -M*ym(j);
modelflux_max(j)$( mod(j) ).. v(j) =l= UpperLimits(j);

*** Constraints on flux of database reactions
* These constraints ensure that y is set when a reaction is active
databaseflux_min(j)$( database(j) and not mod(j) ).. v(j) =g= LowerLimits(j)*yd(j);
databaseflux_max(j)$( database(j) and not mod(j) ).. v(j) =l= UpperLimits(j)*yd(j);

*** Constraints on the production of the metabolite in question by reaction
* Ensure that the binary variable w is only one in the case where at least
* epsilon units are produced.
prodcons_min(i,j)$( blocked(i) and S(i,j) ne 0 ).. S(i,j)*v(j) =g= delta-M*(1-w(i,j));
prodcons_max(i,j)$( blocked(i) and S(i,j) ne 0 ).. S(i,j)*v(j) =l= M*w(i,j);

*** Constraint on the production of the metabolite in question in general
* Ensure that the metabolite in question is produced by at least one
* reaction.
binarycons(i)$( blocked(i) ).. sum(j$((S(i,j) ne 0 and rev(j) ) or (S(i,j) gt 0 and not rev(j) )), w(i,j)) =g= 1;

*** Minimize the number of model reaction modifications
* Minimize the sum of added database reactions and reversed model reactions.
obj .. z =e= sum(j, yd(j)) + sum(j, ym(j));


Model gapfill /
        massbalance
        massbalance_comp
        modelfluxrev_min
        modelfluxirrev_min
        modelflux_max
        databaseflux_min
        databaseflux_max
        prodcons_min
        prodcons_max
        binarycons
        obj
/;

*** Fix the reactions both in model and database to flux 0 ??
v.fx(j)$( mod(j) and database(j) ) = 0;
v.fx(j)$( not mod(j) and not database(j) ) = 0;

solve gapfill using mip minimizing z;
