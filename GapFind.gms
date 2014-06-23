*** GapFind
* This model can identify the root and downstream blocked metabolites in a
* metabolic reconstruction. The binary variable xp will indicate whether a metabolite
* can be produced in the metabolic network.
*
* This implements the GapFind model described by Kumar, Vinay Satish, Madhukar S. Dasika,
* and Costas D. Maranas. "Optimization based automated curation of metabolic reconstructions."
* BMC bioinformatics 8.1 (2007): 212.


Set i metabolites in S (m) /
$include "metabolites.txt"
/;

Set j reactions in S (n) /
$include "rxnnames.txt"
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

Parameters
        UpperLimits(j)          maximum value flux can take
        LowerLimits(j)          minimum value flux can take
        Vmax    /100/
        epsilon /0.001/
        M       /100/
;

Parameter S(i,j) stoichiometric coefficients (S matrix) /
$include "mat.txt"
/;

*** Set lower and upper bounds on fluxes accordingly
LowerLimits(j) = 0;
LowerLimits(j)$rev(j) = -Vmax;
UpperLimits(j) = Vmax;

Variables
        v(j)			flux value through reaction
        z			objective
;	

Binary variables
        xp(i)                   whether metabolite i can be produced by the network
        w(i,j)                  whether reaction j that produces metabolite i is active
;

Equations
        massbalance(i)          mass balance for cytosolic metabolites
        massbalance_comp(i)     mass balance for metabolites in internal compartments
        prodconsirrev_min(i,j)  minimum production binary constraint (irreversible)
        prodconsirrev_max(i,j)  maximum production binary constraint (irreversible)
        prodconsrev_min(i,j)    minimum production binary constraint (reversible)
        prodconsrev_max(i,j)    maximum production binary constrains (reversible)
        binarycons(i)           metabolite production binary constraint
        obj                     objective function
;

*** The mass balance constraints
* There is an implicit metabolite sink for cytosolic metabolites because of
* the diluting effect from cell division, export from diffusion, etc., which means
* that the production of these metabolites only has to be positive. For metabolites
* in internal compartments the uptake and secretion has to even out.
massbalance(i)$(cytosol(i) and not extracellular(i)).. sum(j$S(i,j),S(i,j)*v(j)) =g= 0;
massbalance_comp(i)$(not cytosol(i) and not extracellular(i)).. sum(j$S(i,j),S(i,j)*v(j)) =e= 0;

*** Constraints on production of metabolites in reaction
* Ensure that the binary variable w is only one in the case where at least
* epsilon units are produced.
prodconsirrev_min(i,j)$( (S(i,j) gt 0 and not rev(j)) ).. S(i,j)*v(j) =g= epsilon*w(i,j);
prodconsirrev_max(i,j)$( (S(i,j) gt 0 and not rev(j)) ).. S(i,j)*v(j) =l= M*w(i,j);
prodconsrev_min(i,j)$( (S(i,j) ne 0 and rev(j)) ).. S(i,j)*v(j) =g= epsilon-M*(1-w(i,j));
prodconsrev_max(i,j)$( (S(i,j) ne 0 and rev(j)) ).. S(i,j)*v(j) =l= M*w(i,j);

*** Constraint on production of metabolites in general
* Ensure that the binary variable xp is only one in the case where at least
* one of w for the metabolite in question is one.
binarycons(i).. sum(j$((S(i,j) ne 0 and rev(j) ) or (S(i,j) gt 0 and not rev(j) )), w(i,j)) =g= xp(i);

*** Maximize the number of metabolites that can be produced
obj.. z =e= sum(i,xp(i));


Model gapfind /
        massbalance
        massbalance_comp
        prodconsirrev_min
        prodconsirrev_max
        prodconsrev_min
        prodconsrev_max
        binarycons
        obj
/;

v.lo(j) = LowerLimits(j);
v.up(j) = UpperLimits(j);
	
solve gapfind using mip maximizing z;
