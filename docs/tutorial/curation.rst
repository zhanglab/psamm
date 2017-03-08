
Model Curation
==============

This tutorial will go over how to utilize the curation functions in PSAMM
to correct common errors and ensure that metabolic reconstructions are
accurate representations of the metabolism of an organism.

.. contents::
   :depth: 1
   :local:


{{{Reference to the installation tutorial}}}



Common Errors in Metabolic Reconstructions
------------------------------------------

Many types of errors can be introduced into metabolic models. Some errors
can be introduced during manual editing of model files while others can result
from inconsistent representations of the biology of the system. Various features
in PSAMM are designed ot help identify and fix these problems to ensure that
the reconstruction does not contain these kinds of errors.

Some errors cannot be easily identified without extensive manual inspection of
the model data files. These PSAMM functions are designed to help identify these
errors and make the correction process easier.

PSAMM Warnings
--------------

The most basic way to identify possible errors in a model will be through
reading the warning messages printed out by PSAMM when any functions are run
on a model. These warning messages can be an easy way to identify if something
in the reconstruction is not set up the way that was intended. The following
are examples of the types of warnings that PSAMM will provide and what kinds of
errors they might indicate.

The first type of warning that PSAMM can provide is a waning that there is a
compound that is in a reaction but is not defined in the compound information
of the model. While PSAMM doesn't necessarily know if this is an error,
these warning can help identify compound ids in the reconstruction that
may have typos in them or that need to be defined in the compounds data
for the reconstruction. For example in the warning below it would appear that
the compound id for ATP had been mistyped and included two extra t's in it.
These types of errors can make reactions in a model inconsistent and may lead
to incorrect conclusions from the model if they are not corrected.

.. code-block:: shell

    WARNING: The compound cpd_atttp was not defined in the list of compounds

The second type of warning will similarly help identify if there was an error
introduced in one of the reconstruction's reactions. This warning will indicate
that there is a compound present in the reconstruction that has a compartment that
is not defined elsewhere in the model. In the example below a compound was added in
a reaction as being in the compartment 'X'. Since this compartment was not used in the
model the reaction involving this instance of the compound would become flux inconsistent.

.. code-block:: shell

    WARNING: The compartment X was not defined in the list of compartments.

The third and fourth types of warnings can be useful in identify that the exchange file
is set up correctly for the reconstruction. These two kinds of errors will help identify if
there are compounds that are present in the extracellular compartment but do not have a corresponding
exchange reaction in the boundary conditions. This can be problematic for some models that require certain
sinks for overproduced compounds in the boundary. The other kind of warning will indicate if there are
compounds in the exchange reactions that cannot be utilized by any reactions in the model. This could
indicate that a transport reaction is missing from the model or that the compound could be removed
from the exchange file.

.. code-block:: shell

    WARNING: The compound cpd_chitob was in the extracellular compartment but not defined in the medium
    WARNING: The compound cpd_etoh was defined in the medium but is not in the extracellular compartment

Reaction Consitency in PSAMM
----------------------------

The previous examples of warning messages produced by PSAMM can be helpful as a first step in identifying
possible errors in a model but there are various other types of errors that may be present in models that
specific PSAMM functions can help identify. The first kind of errors are ones related to the balancing of
reactions in model. It is important that metabolic models be balanced in terms of elements, charge, and
stoichiometry. PSAMM has three functions available to identify reactions that are not balanced
in these properties which can help correct them and lead to more accurate and true representations of
metabolism.

Stoichiometric Checking
~~~~~~~~~~~~~~~~~~~~~~~

PSAMM's masscheck tool can be used to check if the reactions in the model are
stoichiometrically consistent and the compounds that are causing the imbalance.
This can be useful when curating the model
because it can assist in easily identify missing compounds in reactions.
A common problem that can be identified using this tool is a loss of
hydrogen atoms during a metabolic reaction. This can occur due to modeling
choices or incomplete reaction equations but is generally easy to identify
using masscheck.

To report on the compounds that are not balanced use the following masscheck
command:

.. code-block:: shell

    (psamm-env) $ psamm-model masscheck

This command will produce an output like the following:

.. code-block:: shell

    ...
    accoa_c	1.0	Acetyl-CoA
    acald_e	1.0	Acetaldehyde
    acald_c	1.0	Acetaldehyde
    h_e	0.0	H
    h_c	0.0	H
    INFO: Consistent compounds: 73/75

The ``masscheck`` command will first try to assign a positive mass to all of
the compounds in the model while balancing the masses such that the left-hand
side and right-hand side add up in every model reaction. All the compound
masses are reported, and the compounds that have been assigned a zero value for
the mass are the ones causing imbalances.

In certain cases a metabolic model can contain compounds that represent electrons,
photons, or some other artificial compound. These compounds can cause problems with
the stoichiometric balance of a reaction because of their unique functions. In order
to deal with this an additional property can be added to the compound entry that
will designate it as a compound with zero mass. This designation will tell PSAMM
to consider these compounds to have no mass during the stoichiometric checking which
will prevent them from causing imbalances in the reactions. An example of how to add
that property to a compound entry can be seen below:

.. code-block:: yaml

    - id: phot
      name: Photon
      zeromass: yes

To report on the specific reactions that may be causing the imbalance, the
following command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-model masscheck --type=reaction
    ...
    FRUKIN	1.0	|Fructose[c]| + |ATP[c]| => |D-Fructose-6-phosphate[c]| + |ADP[c]| + |H[c]|
    INFO: Consistent reactions: 100/101

This check is performed similarly to the compound check. In addition, mass
residual values are introduced for each metabolic reaction in the network.
These mass residuals are then minimized and any reactions that result in a
non-zero mass residual value after minimization are reported as being
stoichiometrically inconsistent. A non-zero residual value after minimization
tells you that the reaction in question may be unbalanced and missing
some mass from it.

Sometimes the residue minimization problem may have multiple solutions. In
these cases the residue value may be reallocated among a few connected
reactions. In this example the unbalanced reaction is the MANNIDEH reaction::

    MANNIDEH    |manni[c]| + |nad_c[c]| => |fru_c[c]| + |nadh_c[c]|

In this reaction equation the right hand side is missing a proton. However
minimization problem can result in the residue being placed on either the
`fru_c` or the `nadh_c` compounds in an attempt to balance the reaction.
Because `nadh_c` occurs in thirteen other reactions in the network, the
program has already determined that that compound is stoichiometrically
consistent. On the other hand `fru_c` only occurs one other time. Since
this compound is less connected the minimization problem will assign the
non-zero residual to this compound. This process results in the FRUKIN reaction
which contains this compound as being identified as being stoichiometrically
inconsistent.

In these cases you will need to manually check the reaction and then use
the ``--checked`` option for the ``masscheck`` command to force the non-zero
residual to be placed on a different reaction. This will rerun the consistency
check and force the residual to be placed on a different reaction. To do this
we would run the following command.

.. code-block:: shell

    (psamm-env) $ psamm-model masscheck --type=reaction --checked FRUKIN

Now, the output should report the `MANNIDEH` reaction and it can be seen that
the reaction equation of `MANNIDEH` is specified incorrectly. It appears that a
hydrogen compound was left out of the reaction for `MANNIDEH`. This would be an
easy problem to correct by simply adding in a hydrogen compound to correct the
lost atom in the equation.

The stoichiometric consistency checking allows for the easy identification of
stoichiometrically inconstent compounds while providing a more targeted subset
of reactions to check to fix the problem. This allows you to quickly identify
problematic reactions rather than having to manually go through the whole
reaction database in an attempt to find the problem.

In some cases there are reactions that are going to be inherently unbalanced
and might cause problems with using these methods. If you know that this is the
case for a specific reaction they can specify that the reaction be excluded
from the mass check so that the rest of the network can be analyzed. To do this
the ``--exclude`` option can be used. For example if you wanted to exclude the
reaction `FRUKIN` from the mass check they could use the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model masscheck --exclude FRUKIN

This exclude option can be helpful in removing inherently unbalanced reactions
like macromolecule synthesis reations or incomplete reactions that would be
identified as being stoichiometrically inconsistent. It is also possible to
create a file that lists multiple reactions to exclude. Put each reaction
identifier on a separate line in the file and refer to the file be prefixing
the file name with a ``@``:

.. code-block:: shell

    (psamm-env) $ psamm-model masscheck --exclude @excluded_reactions.txt

Before we fix the model with the correction to the `MANNIDEH` reaction, let us
first check the model for formula inconsistencies to show how this can also be
used in conjunction with mass checking and other methods to correct model
inconsistencies.

Formula Consistency Checking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Formula checking will check that each reaction in the model is balanced
with respect to the chemical formulas of each compound. To check the model
for formula consistencies run the formula check command:

.. code-block:: shell

    (psamm-env) $ psamm-model formulacheck

The output should appear as follows:

.. code-block:: shell

    INFO: Model: Ecoli_core_model
    INFO: Model Git version: 9812080
    MANNIDEH	C27H40N7O20P2	C27H39N7O20P2		H
    Biomass_Ecoli_core_w_GAM	C1088.0232H1471.1810N446.7617O1236.7018P240.5298S3.7478	C1045.4677H1395.2089N441.3089O1189.0281P236.8511S3.7478		C42.5555H75.9721N5.4528O47.6737P3.6787
    INFO: Unbalanced reactions: 2/80
    INFO: Unchecked reactions due to missing formula: 0/80

In this case two reactions are identified in the model as being unbalanced.
The biomass objective function, `Biomass_Ecoli_core_w_GAM`, and the
reaction that was previously identified through masscheck as being
unbalanced, `MANNIDEH`. In the case of the objective function this is
imbalanced due to the formulation of the objective function. The reaction
functions as a sink for the compounds required for growth and only outputs
depleted energy compounds. This leads to it being inherently formula
imbalanced but it is a necessary feature of the model. The other reaction
is `MANNIDEH`. It can be seen that the total number of atoms on each side
does not match up. PSAMM also outputs what atoms would be needed to balance
the reaction on both sides. In this case there is a missing hydrogen atom
on the right side of the equation. This can be easily rectified by adding
in the missing hydrogen. To do this correction in this tutorial, you
can copy a fixed version of the mannitol pathway from the additional files
folder using the following command:

.. code-block:: shell

    (psamm-env) $ cp ../additional_files/mannitol_pathway_v2.yaml mannitol_pathway.yaml

Once that problem with the new reaction is fixed the model will pass both the
formula check and mass check.

Now this fix can be added to the Git repository so that the latest version
of the model will be the fixed version. To do this the following commands
can be used.

.. code-block:: shell

    (psamm-env) $ git add mannitol_pathway.yaml
    (psamm-env) $ git commit -m 'Fixed mass and formula inconsistencies in Mannitol pathway'

Charge Consistency Checking
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The charge consistency function is similar to the formula consistency function
but instead of using the chemical formulas for the compounds, PSAMM
will use the assigned
charges that are designated in the compounds file and check that these
charges are balanced on both sides of the reaction.

To run a charge consistency check on the model use the chargecheck command:

.. code-block:: shell

    (psamm-env) $ psamm-model chargecheck

This `E. coli` SBML model does not contain charge information for the
compounds. A sample output is provided below to show what the results would
look like for a charge imbalanced model. The output from the charge check will
display any reactions that are charge imbalanced and show what the imbalance is
and then show the reaction equation. This can be used to quickly check for any
missed inconsistencies and identify reactions and compounds that should be
looked at more closely to confirm their correctness.

.. code-block:: shell

    ...
    rxn12510	1.0	|ATP[c]| + |Pantothenate[c]| => |4-phosphopantothenate[c]| + |H+[c]| + |ADP[c]|
    rxn12825	4.0	|hemeO[c]| + |H2O[c]| => |Heme[c]| + (4) |H+[c]|
    rxn13643	1.0	|ADP-glucose[c]| => |Glycogen[c]| + |H+[c]| + |ADP[c]|
    rxn13710	6.0	(5) |D-Glucose[c]| + (4) |ATP[c]| => |Glycogen[c]| + (4) |H+[c]| + (4) |Phosphate[c]| + (4) |H2O[c]| + |ADP[c]|
    INFO: Unbalanced reactions: 94/1093
    INFO: Unchecked reactions due to missing charge: 0/1093

Flux Consistency Checking
~~~~~~~~~~~~~~~~~~~~~~~~~

The flux consistency checking function can be used to identify reactions that
cannot carry flux in the model. This tool can be used as a curation tool as
well as an analysis tool. In this tutorial it will be highlighted for the
curation aspects and later its use in flux analysis will be demonstrated.

To run a flux consistency check on the model use the ``fluxcheck`` command:

.. code-block:: shell

    (psamm-env) $ psamm-model fluxcheck --unrestricted

The unrestricted option with the command will tell PSAMM to
remove any limits on the exchange reactions. This will tell you which
reactions in the model can carry flux if the model is given all compounds in
the media freely. This can be helpful for identifying which reactions may not
be linked to other parts of the metabolism and can be helpful in identifying
gaps in the model. In this case it can be seen that no reactions were identified
as being inconsistent.

In some situations there are pathways that might be
modeled but not necessarily connected to the other aspects of metabolism.
A common occurance of this is with vitamin biosynthesis pathways that are
not incorporated into the biomass in the model. Fluxcheck will identify
these as being flux inconsistent but the modeler will need to identify if this
is due to incomplete information on the pathways or if it is due to some
error in the formulation of the reactions.

PSAMM will tell you how many exchange reactions cannot be used as well
as how many internal model reactions cannot carry flux. PSAMM will also
list the reactions and the equations for the reactions to make curation of
these reactions easier.

Above the fluxcheck command was used with the --unrestricted option which
allowed the exchange reactions to all be active. This command can also be
used to see what reactions cannot carry flux when specific media are
supplied. To run this command on the network with the media that is
specified in the media file run the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model fluxcheck
    INFO: Model: Ecoli_core_model
    INFO: Model Git version: 9812080
    INFO: Using flux bounds to determine consistency.
    ...
    EX_fru_e	|D-Fructose[e]| <=>
    EX_fum_e	|Fumarate[e]| <=>
    EX_glc_e	|D-Glucose[e]| <=>
    EX_gln_L_e	|L-Glutamine[e]| <=>
    EX_mal_L_e	|L-Malate[e]| <=>
    FRUpts2	|D-Fructose[e]| + |Phosphoenolpyruvate[c]| => |D-Fructose-6-phosphate[c]| + |Pyruvate[c]|
    FUMt2_2	(2) |H[e]| + |Fumarate[e]| => (2) |H[c]| + |Fumarate[c]|
    GLCpts	|Phosphoenolpyruvate[c]| + |D-Glucose[e]| => |Pyruvate[c]| + |D-Glucose-6-phosphate[c]|
    GLNabc	|ATP[c]| + |L-Glutamine[e]| + |H2O[c]| => |L-Glutamine[c]| + |ADP[c]| + |H[c]| + |Phosphate[c]|
    MALt2_2	|L-Malate[e]| + (2) |H[e]| => |L-Malate[c]| + (2) |H[c]|
    INFO: Model has 5/80 inconsistent internal reactions (0 disabled by user)
    INFO: Model has 5/21 inconsistent exchange reactions (0 disabled by user)

In this case it can be seen that there are various exchange reactions
blocked as well as various internal reactions related to other carbon
metabolic pathways. The current model should only be supplying mannitol
as a carbon source and this would mean that these other carbon pathways
would be blocked in this condition. In this way, you can use the ``fluxcheck``
command to see what reactions are specific to certain metabolic pathways and
environmental conditions.



Gap Identification in PSAMM
---------------------------

In addition to inconsistencies found within individual reactions there can also be
global inconsistencies for the reactions within a metabolic network. These include
metabolites that can be produced but not consumed, ones that can be consumed by reactions
but are not produced, and reactions that cannot carry flux in a model. PSAMM includes various
functions for the identification of these features in a network including the functions ``gapcheck`` and
``fluxcheck``. Additionally the functions ``gapfill`` and ``fastgapfill`` can be used to help
fill these gaps that are present through the introduction of additional reactions into the network.

Gapcheck in PSAMM
~~~~~~~~~~~~~~~~~

The ``gapcheck`` function in psamm can be used to identify dead end metabolites in a metabolic network.
These dead end metabolites are compounds in the metabolic model that can either be produced but not consumed
or ones that can be consumed but not produced. Reactions that contain these compounds cannot carry flux within
a model and are often the result of knowledge gaps in our understanding of metabolic networks.

The ``gapcheck`` function allows the use of three methods for the identification of these dead end metabolites
within a metabolic network. These are the ``prodcheck``, ``concheck``, and ``gapfind`` methods.

The ``prodcheck`` method is the most straightforward of these methods and can be used to identify any
compounds that cannot be produced in the metabolic network. Will iterate through the reactions in a network
and maximize each one. If the reaction can carry a flux then the metabolites involved in the reaction
are not considered to be blocked.

To use this function the following command can be run:

.. code-block:: shell

    (psamm-env) $ psamm-model gapcheck --method prodcheck

The function will produce output like the following that lists out any metabolites in the
model that cannot be produced in this condition:

.. code-block:: shell


This result indicates that the following metabolites currently cannot be produced in the model.
This only tells part of the story though as this function was run with the defined media that
was set for the model. As a result there are gaps identified like, '', that can be produced in other
conditions. To do a global check using this function on the model without restrictions on the media
the following command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-model gapcheck --method prodcheck --unrestricted-exchange

The unrestricted tag in this function will allow temporarily set all of the exchange
reaction bounds to be -1000 to 1000 allowing all nutrients to be either taken up or produced.
gap-checking in this condition will allow for the identification of gaps that are not media dependent
and may instead be the result of incomplete pathways and knowledge gaps.


The second method implemented in the ``gapcheck`` function is the ``sinkcheck`` method. This method is similar to
``prodcheck`` but is implemented in a way where the flux through each introduced sink for a compound is maximized.
This ensures that the metabolite can be produced in excess from the network for it to not be considered a dead end
metabolite.

.. code-block:: shell

    (psamm-env) $ psamm-model gapcheck --method sinkcheck --unrestricted-exchange


The last method implemented in the ``gapcheck`` function is the ``gapfind`` method. This method is an implementation
of a previously published method to identify gaps in metabolic networks "Include CITATION HERE". This method will use
a network based optimization to identify metabolites with no production pathways present.

.. code-block:: shell

    (psamm-env) $ psamm-model gapcheck --method gapfind --unrestricted-exchange




Search Functions in PSAMM
-------------------------

``psamm-model`` includes a search function that can be used to search the model
information for specific compounds or reactions. To do this the search function
can be used. This can be used for various search methods. For example to search
for the compound named fructose the following command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-model search compound --name 'Fructose'
    INFO: Model: Ecoli_core_model
    INFO: Model Git version: db22229
    id: fru_c
    formula: C6H12O6
    name: Fructose
    Defined in ./compounds.yaml:?:?

To do the same search but instead use the compound ID the following command can
be used:

.. code-block:: shell

    (psamm-env) $ psamm-model search compound --id 'fru_c'

These searches will result in a printout of the relevant information contained
within the model about these compounds. In a similar way reactions can also be
searched. For example to search for a reaction by a specific ID the following
command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-model search reaction --id 'FRUKIN'

Or to search for all reactions that include a specific compound the following
command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-model search reaction --compound 'manni[c]'