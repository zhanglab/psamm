Automated Comparison of Models Using ModelMapping
=================================================

This tutorial will go over how to use the ``modelmapping`` function
in `PSAMM`. This function can be used to automatically compare different models
based on the profiles of metabolites and reactions.

.. contents::
   :depth: 1
   :local:

Materials
---------

For information on how to install `PSAMM` and the associated requirements, as well
how to download the materials required for this tutorial, you can reference the
Installation and Materials section of the tutorial.

For this part of tutorial, we will be mapping a modified version of the Shewanella
core metabolic model to a modified version of the *E. coli* core metabolic model
that has been used in the other sections of the tutorial. To access these models
and other additional files you will need to go into the tutorial-part-5 folder
located in the psamm-tutorial folder.

.. code-block:: shell
    
    (psamm-env) $ cd <PATH>/tutorial-part-5/

Once in this folder, you should see several folders. E_coli_core and Shewanella_core
are the two folders that store the corresponding metabolic models, while the
folder called `prerun` is where pre-run modelmapping result is located.

.. code-block:: shell
    
    (psamm-env) $ ls

    additional_files  E_coli_core  prerun  Shewanella_core

Sub-commands in ``modelmapping`` function
-----------------------------------------

The ``modelmapping`` command in PSAMM have three sub-commands, ``map``,
``manual_curation`` and ``translate_id``. The list of available sub-commands
and their meaning can be shown when you run the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model modelmapping -h

    usage: psamm-model modelmapping [-h] {map,manual_curation,translate_id} ...

    Generate a common model for two distinct metabolic models.

    optional arguments:
      -h, --help            show this help message and exit

    Tools:
        {map,manual_curation,translate_id}
            map                 Bayesian model mapping
            manual_curation     Interactive mode to check mapping result
            translate_id        Translate ids based on mapping result
            
To run a specific sub-command, e.g. Bayesian model mapping, just type in the
function together with the sub-command ``modelmapping map``. Each sub-command
also have a related help document, which can be accessed by setting ``-h``
parameter via command line.

Automatically map metabolic models
----------------------------------
It is a common situation that different models use different naming systems for
their compounds and reactions. For example, for the same compound `oxygen`, in a
KEGG-based model it's called `C00007`, while in a BiGG-based model
it's called `o2`. This kind of inconsistency will make the direct comparison
between different models nearly impossible.

To address this problem, a `ModelMapping` pipeline was developed to map the
compounds and reactions in one model to another using Bayesian probability on
the similarity of several properties. For compound mapping, the properties of
id, name, charge and formula are always compared, and "KEGG ID" is an additional
property to consider if both models have that information for their compounds.
For reaction mapping, the properties of id, name and equation are compared,
while "Gene association" is an add on. 

+-----------------------+
|Properties of compounds|
+=======================+
|ID                     |
+-----------------------+
|Name                   |
+-----------------------+
|Charge                 |
+-----------------------+
|Formula                |
+-----------------------+
|*KEGG ID*              |
+-----------------------+

+-----------------------+
|Properties of reactions|
+=======================+
|ID                     |
+-----------------------+
|Name                   |
+-----------------------+
|Equation               |
+-----------------------+
|*Gene association*     |
+-----------------------+


Basic use of the ``modelmapping map`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``modelmapping map`` command in PSAMM is used to perform the mapping from
one model to another. The basic command looks like the following:

.. code-block:: shell

    (psamm-env) $ psamm-model --model E_coli_core modelmapping map \
    --dest-model Shewanella_core \
    -o modelmapping

The query model is set via the ``--model`` parameter immediately following
``psamm-model``, just like all the other commands in SPAMM. The targeting model
is set via the ``--dest-model``. In this case, the command will map the compounds
and reactions from the *E. coli* core model to the ones in the Shewanella core
model.

This command will create the folder `modelmapping`, then put the two output
files into that folder. The file `bayes_compounds_best.tsv` stores every
compound in the query model and its best mapping results in the target model,
while the file `bayes_reactions_best.tsv` stores the mapping results for
reactions.

The output files are tab-delimited tables. The first column stores the compound or
reaction ids in the query model, the second column stores the best mapping ids in
the target model, the third column stores the Bayesian probability of this
mapping pair. The rest columns show the score of mapping pairs on a specific
property.

.. code-block:: shell

    (psamm-env) $ head modelmapping/bayes_compounds_best.tsv

    e1      e2      p       p_id    p_name  p_charge        p_formula       p_kegg
    10fthf  cpd_10fthf      0.9999472032623282      0.0005017605633802816   0.7516733067729083      0.0040118652717529985   0.6446583143507973      0.0008362676056338028
    12dgr120        Biomass 9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_12dgr       9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_26dap-LL    9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_26dap-M     9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_2ahhmp      9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_2aobut      9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_2dmmq7      9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028
    12dgr120        cpd_2ohph       9.668708566791001e-05   0.0005017605633802816   0.00033473048314771204  0.0040118652717529985   8.372450922181071e-05   0.0008362676056338028

Additional mapping options in ``modelmapping map``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The KEGG ids of the compounds can be compared by setting the parameter
``--map-compound-kegg``.

For reaction comparison, ``--map-reaction-gene`` can
be set to compare the gene association information. If the gene ids in the query
and the target model are not directly comparable, the parameter
``--gene-map-file`` can be set to an additional file listing the mapping of
gene ids.

.. code-block:: shell

    (psamm-nev) $ psamm-model --model E_coli_core modelmapping map \
    --dest-model Shewanella_core \
    --map-compound-kegg \
    --map-reaction-gene \
    --gene-map-file additional_files/gene_map.tsv \
    -o modelmapping

Additional output options in ``modelmapping map``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the command will output the best mapping result for every compound
and reaction in the query model. However, most of them may not be true. The
command can be set to output only the mapping results with a probability higher
than the given threshold. This will significantly reduce the size of output,
but the selection of thresholds is pretty arbitrary, a better way of filtering
the mapping result will be the `interactive manual curation`_ that will be
discussed in the next section.

.. _`interactive manual curation`: `Manually curate the mapping result using interactive interface`_


.. code-block:: shell

    (psamm-env) $  psamm-model --model E_coli_core modelmapping map \
    --dest-model Shewanella_core \
    -o modelmapping/ \
    --threshold-compound 0.01 \
    --threshold-reaction 1e-4

Manually curate the mapping result using using interactive interface
-------------------------------------------------------------------------

The command ``modelmapping manual_curation`` will enter an interactive
interface to assist the manual curation of the mapping result.

.. code-block:: shell

    (psamm-env) $ psamm-model --model E_coli_core/ modelmapping manual_curation \
    --dest-model Shewanella_core/ \
    --compound-map modelmapping/bayes_compounds_best.tsv \
    --reaction-map modelmapping/bayes_reactions_best.tsv \
    --curated-compound-map modelmapping/curated_compound_map.tsv \
    --curated-reaction-map modelmapping/curated_reaction_map.tsv

The parameters ``--compound-map`` and ``--reaction-map`` refer to the output
files from the ``modelmapping map`` command, while ``--curated-compound-map``
and ``--curated-reaction-map`` refer to the files that store the true mapping
results after manual check. Besides the curated mapping files, the command will
also store false mappings into `.false` files, and compounds and reactions
to be ignored can be stored in `.ignore` files. For example, if the
``--curated-compound-map`` is set to `curated_compound_mapping.tsv`, then the
false mappings will be stored
in `curated_compound_mapping.tsv.false`, and the pairs to be ignored
should be stored in `curated_compound_mapping.tsv.ignore`.

.. code-block:: shell

    Below are the compound mapping involved:

    p            0.999947
    p_id         0.000502
    p_name       0.751673
    p_charge     0.004012
    p_formula    0.644658
    p_kegg       0.000836
    Name: cpd_atp, dtype: float64


    id: atp
    charge: -4
    formula: C10H12N5O13P3
    name: ATP
    Defined in E_coli_core/compounds.yaml


    id: cpd_atp
    charge: -4
    formula: C10H12N5O13P3
    name: ATP
    Defined in Shewanella_core/compounds.tsv:170


    True compound match? (y/n/ignore/save/stop, type ignore to ignore this compound in future, type save to save current progress, type stop to save and exit):

.. code-block:: shell

    Here is the reaction mapping:
    p             1.000000
    p_id          1.000000
    p_name        0.001044
    p_equation    0.895363
    p_genes       0.001305
    Name: (HCO3E, HCO3E), dtype: float64


    id: HCO3E
    equation: co2[c] + h2o[c] <=> h[c] + hco3[c]
    equation (compound names): CO2[c] + H2O[c] <=> H+[c] + Bicarbonate[c]
    ec: 4.2.1.1
    genes: NP_414668.1
    name: HCO3 equilibration reaction
    subsystem: Unassigned
    Defined in E_coli_core/reactions.yaml


    id: HCO3E
    equation: cpd_co2[c] + cpd_h2o[c] <=> cpd_h[c] + cpd_hco3[c]
    equation (compound names): CO2[c] + H2O[c] <=> H+[c] + Bicarbonate[c]
    ec: 4.2.1.1
    genes: WP_020912789.1
    name: carbonate dehydratase (HCO3 equilibration reaction)
    subsystem: Unassigned
    Defined in Shewanella_core/reactions_core.yaml


    These two reactions have the following curated compound pairs:

    co2[c] cpd_co2[c] 0.9999472032623282
    h2o[c] cpd_h2o[c] 0.9999472032623282
    h[c] cpd_h[c] 0.9999472032623282
    hco3[c] cpd_hco3[c] 0.9999472032623282


    4 compounds in HCO3E
    4 compounds in HCO3E
    4 curated compound pairs

    True reaction match? (y/n/ignore/save/stop, type ignore to ignore this reaction in future, type save to save current progress, type stop to save and exit):

In the interface,
mapping pairs and their related information will be shown one-by-one. The
current pair can be marked as true pair by type in "y", while false pair can
be marked by type in "n". The current progress can be saved by type in "save",
and type in "stop" will save then exit the interface. It's totally OK to type
in a typo in the interface, the program will ignore that input and ask again.

.. code-block:: shell

    True compound match? (y/n/ignore/save/stop, type ignore to ignore this compound in future, type save to save current progress, type stop to save and exit): saave
    True compound match? (y/n/ignore/save/stop, type ignore to ignore this compound in future, type save to save current progress, type stop to save and exit):

Some of the compounds or reactions in the query model
don't actually have a true mapping in the target model, it will be very
annoying to type in "n" again and again for those false mapping pairs. In such
case, type in "ignore" can ignore the mapping pairs of this 
compound or reaction in the future.

If the curated files already exist, the command will consider them
as the previous progress, then append new curation results. So it's totally safe
to do part of curation at one time, save and exit, then resume the progress
at another time.

Translate ids in the query model based on curated mapping result
----------------------------------------------------------------

A common application of model mapping is to translate the compound and reaction
ids in the query model to the ones used in the target model, therefore the two
models become directly comparable. Once the manual curation is finished, the
translation can be done by the ``modelmapping translate_id`` command.

.. code-block:: shell

    (psamm-env) $ psamm-model --model E_coli_core/ modelmapping translate_id \
    --compound-map modelmapping/curated_compound_map.tsv \
    --reaction-map modelmapping/curated_reaction_map.tsv \
    -o modelmapping/translated

A yaml formatted model will be stored in the `modelmapping/translated` folder.
