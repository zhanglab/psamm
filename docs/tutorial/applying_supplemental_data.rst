Supplementing your model with other data
========================================

This tutorial will go over using supplementary data from outside sources, such
as the utilization of gene expression data to subset a model into "active"
reactions, the use of a template model to generate a draft model based on
reciprocal best hits, or the application of thermodynamic properties to further
constrain modeling results through the gimme, psammotate, and TMFA.

.. contents::
   :depth: 1
   :local:

Materials
---------

The materials used in the part of the tutorial can be found in the `tutorial-part-X`
directory in the psamm-tutorial repository. The previously used E coli model will
be reused here in order to demonstrate the use of the psammotate and gimme functions,
which will be located in the ``tutorial-part-X/E_coli_yaml/`` directory. There
will also be a ``tutorial-part-X/gimme``, ``tutorial-part-X/psammotate``, and
``tutorial-part-X/TMFA`` directories that contain the required materials for each
of the following functions.


.. code-block:: shell

   (psamm-env) $ cd <PATH>/tutorial-part-X/E_coli_yaml


Subsetting a Metabolic Models Using Expression Data with Gimme
--------------------------------------------------------------

This tutorial will go over using the ``gimme`` function to subset a metabolic
model and create a new metabolic model based on gene expression data.


Constructing a Transcriptome File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Gimme algorithm functionally subsets a model based on expression data provided
using the ``--transcriptome-file`` argument. This file should contain two columns:
one with each gene within a model and one with a measure of the expression of the
gene, such as transcripts per million (TPM) or Reads Per Kilobase of transcript,
per Million mapped reads (RPKM). For the purposes of this tutorial, expression
data has been mocked up by randomly generating expression numbers for genes in
the E coli model.

Basic use of the ``Gimme`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run Gimme, you must first specify the directory of the two column transcriptome
file and specify a threshold value that defines at what threshold below which
genes will be dropped from the new model. This produces a new model that drops
any reactions coded by genes that do not meet this threshold. Any reactions assoociated
with multiple genes must meet the threshold across all. Running the command as follows
will print out a list of the internal reactions and a list of the external reaction
that meet the thresholds.

.. code-block:: shell

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100

As you increase the expression threshold, you will see fewer reactions retained
within the model, which is shown below by using the word count function to count
the number of lines returned by each iteration of gimme:

.. code-block:: shell

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 | wc -l
      86

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 500 | wc -l
      74

This application of the gimme algorithm sets the objective biomass to 100%. That
is to say that if a reaction is below the required threshold but is required to
maintain biomass, it is kept. You can lower the required biomass with the
``--biomass-threshold`` flag, which will subsequently drop increasing numbers of
reactions as they are no longer necessary to meet the biomass threshold. This is
shown below by using the word count function to count the number of lines returned
by each iteration of gimme:

.. code-block:: shell

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0.93757758084 | wc -l # Maximum for this model
      86

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0 | wc -l # No Biomass Threshold
      82

In addition to simply writing out a list of reactions that satisfy the expression
and biomass thresholds, by specifying ``--export-model``, you can redirect the
output to create an entirely new model.

.. code-block:: shell

    (psamm-env) $ mkdir gimme_out

    (psamm-env) $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --export-model ./gimme_out/


Generating a Draft Model From a Template Using Psammotate
---------------------------------------------------------

This tutorial will go over using the ``psammotate`` function to generate draft
models based on a reciprocal best hits file between the two models that
provides gene associations based on mapping the genes from a reference file
onto the genes of a draft model.

.. contents::
   :depth: 1
   :local:

Materials
~~~~~~~~~

The materials used in this part of the tutorial can be found in the `tutorial-part-7`
directory in the psamm-tutorial repository. This directory contains a file called
``gene_associations.tsv`` which contains a two column reciprocal best hits mapping,
mapping the genes in the E coli model model to mock genes from a mock organism
(the mock organism gene names are formatted as "imaginary<Integer>" and have been
randomly generated).

Format of the Reciprocal Best Hits File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The psammotate program requires a reciprocal best hits file. This is essentially
a file that must have two columns (among other potential information):
(1) a list of genes from the organism you are drafting a model for
(2) genes from the reference organism that are mapped to (i.e. share a row with)
    genes from the draft organism based on some annotation
This will allow you to create a model based on the curations of the reference
organism and the annotations of the draft organism based on the gene associations.
These columns need not be in any particular location within a table, as you will
specify the index of the columns for the target and template genes.

If you do not have a gene association for every gene, the genes from the template
model are retained by default. these lines may be simply left blank.

Basic use of the ``psammotate`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run psammotate, you must specify the file containing the gene mapping between the
template and the target model. Additionally, you must specify which columns contain
the genes from the template model and which contain the genes from the target,
or draft model, genes. This will by default generate a new reactions file called
``homolo_reactions.yaml`` in the current directory, that is formatted as a
psamm reactions file and contains the new gene mappings from the draft model.

.. code-block:: shell

    (psamm-env) $ psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2

The output file, ``homolo_reactions.yaml`` contains all of the reactions that
were mapped with new gene annotations. Remember that if there is not gene
annotation in ``gene_associations.tsv`` for a reference gene, it is kept by
default with the gene name of "None". This can also be seen in the standard output:

.. code-block:: shell

    ReactionID	Original_Genes	Translated_Genes	In_final_model
    ACALD	b0351 or b1241	imaginary7180 or imaginary2425	True
    ACALDt	s0001	imaginary1481	True
    ACKr	b3115 or b2296 or b1849	imaginary7287 or imaginary956 or imaginary1755	True
    ACONTa	b0118 or b1276	imaginary4907 or imaginary2569	True
    ACONTb	b0118 or b1276	imaginary4907 or imaginary2569	True
    ACt2r	None	None	True

In this output, the first column is the reaction name, the second is the template
gene name, the third is the target gene name, and the last column indicates if the
gene was imported into ``homolo_reactions.yaml`` (True) or dropped from the model
(False). If the reactions not mapped to should be dropped, use the --ignore-na
option (Note, we cannot overwrite homolo_reactions.yaml, so lets remove it first):

.. code-block:: shell

    (psamm-env) $ rm homolo_reactions.yaml

    (psamm-env) $ psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2 --ignore-na

Note the difference in the output, where the reaction ACt2r is now false and has
not been imported into the new draft model:

.. code-block:: shell

    ReactionID	Original_Genes	Translated_Genes	In_final_model
    ACALD	b0351 or b1241	imaginary7180 or imaginary2425	True
    ACALDt	s0001	imaginary1481	True
    ACKr	b3115 or b2296 or b1849	imaginary7287 or imaginary956 or imaginary1755	True
    ACONTa	b0118 or b1276	imaginary4907 or imaginary2569	True
    ACONTb	b0118 or b1276	imaginary4907 or imaginary2569	True
    ACt2r	None	None	False


Output options
~~~~~~~~~~~~~~

There are several options for output file names/directories besides the default
as well. If you would prefer to not use homolo_reactions.yaml, you can specify your
own prefix using ``--output``, as shown below:

.. code-block:: shell

    (psamm-env) $ psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2 --output draft_reactions

Which will output the ``draft_reactions.yaml`` file instead of the ``homolo_reactions.yaml`` file.


Thermodynamics-Based Metabolic Flux Analysis
---------------------------------------------
The TMFA function in psamm is an implementation of the TMFA algorithm as
detailed in [Henry07]_. This method incoperates additional thermodynamic
constraints into the flux balance framework, allowing for the simulation
of growth, while accounting for the thermodynamic feasibility of the
metabolic reactions. Like the other two methods in this part of the tutorial,
TMFA requires additional data to be prepared beforehand. For details on all
of these input files, see the command line interface section related to the
TMFA command :ref:`commands-tmfa`.

For this tutorial example TMFA data has been provided based based on the
available data from another *E. coli* model in [Henry07]_. Since multiple
files are required to run TMFA, the ``tmfa`` command has been set up to use a
central 'config.yaml' file. This file is then used to specify the relative
paths (from where you are running the program) to the various input files.
This config file is specified through providing the path to the file through
the ``--config`` command line argument.

.. code-block:: shell

  (psamm-env) $ psamm-model tmfa --config ./config.yaml ....

This option allows for tmfa to be set up and run without having to specify
paths to multiple files on the command line every time.


Basic TMFA Input Options
~~~~~~~~~~~~~~~~~~~~~~~~
The ``tmfa`` command contains a few options that can be specified through the
command line to designate things like biomass thresholds and the temperature
that the simulation will be run at.

The first of these options is the ``--threshold`` option. This can be used to
specify a value that the biomass flux will be fixed at during the ``tmfa``
simulations. For example to run a tmfa simulation where the biomass flux is
fixed at 0.5, you can use the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config ./config.yaml --threshold 0.1 simulation

The next option that can be specified is the temperature that will be used for
the simulation. Since temperature is a component of the calculation of the
gibbs free energy of reactions, this parameter can affect the thermodynamics in
the model. The temperature is given in Celsius.

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config ./config.yaml --threshold 0.1 --temp 15 simulation


The next option is related to the use of the error estimates in for the gibbs
free energy of reaction values. For most prediction methods there will be some
uncertainty in the estimation of the gibbs free energy values. This uncertainty
can be incopertated into the ``tmfa`` simulation directly through using the
``--err`` option.

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config ./config.yaml --err simulation


The last general option for the ``tmfa`` command is the ``--hamilton`` option.
This option allows the user to run TMFA with a slightly modified version of the
algorithm that makes all reactions reversible, and only constraints the
reversibility based on thermodynamics. This method is further detailed in the
paper [Hamilton13]_. To run the ``tmfa`` command using this option you can use
the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config config.yaml --hamilton simulation


The TMFA command then contains two sub-commands that can be used for debugging,
``util``, and for running simulations, ``simulation``. To access these subcommands
you can run the ``tmfa`` command like so:

.. code-block:: shell

   (psamm-env) $ psamm-model tmfa --config ./config.yaml util ....

   or

   (psamm-env) $ psamm-model tmfa --config ./config.yaml simulation ....

TMFA util functions
~~~~~~~~~~~~~~~~~~~

The TMFA Utility functions can be accessed through the ``util`` subcommand of
the ``tmfa`` comand. The command contains two utility functions, one is to
generate a template configuration file that can be used when setting up new
models to run TMFA. The option ``--generate-config`` can be used to generate a
template configuration file called example-config.yaml.

.. code-block:: shell

   (psamm-env) $ psamm-model tmfa --config ./config.yaml util --generate-config

The other utility function that is provided is the ``--random-addition``
function. This function can be used to randomly add thermodynamic constratins
to the reactions in the model, and test if the biomass falls below a set
threshold. This process can be used to test out the Gibbs free energy
constraints for a model that is not producing biomass, to see what thermodynamic
constraints might be causing problems.

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config ./config.yaml --threshold 0.1 util --random-addition


Running Growth Simulations with TMFA
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TMFA simulations can be run in two ways. By default the simulation will be run
and will produce Flux Variability-like results that provide upper and lower
bounds for the variables in the TMFA problem. This type of simulation can be
run as follows to simulation growth at maximum biomss production with applying
thermoynamic constraints:

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config ./config.yaml simulation

This command will produce output like the following showing the variable type,
name of the variable, the lower bound, and then upper bound:

.. code-block:: shell

    Flux	CS	0.18355092011356122	0.18355092011362384
    DGR	CS	-42.528916818320724	-1.0000000010279564e-06
    Zi	CS	1.0	1.0
    ....

In this example the variables associated with the Citrate Synthase reaction (CS),
are shown. This simulation shows that the model does use the citrate synthase
reaction and that this reaction is thermodynamically feasible (as indicated by
the negative gibbs free energy value). The 'Flux' range shows
the possible uppper and lower bound of the flux
values for the reactions in the simulation, the 'DGR' range shows the possible
range of gibbs free energy of reaction values for the reaction, and lastly the
'Zi' variable shows the binary constraint variable that is used to constrain
reactions to be on or off based on the thermodynamics. For the 'Zi' variable a 1
indicates that the reaction can carry flux, while 0 indicates that it cannot.

Further down the results, after the reactions have been printed out, the compound
concentrations will be printed out. Similarly they show the compound ID and the
lower and upper bounds of the compound conecentrations. The concentrations are
printed as molar values. Due to the small size of the *E. coli* core model,
most of the metabolites are largely unconstrained and able to vary between the
lower and upper bounds. But a few of the central metabolites do end up being
constrained. For example in this model, many of the metabolites in the central
carbon metabolism are constrained to some extent.

.. code-block:: shell

    CONC	icit_c[c]	9.999e-06	0.007
    CONC	co2_c[c]	9.999e-06	9.999e-05
    CONC	2pg_c[c]	9.999e-06	0.0127


Some additional TMFA simulation options are provided in addition to the
default FVA-like option. The first of these options runs a single FBA-like
``tmfa`` simulation that just provides one solution to the problem without
simulating the variability of the variables. This type of simulation can be run
with the following command:

.. code-block:: shell

    (psamm-env) $ psamm-model tmfa --config config.yaml simulation --single-solution fba

    or

    psamm-model tmfa --config config.yaml simulation --single-solution l1min

These commands will produce just single values for each reactions or compound
instead of providing a range of values. These functions are useful for testing
and debugging, but will miss some of the inherent variability in the simulation.

.. code-block:: shell

    Flux	NH4t_forward	0.9276730532906614
    Flux	NH4t_reverse	0.0
    Flux	O2t_forward	0.0
    Flux	O2t_reverse	0.0
    Flux	PDH	0.0
    Flux	PFK	9.830773843510082
    Flux	PFL	18.235468131213306

The last type of simulation function provided in the ``tmfa`` command is
the randomsparse functions. These commands work in the same way as the
``randomsparse`` function in `PSAMM` and can be used to either do deletions
based on genes or based on reactions.

.. code-block:: shell

    (psamm-model) $ psamm-model tmfa --config config.yaml simulation --randomsparse

    or

    (psamm-model) $ psamm-model tmfa --config config.yaml simulation --randomsparse_genes


Overall the ``tmfa`` function can be used to explore a variety of metabolic
features and provide a way to futher explore the relationships between
metabolic reactions through their thermodynamics.
