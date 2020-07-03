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
-----------------------------------

The Gimme algorithm functionally subsets a model based on expression data provided
using the ``--transcriptome-file`` argument. This file should contain two columns:
one with each gene within a model and one with a measure of the expression of the
gene, such as transcripts per million (TPM) or Reads Per Kilobase of transcript,
per Million mapped reads (RPKM). For the purposes of this tutorial, expression
data has been mocked up by randomly generating expression numbers for genes in
the E coli model.

Basic use of the ``Gimme`` command
-----------------------------------

To run Gimme, you must first specify the directory of the two column transcriptome
file and specify a threshold value that defines at what threshold below which
genes will be dropped from the new model. This produces a new model that drops
any reactions coded by genes that do not meet this threshold. Any reactions assoociated
with multiple genes must meet the threshold across all. Running the command as follows
will print out a list of the internal reactions and a list of the external reaction
that meet the thresholds.

.. code-block:: shell
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100

As you increase the expression threshold, you will see fewer reactions retained
within the model, which is shown below by using the word count function to count
the number of lines returned by each iteration of gimme:

.. code-block:: shell
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 | wc -l
      86
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 500 | wc -l
      74

This application of the gimme algorithm sets the objective biomass to 100%. That
is to say that if a reaction is below the required threshold but is required to
maintain biomass, it is kept. You can lower the required biomass with the
``--biomass-threshold`` flag, which will subsequently drop increasing numbers of
reactions as they are no longer necessary to meet the biomass threshold. This is
shown below by using the word count function to count the number of lines returned
by each iteration of gimme:

.. code-block:: shell
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0.93757758084 | wc -l # Maximum for this model
      86
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0 | wc -l # No Biomass Threshold
      82

In addition to simply writing out a list of reactions that satisfy the expression
and biomass thresholds, by specifying ``--export-model``, you can redirect the
output to create an entirely new model.

.. code-block:: shell
    $ mkdir gimme_out
    $ psamm-model gimme --transcriptome-file ../gimme/transcriptome.tsv --expression-threshold 100 --export-model ./gimme_out/


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
---------

The materials used in this part of the tutorial can be found in the `tutorial-part-7`
directory in the psamm-tutorial repository. This directory contains a file called
``gene_associations.tsv`` which contains a two column reciprocal best hits mapping,
mapping the genes in the E coli model model to mock genes from a mock organism
(the mock organism gene names are formatted as "imaginary<Integer>" and have been
randomly generated).

Format of the Reciprocal Best Hits File
----------------------------------------

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
----------------------------------------

To run psammotate, you must specify the file containing the gene mapping between the
template and the target model. Additionally, you must specify which columns contain
the genes from the template model and which contain the genes from the target,
or draft model, genes. This will by default generate a new reactions file called
``homolo_reactions.yaml`` in the current directory, that is formatted as a
psamm reactions file and contains the new gene mappings from the draft model.

.. code-block:: shell

    $ psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2

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

    $ rm homolo_reactions.yaml
    $ psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2 --ignore-na

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
---------------

There are several options for output file names/directories besides the default
as well. If you would prefer to not use homolo_reactions.yaml, you can specify your
own prefix using ``--output``, as shown below:

.. code-block:: shell

    $  psamm-model psammotate --rbh ../psammotate/gene_associations.tsv --template 1 --target 2 --output draft_reactions

Which will output the ``draft_reactions.yaml`` file instead of the ``homolo_reactions.yaml`` file.
