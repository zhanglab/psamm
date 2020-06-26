Subsetting a Metabolic Models Using Expression Data
====================================================

This tutorial will go over using the ``gimme`` function to subset a metabolic
model and create a new metabolic model based on gene expression data.

.. contents::
   :depth: 1
   :local:

Materials
---------

The materials used in the part of the tutorial can be found in the `tutorial-part-6`
directory in the psamm-tutorial repository. This directory contains a file called
``transcriptome.tsv`` that contains associations between all genes in the E coli model
provided for the earlier tutorials, which will be provided again for this tutorial.
This expression data has been randomly generated from 1-1000 for all genes in the
model. The model file has been resupplied in the ``E_coli_yaml/`` directory and
the expression file is in the additional_files directory.

.. code-block:: shell

   (psamm-env) $ cd <PATH>/tutorial-part-6/E_coli_yaml


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
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 100

As you increase the expression threshold, you will see fewer reactions retained
within the model, which is shown below by using the word count function to count
the number of lines returned by each iteration of gimme:

.. code-block:: shell
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 100 | wc -l
      86
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 500 | wc -l
      74

This application of the gimme algorithm sets the objective biomass to 100%. That
is to say that if a reaction is below the required threshold but is required to
maintain biomass, it is kept. You can lower the required biomass with the
``--biomass-threshold`` flag, which will subsequently drop increasing numbers of
reactions as they are no longer necessary to meet the biomass threshold. This is
shown below by using the word count function to count the number of lines returned
by each iteration of gimme:

.. code-block:: shell
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0.93757758084 | wc -l # Maximum for this model
      86
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 100 --biomass-threshold 0 | wc -l # No Biomass Threshold
      82

In addition to simply writing out a list of reactions that satisfy the expression
and biomass thresholds, by specifying ``--export-model``, you can redirect the
output to create an entirely new model.

.. code-block:: shell
    $ mkdir gimme_out
    $ psamm-model gimme --transcriptome-file ../additional_files/transcriptome.tsv --expression-threshold 100 --export-model ./gimme_out/
