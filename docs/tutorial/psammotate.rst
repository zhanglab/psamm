Generating a Draft Model From a Template
=========================================

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

    $ psamm-model psammotate --rbh ../additional_files/gene_associations.tsv --template 1 --target 2

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
    $ psamm-model psammotate --rbh ../additional_files/gene_associations.tsv --template 1 --target 2 --ignore-na

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

    $  psamm-model psammotate --rbh ../additional_files/gene_associations.tsv --template 1 --target 2 --output draft_reactions

Which will output the ``draft_reactions.yaml`` file instead of the ``homolo_reactions.yaml`` file.
