Generating Metabolic Models Based On Functional Annotations
=========================================================================

This tutorial will go over how to use the ``psamm-generate-model`` functions in
`PSAMM`. These functions are designed to allow users to generate draft models
based on basic annotation information. Within ``psamm-generate-model``, there
are three functions:
- `generate-database`: generates a database of metabolic reactions
- `generate-biomass`: generates a draft biomass function
- `generate-transporters`: generates a draft database of transporters

.. contents::
   :depth: 1
   :local:

Materials
---------

For information on how to install `PSAMM` and the associated requirements, as
well how to download the materials required for this tutorial, you can
reference the Installation and Materials section of the tutorial.

In addition to the basic installation of `PSAMM`, these functions utilize an
API to interest with the Kyoto Encyclopedia of Genes and Genomes (KEGG) to
download model information. In order to interact with this API, the `biopython`
package must be installed. The `libchebipy` package, which allows the user to
query curated compound information from chebi.

To use this function, the user will have to generate an annotation file that
contains gene ids in the first column and another column containing an
annotation in one of other columns, which may be specified. For convenience
of the user, the default format for the input of ``psamm-generate-model``
is the default output format for the eggnog mapper annotation software. This
file will be required for `generate-database` and `generate-transporters`

Additional information required are the genome, proteome, and a gff file
for the organism of interest. These files will be required for the
`generate-biomass` function. It is recommended to use prodigal to generate
the proteome and gff file.

For the purposes of this tutorial, these files have been provided for an
organism of interest within the class Mollicutes.

.. code-block:: shell

   (psamm-env) $ cd <PATH>/tutorial-part-7/


Basic use of the ``generate-database`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As noted above, the `generate-database` command uses an annotation file to
query the KEGG database and download reaction data. You will need an
active internet connection for this. There are three types of anntotation that
can be used in order to generate the model (i.e. Reaction number - R00001;
Enzyme commission number - 1.1.1.1; and Kegg Orthology - K00001, specified by
R, EC, and KO, respectively). Reaction number directly queries the Kegg
reaction. Enzyme commission numbers query all reacion numbers associated with
that EC. Kegg orthology numbers query all reaction numbers associated with
that KO.

The most basic usage of this is the following (showing an example for R, EC,
and KO):

.. code-block:: shell

   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type R \
                   --out mollicute_model_R
   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type EC \
                   --out mollicute_model_EC
   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type KO \
                   --out mollicute_model_KO

The output of these commands will be a model compatible with the psamm
program, with a model.yaml, compounds.yaml, and reactions.yaml. Additional
files will be created to ease curation of this draft model, specifically
a log file that tracks compounds and reactions that have probably incorrect
formulations, along with addition compounds and reactions files called
compounds_generic.yaml and reactions_generic.yaml.

Generation of the initial model may take some time (expect 10-20 minutes for
a bacterial genome). If the user is interested in tracking the progress of the
program, it can be run in verbose mode.

.. code-block:: shell

   (psamm-env) $ rm -r mollicute_model_R
   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type R \
                   --out mollicute_model_R --verbose

The above scripts generate the model with a default compartment of "C". If
a different initial compartment is required, use the `--default_compartment`
flag. The following will assign "cyt" as the default compartment

.. code-block:: shell

   (psamm-env) $ rm -r mollicute_model_R
   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type R \
                   --out mollicute_model --verbose --default_compartment cyt

The models output above use the default compound formulation in Kegg based
on the first CHEBI assignment in the KEGG dblinks. In order to generate a
more realistic model, the ability to adjust compound formula based on
protonation states at a typical default pH (7.3) has been added using the
Rhea database. This option is specified using the `--rhea` flag.

.. code-block:: shell

   (psamm-env) $ rm -r mollicute_model_R
   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations mollicute_eggnog.tsv --type R \
                   --out mollicute_model --verbose --rhea

If the user has a custom formulated annotation table, this may also be used
to generate the model. In this case, the gene should be the first column
in the table and the `--col` option can be used to specify the index of the
column in the table specified with `--annotations`

.. code-block:: shell

   (psamm-env) $ psamm-generate-model generate-database \
                   --annotations custom.tsv --type R \
                   --out custom_model --verbose --rhea --col 2


Basic use of the ``generate-transporters`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the database of metabolic reactions, another important
component of metabolic models is the presence of transporters. These
transporters are also predicted in the default eggnog function based on
the classification from the Transporter Classification Database (TCDB). This
Function can be run after the `generate-database` command and generates a new
transporters.yaml file and transporter_log.tsv file. The basic usage is below:

.. code-block:: shell

(psamm-env) $ psamm-generate-model generate-database \
                --annotations mollicute_eggnog.tsv \
                --model mollicute_model

The default compartments for this basic usage are "c" for internal compartment
and "e" for external compartment, but these can be changed with
`--compartment_in` and `--compartment_out`, as in the following:

.. code-block:: shell

(psamm-env) $ psamm-generate-model generate-database \
                --annotations mollicute_eggnog.tsv \
                --model mollicute_model --compartment_in cyt \
                --compartment_out ext

If a custom annotation table is provided, it is handled similarly to in
`generate-database`, where `--col` specifies the index of the column
of the TCDB id.

(psamm-env) $ psamm-generate-model generate-database \
                --annotations custom_transport.tsv \
                --model custom_model

It is also worth noting that the substrate and family information for these
transporters are included in PSAMM as external files; however, if you would
like to use custom annotation tables, these can be provided with
`--db_substrates` and `--db_families`.


Basic use of the ``generate-biomass`` command
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
