
Importing, Exporting, and working with Models with PSAMM
========================================================

This part of the tutorial will focus on how to use PSAMM to convert
files between the YAML format and other popular formats. An additional
description of the YAML model format and its features is also provided
here.

.. contents::
   :depth: 1
   :local:

Import Functions in PSAMM
_________________________

{{{Reference to the installation tutorial}}}

Importing Existing Models (psamm-import)
________________________________________

In order to work with a metabolic model in PSAMM the model must be converted
to the PSAMM-specific YAML format. This format allows for a human readable
representation of the model components and allows for enhanced customization
with respect to the organization of the metabolic model. This enhanced
organization will allow for a more direct interaction with the metabolic
model and make the model more accessible to both the modeler and experimental
biologists.

Import Formats
~~~~~~~~~~~~~~

``psamm-import`` supports the import of models in various formats. For the SBML
format, it supports the COBRA-compliant SBML specifications, the FBC
specifications, and the basic SBML specifications in levels 1, 2, and 3;
for the JSON format, it supports the import of JSON files directly from the
`BiGG`_ database or from locally downloaded versions;
the support for importing from Excel file is model specific and are available
for 17 published models. There is also a generic Excel import for models
produced by the ModelSEED pipeline. To see a list of these models or model
formats that are supported, use the command:

.. _BiGG: http://bigg.ucsd.edu

.. code-block:: shell

    (psamm-env) $ psamm-import list

In the output, you will see a list of specific Excel models that are supported
by ``psamm-import`` as well as the different SBML parsers that are available in
PSAMM:

.. code-block:: shell

    Generic importers:
    json          COBRA JSON
    modelseed     ModelSEED model (Excel format)
    sbml          SBML model (non-strict)
    sbml-strict   SBML model (strict)

    Model-specific importers:
    icce806       Cyanothece sp. ATCC 51142 iCce806 (Excel format), Vu et al., 2012
    ecoli_textbook  Escerichia coli Textbook (core) model (Excel format), Orth et al., 2010
    ijo1366       Escerichia coli iJO1366 (Excel format), Orth et al., 2011
    gsmn-tb       Mycobacterium tuberculosis GSMN-TB (Excel format), Beste et al., 2007
    inj661        Mycobacterium tuberculosis iNJ661 (Excel format), Jamshidi et al., 2007
    inj661m       Mycobacterium tuberculosis iNJ661m (Excel format), Fang et al., 2010
    inj661v       Mycobacterium tuberculosis iNJ661v (Excel format), Fang et al., 2010
    ijn746        Pseudomonas putida iJN746 (Excel format), Nogales et al., 2011
    ijp815        Pseudomonas putida iJP815 (Excel format), Puchalka et al., 2008
    stm_v1.0      Salmonella enterica STM_v1.0 (Excel format), Thiele et al., 2011
    ima945        Salmonella enterica iMA945 (Excel format), AbuOun et al., 2009
    irr1083       Salmonella enterica iRR1083 (Excel format), Raghunathan et al., 2009
    ios217_672    Shewanella denitrificans OS217 iOS217_672 (Excel format), Ong et al., 2014
    imr1_799      Shewanella oneidensis MR-1 iMR1_799 (Excel format), Ong et al., 2014
    imr4_812      Shewanella sp. MR-4 iMR4_812 (Excel format), Ong et al., 2014
    iw3181_789    Shewanella sp. W3-18-1 iW3181_789 (Excel format), Ong et al., 2014
    isyn731       Synechocystis sp. PCC 6803 iSyn731 (Excel format), Saha et al., 2012

Now the model can be imported using the ``psamm-import`` functions. Return to
the ``psamm-tutorial`` folder if you have left it using the following command:

.. code-block:: shell

    (psamm-env) $ cd <PATH>/psamm-tutorial

Importing an SBML Model
~~~~~~~~~~~~~~~~~~~~~~~

In this tutorial, we will use the `E. coli` textbook core model [1]_ as an
example to demonstrate the functions in PSAMM. The model should be imported
from the SBML model. To import the ``E_coli_core.xml`` model to YAML format run
the following command:

.. code-block:: shell

    (psamm-env) $ psamm-import sbml --source E_coli_sbml/ecoli_core_model.xml --dest E_coli_yaml

This will convert the SBML file in the ``E_coli_sbml`` directory into the YAML
format that will be stored in the ``E_coli_yaml/`` directory. The output will
give the basic statistics of the model and should look like this:

.. code-block:: shell

    ...
    INFO: Detected biomass reaction: R_Biomass_Ecoli_core_w_GAM
    INFO: Removing compound prefix 'M_'
    INFO: Removing reaction prefix 'R_'
    INFO: Removing compartment prefix 'C_'
    Model: Ecoli_core_model
    - Biomass reaction: Biomass_Ecoli_core_w_GAM
    - Compounds: 72
    - Reactions: 95
    - Genes: 137
    INFO: e is extracellular compartment
    INFO: Using default flux limit of 1000.0
    INFO: Converting exchange reactions to medium definition

``psamm-import`` will produce some warnings if there are any aspects of the
model that are going to be changed during import. In this case the warnings are
notifying you that the species with a ``_b`` suffix have been converted to a
boundary condition. You should also see information on whether the biomass
reaction was identified, as well as some basic information on the model name,
size and the default flux settings.

Importing an Excel Model
~~~~~~~~~~~~~~~~~~~~~~~~

The process of importing an Excel model is the same as importing an SBML model
except that you will need to specify the specific model name in the command.
The list of supported models can be seen using the list function above. An
example of an Excel model import is below:

.. code-block:: shell

    (psamm-env) $ psamm-import ecoli_textbook --source E_coli_excel/ecoli_core_model.xls --dest converted_excel_model

This will produce a YAML version of the Excel model in the
``converted_excel_model/`` directory.

Since the Excel models are not in a standardized format these parsers need to
be developed on a model-by-model basis in order to parse all of the relevant
information out of the model. Future support may be added for more Excel-based
models as the parsers are developed.

Importing a JSON Model
~~~~~~~~~~~~~~~~~~~~~~

``psamm-import`` also supports the conversion of JSON format models that follows
the conventions in COBRApy. If the JSON model is stored locally, it can be
converted with the following command:

.. code-block:: shell

    (psamm-env) $ psamm-import json --source E_coli_json/e_coli_core.json --dest converted_json_model/

Alternatively, an extension of the JSON importer has been provided,
``psamm-import-bigg``, which can be applied to convert JSON models from `BiGG`_
database. To see the list of available models on the BiGG database the
following command can be used:

.. code-block:: shell

    (psamm-env) $ psamm-import-bigg list

This will show the available models as well as their names. You can then
import any of these models to YAML format. For example, using the following
command to import the `E. coli` iJO1366 [2]_ model from the BiGG database:

.. code-block:: shell

    (psamm-env) $ psamm-import-bigg iJO1366 --dest converted_json_model_bigg/

.. note::
    To use ``psamm-import-bigg`` you must have internet access so
    that the models can be downloaded from the online BiGG database.

YAML Format and Model Organization
__________________________________

The PSAMM YAML format stores individual models under a designated directory,
in which there will be a number of files that stores the information of the
model and specifies the simulation conditions. The entry point of the YAML
model is a file named ``model.yaml``, which points to additional files that
store the information of the model components, including compounds, reactions,
flux limits, medium conditions, etc. While we recommend that you use the name
``model.yaml`` for the central reference file, the file names for the included
files are flexible and can be customized as you prefer. In this tutorial, we
simply used the names: ``compounds.yaml``, ``reactions.yaml``, ``limits.yaml``,
and ``medium.yaml`` for the included files.

First change directory into ``E_coli_yaml``:

.. code-block:: shell

    (psamm-env) $ cd E_coli_yaml/

The directory contains the main ``model.yaml`` file as well as the included
files:

.. code-block:: shell

    (psamm-env) $ ls
    compounds.yaml
    limits.yaml
    medium.yaml
    model.yaml
    reactions.yaml

These files can be opened using any standard text editor. We highly recommend
using an editor that includes syntax highlighting for the YAML language (we
recommend the Atom_ editor which includes built-in support for YAML and is
available for macOS, Linux and Windows). You can also use a command like
``less`` to quickly inspect the files:

.. _Atom: https://atom.io/

.. code-block:: shell

    (psamm-env) $ less <file_name>.yaml

The central file in this organization is the ``model.yaml`` file. The following
is an example of the ``model.yaml`` file that is obtained from the import of
the `E. coli` textbook model. The ``model.yaml`` file for this imported SBML
model should look like the following:

.. code-block:: yaml

    name: Ecoli_core_model
    biomass: Biomass_Ecoli_core_w_GAM
    default_flux_limit: 1000.0
    compounds:
    - include: compounds.yaml
    reactions:
    - include: reactions.yaml
    media:
    - include: medium.yaml
    limits:
    - include: limits.yaml

The ``model.yaml`` file defines the basic components of a metabolic model,
including the model name (`Ecoli_core_model`), the biomass function
(`Biomass_Ecoli_core_w_GAM`), the compound files (``compounds.yaml``), the
reaction files (``reactions.yaml``), the flux boundaries (``limits.yaml``), and
the medium conditions (``medium.yaml``). The additional files are defined using
include functions. This organization allows you to easily change
aspects of the model like the exchange reactions by simply referencing a
different media file in the central ``model.yaml`` definition.

This format can also be used to include multiple files in the list of
reactions and compounds. This feature can be useful, for example, if you
want to name different reaction files based on the subsystem designations or
cellular compartments, or if you want to separate the temporary reactions
that are used to fill reaction gaps from the main model. An example of how you
could designate multiple reaction files is found below. This file can be found
in the additional files folder in the file ``complex_model.yaml``.

.. code-block:: yaml

    name: Ecoli_core_model
    biomass: Biomass_Ecoli_core_w_GAM
    default_flux_limit: 1000.0
    model:
    - include: core_model_definition.tsv
    compounds:
    - include: compounds.yaml
    reactions:
    - include: reactions/cytoplasm.yaml
    - include: reactions/periplasm.yaml
    - include: reactions/transporters.yaml
    - include: reactions/extracellular.yaml
    media:
    - include: medium.yaml
    limits:
    - include: limits.yaml


As can be see here the modeler chose to distribute their reaction database
files into different files representing various cellular compartments and roles.
This organization can be customized to suit your preferred workflow.

There are also situations where you may wish to designate only a subset
of the reaction database in a metabolic simulation. In these situations it is
possible to use a model definition file to identify which subset of reactions
will be used from the larger database. The model definition file is simply a
list of reaction IDs that will be included in the simulation.

An example of how to include a model definition file can be found below.

.. code-block:: yaml

    name: Ecoli_core_model
    biomass: Biomass_Ecoli_core_w_GAM
    default_flux_limit: 1000.0
    model:
    - include: subset.tsv
    compounds:
    - include: compounds.yaml
    reactions:
    - include: reactions.yaml
    media:
    - include: medium.yaml
    limits:
    - include: limits.yaml

.. note::
    When the model definition file is not identified, PSAMM will include
    the entire reaction database in the model. However, when it is identified,
    PSAMM will only include the reactions that are listed in the model
    definition file in the model. This design can be useful when you want to
    make targeted tests on a subset of the model or when you want to include a
    larger database for use with the gap filling functions.

Reactions
~~~~~~~~~

The ``reactions.yaml`` file is where the reaction information is stored in the
model. A sample of this kind of file can be seen below:

.. code-block:: yaml

    - id: ACALDt
      name: acetaldehyde reversible transport
      genes: s0001
      equation: '|acald_e[e]| <=> |acald_c[c]|'
      subsystem: Transport, Extracellular
    - id: ACKr
      name: acetate kinase
      genes: b3115 or b2296 or b1849
      equation: '|ac_c[c]| + |atp_c[c]| <=> |actp_c[c]| + |adp_c[c]|'
      subsystem: Pyruvate Metabolism

Each reaction entry is designated with the reaction ID first. Then the various
properties of the reaction can be listed below it. The required properties for
a reaction are ID and equation. Along with these required attributes others
can be included as needed in a specific project. These can include but are not
limited to EC numbers, subsystems, names, and genes associated with the
reaction. For example, in a collaborative reconstruction you may want to
include a field named ``authors`` to identify which authors have contributed to
the curation of a particular reaction.

Reaction equations can be formatted in multiple ways to allow for more
flexibility during the modeling process. The reactions can be formatted in a
string format based on the ModelSEED reaction format. In this representation
individual compounds in the reaction are represented as compound IDs followed by
the cellular compartment in brackets, bordered on both sides by single pipes.
For example if a hydrogen compound in the cytosol was going to be in an equation
it could be represented as follows:

.. code-block:: shell

    |Hydr[cytosol]|

These individual compounds can be assigned stoichiometric coefficients by
adding a number in parentheses before the compound. For example if a reaction
contained two hydrogens it could appear as follows:

.. code-block:: shell

    (2) |Hydr[cytosol]|

These individual components are separated by + signs in the reaction string. The
separation of the reactants and products is through the use of an equal sign
with greater than or less than signs designating directionality. These could
include => or <= for reactions that can only progress in one direction or <=>
for reactions that can progress in both directions. An example of a correctly
formatted reaction could be as follows:

.. code-block:: shell

    '|ac_c[c]| + |atp_c[c]| <=> |actp_c[c]| + |adp_c[c]|'

For longer reactions the YAML format
provides a way to list each reaction component on a single line. For example a
reaction could be represented as follows:

.. code-block:: yaml

    - id: ACKr
      name: acetate kinase
      equation:
        compartment: c
        reversible: yes
        left:
          - id: ac_c
            value: 1
          - id: atp_c
            value: 1
        right:
          - id: actp_c
            value: 1
          - id: adp_c
            value: 1
      subsystem: Pyruvate Metabolism

This line based format can be especially helpful when dealing with larger
equations like biomass reactions where there can be dozens of components in
a single reaction.

Gene associations for the reactions in a model can also be included in the
reaction definitions so that gene essentiality experiments can be performed
with the model. These genes associations are included by adding the ``genes``
property to the reaction like follows:

.. code-block:: yaml

    - id: ACALDt
      name: acetaldehyde reversible transport
      equation: '|acald_e[e]| <=> |acald_c[c]|'
      subsystem: Transport, Extracellular
      genes: gene_0001


More complex gene associations can also be included by using logical and/or
statements in the genes property. When performing gene essentiality simulations
this logic will be taken into account. Some examples of using this logic with
the genes property can be seen below:

.. code-block:: yaml

    genes: gene_0001 or gene_0002

    genes: gene_0003 and gene_0004

    genes: gene_0003 and gene_0004 or gene_0005 and gene_0006

    genes: gene_0001 and (gene_0002 or gene_0003)


Compounds
~~~~~~~~~

The ``compounds.yaml`` file is organized in the same way as the
``reactions.yaml``. An example can be seen below.

.. code-block:: yaml

    - id: 13dpg_c
      name: 3-Phospho-D-glyceroyl-phosphate
      formula: C3H4O10P2
    - id: 2pg_c
      name: D-Glycerate-2-phosphate
      formula: C3H4O7P
    - id: 3pg_c
      name: 3-Phospho-D-glycerate
      formula: C3H4O7P

The compound entries begin with a compound ID which is then followed by the
compound properties. These properties can include a name, chemical formula,
and charge of the compound.

Limits
~~~~~~

The limits file is used to designate reaction flux limits when it is different
from the defaults in PSAMM. By default, PSAMM would assign the lower and
upper bounds to reactions based on their reversibility, i.e. the boundary of
reversible reactions are :math:`-1000 \leq v_j \leq 1000`, and the boundary for
irreversible reactions are :math:`0 \leq v_j \leq 1000`. Therefore, the
``limits.yaml`` file will consist of only the reaction boundaries that are
different from these default values. For example, if you want to force flux
through an artificial reaction like the ATP maintenance reaction `ATPM` you
can add in a lower limit for the reaction in the limits file like this:

.. code-block:: yaml

    - reaction: ATPM
      lower: 8.39

Each entry in the limits file includes a reaction ID followed by upper and
lower limits. Note that when a model is imported only the non-default flux
limits are explicitly stated, so some of the imported models will not contain
a predefined limits file. In the `E. coli` core model, only one reaction has a
non-default limit. This reaction is an ATP maintenance reaction and the
modelers chose to force a certain level of flux through it to simulate the
general energy cost of cellular maintenance.

Medium
~~~~~~

The medium file is where you can designate the boundary conditions for the
model. The compartment of the medium compounds can be designated using the
``compartment`` tag, and if omitted, the extracellular compartment (`e`) will
be assumed. An example of the medium file can be seen below.

.. code-block:: yaml

    name: Default medium
    compounds:
    - id: ac_e
      reaction: EX_ac_e
      lower: 0.0
    - id: acald_e
      reaction: EX_acald_e
      lower: 0.0
    - id: akg_e
      reaction: EX_akg_e
      lower: 0.0
    - id: co2_e
      reaction: EX_co2_e


Each entry starts with the ID of the boundary compound and followed by lines
that defines the lower and upper limits of the compound flux. Internally,
PSAMM will translate these boundary compounds into exchange reactions in
metabolic models. Additional properties can be designated for the exchange
reactions including an ID for the reaction, the compartment for the reaction,
and lower and upper flux bounds for the reaction. In the same way that only
non-standard limits need to be specified in the limits file, only non-standard
exchange limits need to be specified in the media file.


Model Format Customization
~~~~~~~~~~~~~~~~~~~~~~~~~~

The YAML model format is highly customizable to suit your preferences.
File names can be changed according to your own design. These customizations
are all allowed by PSAMM as long as the central ``model.yaml`` file is also
updated to reflect the different file names referred. While all the file names
can be changed it is recommended that the central ``model.yaml`` file name does
not change. PSAMM will automatically detect and read the information from the
file if it is named ``model.yaml``. If you *do* wish to also alter the name of
this file you can do so but whenever any PSAMM commands are run you will need
to specify the path of your model file using the ``--model`` option. For
example, to run FBA with a different central model file named
``ecoli_model.yaml``, you could run the command like this:

.. code-block:: shell

    (psamm-env) $ psamm-model --model ecoli_model.yaml fba

