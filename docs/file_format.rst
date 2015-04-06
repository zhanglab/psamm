
Model file format
=================

The primary model definition file is the ``model.yaml`` file. When creating a
new model this file should be placed in a new directory. The following can be
used as a template:

.. code-block:: yaml

    ---
    name: Escherichia coli test model
    biomass: Biomass
    compounds:
      - include: ../path/to/ModelSEED_cpds.tsv
        format: modelseed
    reactions:
      - include: reactions/reactions.tsv
      - include: reactions/biomass.yaml
    media:
      - include: medium.yaml

Biomass
-------

The optional ``biomass`` key specifies the default reaction to use for
various analyses (e.g. FBA, FVA, etc.)

Compounds
---------

The optional ``compounds`` key is a list of compound information. For some
of the model checks the compound information is required. This section can also
include external files that contain compound information. If the file is a
ModelSEED compound table, the ``format`` key must be set to ``modelseed``. If
the file is a YAML file, the file should have a ``.yaml`` extension. The
following fragment is an example of a YAML formatted compound file:

.. code-block:: yaml

    - id: ac
      name: Acetate
      formula: C2H3O2
      charge: -1

    - id: acac
      name: Acetoacetate
      # ...

The following compound properties are recognized:

========  =======  ================================================
Property  Type     Description
========  =======  ================================================
id        string   Compound ID (*required*)
name      string   Name of compound
formula   string   Compound formula (e.g. C6H12O6)
charge    integer  Formal charge
kegg      string   KEGG ID (reference to compound in KEGG database)
cas       string   CAS number
========  =======  ================================================

Reactions
---------

The key ``reactions`` specifies a list of files that will be used to define
the reactions in the model. The reaction files can be formatted as either
tab-separated (``.tsv``) or YAML files (``.yaml``). The TSV file may be
adequate for most of the reaction definitions while certain particularly
complex reactions (e.g. biomass reaction) may be specified using a YAML file.

The TSV format is a tab-separated table where each row contains the reaction ID
in addition to other data columns. The header must specify the type of each
column. The column ``equation`` will be parsed as ModelSEED reaction equations.

::

    id      equation
    ADE2t   |ade[e]| + |h[e]| <=> |ade[c]| + |ade[c]|
    ADK1    |amp| + |atp| <=> (2) |adp|

Any ``.yaml`` or ``.yml`` file in the ``reactions`` specification will be
parsed as a reaction definition file but in YAML format. This format is
particularly useful for very long reactions containing many different compounds
(e.g. the biomass reaction). It also allows adding more annotations because of
the structured nature of the YAML format. The following snippet is an example
of a YAML reaction file:

.. code-block:: yaml

    # Biomass composition
    - id: Biomass
      equation:
      reversible: no
      left:
        - id: cpd00032 # Oxaloacetate
          value: 1
        - id: cpd00022 # Acetyl-CoA
          value: 1
        - id: cpd00035 # L-Alanine
          value: 0.02
        # ...
      right:
        - id: Biomass
          value: 1
        # ...

Reactions in YAML files can also be defined using ModelSEED formatted reaction
equations. The ``|`` is a special character in YAML so the reaction equations
have to be quoted with ``'`` or, alternatively, using the ``>`` for a multiline
quote:

.. code-block:: yaml

    - id: ADE2t
      equation: >
        |ade[e]| + |h[e]| <=>
        |ade[c]| + |h[c]|
    - id: ADK1
      equation: '|amp| + |atp| <=> (2) |adp|'

The following reaction properties are recognized:

========  ===============  ==========================================
Property  Type             Description
========  ===============  ==========================================
id        string           Reaction ID (*required*)
name      string           Name of reaction
equation  string or dict   Reaction equation formula
ec        string           EC number
genes     list of strings  List of genes associated with the reaction
========  ===============  ==========================================

Media
-----

The optional ``media`` key provides a way of defining the medium (boundary
conditions) for the model. The medium is defined by a set of compounds that are
able enter or leave the model system. The following fragment is an example of
the ``medium.yaml`` file:

.. code-block:: yaml

    compartment: e  # default compartment
    compounds:
      - id: ac      # Acetate
      - id: co2
      - id: o2
      - id: glcC    # D-Glucose with uptake limit of 10
        lower: -10
      - id: compound_x
        compartment: c
        lower: 0    # Provide a sink for compound_x
      # ...

When a medium file is specified, the corresponding exchange reactions are
automatically added. For example, if the compounds ``o2`` in compartment ``e``
is in the medium, the exchange reaction ``EX_o2_e`` is added to the model.

Reaction flux limits
--------------------

The optional ``limits`` property lists the files that are to be combined and
applied as the reaction flux limits. This can be used to limit certain
reactions in the model. These files are currently parsed as tables containing
the reaction ID, the lower limit and (optionally) the upper limit. The
following fragment is an example of a limits file::

    # Make ADE2t irreversible by imposing a lower bound of 0
    ADE2t    0
    # Only allow limited flux on ADK1
    ADK1     -10    10
