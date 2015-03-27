Metabolic modelling tools
=========================

Tools related to metabolic modelling, reconstruction, data parsing and formatting.

Install
-------

The Python module can be installed using `pip`. It is also possible to install a
specific tag or branch by appending the name to the Git URL (e.g. append `@v0.1` to
get the tag `v0.1`). This will typically require `root` permission.

**Note that currently only the Python module and the `model.py` script will be installed
using this method.**

``` shell
$ pip install git+ssh://git@github.com/zhanglab/model_script.git
```

Another option is to use `virtualenv`. First set up a new environment in the project directory,
and activate it.

``` shell
$ virtualenv env
$ . env/bin/activate
```

Now the Python module can be installed in the virtual environment using the `pip` command without
requiring `root` permissions. When returning to the project, simply reactivate the environment by
running the second command.

The Cplex Python bindings will have to be installed manually in the virtual environment as well.
This should be done after activating the virtual environment.

1. Go to the Cplex install directory: `cd /path/to/Cplex`
2. Go to the appropriate subdirectory based on your platform: `cd cplex/python/<platform>`
3. Run `./setup.py install`

Dependencies
------------

- Linear programming solver (_Cplex_, _QSopt_ex_)
- PyYAML (for reading the native model format)
- Numpy (optional; model matrix can be exported to Numpy matrix if available)

There are no strictly required dependencies but most analyses rely on a linear
programming (LP) solver. The LP solver _Cplex_ is the preferred solver. The rational
solver _QSopt_ex_ is also supported (through _python-qsoptex_).

Parsing annotations
-------------------

This script is run to parse the annotation in order to build reaction, compound,
and gene files.

``` shell
$ ~/model_script/make_rxn_table.pl Suden_network \
	ModelSEED_rxns_edited.tsv ModelSEED_cpds_edited.tsv
$ ~/model_script/make_rxn_table.pl Suaut_network \
	ModelSEED_rxns_edited.tsv ModelSEED_cpds_edited.tsv
$ ~/model_script/make_rxn_table.pl Cjr_network \
	ModelSEED_rxns_edited.tsv ModelSEED_cpds_edited.tsv
```

Model format
------------

The primary model definition file is the `model.yaml` file. When creating a new
model this file should be placed in a new directory. The following can be used
as a template:

``` yaml
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
```

### Biomass

The optional **`biomass`** key specifies the default reaction to use for
various analyses (e.g. FBA, FVA, etc.)

### Compounds

The optional **`compounds`** key is a list of compound information. For some
of the model checks the compound information is required. This section can also
include external files that contain compound information. If the file is a
ModelSEED compound table, the `format` key must be set to `modelseed`. If the
file is a YAML file, the file should have a `.yaml` extension. The following
fragment is an example of a YAML formatted compound file:

``` yaml
- id: ac
  name: Acetate
  formula: C2H3O2
  charge: -1

- id: acac
  name: Acetoacetate
  # ...
```

The following compound properties are recognized:

| Property | Type | Description |
|----------|------|------|
| id       | string | Compound ID (_required_) |
| name     | string | Name of compound |
| formula  | string | Compound formula (e.g. C6H12O6) |
| charge   | integer | Formal charge |
| kegg     | string | KEGG ID (reference to compound in KEGG database) |
| cas      | string | CAS number |

### Reactions

The key **`reactions`** specifies a list of files that will be used to define
the reactions in the model. The reaction files can be formatted as either
tab-separated (`.tsv`) or YAML files (`.yaml`). The TSV file may be adequate
for most of the reaction definitions while certain particularly complex
reactions (e.g. biomass reaction) may be specified using a YAML file.

The TSV format is a tab-separated table where each row contains the reaction ID
in addition to other data columns. The header must specify the type of each
column. The column `equation` will be parsed as ModelSEED reaction equations.

```
id      equation
ADE2t   |ade[e]| + |h[e]| <=> |ade[c]| + |ade[c]|
ADK1    |amp| + |atp| <=> (2) |adp|
```

Any `.yaml` or `.yml` file in the `reactions` specification will be parsed as a
reaction definition file but in YAML format. This format is particularly useful
for very long reactions containing many different compounds (e.g. the biomass
reaction). It also allows adding more annotations because of the structured
nature of the YAML format. The following snippet is an example of a YAML
reaction file:

``` yaml
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
```

Reactions in YAML files can also be defined using ModelSEED formatted reaction
equations. The `|` is a special character in YAML so the reaction equations
have to be quoted with `'` or, alternatively, using the `>` for a multiline
quote:

``` yaml
- id: ADE2t
  equation: >
    |ade[e]| + |h[e]| <=>
    |ade[c]| + |h[c]|
- id: ADK1
  equation: '|amp| + |atp| <=> (2) |adp|'
```

The following reaction properties are recognized:

| Property | Type | Description |
|----------|------|------|
| id       | string | Reaction ID (_required_) |
| name     | string | Name of reaction |
| equation | string or dict | Reaction equation formula |
| ec       | string | EC number |
| genes    | list of strings | List of genes associated with the reaction |

### Media

The optional **`media`** key provides a way of defining the medium (boundary
conditions) for the model. The medium is defined by a set of compounds that are
able enter or leave the model system. The following fragment is an example of
the `medium.yaml` file:

``` yaml
- compartment: e
  compounds:
	- id: ac    # Acetate
	- id: co2
	- id: o2
	- id: glcC  # D-Glucose with uptake limit of 10
	  lower: -10
- compartment: c
  compounds:
	- id: compound_x
	  lower: 0  # Provide a sink for compound_x
	# ...
```

When a medium file is specified, the corresponding exchange reactions are
automatically added. For example, if the compounds `o2` in compartment `e` is
in the medium, the exchange reaction `EX_o2_e` is added to the model.

### Reaction flux limits

The optional **`limits`** property lists the files that are to be combined and
applied as the reaction flux limits. This can be used to limit certain
reactions in the model. These files are currently parsed as tables containing
the reaction ID, the lower limit and (optionally) the upper limit. The
following fragment is an example of a limits file:

```
# Make ADE2t irreversible by imposing a lower bound of 0
ADE2t    0
# Only allow limited flux on ADK1
ADK1     -10    10
```

Model tools
-----------

The tools that can be applied to metabolic models are run through the
`model.py` program. To see the full help text of the program use

``` shell
$ model.py --help
```

This program allows you to specify a metabolic model and a command to apply to
the given model. The available commands can be seen using the help command
given above, and are also described in more details below.

To run the program with a model, use

``` shell
$ model.py --model model.yaml command [...]
```

In most cases you will probably be running the command from the same directory
as where the `model.yaml` file is located, and in that case you can simply run

``` shell
$ model.py command [...]
```

To see the help text of a command use

``` shell
$ model.py command --help
```

### Flux balance analysis (`fba`)

This command will first try to maximize the flux of the biomass reaction
defined in the model. It is also possible to provide a different reaction on
the command line to maximize.

By default, this is followed by running the flux balance analysis with
thermodynamic constraints (tFBA) in order to remove internal flux cycles. The
results of these two analyses is presented in a two-column table along with the
reaction IDs.

If the parameter `--no-tfba` is given, the second column instead represents a
flux minimization in which the FBA maximum is fixed while the sum of the fluxes
is minimized. This will often eliminate loops as well.

To run FBA use:

``` shell
$ model.py fba
```

or with a specific reaction:

``` shell
$ model.py fba ATPM
```

### Robustness (`robustness`)

Given a reaction to maximize and a reaction to vary, the robustness analysis
will run flux balance analysis and flux minimization while fixing the reaction
to vary at each iteration. The reaction will be fixed at a given number of
steps between the minimum and maximum flux value specified in the model.

``` shell
$ model.py robustness \
	--steps 200 --minimum -20 --maximum 160 EX_Oxygen
```

In the example above, the biomass reaction will be maximized while the
`EX_Oxygen` (oxygen exchange) reaction is fixed at a certain flux in each
iteration. The fixed flux will vary between the minimum and maximum flux. The
number of iterations can be set using `--steps`. In each iteration, all
reactions and the corresponding fluxes will be shown in a table, as well as
the value of the fixed flux. If the fixed flux results in an infeasible model,
no output will be shown for that iteration.

### Random sparse network (`randomsparse`)

Delete reactions randomly until the flux of the biomass reaction falls below
the threshold. Keep deleting reactions until no more reactions can be deleted.
This can also be applied to other reactions than the biomass reaction by
specifying the reaction explicitly.

``` shell
$ model.py randomsparse 0.95
```

When the given reaction is the biomass reaction, this results in a smaller
model which is still producing biomass within the tolerance given by the
threshold. Aggregating the results from multiple random sparse networks allows
classifying reactions as essential, semi-essential or non-essential.

### Mass consistency check (`masscheck`)

A model or reaction database can be checked for mass inconsistencies. The basic
idea is that we should be able to assign a positive mass to each compound in the
model and have each reaction be balanced with respect to these mass assignments.
If it can be shown that assigning the masses is impossible, we have discovered
an inconsistency.

Some variants of this idea is implemented in the `metnet.massconsistency`
module. The mass consistency check can be run using

``` shell
$ model.py masscheck
```

This will first try to assign a positive mass to as many compounds as possible.
This will indicate whether or not the model is consistent but in case it is
_not_ consistent it is often hard to figure out how to fix the model from this
list of masses.

Next, a different check is run where the residual mass is minimized for all
reactions in the model. This will often give a better idea of which reactions
need fixing.

### Formula consistency check (`formulacheck`)

Similarly, a model or reaction database can be checked for formula
inconsistencies when the chemical formulae of the compounds in the model are
known.

``` shell
$ model.py formulacheck
```

For each inconsistent reaction, the reaction identifier will be printed
followed by the elements ("atoms") in, respectively, the left- and right-hand
side of the reaction, followed by the elements needed to balance the left- and
right-hand side, respectively.

### GapFind/GapFill (`gapfill`)

The GapFind algorithms can be used to identify the compounds that are needed by reactions
in the model but cannot be produced in the model. The GapFill algorithm will extend the
model with reactions from the parent database and try to find a minimal subset that allows
all blocked compounds to be produced. This command will run GapFind to identify the
blocked compounds and then uses GapFill to try to reconstruct a model that allows these
compounds to be produced.

These algorithms are defined in terms of MILP problems and are therefore (particularly
GapFill) computationally expensive to run for larger models.

``` shell
$ model.py gapfill
```

### FastGapFill (`fastgapfill`)

The FastGapFill algorithm tries to reconstruct a flux consistent model (i.e. a model
where every reaction takes a non-zero flux for at least one solutions). This is done by
extending the model with reactions from the parent database and trying to find a minimal
subset that is flux consistent. The solution is approximate.

The database reactions can be assigned a weight (or "cost") using the `--penalty` option.
These weights are taken into account when determining the minimal solution.

``` shell
$ model.py fastgapfill --penalty penalty.tsv
```

### Search (`search`)

This command can be used to search in a database for compounds or reactions. To search
for a compound use

``` shell
$ model.py search compound [...]
```

Use the `--name` option to search for a compound with a specific name or use the
`--id` option to search for a compound with a specific identifier.

To search for a reaction use

``` shell
$ model.py search reaction [...]
```

Use the `--id` option to search for a reaction with a specific identifier. The `--compound`
option can be used to search for reactions that include a specific compound. If more that
one compound identifier is given (comma-separated) this will find reactions that include
all of the given compounds.

### Console (`console`)

This command will start a Python session where the database, model and compounds have
been loaded into the corresponding Python object representation.

``` shell
$ model.py console
```

GAMS models
-----------

In order to use the GAMS models it is necessary to convert the model into a format
that can be loaded by the GAMS models. This can be done using

``` shell
$ ~/model_script/export_gapfind_parameters.py --database database.tsv \
	model.tsv
```

This script will produce the necessary output files that will be read by the GAMS solver.
These files include

* `rxnnames.txt`: a list of all reaction ids
* `rev.txt`: a list of reversible reactions
* `metabolites.txt`: a list of all metabolites
* `cytosol_metabolites.txt`: a list of metabolites found in the cytosol
* `extracellular_metabolites.txt`: a list of metabolites found in the extracellular space
* `modelrxn.txt`: a list of reactions defined by the model
* `databaserxn.txt`: a list of reactions in the database

### GapFind/GapFill

After exporting the GAMS parameter files the GAMS model can be run using

``` shell
$ gams ~/model_script/GAMS/GapFind.gms ps=0
```

The output of the GapFind model can be parsed using

``` shell
$ ~/model_script/parse_gapfind_noproduction.py GapFind.lst
```

This will show a detailed list of the compounds that could not be produced in
the model, and will also write a list to `blocked.txt`.

Next the GapFill problem can be solved which will try to resolve the model
by adding reactions from the database or by reversing irreversible reactions
in the model. Since we may need some transport and exchange reactions to resolve
the model that were not defined in the database, we can add some additional
reactions to the database. `--transport` will allow a reaction for every
compound in the cytosol that brings that compound into the cytosol from the
extracellular space. Similarly, `--exchange` will allow a pseudoreaction for
every compound in the extracellular space which makes the compound readily
available. Lastly, we can add any number of additional database files.

``` shell
$ ~/model_script/export_gapfind_parameters.py --exchange --transport \
	--database database.tsv \
    --database custom_database.tsv model.tsv
```

This command will recreate the output files with the additional reactions
in the database. Now GAMS can be used to solve the GapFill problem.

``` shell
$ gams ~/model_script/GAMS/GapFill.gms ps=0
```

The output will provide a list of added exchange reactions and/or reversed
reactions that resolved the model. We can run the following command to obtain the list
of reactions that had to be added from the database.

``` shell
$ ~/model_script/parse_gapfill.py GapFill.lst
```

### Flux balance analysis

First we need a file listing the exchange reactions that are defined in the
model (`exchangerxn.txt`). Secondly the lower limits of the exchange reactions
are listed in a separate file, `exchangelimit.txt` (tab-separated; reaction id
followed by a number). This file defines the maximum availability of the
compounds in each exchange reaction. By default these reactions will have a
lower flux limit of zero indicating the the reaction does not make any
compounds available.

When running the flux balance analysis the flux of the reaction `Biomass`
will be maximized. This reaction should correspond to the actual biomass
reaction, or alternatively the reaction id can be changed in the GAMS file.

``` shell
$ gams ~/model_script/GAMS/FluxBalance.gms ps=0
```

The resulting flux values can be parsed out using

``` shell
$ ~/model_script/parse_fluxbalance.py FluxBalance.lst rxn_list
```

Test suite
----------

The python modules have test suites that allows us to automatically test
various aspects of the module implementation. It is important to make sure that
all tests run without failure _before_ committing changes to any of the
modules. The test suite is run by changing to the project directory and running

``` shell
$ ./setup.py test
```

Adding or improving tests for python modules is highly encouraged. A test suite for
a new module should be created in `tests/test_<modulename>.py`. These test suites
use the built-in `unittest` module.

In addition, some modules have documentation that can be tested using the `doctest`
module. These test suites should also run without failure before any commits. They
can be run from the `model_script` directory by specifying the particular module
(e.g the `affine` module in `expression`) using

``` shell
$ python -m metnet.expression.affine -v
```
