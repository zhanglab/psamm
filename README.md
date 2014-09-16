model_script
=================

Scripts related to metabolic reconstruction, data parsing and formating, preparation for GAMS simulation

#-parse the annotation to build rxn, cpd, and gene files
$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suden/Suden_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suaut/Suaut_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Cjr/Cjr_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

GapFind/GapFill
---------------

The scripts in this repository can be used to convert a definition of reactions
of an organism model into a problem that can be solved by GAMS. The organism
definition is simply a list of reaction ids, identifying which reactions from
a larger database that constitute the model. The reaction database is given as
a table of reaction ids and equations. The reaction database can be generated
from multiple sources using the appropriate `XX_format_conversion.py` script,
for example:

``` shell
$ ~/model_script/ModelSEED_format_conversion.py ModelSEED_rxns.tsv \
	ModelSEED_cpds.tsv > ModelSEED_database.tsv
```

The first problem is the GapFind problem which allows us to identify the compounds
that cannot be produced but are consumed by the model. The input files for GapFind
need to be prepared by running

``` shell
$ ~/model_script/export_gapfind_parameters.py --database database.tsv \
	rxn_list
```

In this case the `database.tsv` file could be `ModelSEED_database.tsv` as
prepared previously, and `rxn_list` would be a file containing each reaction id
on a separate line for each reaction included in the model. This script will
produce the necessary output files that will be read by the GAMS solver.

These files include

* `rxnnames.txt`: a list of all reaction ids
* `rev.txt`: a list of reversible reactions
* `metabolites.txt`: a list of all metabolites
* `cytosol_metabolites.txt`: a list of metabolites found in the cytosol
* `extracellular_metabolites.txt`: a list of metabolites found in the extracellular space
* `modelrxn.txt`: a list of reactions defined by the model
* `databaserxn.txt`: a list of reactions in the database

Now GAMS can be used to solve the GapFind problem using the list of parameters
that was just produced.

``` shell
$ gams ~/model_script/GAMS/GapFind.gms ps=0
```

A script is used to parse the output from running GapFind.

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
	--database ModelSEED_database.tsv \
        --database custom_database.tsv rxn_list
```

This command will recreate the output files with the additional reactions
in the database. Now GAMS can be used to solve the GapFill problem.

``` shell
$ gams ~/model_script/GAMS/GapFill.gms ps=0
```

The output will provide a list of added exchange reactions and/or reversed
reactions that resolved the model. We can use a script to obtain the list
of reactions that had to be added from the database.

``` shell
$ ~/model_script/parse_gapfill.py GapFill.lst
```

Flux balance analysis
---------------------

We can use the same files as produced for GapFind and GapFill to actually run
a flux balance analysis on the model. To do this we will first need to define
the boundary conditions of the simulation (i.e. the simulated growth medium).

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

Alternatively, the python implementation of FBA can be used. This will operate
directly on the reaction database and model definition, and we will avoid
having to export the `txt` files entirely. The `exchangelimit.txt` can still be
used to define the exchange reaction limits. This method can be run using

``` shell
$ ~/model_script/run_fluxanalysis.py --database ModelSEED_database.tsv \
        --database custom_database.tsv rxn_list Biomass
```

More advanced analysis and data processing can be done by using the python
module `metnet.fluxanalysis` directly. A demonstration of how to accomplish
this can be seen in the `run_fluxanalysis.py` script.

Robustness analysis
-------------------

Robustness analysis can be run using `run_robustness.py`. The options
that can be given are similar to the ones given to the previously described
programs:

``` shell
$ ~/model_script/run_robustness.py --database ModelSEED_database.tsv \
        --database custom_database.tsv rxn_list Biomass EX_Oxygen
```

In the example above, the `Biomass` reaction will be maximized while the
`EX_Oxygen` (oxygen exchange) reaction is fixed at a certain flux in each
iteration. The fixed flux will vary between the minimum and maximum flux.
The number of iterations can be set using `--steps`. In each iteration,
all reactions and the corresponding fluxes will be printed, as well as
the value of the fixed flux. If the fixed flux results in an infeasible
model, no output will be printed for that iteration.

Mass consistency
----------------

A model or reaction database can be checked for mass inconsistencies. The basic
idea is that we should be able to assign a positive mass to each compound in the
model and have each reaction be balanced with respect to these mass assignments.
If it can be shown that assigning the masses is impossible, we have discovered
an inconsistency.

Some variants of this idea is implemented in the `metnet.massconsistency` module.
The mass consistency check can also be run using

``` shell
$ ~/model_script/run_fluxanalysis.py --database ModelSEED_database.tsv \
        --database custom_database.tsv --compounds ModelSEED_cpds.tsv rxn_list
```

In addition, the chemical formula of compounds can be used to more closely
point out why a reaction is inconsistent. This check can be done using the
following script

``` shell
$ ~/model_script/check_formula_balance.py ModelSEED_database.tsv \
        ModelSEED_cpds.tsv
```

For each inconsistent reaction, the reaction id will be printed followed by
the elements in, respectively, the left- and right-hand side of the reaction,
followed by the elements needed to balance the left- and right-hand side,
respectively.

Test suite
----------

The python modules have test suites that allows us to automatically test various
aspects of the module implemention. It is important to make sure that all tests
run without failure _before_ committing changes to any of the modules. The test
suite is run by changing to the project directory and running

``` shell
$ ./test.py
```

Adding or improving tests for python modules is highly encouraged. A test suite for
a new module should be created in `tests/test_<modulename>.py`. These test suites
use the built-in `unittest` module.

In addition, some modules have documentation that can be tested using the `doctest`
module. These test suites should also run without failure before any commits. They
can be run from the `model_script` directory by specifying the particular module
(e.g the `expression` module) using

``` shell
$ python -m metnet.expression -v
```
