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
