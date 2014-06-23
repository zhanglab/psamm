model_script
=================

Scripts related to metabolic reconstruction, data parsing and formating, preparation for GAMS simulation

#-parse the annotation to build rxn, cpd, and gene files
$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suden/Suden_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv 

$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suaut/Suaut_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

$/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Cjr/Cjr_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

GapFind/GapFill
---------------

The scripts in this repository can be used to convert a table of reactions
of an organism model into a problem that can be solved by GAMS. The first problem
is the GapFind problem which allows us to identify the compounds that cannot be
produced but are consumed by the model.

``` shell
$ ~/model_script/export_gapfind_parameters.py rxn_table.csv
```

This script will produce the necessary output files that will be read by the
GAMS solver. These files include a list of all reaction ids (`rxnnames.txt`);
a list of reversible reactions (`rev.txt`); a list of all metabolites
(`metabolites.txt`); a list of metabolites found in the cytosol
(`cytosol_metabolites.txt`); a list of metabolites found in the extracellular
space (`extracellular_metabolites.txt`); a list of reactions defined by the model
(`modelrxn.txt`); and a list of reactions that should be considered
(`databaserxn.txt`). In the simple case, the last two files will be a list of
all the reactions and an empty list, respectively.

Now GAMS can be used to solve the GapFind problem using the list of parameters
that was just produced.

``` shell
$ gams ~/model_script/GapFind.gms ps=0
```

A script is used to parse the output from running GapFind.

``` shell
$ ~/model_script/parse_gapfind_noproduction.py GapFind.lst cpd_table.csv
```

This will show a detailed list of the compounds that could not be produced in
the model, and will also write a list to `blocked.txt`.

Next the GapFill problem can be solved which will try to resolve the model
by adding reactions from the database or by reversing irreversible reactions
in the model. Since the database (`databaserxn.txt`) we produced previously
was empty this will usually not be able to provide an interesting result at
this time. We can run the export script again and tell it to put some reactions
into the database this time.

``` shell
$ ~/model_script/export_gapfind_parameters.py --exchange --transport rxn_table.csv
```

This command will add exchange and transport reactions to the database.
`--transport` will allow a reaction for every compound in the cytosol that
brings that compound into the cytosol from the extracellular space.
Similarly, `--exchange` will allow a pseudoreaction for every compound in the
extracellular space which makes the compound readily available.

Now GAMS can be used to solve the GapFill problem.

``` shell
$ gams ~/model_script/GapFill.gms ps=0
```

The output will provide a list of added exchange reactions and/or reversed
reactions that resolved the model.
