#!/bin/bash 

    #PBS -l walltime=99:00:00
    #PBS -l nice=10
    #PBS -q default
    #PBS -l nodes=1:ppn=20
#    #PBS -t 1-40
echo "Starting Run"
date

#mpiexec hostname
#cd ../pan/pan_script
/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suden/Suden_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv 

/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Suaut/Suaut_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

/home/ying/model/model_script/make_rxn_table.pl /home/ying/model/Cjr/Cjr_network /home/ying/db/ModelSEED_04032014/ModelSEED_rxns_edited.tsv /home/ying/db/ModelSEED_04032014/ModelSEED_cpds_edited.tsv

echo "END"
date
