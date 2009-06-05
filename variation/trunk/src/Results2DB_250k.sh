#!/bin/bash

if test $# -ne 6
then
	echo "Usage:"
	echo "    Results2DB_250k.sh FILE_PREFIX CALL_METHOD_ID ANALYSIS_METHOD_ID DB_PASSWORD START_PHENOTYPE_ID END_PHENOTYPE_ID"
	echo
	echo "This script calls Results2DB_250k.py to submit results into database."
	echo
	echo "Example:"
	echo " Results2DB_250k.sh /Network/Data/250k/tmp-bvilhjal/marg_results/Marg_newDataset_ 17 5 secret 1 180"
	exit
fi

file_prefix=$1
call_method_id=$2
analysis_method_id=$3
db_password=$4
start_phenotype_id=$5
end_phenotype_id=$6

for((i=$start_phenotype_id;i<=$end_phenotype_id;i++))
	do echo $i;
	~/script/variation/src/Results2DB_250k.py -a $call_method_id -e $i -i $file_prefix$i\_*[^s][^r].tsv -l $analysis_method_id -u yh -p $db_password -c
done

