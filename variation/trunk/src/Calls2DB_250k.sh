#!/bin/sh

if test $# -lt 1
then
	echo "Usage:"
	echo "    Calls2DB_250k.sh input_dir"
	echo
	echo " This program wraps around Calls2DB_250k.py."
	echo
	exit
fi
input_dir=$1

for filename in `ls $input_dir`
do
	input_fname=$input_dir/$filename
	while (true)
	do
		echo -n "Enter the strain name for $input_fname: "
		read strain_name
		echo ~/script/variation/src/Calls2DB_250k.py -i $input_fname -n $strain_name  -e "$input_fname" -c
		~/script/variation/src/Calls2DB_250k.py -i $input_fname -n $strain_name  -e "$input_fname" -c
		echo -n "Press Enter for next file or type any character to retry: "
		read retry
		if test -z $retry
		then
			break
		else
			continue
		fi
	done
done
