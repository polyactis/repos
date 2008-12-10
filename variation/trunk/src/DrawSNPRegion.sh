#!/bin/sh
if test $# -lt 1
then
	echo "Usage:"
	echo "    DrawSNPRegion.sh INPUTFILENAME [STEP] [PlotTypeShortName] [COMMIT]"
	echo
	echo "This script calls DrawSNPRegion.py to roll over blocks of entries in input file (to avoid memory blowup)."
	echo "	INPUTFILENAME is gonna fed to -i of DrawSNPRegion.py."
	echo "	STEP is optional, it's the number of lines included in each block. default is 20."
	echo "	PlotTypeShortName is optional, used if the program is gonna commit database transaction."
	echo "	COMMIT is optional, specify 1 to indicate whether the program is gonna commit db transaction."
	echo
	echo "	Output goes to /Network/Data/250k/tmp-yh/snp_region if no db commit."
	echo 
	echo "Example:"
	echo " DrawSNPRegion.sh /Network/Data/250k/tmp-yh/tmp/Flowering_Emma_finished_send_to_yu_2.csv 40 SuziHandPick20081209 1"
	exit
fi

fname=$1
n=`wc -l "$fname"|awk '{print $1}'`	#no of lines in total
content_n=`echo $n-1|bc`	#number of lines excluding the header
if test -n "$2" ; then
	step=$2
else
	step=20
fi
if  test -n "$3" ; then
	plot_type_short_name="-y $3"
else
	plot_type_short_name=""
fi
if [ $4 = "1" ]; then
	commit="-c"
else
	commit=""
fi

m=`echo $n/$step+1|bc`	#how many blocks
echo "working on $fname"

for ((i=0;i<m;i++)); do
	echo $i
	bottom_index=`echo "($i+1)*$step"|bc`	#how many lines far down
	tmp_fname=/tmp/$i\.csv	#temporary file to store the block
	head -n 1 "$fname" >$tmp_fname	#generate a header first
	echo "tail -n $content_n $fname|head -n $bottom_index -|tail -n $step >>$tmp_fname"	#1. remove the header, 2. take everything above the bottom_index, 3 tail the step block 
	tail -n $content_n "$fname"|head -n $bottom_index -|tail -n $step >>$tmp_fname
	#sleep 3s
	~/script/variation/src/DrawSNPRegion.py -i $tmp_fname -I /Network/Data/250k/tmp-yh/call_method_17.tsv -N /Network/Data/250k/tmp-yh/phenotype.tsv -l 28 -o /Network/Data/250k/tmp-yh/snp_region -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -u yh -s -p yh324 $plot_type_short_name $commit
done
    
