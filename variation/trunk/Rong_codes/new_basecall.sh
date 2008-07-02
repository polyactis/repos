# flename: new_basecall.sh
# content: script to call genotypes for array with ID = $1 and cutoff 
#  output: tsv format

#array_id
i=$1
#cutoff for confidence measure, i.e., the second column in $i_genocall.txt
c=$2
output_dir=$3

echo -n "Array " $i
dir_prefix="/Network/Data/250k/db/raw_data/"
call_output="$output_dir/call_output_$c"
other_output="$output_dir/other_output"
mkdirhier $call_output
mkdirhier $other_output
rm raw_data.cel
ln -s $dir_prefix$i\_raw_data.cel raw_data.cel
R CMD BATCH new_basecall.R 
cp mprobe_mean.rda $other_output/${i}_mprobe_mean.rda
cp mprobe_norm.rda $other_output/${i}_mprobe_norm.rda		#after quantile-normalization
cp genocalls.txt $other_output/${i}_genocalls.txt		#called genotypes with probability
echo ""

cp genocalls.txt calls.txt
./process-calls $i $c
cp calls.tsv $call_output/${i}_call.tsv	#base calls in DB tsv format

# output the set of SNPs with at least 3 strains for AA and BB in the training set
#./process-calls-snpset $i $c
#cp calls.tsv ${i}_call.tsv

echo ""
