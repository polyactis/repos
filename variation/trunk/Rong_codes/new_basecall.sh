# flename: new_basecall.sh
# content: script to call genotypes for array with ID = $1 and cutoff 
#  output: tsv format

#array_id
i=$1
#cutoff for confidence measure, i.e., the second column in $i_genocall.txt
c=$2

cp ${i}_raw_data.cel raw_data.cel
R CMD BATCH new_basecall.R 
cp mprobe.norm.rda ${i}_mprobe.norm.rda		#after quantile-normalization
cp genocalls.txt ${i}_genocalls.txt		#called genotypes with probability

cp genocalls.txt calls.txt
./process-calls $i $c
cp calls.tsv ${i}_call.tsv	#base calls in DB tsv format

# output the set of SNPs with at least 3 strains for AA and BB in the training set
#./process-calls-snpset $i $c
#cp calls.tsv ${i}_call.tsv

