# base calling for arrays

c=0.85	#cutoff probability 
start_array_id=$1
end_array_id=$2
output_dir=$3
for(( i=$start_array_id; i<=$end_array_id; i++ ));
do ./new_basecall.sh $i $c $output_dir;
done

