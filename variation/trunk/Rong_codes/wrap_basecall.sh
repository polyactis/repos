# base calling for arrays

c=0.85 //cutoff 

for((i = 164; i <= 386; ++i)); do
	./new_basecall.sh $i $c
done

