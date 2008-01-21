#!/bin/bash
cd lindna2ps

for i in $(ls); do
	echo item: $i
	pstopnm $i
 
done
