#!/bin/bash
for i in $(ls 2lindna); do
	echo item: $i
	lindna -inputfile 2lindna/$i -graphout cps -ruler -blocktype filled -textheight 2
	mv lindna.ps lindna2ps/$i.ps
done
