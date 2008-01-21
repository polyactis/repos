#!/bin/bash
for i in $(ls gp); do
	echo item: $i
	/home/gtkusr/script/pfam2lindna.py gp/$i
	mv outfile.2lindna 2lindna/$i.2lindna
done
