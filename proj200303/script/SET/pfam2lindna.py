#!/usr/bin/python

import lindnaface
import sys

if len(sys.argv)>1:
	infile=sys.argv[1]
else:
	infile=raw_input("please enter the Pfam file:\n")

inf=open(infile,'r')

outf=open('outfile.2lindna','w')		
	
gl=lindnaface.infile(inf)

m=lindnaface.mark(gl)

for i in range(len(m)-1):
	lindnaface.sort(gl,m[i]+1,m[i+1])

lindnaface.outfile(outf,gl)

inf.close()
outf.close()

