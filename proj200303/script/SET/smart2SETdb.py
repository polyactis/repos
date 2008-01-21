#!/usr/bin/python

import tran,sys

if len(sys.argv)>1:
	infile=sys.argv[1]
else:
	infile=raw_input("Enter the sequence file in fasta format:\n")

inf=open(infile,'r')

outf=open('outfile.smart','w')

tran.smart2SETdb(inf,outf)

inf.close()
outf.close()
