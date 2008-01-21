#!/usr/bin/python

import tran,sys

if len(sys.argv)>1:
	infile=sys.argv[1]
else:
	infile=raw_input("Enter the sequence file in chromdb's fasta format:\n")

inf=open(infile,'r')

outf=open('outfile.chromdbnr','w')

tran.chromdbnt2SETdb(inf,outf)

inf.close()
outf.close()
