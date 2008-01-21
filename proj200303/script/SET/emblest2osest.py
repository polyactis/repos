#!/usr/bin/python

import tran,sys


if len(sys.argv)>1:
	infname=sys.argv[1]
	if len(sys.argv)>2:
		outfname=sys.argv[2]
	else:
		outfname='outfile.2osest'
else:
	infname=raw_input("please enter the files containing paccs:\n")
	outfname='outfile.2osest'


inf=open(infname,'r')

outf=open(outfname,'w')

tran.emblest2osest(inf,outf)

inf.close()
outf.close()
