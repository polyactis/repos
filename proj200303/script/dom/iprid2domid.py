#!/usr/bin/python

import re,sys

lmaps={'Nuclear protein SET':'SET',
	'Zn-finger, C2H2 subtype':'zf-C2H2',
	'Zn-finger, C2H2 type':'zf-C2H2',
	'PWWP domain':'PWWP',
	'Zn-finger, RING':'PHD',
	'Zn-finger-like, PHD finger':'PHD',
	'SET-related region':'Post-SET',
	'Nuclear protein Zn2\+-binding':'Pre-SET',
	'High mobility group proteins HMG-I and HMG-Y':'AT_hook',
	'HMG-I and HMG-Y DNA-binding domain \(A\+T-hook\)':'AT_hook',
	'Zn-finger, MYND type':'zf-MYND',
	'Protein of unknown function DUF260':'DUF260',
	'Nuclear protein G9a':'YDG_SRA'}

inf=open(sys.argv[1],'r')
outf=open(sys.argv[2],'w')

tmp=inf.readline()

while(tmp):
	for i in range(len(lmaps)):
		long=re.compile(lmaps.keys()[i])
		
		tmp2=long.sub(lmaps.values()[i],tmp)
		if tmp2!=tmp:
			tmp=tmp2
			break
	outf.write(tmp)
	tmp=inf.readline()

inf.close()
outf.close()
