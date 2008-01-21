#!/usr/bin/python

import sys,re

class extract:

	def __init__(self,infname,outfname):

		self.inf=open(infname,'r')
		self.outf=open(outfname,'w')
		
	def run(self):
	
		counter=0
		iswanted=0
		record=''
		
		line=self.inf.readline()

		organism=re.compile(r'^Organism')
		species=re.compile(r'Saccharum')
		tail=re.compile(r'^\|\|')

		while line:
		
			record=record+line
			
			if organism.match(line) and species.search(line):
				iswanted=1
				counter=counter+1
			
			if tail.match(line): 
				if	iswanted==1:
					self.outf.write(record)
					record=''
					iswanted=0
				else:
					record=''
					iswanted=0

			line=self.inf.readline()


if __name__=='__main__':

	instance=extract(sys.argv[1],sys.argv[2])
	instance.run()

