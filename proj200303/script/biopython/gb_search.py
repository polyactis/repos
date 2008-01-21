#!/usr/bin/python

import sys,re
from Bio import GenBank
import getopt

def usage():
	print "Usage: gb_search.py -i inputfile  [OPTION]"
	print "inputfile will be one-columned or two-columned"
	print "output to the stdout"
	print "\t-r\t--range=RANGE\tfetch records from specified range eg. 2-20 \n\t\t\t\tif not given, the whole range"
	print "\t-h\t--help\t\tdisplay this help and exit"
	print "written by hy"
	
class gb_search:

	def __init__(self, inf):
		self.inf = inf
		self.accs=[]
		
	def readaccs(self):
		
		tmp=self.inf.readline()
		anti_content=re.compile(r'\s')

		while tmp:
			tmp=anti_content.sub('',tmp)
			self.accs.append(tmp)
			tmp=self.inf.readline()
	
	def reademblaccs(self):
		tmp=self.inf.readline()
		anti_content=re.compile(r'\s')
		
		while tmp:
			tmp=tmp.split()[1]
			tmp=tmp.split('|')
			self.accs=self.accs+tmp
			tmp=self.inf.readline()
		
			
	def search(self,start,end):
	
		ncbi_dict = GenBank.NCBIDictionary()
		j=1
		if start<1:
			start =1	
		if end >len(self.accs) or end == 0:
			end = len(self.accs)
			
		for k in range(start-1,end):
			sys.stderr.write("No " + repr(j) + ": " + self.accs[k] +'\n')
			j=j+1
			gi_list = GenBank.search_for(self.accs[k],database='protein')
		
			for i in range(0,len(gi_list)):
				try:
					gb_record = ncbi_dict[gi_list[i]]
					sys.stdout.write('>'+self.accs[k]+'\n')
					sys.stdout.write( gb_record)
				except:
					sys.stderr.write( self.accs[k] + " fetching error \n")
				

if __name__=='__main__':
	try:
		opts,args = getopt.getopt(sys.argv[1:],"o:i:r:h",["input=","output=","range=", "help"])
	except:
		print "GetOpt Error"
		sys.exit(2)

	start = 0
	end = 0 
	
	for o, a in opts:
		if o in ("-r", "--range"):
			start=int(a.split("-")[0])
			end=int(a.split("-")[1])
		if o in ("-h", "--help"):
			usage()
			sys.exit(2)
		if o in ("-i", "--input"):
			inf = open(a, 'r')
			
	try:
		instance=gb_search(inf)
		instance.reademblaccs()
		instance.search(start,end)
	except:
		print "Error Encounted"
		usage()

