#!/usr/bin/python

import os,sys,re,cStringIO
from Bio import GenBank
from Bio.WWW import NCBI
import getopt

def usage():
	print "Usage: entrez.py <-i> inputfile <-o> outputfile [OPTION]"
	print "inputfile will be one-columned or two-columned"
	print "\t-p\t--proteinfile=FILE1\tfile to hold protein sequences"
	print "\t-c\t--CDSfile=FILE2\tfile to hold CDS sequences"
	print "\t-r\t--range=RANGE\tfetch records from specified range eg. 2-20"
	print "\t-h\t--help\t\tdisplay this help and exit"
	print "written by hy"
	
class gb_search:

	def __init__(self, inf,proteinfile,CDSfile):
		self.proteinfile=proteinfile
		self.CDSfile=CDSfile
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
	
		search_command = 'Search'
		search_database = 'Protein'
		return_format = 'GenPept'
		Protein_gi_line = re.compile(r'Protein\<\/a\>')
		CDS_gi_line = re.compile(r'CDS\<\/a\>')
		pattern_gi = re.compile(r'val=(\d*)')
		
		j=1
		if start<1:
			start =1	
		if end >len(self.accs) or end == 0:
			end = len(self.accs)
			
		for k in range(start-1,end):
			sys.stderr.write("No " + repr(j) + ": " + self.accs[k] +'\n')
			j=j+1
			try :
				result_handle = NCBI.query(search_command, search_database, term=self.accs[k], doptcmdl = return_format)
				tmp=result_handle.readline()
				sys.stdout.write(tmp)
				while tmp:
					sys.stdout.write(tmp)
					if Protein_gi_line.search(tmp):
						protein_gi=pattern_gi.search(tmp).group(1)
						sys.stderr.write( "fetching "+protein_gi + '\n')
						protein_handle=NCBI.query('Text',search_database, uid=protein_gi,dopt=return_format)
						proteinfile.write(protein_handle.read())
					
					if CDS_gi_line.search(tmp):
						CDS_gi=pattern_gi.search(tmp).group(1)
						#sys.stderr.write( "fetching " + CDS_gi + '\n')
						#CDS_handle=NCBI.query('Text','nucleotide', uid=CDS_gi,dopt='GenBank')
						#CDSfile.write(CDS_handle.read())
					
					tmp = result_handle.readline()
					
				
			except:
				print "Entrez query error"
			
			
				

if __name__=='__main__':
	try:
		opts,args = getopt.getopt(sys.argv[1:],"c:p:o:i:r:h",["input=","output=","range=","proteinfile=","CDSfile=", "help"])
	except:
		print "GetOpt Error"
		sys.exit(2)

	start = 0
	end = 0 
	
	proteinfile=open('entrez.proteinfile','w')
	CDSfile=open('entrez.CDSfile','w')
	
		
	for o, a in opts:
		if o in ("-r", "--range"):
			start=int(a.split("-")[0])
			end=int(a.split("-")[1])
		if o in ("-h", "--help"):
			usage()
			sys.exit(2)
		if o in ("-i", "--input"):
			inf = open(a, 'r')
		if o in ("-p", "--proteinfile"):
			proteinfile=open(a,'w')
		if o in ("-c", "--CDSfile"):
			CDSfile=open(a,'w')
			
	#try:
	instance=gb_search(inf,proteinfile,CDSfile)
		
	instance.reademblaccs()
	instance.search(start,end)
	#except:
	#	print "Error Encounted"
	#	usage()

