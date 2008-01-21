#!/usr/bin/python

from Bio.WWW import ExPASy
import sys,re
from Bio import File
import getopt

ids = ['O23729', 'O23730', 'O23731']

def usage():
	print "Usage: sprot_search.py <-i> inputfile <-o> outputfile [options]"
	print "\t-i\t--input\t\tspecifiy the input file"
	print "\t-o\t--output\t\tspecify the output file"
	print "\t-h\t--help\t\tshow what you have seen"

def readaccs(inf):
	
	tmp=inf.readline()
	anti_content=re.compile(r'\s')
	accs=[]
	
	while tmp:
		tmp=anti_content.sub('',tmp)
		accs.append(tmp)
		tmp=inf.readline()
	
	return accs
	
def sprot_search(accs, outf=sys.stdout):

	all_results = ''
	i=0
	
	for acc in accs:
	
		outf.write('>'+acc+'\n')
		sys.stderr.write( acc + "\tNo. " + repr(i) +'\n')
		try:
			results=ExPASy.get_sprot_raw(acc)
			
			outf.write(results.read())
		

		except IOError:
			sys.stderr.write( "I/O error: No results \n")

		i=i+1
	

if __name__ == '__main__':

	try :
		opts, args =getopt.getopt(sys.argv[1:],"i:o:h", ["help","input=","output="])
	except getopt.GetoptError:
		print "GetOption Error"
		sys.exit(2)
	for o,a in opts:
		if o in ("-o", "--output"):
			outf=open(a,'w')

		if o in ("-i", "--input"):
			inf = open(a,'r')

		if o in ("-h", "--help"):
			usage()
			sys.exit(3)
	
	
	try:
	
		accs=readaccs(inf)
		sprot_search(accs,outf=outf)
	except:
		print "Error Encountered"
		usage()
		sys.exit(4)
