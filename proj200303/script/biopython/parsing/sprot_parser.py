#!/usr/bin/python

import SProt
from Bio import File
import sys,getopt,string

def usage():
	print "Usage sprot_parser.py -i inputfile  [OPTION]"
	print "\t-u\t--updateaccs=FILE1\tfile to contain updated accs"
	print "\t-e\t--emblaccs=FILE2\tfile to contain referenced EMBL accs"
	print "\t-h\t--help\t\tdisplay this help information and exit"

try:
	opts, args =getopt.getopt(sys.argv[1:], "i:o:e:u:h", ["help","updateaccs=","emblaccs=","output=","input="])
except getopt.GetoptError:
	print "GetOption Error"
	sys.exit(2)

for o , a in opts:
	if o in ("-i", "--input"):
		s_handle = open(a, 'r')
	if o in ("-h","--help"):
		usage()
		sys.exit(3)
		pass
	if o in ("-u","--updateaccs"):
		updateaccs=open(a,'w')
	if o in ("-e","--emblaccs"):
		emblaccs = open( a,'w')

#s_handle = open(sys.argv[1], 'r')
s_parser = SProt.RecordParser()
try:
	s_iterator = SProt.Iterator(s_handle, s_parser)
except:
	print "Error encountered when making s_iterator"
	usage()
	sys.exit(4)

args=string.join(sys.argv[1:])
if  args.find('--updateaccs') == -1:
	updateaccs = open("sprot_parser.updateaccs", 'w')

if args.find('--emblaccs') == -1:
	emblaccs = open("sprot_parser.emblaccs", 'w')

while 1:
	embl_accs = ()
	cur_record = s_iterator.next()
	
	if cur_record is None:
		break
		
	updateaccs.write(cur_record.fasta + '\t')
	emblaccs.write(cur_record.fasta + '\t')
	
#	print dir(cur_record)
	for cref in cur_record.cross_references :
		if cref[0] == 'EMBL':
			embl_accs = embl_accs + cref[1:len(cref)-1]
	emblaccs.write('|'.join(embl_accs) + '\n')
	 

	sys.stdout.write(cur_record.fasta+'\n')
	updateaccs.write( '|'.join(cur_record.accessions) + '\n')

