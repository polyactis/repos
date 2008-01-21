#!/usr/bin/python

import SProt
import sys

def db_sp_make(inf):
	
	s_parser = SProt.RecordParser()
	s_iterator = SProt.Iterator(inf,s_parser)

	while 1:
		record = s_iterator.next()

		if record is None:
			break

		for acc in record.accessions:
			sys.stdout.write('%s\t%s\t%s\n' % ( record.fasta,acc,record.sequence))

		#print dir(record)

if __name__ == '__main__':
	inf = open(sys.argv[1],'r')

	db_sp_make(inf)
