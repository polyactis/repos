#!/usr/bin/python

from Bio import GenBank
import sys

gb_handle = open(sys.argv[1], 'r')
parser = GenBank.RecordParser()
gb_iterator = GenBank.Iterator(gb_handle, parser)
i=1


"""

while 1:
	try:
		cur_record = gb_iterator.next()
		if cur_record is None:
			break
		# now do something with the record
		#print dir(cur_record.seq)
		#print cur_record.seq.data
		sys.stdout.write('No ' + repr(i) + ':\t')
		i=i+1
		sys.stdout.write('|'.join(cur_record.accession) + '\n')
		print dir(cur_record)
		
		for feature in cur_record.features:
			if feature.key == 'CDS':
				print dir(feature)
				print dir(feature.location)
				#print dir(feature.qualifiers)
				for qualifier in feature.qualifiers:
					pass
					if qualifier.has_key('protein_id'):
						print qualifier('protein_id')
					if qualifier.key=='protein_id":
						print qualifier.value
					if qualifier.key=='translation':
						print qualifer.value
					if qualifier.key=='gene':	
						print qualifer.value
	except:
		print "Error when processing the data: No %d" % i
		i=i+1


"""
