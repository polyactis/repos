#!/usr/bin/env python
"""
Usage: MarkSNPSelected.py [OPTIONS] -i INPUT_FILE -o OUTPUT_TABLE

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input file
	-o ...,	output table, 'snp_locus'(default)
	-c,	commit
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	MarkSNPSelected.py -i data/justin_data_y_b.csv -r

Description:
	

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import psycopg, sys, getopt, csv, re
from codense.common import db_connect
from sets import Set

class MarkSNPSelected:
	"""
	2007-07-09
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_fname=None, \
		output_table='snp_locus', commit=0, debug=0, report=0):
		"""
		2007-07-09
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_table = output_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def readSNPMarkers(self, input_fname):
		sys.stderr.write("Getting SNP markers...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		snp_acc_ls = header[2:]
		del reader
		sys.stderr.write("Done.\n")
		return snp_acc_ls
	
	def markSelected(self, curs, output_table, snp_acc_ls):
		sys.stderr.write("Getting SNP markers...")
		snp_acc_set = Set(snp_acc_ls)
		curs.execute("DECLARE crs CURSOR FOR select id,acc from %s"%(output_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				snp_id, acc = row
				if acc in snp_acc_set:
					curs.execute("update %s set selected=1 where id=%s"%(output_table, snp_id))
				else:
					curs.execute("update %s set selected=0 where id=%s"%(output_table, snp_id))
				counter += 1
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done.\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		snp_acc_ls = self.readSNPMarkers(self.input_fname)
		self.markSelected(curs, self.output_table, snp_acc_ls)
		if self.commit:
			curs.execute("end")
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'dl324b-1'
	dbname = 'yhdb'
	schema = 'dbsnp'
	input_fname = None
	output_table = 'snp_locus'
	commit = 0
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if input_fname and output_table and hostname and dbname and schema:
		instance = dbSNP2data(hostname, dbname, schema, input_fname, output_table, \
			commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)