#!/usr/bin/env python
"""
Usage: data2dbSNP.py [OPTIONS] -i INPUT_FILE -o OUTPUT_TABLE

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	INPUT_FILE
	-o ...,	SNP output table
	-s ...,	strain_info table, 'strain_info'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-g ...	the two-letter abbreviation of organism, at(default)
		only for type =1
	-y ...	parser type(1, default)
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-c, --commit	commit this database transaction
	-h, --help	show this help

Examples:
	data2dbSNP.py -i UCHIC_HO_011107.csv -o justin_data -r -c

Description:
	put SNP data into database schema
	
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
from codense.common import db_connect, org_short2long, org2tax_id

class data2dbSNP:
	"""
	2007-02-17
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_fname=None, \
		output_table=None, strain_info_table='strain_info', snp_locus_table='snp_locus', \
		organism='hs', type=1, debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_table = output_table
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		self.tax_id = org2tax_id(org_short2long(organism))
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
		
		self.snp_acc_category_pattern = re.compile("([a-zA-Z]*[\-]*[a-zA-Z]+)[\-_ ]*[\w]+")
		#- (followed by number, 2nd - followed by alphanumeric), _ or space is separator
	
	def get_strain_acc2id(self, curs, strain_info_table, tax_id):
		sys.stderr.write("Getting mapping data from %s ..."%strain_info_table)
		strain_acc2id = {}
		curs.execute("DECLARE crs CURSOR FOR select acc, id from %s where tax_id=%s"%(strain_info_table, tax_id))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				acc, id = row
				strain_acc2id[acc] = id
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return strain_acc2id
	
	def get_snp_acc2id(self, curs, snp_locus_table, tax_id):
		snp_acc2id = self.get_strain_acc2id(curs, snp_locus_table, tax_id)
		return snp_acc2id
	
	def submit_to_strain_info_table(self, curs, strain_acc, category, new_id, tax_id, strain_info_table):
		curs.execute("insert into %s(acc, category, id, tax_id) values('%s', '%s', %s, %s)"%\
			(strain_info_table, strain_acc, category, new_id, tax_id))
	
	def expand_snp_locus_table(self, curs, snp_locus_table, snp_acc_list, snp_acc2id, tax_id):
		sys.stderr.write("Expanding snp_locus_table ...")
		for snp_acc in snp_acc_list:
			if snp_acc not in snp_acc2id:
				new_id = len(snp_acc2id)+1
				curs.execute("insert into %s(acc, id, tax_id) values('%s', %s, %s)"%(snp_locus_table, snp_acc, new_id, tax_id))
				snp_acc2id[snp_acc] = new_id
		sys.stderr.write("Done.\n")
		return snp_acc2id
	
	def submit_to_snp_table(self, curs, strain_id, snp_index2acc, snp_acc2id, call_list, output_table):
		snp_counter = 0
		for i in range(len(call_list)):
			call = call_list[i]
			if call!='NA' and call!='0':	#0 or NA means Not Available
				snp_id = snp_acc2id[snp_index2acc[i]]
				curs.execute("insert into %s(strain_id, snp_id, call) values(%s, %s, '%s')"%(output_table, strain_id, snp_id, call))
				snp_counter += 1
		return snp_counter
	
	def parse_file(self, input_fname, curs, output_table, tax_id, snp_acc_category_pattern, strain_acc2id, snp_acc2id, snp_locus_table):
		sys.stderr.write("Starting to parse %s ...\n"%input_fname)
		reader = csv.reader(open(input_fname), delimiter='\t')
		snp_acc_list = reader.next()[1:]
		snp_index2acc = dict(zip(range(len(snp_acc_list)), snp_acc_list))
		snp_acc2id = self.expand_snp_locus_table(curs, snp_locus_table, snp_acc_list, snp_acc2id, tax_id)
		counter = 0
		snp_counter = 0
		for row in reader:
			strain_acc = row[0]
			call_list = row[1:]
			if strain_acc not in strain_acc2id:
				new_id = len(strain_acc2id)+1
				try:
					category = snp_acc_category_pattern.match(strain_acc).groups()[0]
				except:
					print row
					print counter
					print strain_acc
				self.submit_to_strain_info_table(curs, strain_acc, category, new_id, tax_id, strain_info_table)
				strain_acc2id[strain_acc] = new_id
			snp_counter += self.submit_to_snp_table(curs, new_id, snp_index2acc, snp_acc2id, call_list, output_table)
			counter += 1
			if self.report and counter%100==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*20, counter, snp_counter))
		if self.report:
			sys.stderr.write('%s%s\t%s'%('\x08'*20, counter, snp_counter))
		del reader
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2007-02-17
		
		-db_connect()
		-get_strain_acc2id()
		-get_snp_acc2id()
		-parse_file()
			-expand_snp_locus_table()
			-submit_to_strain_info_table()
			-submit_to_snp_table()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		strain_acc2id = self.get_strain_acc2id(curs, self.strain_info_table, self.tax_id)
		snp_acc2id = self.get_snp_acc2id(curs, self.snp_locus_table, self.tax_id)
		self.parse_file(self.input_fname, curs, self.output_table, self.tax_id, self.snp_acc_category_pattern, strain_acc2id, snp_acc2id, self.snp_locus_table)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:g:y:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'dl324b-1'
	dbname = 'yhdb'
	schema = 'dbsnp'
	input_fname = None
	output_table = None
	strain_info_table = 'strain_info'
	snp_locus_table = 'snp_locus'
	organism = 'at'
	type = 1
	debug = 0
	report = 0
	commit = 0
	
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
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if input_fname and output_table and hostname and dbname and schema:
		instance = data2dbSNP(hostname, dbname, schema, input_fname, output_table, \
			strain_info_table, snp_locus_table, organism, type, debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
