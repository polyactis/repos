#!/usr/bin/env python
"""
Usage: dbSNP2data.py [OPTIONS] -i INPUT_TABLE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input table
	-o ...,	output file
	-s ...,	strain_info table, 'strain_info'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-y,	need heterozygous call
	-w,	with header line
	-a,	use alphabet to represent nucleotide, not number
	-b, --debug	just test running the program, no daytime restriction
	-r, --report	report the progress (time before each query)
	-h, --help	show this help

Examples:
	dbSNP2data.py -i justin_data -o justin_data.csv -r

Description:
	output SNP data from database schema
	
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
from common import nt2number, number2nt
import Numeric as num

class dbSNP2data:
	"""
	2007-02-19
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_table=None, \
		output_fname=None, strain_info_table='strain_info', snp_locus_table='snp_locus', \
		organism='hs', need_heterozygous_call=0, with_header_line=0, nt_alphabet=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_fname = output_fname
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		#self.tax_id = org2tax_id(org_short2long(organism))
		self.need_heterozygous_call = int(need_heterozygous_call)
		self.with_header_line = int(with_header_line)
		self.nt_alphabet = int(nt_alphabet)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_snp_id2index(self, curs, input_table):
		sys.stderr.write("Getting snp_id2index ...")
		snp_id2index = {}
		snp_id_list = []
		curs.execute("select distinct snp_id from %s order by snp_id"%(input_table))
		rows = curs.fetchall()
		for row in rows:
			snp_id = row[0]
			snp_id_list.append(snp_id)
			snp_id2index[snp_id] = len(snp_id2index)
		sys.stderr.write("Done.\n")
		return snp_id2index, snp_id_list
	
	def get_strain_id2index(self, curs, input_table):
		sys.stderr.write("Getting strain_id2index ...")
		strain_id2index = {}
		strain_id_list = []
		curs.execute("select distinct strain_id from %s order by strain_id"%(input_table))
		rows = curs.fetchall()
		for row in rows:
			strain_id = row[0]
			strain_id_list.append(strain_id)
			strain_id2index[strain_id] = len(strain_id2index)
		sys.stderr.write("Done.\n")
		return strain_id2index, strain_id_list
	
	def get_strain_id_info(self, curs, strain_id_list, strain_info_table):
		sys.stderr.write("Getting strain_id_info ...")
		strain_id2acc = {}
		strain_id2category = {}
		for strain_id in strain_id_list:
			curs.execute("select acc, category from %s where id=%s"%(strain_info_table, strain_id))
			rows = curs.fetchall()
			acc, category = rows[0]
			strain_id2acc[strain_id] = acc
			strain_id2category[strain_id] = category
		sys.stderr.write("Done.\n")
		return strain_id2acc, strain_id2category
	
	def get_snp_id_info(self, curs, snp_id_list, snp_locus_table):
		sys.stderr.write("Getting snp_id_info ...")
		snp_id2acc = {}
		for snp_id in snp_id_list:
			curs.execute("select acc from %s where id=%s"%(snp_locus_table, snp_id))
			rows = curs.fetchall()
			acc = rows[0][0]
			snp_id2acc[snp_id] = acc
		sys.stderr.write("Done.\n")
		return snp_id2acc
	
	def get_data_matrix(self, curs, strain_id2index, snp_id2index, nt2number, input_table, need_heterozygous_call):
		sys.stderr.write("Getting data_matrix ...\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)])
		curs.execute("DECLARE crs CURSOR FOR select strain_id, snp_id, call from %s"%(input_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				strain_id, snp_id, call = row
				call_number = nt2number[call]
				if need_heterozygous_call:
					data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
				elif call_number<=4:	#single letter or NA
					data_matrix[strain_id2index[strain_id], snp_id2index[snp_id]] = call_number
				counter += 1
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def write_data_matrix(self, data_matrix, output_fname, strain_id_list, snp_id_list, snp_id2acc, with_header_line, nt_alphabet, strain_id2acc=None, strain_id2category=None):
		"""
		2007-02-19
			if strain_id2acc is available, translate strain_id into strain_acc,
			if strain_id2category is available, add 'category'
		2007-02-25
			if one strain's SNP row is all NA, it'll be skipped
		"""
		sys.stderr.write("Writing data_matrix ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		if with_header_line:
			header_row = ['strain']
			if strain_id2category:
				header_row.append('category')
			for snp_id in snp_id_list:
				header_row.append(snp_id2acc[snp_id])
			writer.writerow(header_row)
		for i in range(len(data_matrix)):
			if strain_id2acc:
				new_row = [strain_id2acc[strain_id_list[i]]]
			else:
				new_row = [strain_id_list[i]]
			if strain_id2category:
				new_row.append(strain_id2category[strain_id_list[i]])
			if data_matrix[i]!=0:	#rows with all NAs(after heterozygous calles are removed) skipped
				for j in data_matrix[i]:
					if nt_alphabet:
						j = number2nt[j]
					new_row.append(j)
				writer.writerow(new_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def sort_file(self, ofname):
		"""
		2007-02-19
		1. sort file1 into file2
		2. move file2 to file1
		"""
		sys.stderr.write("\tSorting %s"%ofname)
		commandline = 'sort %s > %s.post'%(ofname, ofname)
		exit_code = system_call(commandline)
		commandline = 'mv %s.post %s'%(ofname, ofname)
		exit_code = system_call(commandline)
		sys.stderr.write(".\n")
		
	def run(self):
		"""
		2007-02-19
			--db_connect
			--get_snp_id2index()
			--get_strain_id2index()
			--get_strain_id_info()
			--get_snp_id_info()
			--get_data_matrix()
			--write_data_matrix()
			#--sort_file()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		snp_id2index, snp_id_list = self.get_snp_id2index(curs, self.input_table)
		strain_id2index, strain_id_list = self.get_strain_id2index(curs, self.input_table)
		
		strain_id2acc, strain_id2category = self.get_strain_id_info(curs, strain_id_list, self.strain_info_table)
		snp_id2acc = self.get_snp_id_info(curs, snp_id_list, self.snp_locus_table)
		data_matrix = self.get_data_matrix(curs, strain_id2index, snp_id2index, nt2number, self.input_table, self.need_heterozygous_call)
		self.write_data_matrix(data_matrix, self.output_fname, strain_id_list, snp_id_list, snp_id2acc, self.with_header_line,\
			self.nt_alphabet, strain_id2acc, strain_id2category)
		#self.sort_file(self.output_fname)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:g:ywabrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'dl324b-1'
	dbname = 'yhdb'
	schema = 'dbsnp'
	input_table = None
	output_fname = None
	strain_info_table = 'strain_info'
	snp_locus_table = 'snp_locus'
	organism = 'at'
	need_heterozygous_call = 0
	with_header_line = 0
	nt_alphabet = 0
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
			input_table = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-y",):
			need_heterozygous_call = 1
		elif opt in ("-w",):
			with_header_line = 1
		elif opt in ("-a",):
			nt_alphabet = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_table and output_fname and hostname and dbname and schema:
		instance = dbSNP2data(hostname, dbname, schema, input_table, output_fname, \
			strain_info_table, snp_locus_table, organism, need_heterozygous_call, with_header_line, nt_alphabet, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)