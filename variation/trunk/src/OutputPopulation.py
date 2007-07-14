#!/usr/bin/env python
"""
Usage: OutputPopulation.py [OPTIONS] -o -p -f

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input table, 'calls'(default)
	-o ...,	output file
	-s ...,	strain_info table, 'ecotype'(default)
	-n ...,	snp_locus_table, 'snps'(default)
	-p ...,	population table
	-f ...,	snpacc_fname, the header of this data matrix file specifies which snps are needed.
	-m ...,	min_no_of_strains_per_pop, 5(default)
	-t ...,	output type, 1(each population is a matrix file, default), 2(RMES)
	-w,	with header line (for output type=1)
	-a,	use alphabet to represent nucleotide, not number
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	OutputPopulation.py  -o /tmp/pop.txt -p populations_50 -f /tmp/chicago.data -t 2
	
	OutputPopulation.py  -o /tmp/pop.txt -p populations_50 -f /tmp/chicago.data.y.filtered.bad  -t 1 -w
	
Description:
	output SNP data in units of population
	
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
from common import nt2number, number2nt
import Numeric as num
from sets import Set

class OutputPopulation:
	"""
	2007-07-12
	"""
	def __init__(self, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', input_table='calls', \
		output_fname=None, strain_info_table='strain_info', snp_locus_table='snp_locus', \
		population_table=None, snpacc_fname='', min_no_of_strains_per_pop=5, output_type=1, \
		with_header_line=0, nt_alphabet=0,\
		debug=0, report=0):
		"""
		2007-07-12
		2007-07-13
			add output_type, need_heterozygous_call, with_header_line, nt_alphabet
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_fname = output_fname
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		self.population_table = population_table
		self.snpacc_fname = snpacc_fname
		self.min_no_of_strains_per_pop = int(min_no_of_strains_per_pop)
		self.output_type = int(output_type)
		self.with_header_line = int(with_header_line)
		self.nt_alphabet = int(nt_alphabet)
		self.debug = int(debug)
		self.report = int(report)
		
		self.OutputPop_dict = {1: self.OutputPopMatrixFormat,
			2: self.OutputPopRMES}
	
	def get_snp_struc(self, curs, snpacc_fname, snp_locus_table):
		"""
		2007-07-12
		"""
		sys.stderr.write("Getting snp_id2index ...")
		
		snp_id2index = {}
		snp_id_list = []
		snp_id2acc = {}
		reader = csv.reader(open(snpacc_fname), delimiter='\t')
		header = reader.next()
		snp_acc_list = header[2:]
		for snp_acc in snp_acc_list:
			curs.execute("select id, snpid from %s where snpid='%s'"%(snp_locus_table, snp_acc))
			rows = curs.fetchall()
			id, snpid = rows[0]
			snp_id_list.append(id)
			snp_id2index[id] = len(snp_id2index)
			snp_id2acc[id] = snpid
		sys.stderr.write("Done.\n")
		return snp_id2index, snp_id_list, snp_acc_list, snp_id2acc
	
	def get_popid2ecotypeid_ls(self, curs, population_table):
		"""
		2007-07-12
		"""
		sys.stderr.write("Getting popid2ecotypeid_ls ...")
		curs.execute("select popid, ecotypeid from %s"%population_table)
		rows = curs.fetchall()
		popid2ecotypeid_ls = {}
		for row in rows:
			popid, ecotypeid = row
			if popid not in popid2ecotypeid_ls:
				popid2ecotypeid_ls[popid] = []
			popid2ecotypeid_ls[popid].append(ecotypeid)
		sys.stderr.write("Done.\n")
		return popid2ecotypeid_ls
	
	def get_strain_id2index(self, popid2ecotypeid_ls, min_no_of_strains_per_pop):
		sys.stderr.write("Getting strain_id2index ...")
		strain_id_list = []
		popid_ls = popid2ecotypeid_ls.keys()
		for popid in popid_ls:
			ecotypeid_ls = popid2ecotypeid_ls[popid]
			if len(ecotypeid_ls)<min_no_of_strains_per_pop:
				del popid2ecotypeid_ls[popid]
			else:
				strain_id_list += ecotypeid_ls
		strain_id_list.sort()
		strain_id2index = dict(zip(strain_id_list, range(len(strain_id_list))))
		sys.stderr.write("Done.\n")
		return strain_id2index, strain_id_list
	
	def OutputPopMatrixFormat(self, data_matrix, popid2ecotypeid_ls, strain_id2index, output_fname_prefix, strain_id2acc, strain_id2category, snp_acc_list, with_header_line, nt_alphabet):
		"""
		2007-07-13
		"""
		sys.stderr.write("Outputting population data in Matrix format (lots of files)...\n")
		header = ['strain', 'category'] + snp_acc_list
		no_of_rows, no_of_cols = data_matrix.shape
		for popid, ecotypeid_ls in popid2ecotypeid_ls.iteritems():
			sys.stderr.write("\tPopulation %s"%popid)
			writer = csv.writer(open('%s.%s'%(output_fname, popid), 'w'), delimiter='\t')
			if with_header_line:
				writer.writerow(header)
			for ecotypeid in ecotypeid_ls:
				strain_acc = strain_id2acc[ecotypeid]
				strain_category = strain_id2category[ecotypeid]
				row = [strain_acc, strain_category]
				i = strain_id2index[ecotypeid]
				for j in range(no_of_cols):
					if nt_alphabet:
						row.append(number2nt[data_matrix[i][j]])
					else:
						row.append(data_matrix[i][j])
				writer.writerow(row)
			del writer
			sys.stderr.write(".\n")
		sys.stderr.write("Done.\n")
		
	
	def OutputPopRMES(self, data_matrix, popid2ecotypeid_ls, strain_id2index, output_fname, strain_id2acc=None, strain_id2category=None, snp_acc_list=None, with_header_line=0, nt_alphabet=0):
		"""
		2007-07-12
			format for RMES
			homozygous = 0
			heterozygous = 1
			NA = other integer -99
		"""
		sys.stderr.write("Outputting population data in RMES format ...")
		no_of_strains, no_of_snps = data_matrix.shape
		writer = csv.writer(open(output_fname, 'w'), delimiter=' ')
		popid_ls = []
		writer.writerow([len(popid2ecotypeid_ls)])
		for popid, ecotypeid_ls in popid2ecotypeid_ls.iteritems():
			writer.writerow(['population%s'%popid])
			writer.writerow([len(ecotypeid_ls)])
			writer.writerow([no_of_snps])
			popid_ls.append(popid)
		for popid in popid_ls:
			ecotypeid_ls = popid2ecotypeid_ls[popid]
			for ecotypeid in ecotypeid_ls:
				data_row = data_matrix[strain_id2index[ecotypeid],]
				new_data_row = []
				for data_point in data_row:
					if data_point>4:
						new_data_row.append(1)
					elif data_point==0:
						new_data_row.append(-99)
					else:
						new_data_row.append(0)
				writer.writerow(new_data_row)
		del writer
		sys.stderr.write("Done.\n")
	
	
	def run(self):
		"""
		2007-07-12
		"""
		from dbSNP2data import dbSNP2data
		dbSNP2data_instance = dbSNP2data()
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		
		snp_id2index, snp_id_list, snp_acc_list, snp_id2acc = self.get_snp_struc(curs, self.snpacc_fname, self.snp_locus_table)
		popid2ecotypeid_ls = self.get_popid2ecotypeid_ls(curs, self.population_table)
		
		strain_id2index, strain_id_list = self.get_strain_id2index(popid2ecotypeid_ls, self.min_no_of_strains_per_pop)
		
		strain_id2acc, strain_id2category = dbSNP2data_instance.get_strain_id_info_m(curs, strain_id_list, self.strain_info_table)
		data_matrix = dbSNP2data_instance.get_data_matrix_m(curs, strain_id2index, snp_id2index, nt2number, self.input_table, need_heterozygous_call=1)
		
		self.OutputPop_dict[self.output_type](data_matrix, popid2ecotypeid_ls, strain_id2index, self.output_fname, strain_id2acc,\
			 strain_id2category, snp_acc_list, self.with_header_line, self.nt_alphabet)
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:p:f:m:t:wabrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	input_table = 'calls'
	output_fname = None
	strain_info_table = 'ecotype'
	snp_locus_table = 'snps'
	population_table = ''
	snpacc_fname = ''
	min_no_of_strains_per_pop = 5
	output_type = 1
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
		elif opt in ("-p",):
			population_table = arg
		elif opt in ("-f",):
			snpacc_fname = arg
		elif opt in ("-m",):
			min_no_of_strains_per_pop = int(arg)
		elif opt in ("-t",):
			output_type = int(arg)
		elif opt in ("-w",):
			with_header_line = 1
		elif opt in ("-a",):
			nt_alphabet = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_table and output_fname and hostname and dbname and schema and population_table and snpacc_fname:
		instance = OutputPopulation(hostname, dbname, schema, input_table, output_fname, \
			strain_info_table, snp_locus_table, population_table, snpacc_fname, \
			min_no_of_strains_per_pop, output_type, with_header_line, nt_alphabet, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)