#!/usr/bin/env python
"""
Usage: CmpAccession2Ecotype.py [OPTIONS] -p XXX -o OUTPUT_FNAME -j XXX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete'/'at.ecotype2accession'
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-l ...,	calls table, 'stock20071008.calls'(default)
	-o ...,	output_fname
	-j ...,	sub_justin_output_fname, (output the sub_data_matrix extracted from input_fname)
	
	- ...,	what type of difference should be outputted, 5(default, call_ineq_call)
	-f ...,	output_fname which stores the strain and snp where difference happens
	-c,	comparison_only, output_fname already contains the data. no extraction
		of SNP_matrix_2010.
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/CmpAccession2Ecotype.py -p at.ecotype2accession -o 2010_with_149snps_ecotype2accession.csv -j 149snp_with_2010_ecotype2accession.csv -r
	
	./src/CmpAccession2Ecotype.py -p at.accession2ecotype_complete -o 2010_with_149snps_accession2ecotype_complete.csv -j 149snp_with_2010_accession2ecotype_complete.csv -r

Description:
	program to compare the 149snp calls based on the common strains inn 2010 pcr data and Justin's sequenom data.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/transfac/src')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/transfac/src')))
import psycopg, sys, getopt, csv, re, cStringIO, numpy
from codense.common import db_connect
from common import nt2number, number2nt
import Numeric
from sets import Set
from transfacdb import fasta_block_iterator

class CmpAccession2Ecotype:
	"""
	2007-10-31
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', ecotype_table='ecotype',\
		accession2ecotype_table=None, alignment_table='at.alignment', sequence_table='at.sequence',\
		snp_locus_table='stock20071008.snps', calls_table='stock20071008.calls', output_fname=None,\
		sub_justin_output_fname=None,\
		diff_type_to_be_outputted=5, diff_output_fname=None, comparison_only=0, debug=0, report=0):
		"""
		2007-10-31
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ecotype_table = ecotype_table
		self.accession2ecotype_table = accession2ecotype_table
		self.alignment_table = alignment_table
		self.sequence_table = sequence_table
		self.snp_locus_table = snp_locus_table
		self.calls_table = calls_table
		self.output_fname = output_fname
		self.sub_justin_output_fname = sub_justin_output_fname
		
		self.diff_type_to_be_outputted = int(diff_type_to_be_outputted)
		self.diff_output_fname = diff_output_fname
		self.comparison_only = int(comparison_only)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_nt_number2diff_matrix_index(self, nt2number):
		"""
		2007-10-31
			nucleotide number ranges from -1 to 10.
			the diff_matrix_index ranges from 0 to 11.
		"""
		sys.stderr.write("Getting nt_number2diff_matrix_index ...")
		nt_number2diff_matrix_index = {}
		for nt, number in nt2number.iteritems():
			nt_number2diff_matrix_index[number] = number+1	#'-'(deletion) is -1
		sys.stderr.write("Done.\n")
		return nt_number2diff_matrix_index
	
	def setup_SNP_dstruc(self, curs, snp_locus_table):
		sys.stderr.write("Setting up SNP data structure ...")
		SNPpos2col_index = {}
		snpid2col_index = {}
		snp_acc_ls = []
		curs.execute("select id, snpid, chromosome, position from %s order by chromosome, position"%(snp_locus_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, snp_acc, chromosome, position = row
			snpid2col_index[snpid] = len(snpid2col_index)
			SNPpos2col_index[(chromosome, position)] = snpid2col_index[snpid]
			snp_acc_ls.append(snp_acc)
		sys.stderr.write("Done.\n")
		return SNPpos2col_index, snpid2col_index, snp_acc_ls
	
	def setup_accession_ecotype_dstruc(self, curs, accession2ecotype_table, ecotype_table):
		"""
		2007-10-31
		"""
		sys.stderr.write("Setting up accession or ecotype data structure ...")
		ecotype_id2accession_id = {}	#not accession2ecotype because multiple ecotype ids correspond to same accession_id
		ecotype_id2row_index = {}
		ecotype_id_ls = []
		ecotype_id2nativename = {}
		accession_id2row_index = {}
		accession_id_ls = []
		curs.execute("select ea.ecotype_id, e.nativename, ea.accession_id from %s ea, %s e where e.id=ea.ecotype_id order by nativename"%(accession2ecotype_table, ecotype_table))
		rows = curs.fetchall()
		for row in rows:
			ecotype_id, nativename, accession_id = row
			ecotype_id2accession_id[ecotype_id] = accession_id
			ecotype_id2row_index[ecotype_id] = len(ecotype_id2row_index)
			ecotype_id2nativename[ecotype_id] = nativename
			ecotype_id_ls.append(ecotype_id)
			if accession_id not in accession_id2row_index:
				accession_id2row_index[accession_id] = len(accession_id2row_index)
				accession_id_ls.append(accession_id)
		sys.stderr.write("Done.\n")
		return ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2nativename, ecotype_id_ls, accession_id2row_index, accession_id_ls
	
	def get_alignment_id2positions_to_be_checked_ls(self, curs, alignment_table):
		"""
		2007-11-05
		"""
		sys.stderr.write("Getting alignment_id2positions_to_be_checked_ls ...")
		alignment_id2positions_to_be_checked_ls = {}
		curs.execute("select id, target from %s"%alignment_table)
		rows = curs.fetchall()
		for row in rows:
			alignment_id, target_seq = row
			positions_to_be_checked_ls = []
			for i in range(len(target_seq)):
				if target_seq[i]!='-':
					positions_to_be_checked_ls.append(i)
			alignment_id2positions_to_be_checked_ls[alignment_id] = positions_to_be_checked_ls
		sys.stderr.write("Done.\n")
		return alignment_id2positions_to_be_checked_ls
	
	def get_ecotype_X_snp_matrix(self, curs, ecotype_id2row_index, snpid2col_index, calls_table):
		sys.stderr.write("Getting ecotype_X_snp_matrix ... \n")
		data_matrix = numpy.zeros([len(ecotype_id2row_index), len(snpid2col_index)], numpy.integer)
		curs.execute("select ecotypeid, snpid, call1, callhet from %s"%(calls_table))
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				ecotype_id, snpid, call, callhet = row
				call = call.upper()
				if callhet:
					callhet.upper()	#2007-09-20	just in case
					call = call+callhet
				if ecotype_id in ecotype_id2row_index and snpid in snpid2col_index:
					call_number = nt2number[call]
					data_matrix[ecotype_id2row_index[ecotype_id], snpid2col_index[snpid]] = call_number
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def get_accession_X_snp_matrix(self, curs, accession_id2row_index, SNPpos2col_index, sequence_table, alignment_table, alignment_id2positions_to_be_checked_ls):
		"""
		2007-11-05
			add alignment_id2positions_to_be_checked_ls to skip '-' columns in reference genome to match the real genome position
		"""
		sys.stderr.write("Getting accession_X_snp_matrix ...\n")
		data_matrix = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		curs.execute("select s.accession, s.alignment, s.bases, a.chromosome, a.start from %s s, %s a where s.alignment=a.id"%(sequence_table, alignment_table))
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				accession_id, alignment_id, bases, chromosome, start_pos = row
				if accession_id in accession_id2row_index:
					positions_to_be_checked_ls = alignment_id2positions_to_be_checked_ls[alignment_id]
					for i in range(len(positions_to_be_checked_ls)):	#i is the real index after reference '-' is removed.
						inflated_index = positions_to_be_checked_ls[i]
						SNP_pos_key = (chromosome, start_pos+i)
						if SNP_pos_key in SNPpos2col_index:
							if self.debug:
								import pdb
								pdb.set_trace()
							SNP_index = SNPpos2col_index[SNP_pos_key]
							data_matrix[accession_id2row_index[accession_id], SNP_index] = nt2number[bases[inflated_index]]
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def cmp_two_matricies(self, accession_X_snp_matrix, ecotype_X_snp_matrix, nt_number2diff_matrix_index, ecotype_id2accession_id, ecotype_id2row_index, accession_id2row_index):
		sys.stderr.write("Comparing two matricies ...")
		diff_matrix = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		no_of_snps = ecotype_X_snp_matrix.shape[1]
		for ecotype_id, accession_id in ecotype_id2accession_id.iteritems():
			e_row_index = ecotype_id2row_index[ecotype_id]
			a_row_index = accession_id2row_index[accession_id]
			for i in range(no_of_snps):
				a_nt_diff_index = nt_number2diff_matrix_index[accession_X_snp_matrix[a_row_index, i]]
				e_nt_diff_index = nt_number2diff_matrix_index[ecotype_X_snp_matrix[e_row_index, i]]
				diff_matrix[a_nt_diff_index, e_nt_diff_index] += 1
		sys.stderr.write("Done.\n")
		return diff_matrix
	
	def run(self):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		
		import MySQLdb
		#conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		nt_number2diff_matrix_index = self.get_nt_number2diff_matrix_index(nt2number)
		SNPpos2col_index, snpid2col_index, snp_acc_ls = self.setup_SNP_dstruc(curs, self.snp_locus_table)
		ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2nativename, ecotype_id_ls, accession_id2row_index, accession_id_ls = self.setup_accession_ecotype_dstruc(curs, self.accession2ecotype_table, self.ecotype_table)
		
		ecotype_X_snp_matrix = self.get_ecotype_X_snp_matrix(curs, ecotype_id2row_index, snpid2col_index, self.calls_table)
		if self.sub_justin_output_fname:
			header = ['ecotype_id', 'ecotype_id'] + snp_acc_ls
			FilterStrainSNPMatrix_instance.write_data_matrix(ecotype_X_snp_matrix, self.sub_justin_output_fname, header, ecotype_id_ls, ecotype_id_ls)
		
		alignment_id2positions_to_be_checked_ls = self.get_alignment_id2positions_to_be_checked_ls(curs, self.alignment_table)
		accession_X_snp_matrix = self.get_accession_X_snp_matrix(curs, accession_id2row_index, SNPpos2col_index, self.sequence_table, self.alignment_table, alignment_id2positions_to_be_checked_ls)
		
		if self.output_fname:
			header = ['accession_id', 'accession_id'] + snp_acc_ls
			FilterStrainSNPMatrix_instance.write_data_matrix(accession_X_snp_matrix, self.output_fname, header, accession_id_ls, accession_id_ls)
		diff_matrix = self.cmp_two_matricies(accession_X_snp_matrix, ecotype_X_snp_matrix, nt_number2diff_matrix_index, ecotype_id2accession_id, ecotype_id2row_index, accession_id2row_index)
		print diff_matrix


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:e:p:s:a:n:l:o:j:f:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype_table = 'stock20071008.ecotype'
	accession2ecotype_table = None
	sequence_table = 'at.sequence'
	alignment_table = 'at.alignment'
	snp_locus_table = 'stock20071008.snps'
	calls_table = 'stock20071008.calls'
	output_fname = None
	sub_justin_output_fname = None
	diff_type_to_be_outputted = 5
	diff_output_fname = None
	comparison_only = 0
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
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-p",):
			accession2ecotype_table = arg
		elif opt in ("-s",):
			sequence_table = arg
		elif opt in ("-a",):
			alignment_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-l",):
			calls_table = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-j",):
			sub_justin_output_fname = arg
		elif opt in ("-f",):
			diff_output_fname = arg
		elif opt in ("-c",):
			comparison_only = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if ecotype_table and accession2ecotype_table and sequence_table and alignment_table and snp_locus_table and calls_table and hostname and dbname and schema:
		instance = CmpAccession2Ecotype(hostname, dbname, schema, ecotype_table,\
			accession2ecotype_table, alignment_table, sequence_table,\
			snp_locus_table, calls_table, output_fname,\
			sub_justin_output_fname,\
			diff_type_to_be_outputted, diff_output_fname, comparison_only, debug, report)
		
		instance.run()
	else:
		print __doc__
		sys.exit(2)
