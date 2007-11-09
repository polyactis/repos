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
	-f ...,	latex_output_fname which stores both the summary and detailed difference tables
		specify this latex_output_fname if you want output
	-c,	comparison_only, output_fname and sub_justin_output_fname are already there
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/CmpAccession2Ecotype.py -p at.ecotype2accession -o 2010pcr_with_sequenom_149snps_ecotype2accession.csv -j sequenom_with_strains_matched_to_2010pcr_ecotype2accession.csv -r
	
	./src/CmpAccession2Ecotype.py -p at.accession2ecotype_complete -o 2010pcr_with_sequenom_149snps_accession2ecotype_complete.csv -j  sequenom_with_strains_matched_to_2010pcr_accession2ecotype_complete.csv -r -f 2010pcr_vs_sequenom_diff.tex

Description:
	program to compare the 149snp calls based on the common strains inn 2010 pcr data and Justin's sequenom data.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/transfac/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import psycopg, sys, getopt, csv, re, cStringIO, numpy
from codense.common import db_connect
from common import nt2number, number2nt
import Numeric
from sets import Set
from transfacdb import fasta_block_iterator

class CmpAccession2Ecotype:
	"""
	2007-10-31
	2007-11-06
		ecotype_id is (ecotype_id, duplicate)
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', ecotype_table='ecotype',\
		accession2ecotype_table=None, alignment_table='at.alignment', sequence_table='at.sequence',\
		snp_locus_table='stock20071008.snps', calls_table='stock20071008.calls', output_fname=None,\
		sub_justin_output_fname=None,\
		diff_type_to_be_outputted=5, latex_output_fname='', comparison_only=0, debug=0, report=0):
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
		self.latex_output_fname = latex_output_fname
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
		"""
		2007-11-06
			add snp_index2snp_info_ls
		"""
		sys.stderr.write("Setting up SNP data structure ...")
		SNPpos2col_index = {}
		snpid2col_index = {}
		snp_acc_ls = []
		snp_index2snp_info_ls = {}
		curs.execute("select id, snpid, chromosome, position from %s order by chromosome, position"%(snp_locus_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, snp_acc, chromosome, position = row
			snpid2col_index[snpid] = len(snpid2col_index)
			SNPpos2col_index[(chromosome, position)] = snpid2col_index[snpid]
			snp_acc_ls.append(snp_acc)
			snp_index2snp_info_ls[snpid2col_index[snpid]] = [snp_acc, chromosome, position]
		sys.stderr.write("Done.\n")
		return SNPpos2col_index, snpid2col_index, snp_acc_ls, snp_index2snp_info_ls
	
	def setup_accession_ecotype_dstruc(self, curs, accession2ecotype_table, ecotype_table, calls_table):
		"""
		2007-10-31
		2007-11-06
			ecotype_id = (ecotype_id, duplicate)
			add accession_id2ecotype_id_ls
		"""
		sys.stderr.write("Setting up accession or ecotype data structure ...")
		ecotype_id2accession_id = {}	#not accession2ecotype because multiple ecotype ids correspond to same accession_id
		ecotype_id2row_index = {}
		ecotype_id_ls = []
		ecotype_id2info_ls = {}
		accession_id2row_index = {}
		accession_id2ecotype_id_ls = {}
		accession_id_ls = []
		curs.execute("select distinct ea.ecotype_id, c.duplicate, e.nativename, e.stockparent, ea.accession_id from %s ea, %s e, %s c where e.id=ea.ecotype_id and e.id=c.ecotypeid order by nativename"%(accession2ecotype_table, ecotype_table, calls_table))
		rows = curs.fetchall()
		for row in rows:
			ecotype_id, duplicate, nativename, stockparent, accession_id = row
			ecotype_id = (ecotype_id, duplicate)
			ecotype_id2accession_id[ecotype_id] = accession_id
			ecotype_id2row_index[ecotype_id] = len(ecotype_id2row_index)
			ecotype_id2info_ls[ecotype_id] = [nativename, stockparent]
			ecotype_id_ls.append(ecotype_id)
			if accession_id not in accession_id2row_index:
				accession_id2row_index[accession_id] = len(accession_id2row_index)
				accession_id_ls.append(accession_id)
			if accession_id not in accession_id2ecotype_id_ls:
				accession_id2ecotype_id_ls[accession_id] = []
			accession_id2ecotype_id_ls[accession_id].append(ecotype_id)
		sys.stderr.write("Done.\n")
		return ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2info_ls, ecotype_id_ls, accession_id2row_index, accession_id_ls, accession_id2ecotype_id_ls
	
	def get_alignment_id2positions_to_be_checked_ls(self, curs, alignment_table):
		"""
		2007-11-05
		2007-11-06
			add alignment_id2start
		"""
		sys.stderr.write("Getting alignment_id2positions_to_be_checked_ls ...")
		alignment_id2positions_to_be_checked_ls = {}
		alignment_id2start = {}
		curs.execute("select id, start, target from %s"%alignment_table)
		rows = curs.fetchall()
		for row in rows:
			alignment_id, start, target_seq = row
			positions_to_be_checked_ls = []
			for i in range(len(target_seq)):
				if target_seq[i]!='-':
					positions_to_be_checked_ls.append(i)
			alignment_id2positions_to_be_checked_ls[alignment_id] = positions_to_be_checked_ls
			alignment_id2start[alignment_id] = start
		sys.stderr.write("Done.\n")
		return alignment_id2positions_to_be_checked_ls, alignment_id2start
	
	def get_ecotype_X_snp_matrix(self, curs, ecotype_id2row_index, snpid2col_index, calls_table):
		"""
		2007-11-06
			add data_matrix_touched
			ecotype_id = (ecotype_id, duplicate)
		"""
		sys.stderr.write("Getting ecotype_X_snp_matrix ... \n")
		data_matrix = numpy.zeros([len(ecotype_id2row_index), len(snpid2col_index)], numpy.integer)
		data_matrix_touched = numpy.zeros([len(ecotype_id2row_index), len(snpid2col_index)], numpy.integer)
		curs.execute("select ecotypeid, snpid, call1, callhet, duplicate from %s"%(calls_table))
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				ecotype_id, snpid, call, callhet, duplicate = row
				ecotype_id = (ecotype_id, duplicate)
				call = call.upper()
				if callhet:
					callhet.upper()	#2007-09-20	just in case
					call = call+callhet
				if ecotype_id in ecotype_id2row_index and snpid in snpid2col_index:
					call_number = nt2number[call]
					data_matrix[ecotype_id2row_index[ecotype_id], snpid2col_index[snpid]] = call_number
					data_matrix_touched[ecotype_id2row_index[ecotype_id], snpid2col_index[snpid]] = 1
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return data_matrix, data_matrix_touched
	
	def get_accession_X_snp_matrix(self, curs, accession_id2row_index, SNPpos2col_index, sequence_table, alignment_table, alignment_id2positions_to_be_checked_ls):
		"""
		2007-11-05
			add alignment_id2positions_to_be_checked_ls to skip '-' columns in reference genome to match the real genome position
		2007-11-06
			try to differentiate 'missing'(which is 'N') and 'not sequenced' (which is not in the db) using the data_matrix_touched
			0 means not sequenced, 1 means sequenced
		2007-11-07
			add snp_index2alignment_id
		"""
		sys.stderr.write("Getting accession_X_snp_matrix ...\n")
		data_matrix = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		data_matrix_touched = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		snp_index2alignment_id = {}
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
							data_matrix_touched[accession_id2row_index[accession_id], SNP_index] = 1
							if SNP_index not in snp_index2alignment_id:
								snp_index2alignment_id[SNP_index] = alignment_id
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return data_matrix, data_matrix_touched, snp_index2alignment_id
	
	def cmp_two_matricies(self, accession_X_snp_matrix, accession_X_snp_matrix_touched, ecotype_X_snp_matrix, ecotype_X_snp_matrix_touched, nt_number2diff_matrix_index, ecotype_id2accession_id, ecotype_id2row_index, accession_id2row_index, snp_column=-1, need_diff_details_ls=0):
		"""
		2007-11-06
			add accession_X_snp_matrix_touched
			split diff_matrix into diff_matrix_for_touched and diff_matrix_for_untouched
		2007-11-07
			add snp_column and need_diff_details_ls to deal with 2-snp-column comparison
		"""
		if self.report or self.debug:
			sys.stderr.write("Comparing two matricies ...")
		diff_matrix_touched_accession_vs_touched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_touched_accession_vs_untouched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_untouched_accession_vs_touched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_untouched_accession_vs_untouched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		if snp_column!=-1:
			snp_columns_to_be_checked = [snp_column]
		else:
			no_of_snps = ecotype_X_snp_matrix.shape[1]
			snp_columns_to_be_checked = range(no_of_snps)
		diff_details_ls = []
		for ecotype_id, accession_id in ecotype_id2accession_id.iteritems():
			e_row_index = ecotype_id2row_index[ecotype_id]
			a_row_index = accession_id2row_index[accession_id]
			for i in snp_columns_to_be_checked:
				a_nt_diff_index = nt_number2diff_matrix_index[accession_X_snp_matrix[a_row_index, i]]
				e_nt_diff_index = nt_number2diff_matrix_index[ecotype_X_snp_matrix[e_row_index, i]]
				if accession_X_snp_matrix_touched[a_row_index, i]==1 and ecotype_X_snp_matrix_touched[e_row_index, i]==1:
					diff_matrix_touched_accession_vs_touched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
					if need_diff_details_ls:
						if a_nt_diff_index!=e_nt_diff_index:
							diff_details_ls.append([accession_id, accession_X_snp_matrix[a_row_index, i], ecotype_id, ecotype_X_snp_matrix[e_row_index, i]])
				elif accession_X_snp_matrix_touched[a_row_index, i]==1 and ecotype_X_snp_matrix_touched[e_row_index, i]==0:
					diff_matrix_touched_accession_vs_untouched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
				elif accession_X_snp_matrix_touched[a_row_index, i]==0 and ecotype_X_snp_matrix_touched[e_row_index, i]==1:
					diff_matrix_untouched_accession_vs_touched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
				else:
					diff_matrix_untouched_accession_vs_untouched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return [diff_matrix_touched_accession_vs_touched_ecotype, diff_matrix_touched_accession_vs_untouched_ecotype, diff_matrix_untouched_accession_vs_touched_ecotype, diff_matrix_untouched_accession_vs_untouched_ecotype], diff_details_ls
	
	def wrap_diff_matrix_with_row_col_names(self, diff_matrix):
		if self.report or self.debug:
			sys.stderr.write("Wrapping diff_matrix with row, column names ...")
		row_name_ls = ['-', 'NA', 'A', 'C', 'G', 'T', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT']
		wrapped_diff_matrix = [ [''] + row_name_ls]
		for i in range(diff_matrix.shape[0]):
			wrapped_diff_matrix.append([row_name_ls[i]] + diff_matrix[i,:].tolist())
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return wrapped_diff_matrix
	
	def cmp_two_lists(self, accession_X_snp_ls, accession_X_snp_ls_touched, ecotype_X_snp_ls, ecotype_X_snp_ls_touched, nt_number2diff_matrix_index):
		"""
		2007-11-06
			similar to cmp_two_matricies(), but compare two lists
		"""
		if self.report or self.debug:
			sys.stderr.write("Comparing two lists ...")
		diff_matrix_touched_accession_vs_touched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_touched_accession_vs_untouched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_untouched_accession_vs_touched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_matrix_untouched_accession_vs_untouched_ecotype = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		diff_details_ls = []
		for i in range(len(accession_X_snp_ls)):
			a_nt_diff_index = nt_number2diff_matrix_index[accession_X_snp_ls[i]]
			e_nt_diff_index = nt_number2diff_matrix_index[ecotype_X_snp_ls[i]]
			if accession_X_snp_ls_touched[i]==1 and ecotype_X_snp_ls_touched[i]==1:
				diff_matrix_touched_accession_vs_touched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
				if a_nt_diff_index!=e_nt_diff_index:
					diff_details_ls.append([i, accession_X_snp_ls[i], ecotype_X_snp_ls[i]])
			elif accession_X_snp_ls_touched[i]==1 and ecotype_X_snp_ls_touched[i]==0:
				diff_matrix_touched_accession_vs_untouched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
			elif accession_X_snp_ls_touched[i]==0 and ecotype_X_snp_ls_touched[i]==1:
				diff_matrix_untouched_accession_vs_touched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
			else:
				diff_matrix_untouched_accession_vs_untouched_ecotype[a_nt_diff_index, e_nt_diff_index] += 1
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return [diff_matrix_touched_accession_vs_touched_ecotype, diff_matrix_touched_accession_vs_untouched_ecotype, diff_matrix_untouched_accession_vs_touched_ecotype, diff_matrix_untouched_accession_vs_untouched_ecotype], diff_details_ls
	
	def beautify_diff_details_ls(self, diff_details_ls, snp_index2snp_info_ls, alignment_id2start, snp_index2alignment_id):
		"""
		2007-11-07
			add snp_index2alignment_id
		"""
		if self.report or self.debug:
			sys.stderr.write("Beautifying diff_details_ls ...")
		new_diff_details_ls = []
		for row in diff_details_ls:
			snp_index, pcr_call, sequenom_call = row
			pcr_call = number2nt[pcr_call]
			sequenom_call = number2nt[sequenom_call]
			snp_acc, snp_chr, snp_pos = snp_index2snp_info_ls[snp_index]
			alignment_id = snp_index2alignment_id[snp_index]
			alignment_start = alignment_id2start[alignment_id]
			new_diff_details_ls.append([snp_acc, snp_chr, snp_pos, alignment_id, alignment_start, pcr_call, sequenom_call])
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return new_diff_details_ls
	
	def beautify_snp_diff_details_ls(self, diff_details_ls, ecotype_id2info_ls):
		"""
		2007-11-07
		"""
		if self.report or self.debug:
			sys.stderr.write("Beautifying SNP diff_details_ls ...")
		new_diff_details_ls = []
		diff_details_ls.sort()	#sort based on accession_id
		for row in diff_details_ls:
			accession_id, pcr_call, ecotype_id, sequenom_call = row
			nativename, stockparent = ecotype_id2info_ls[ecotype_id]
			new_diff_details_ls.append([nativename, stockparent, ecotype_id[0], ecotype_id[1], accession_id, number2nt[pcr_call], number2nt[sequenom_call]])
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return new_diff_details_ls
	
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
		SNPpos2col_index, snpid2col_index, snp_acc_ls, snp_index2snp_info_ls = self.setup_SNP_dstruc(curs, self.snp_locus_table)
		ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2info_ls, ecotype_id_ls, accession_id2row_index, accession_id_ls, accession_id2ecotype_id_ls = self.setup_accession_ecotype_dstruc(curs, self.accession2ecotype_table, self.ecotype_table, self.calls_table)
		
		ecotype_X_snp_matrix, ecotype_X_snp_matrix_touched = self.get_ecotype_X_snp_matrix(curs, ecotype_id2row_index, snpid2col_index, self.calls_table)
		if self.sub_justin_output_fname:
			header = ['ecotype_id', 'ecotype_id'] + snp_acc_ls
			FilterStrainSNPMatrix_instance.write_data_matrix(ecotype_X_snp_matrix, self.sub_justin_output_fname, header, ecotype_id_ls, ecotype_id_ls)
		
		alignment_id2positions_to_be_checked_ls, alignment_id2start = self.get_alignment_id2positions_to_be_checked_ls(curs, self.alignment_table)
		accession_X_snp_matrix, accession_X_snp_matrix_touched, snp_index2alignment_id= self.get_accession_X_snp_matrix(curs, accession_id2row_index, SNPpos2col_index, self.sequence_table, self.alignment_table, alignment_id2positions_to_be_checked_ls)
		
		if self.output_fname:
			header = ['accession_id', 'accession_id'] + snp_acc_ls
			FilterStrainSNPMatrix_instance.write_data_matrix(accession_X_snp_matrix, self.output_fname, header, accession_id_ls, accession_id_ls)
		summary_diff_matrix_ls, diff_details_ls = self.cmp_two_matricies(accession_X_snp_matrix, accession_X_snp_matrix_touched, ecotype_X_snp_matrix, ecotype_X_snp_matrix_touched, nt_number2diff_matrix_index, ecotype_id2accession_id, ecotype_id2row_index, accession_id2row_index)
		print "diff_matrix_touched_accession_vs_touched_ecotype"
		print summary_diff_matrix_ls[0]
		print "diff_matrix_touched_accession_vs_untouched_ecotype"
		print summary_diff_matrix_ls[1]
		print "diff_matrix_untouched_accession_vs_touched_ecotype"
		print summary_diff_matrix_ls[2]
		print "diff_matrix_untouched_accession_vs_untouched_ecotype"
		print summary_diff_matrix_ls[3]
		
		summary_diff_matrix_caption_ls = ['PCR-tried vs sequenom-tried', 'PCR-tried vs sequenom-untried', 'PCR-untried vs sequenom-tried', 'PCR-untried vs sequenom-untried']
		
		if self.latex_output_fname:
			outf = open(self.latex_output_fname, 'w')
			outf.write('\\section{2010 PCR versus sequenom. summary} \\label{section_summary}\n')
			for i in range(len(summary_diff_matrix_ls)):
				from pymodule.latex import outputMatrixInLatexTable
				wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(summary_diff_matrix_ls[i])
				table_label = 'table_dm%s'%i
				outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, summary_diff_matrix_caption_ls[i], table_label))
			#Strain-wise comparison
			outf.write('\\section{2010 PCR versus sequenom for each strain} \\label{section_strain_wise}\n')
			accession_id_ls.sort()
			table_no = i
			for accession_id in accession_id_ls:
				ecotype_id_ls = accession_id2ecotype_id_ls[accession_id]
				outf.write('\\subsection{strain %s(accession id=%s)}\n'%(ecotype_id2info_ls[ecotype_id_ls[0]][0], accession_id))
				for ecotype_id in ecotype_id_ls:
					outf.write('\\subsubsection{corresponding ecotype %s(stkparent=%s, ecotype id=%s, duplicate=%s)}\n'%(ecotype_id2info_ls[ecotype_id][0], ecotype_id2info_ls[ecotype_id][1], ecotype_id[0], ecotype_id[1]))
					e_row_index = ecotype_id2row_index[ecotype_id]
					a_row_index = accession_id2row_index[accession_id]
					
					diff_matrix_ls, diff_details_ls= self.cmp_two_lists(accession_X_snp_matrix[a_row_index,:], accession_X_snp_matrix_touched[a_row_index,:], ecotype_X_snp_matrix[e_row_index,:], ecotype_X_snp_matrix_touched[e_row_index,:], nt_number2diff_matrix_index)
					wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(diff_matrix_ls[0])
					table_no += 1
					table_label = 'table_dm%s'%table_no
					caption = 'accession id=%s vs ecotype id=%s, duplicate=%s(nativename=%s, stockparent=%s)'%(accession_id, ecotype_id[0], ecotype_id[1], ecotype_id2info_ls[ecotype_id][0], ecotype_id2info_ls[ecotype_id][1])
					outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, caption, table_label))
					
					if diff_details_ls:
						diff_details_ls = self.beautify_diff_details_ls(diff_details_ls, snp_index2snp_info_ls, alignment_id2start, snp_index2alignment_id)
						table_no += 1
						table_label = 'table_dm%s'%table_no
						caption = 'detailed difference for accession id=%s vs ecotype id=%s, duplicate=%s'%(accession_id, ecotype_id[0], ecotype_id[1])
						outf.write(outputMatrixInLatexTable(diff_details_ls, caption, table_label, header_ls=['snp', 'chromosome', 'position', 'alignment_id', 'alignment_start', 'pcr_call', 'sequenom_call']))
			#SNP-wise comparison
			outf.write('\\section{2010 PCR versus sequenom for each SNP} \\label{section_snp_wise}\n')
			for snp_column in range(accession_X_snp_matrix.shape[1]):
				snp_acc, chromosome, position = snp_index2snp_info_ls[snp_column]
				alignment_id = snp_index2alignment_id[snp_column]
				alignment_start = alignment_id2start[alignment_id]
				outf.write('\\subsection{SNP %s(chrom=%s, pos=%s, alignment id=%s, alignment start=%s)}\n'%(snp_acc, chromosome, position, alignment_id, alignment_start))
				
				diff_matrix_ls, diff_details_ls = self.cmp_two_matricies(accession_X_snp_matrix, accession_X_snp_matrix_touched, ecotype_X_snp_matrix, ecotype_X_snp_matrix_touched, nt_number2diff_matrix_index, ecotype_id2accession_id, ecotype_id2row_index, accession_id2row_index, snp_column=snp_column, need_diff_details_ls=1)
				wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(diff_matrix_ls[0])
				table_no += 1
				table_label = 'table_dm%s'%table_no
				caption = 'SNP %s(chromosome=%s, position=%s, alignment id=%s, alignment start=%s)'%(snp_acc, chromosome, position, alignment_id, alignment_start)
				outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, caption, table_label))
				
				if diff_details_ls:
					diff_details_ls = self.beautify_snp_diff_details_ls(diff_details_ls, ecotype_id2info_ls)
					table_no += 1
					table_label = 'table_dm%s'%table_no
					caption = 'detailed difference for SNP %s'%(snp_acc)
					header_ls = ['nativename', 'stkparent', 'ecotype_id', 'duplicate', 'accession_id', 'pcr_call', 'sequenom_call']
					outf.write(outputMatrixInLatexTable(diff_details_ls, caption, table_label, header_ls))
			del outf

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
	latex_output_fname = None
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
			latex_output_fname = arg
		elif opt in ("-c",):
			comparison_only = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if ecotype_table and accession2ecotype_table and sequence_table and alignment_table and snp_locus_table and calls_table and hostname and dbname and schema:
		instance = CmpAccession2Ecotype(hostname, dbname, schema, ecotype_table,\
			accession2ecotype_table, alignment_table, sequence_table,\
			snp_locus_table, calls_table, output_fname, sub_justin_output_fname,\
			diff_type_to_be_outputted, latex_output_fname, comparison_only, debug, report)
		
		instance.run()
	else:
		print __doc__
		sys.exit(2)
