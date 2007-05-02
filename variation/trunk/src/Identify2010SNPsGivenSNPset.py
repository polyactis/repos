#!/usr/bin/env python
"""
Usage: Identify2010SNPsGivenSNPset.py [OPTIONS] -i INPUT_FNAME -a DATA_DIR -o OUTPUT_FNAME

Option:
	-z ..., --hostname=...	the hostname, dl324b-1(default)
	-d ..., --dbname=...	the database name, yhdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input_fname
	-a ...,	2010 aligned sequence data directory
	-o ...,	output_fname
	-s ...,	strain_info table, 'strain_info'(default)
	-t ..., strain_info_2010_table, 'strain_info_2010'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-j ...,	sub_justin_output_fname, (output the sub_data_matrix extracted from input_fname
	-e ...,	what type of difference should be outputted, 5(default, call_ineq_call)
	-f ...,	output_fname which stores the strain and snp where difference happens
	-c,	comparison_only, output_fname already contains the data. no extraction
		of SNP_matrix_2010.
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Identify2010SNPsGivenSNPset.py -i ../data/justin_data.csv -a ../data/2010/alignments_042006/ -o /tmp/snp_matrix_2010
	
	./src/Identify2010SNPsGivenSNPset.py -i data/justin_data_y.csv -a data/2010/alignments_042006/ -o data/snp_matrix_2010.justin.149snps -c -j data/justin_data_y.2010.96strains.149snp.csv
	
Description:
	Based on the SNPs given by input_fname, check the alignment data from 2010 project
	and figure out the SNP matrix for the 96 strains in 2010 project.
	'del_vs_NA':0,
	'del_vs_call':1,
	'NA_vs_NA':2,
	'NA_vs_call':3,
	'call_vs_NA':4,
	'call_ineq_call':5,
	'call_eq_call':6
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
import psycopg, sys, getopt, csv, re, cStringIO
from codense.common import db_connect
from common import nt2number, number2nt
import Numeric
from sets import Set
from transfacdb import fasta_block_iterator

class Identify2010SNPsGivenSNPset:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', \
		input_fname=None, data_dir_2010=None, output_fname=None, strain_info_table='strain_info', \
		strain_info_2010_table='strain_info_2010', snp_locus_table='snp_locus', sub_justin_output_fname=None,\
		diff_type_to_be_outputted=5, diff_output_fname=None, comparison_only=0, debug=0, report=0):
		"""
		2007-03-29
		2007-04-03
			add sub_justin_output_fname and comparison_only
		2007-05-01
			add diff_type_to_be_outputted, diff_output_fname
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.data_dir_2010 = data_dir_2010
		self.output_fname = output_fname
		self.strain_info_table = strain_info_table
		self.strain_info_2010_table = strain_info_2010_table
		self.snp_locus_table = snp_locus_table
		self.sub_justin_output_fname = sub_justin_output_fname
		self.diff_type_to_be_outputted = int(diff_type_to_be_outputted)
		self.diff_output_fname = diff_output_fname
		self.comparison_only = int(comparison_only)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_SNPpos2index(self, curs, snp_acc_ls, snp_locus_table):
		sys.stderr.write("Getting SNPpos2index ...")
		SNPpos2index = {}
		for i in range(len(snp_acc_ls)):
			snp_acc = snp_acc_ls[i]
			curs.execute("select chromosome, position from %s where acc='%s'"%(snp_locus_table, snp_acc))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			SNPpos2index[(chromosome, position)] = i
		sys.stderr.write("Done.\n")
		return SNPpos2index
	
	def get_positions_to_be_checked_ls(self, block):
		block = cStringIO.StringIO(block)
		header_line = block.readline()
		sequence = ''
		for line in block:
			sequence += line[:-1]
		positions_to_be_checked_ls = []
		for i in range(len(sequence)):
			if sequence[i]!='-':
				positions_to_be_checked_ls.append(i)
		return positions_to_be_checked_ls
	
	def get_align_matrix_from_fname(self, fname):
		"""
		2007-03-29
		"""
		inf = open(fname)
		iter = fasta_block_iterator(inf)
		first_block = iter.next()
		positions_to_be_checked_ls = self.get_positions_to_be_checked_ls(first_block)
		
		align_matrix = []
		abbr_name_ls = []
		for block in iter:
			block = cStringIO.StringIO(block)
			header_line = block.readline()
			abbr_name_ls.append(header_line[1:-1])
			sequence = ''
			for line in block:
				sequence += line[:-1]
			align_one_seq = []
			for i in positions_to_be_checked_ls:
				align_one_seq.append(nt2number[sequence[i]])
			align_matrix.append(align_one_seq)
		del iter
		del inf
		return abbr_name_ls, align_matrix
	
	def get_align_length_from_fname(self, fname):
		"""
		2007-03-29
		"""
		inf = open(fname)
		iter = fasta_block_iterator(inf)
		first_block = iter.next()
		positions_to_be_checked_ls = self.get_positions_to_be_checked_ls(first_block)
		del iter
		del inf
		return len(positions_to_be_checked_ls)
	
	def create_SNP_matrix_2010(self, SNPpos2index, data_dir_2010):
		"""
		2007-03-29
		"""
		sys.stderr.write("Creating SNP_matrix_2010 ...\n")
		dirs = os.listdir(data_dir_2010)
		sys.stderr.write("\tTotally, %d dirs to be processed.\n"%len(dirs))
		first_dir = os.path.join(data_dir_2010, dirs[0])
		first_fname = os.listdir(first_dir)[0]
		pathname = os.path.join(first_dir, first_fname)
		abbr_name_ls, align_matrix = self.get_align_matrix_from_fname(pathname)	#the first file is used to detect how many strains
		no_of_strains = len(abbr_name_ls)
		no_of_snps = len(SNPpos2index)
		SNP_matrix_2010 = Numeric.zeros([no_of_strains, no_of_snps], Numeric.Int)
		for sub_dir in dirs:
			dir_pathname = os.path.join(data_dir_2010, sub_dir)
			files = os.listdir(dir_pathname)
			sys.stderr.write("  Dir:%s has %s files.\n"%(sub_dir, len(files)))
			for f in files:
				sys.stderr.write("\t Dir:%s File: %d/%d:\t%s\n"%(sub_dir, files.index(f)+1,len(files),f))
				file_pathname = os.path.join(dir_pathname, f)
				align_length = self.get_align_length_from_fname(file_pathname)
				chromosome = sub_dir	#the directory name is the chromosome
				start_pos = int(f)	#the file name is the starting position on the chromosome
				align_matrix = None
				for i in range(start_pos, start_pos+align_length):
					SNP_pos_key = (chromosome, i)
					if SNP_pos_key in SNPpos2index:
						if align_matrix==None:
							abbr_name_ls, align_matrix = self.get_align_matrix_from_fname(file_pathname)
						SNP_index = SNPpos2index[SNP_pos_key]
						if self.debug:
							import pdb
							pdb.set_trace()
						for j in range(no_of_strains):
							SNP_matrix_2010[j, SNP_index] = align_matrix[j][i-start_pos]
		sys.stderr.write("Done.\n")
		return abbr_name_ls, SNP_matrix_2010
	
	def get_mapping_info_regarding_strain_acc(self, curs, strain_info_table, strain_info_2010_table, abbr_name_ls):
		"""
		2007-03-29
			strain_acc_ls is ordered by strain_info_table.id
		"""
		sys.stderr.write("Getting mapping info about strain_acc, abbr_name and matrix index ...")
		abbr_name2index = dict(zip(abbr_name_ls, range(len(abbr_name_ls))))
		strain_acc2abbr_name = {}
		strain_acc_ls = []
		strain_acc2index = {}
		curs.execute("select s1.id, s1.acc, s2.acc from %s s1, %s s2 where s1.parental_stock_acc = s2.stock_center_number \
			order by id"%(strain_info_table, strain_info_2010_table))
		rows = curs.fetchall()
		for row in rows:
			id, acc, abbr_name = row
			strain_acc_ls.append(acc)
			strain_acc2abbr_name[acc] = abbr_name
			strain_acc2index[acc] = abbr_name2index[abbr_name]
		sys.stderr.write("Done.\n")
		return strain_acc_ls, strain_acc2abbr_name, strain_acc2index
	
	def shuffle_data_matrix_according_to_strain_acc_ls(self, SNP_matrix, strain_acc_ls, strain_acc2index):
		"""
		2007-03-29
		"""
		sys.stderr.write("Shuffling SNP_matrix ...")
		SNP_matrix_new = Numeric.zeros(SNP_matrix.shape, Numeric.Int)
		for i in range(len(strain_acc_ls)):
			SNP_matrix_new[i,:] = SNP_matrix[strain_acc2index[strain_acc_ls[i]],:]
		sys.stderr.write("Done.\n")
		return SNP_matrix_new
	
	def extract_sub_data_matrix(self, src_strain_acc_list, data_matrix, tg_strain_acc_list):
		"""
		2007-04-01
		"""
		sys.stderr.write("Extracting sub data_matrix ...")
		sub_data_matrix = Numeric.zeros([len(tg_strain_acc_list), data_matrix.shape[1]], Numeric.Int)
		src_strain_acc2index = dict(zip(src_strain_acc_list, range(len(src_strain_acc_list))))
		for i in range(len(tg_strain_acc_list)):
			sub_data_matrix[i,:] = data_matrix[src_strain_acc2index[tg_strain_acc_list[i]],:]
		sys.stderr.write("Done.\n")
		return sub_data_matrix
	
	def compare_two_SNP_matrix(self, SNP_matrix_1, SNP_matrix_2):
		"""
		2007-04-01
			SNP_matrix_1 is 2010 matrix, containing -1(deletion)
			SNP_matrix_2 is justin_data, >=0 (no deletion)
		"""
		sys.stderr.write("Comparing two SNP-matrix ...\n")
		no_of_rows, no_of_cols = SNP_matrix_1.shape
		diff_matrix = Numeric.zeros([no_of_rows, no_of_cols], Numeric.Int)
		diff_tag2counter = {}
		diff_tag_dict = {
			'del_vs_NA':0,
			'del_vs_call':1,
			'NA_vs_NA':2,
			'NA_vs_call':3,
			'call_vs_NA':4,
			'call_ineq_call':5,
			'call_eq_call':6}
		for tag in diff_tag_dict:
			diff_tag2counter[tag] = 0
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				if SNP_matrix_1[i,j]>4:	#2007-05-01
					SNP_matrix_1[i,j]=0
				if SNP_matrix_1[i,j]==-1:
					if SNP_matrix_2[i,j]==0:
						tag='del_vs_NA'
					else:
						tag='del_vs_call'
				elif SNP_matrix_1[i,j]==0:
					if SNP_matrix_2[i,j]==0:
						tag='NA_vs_NA'
					else:
						tag='NA_vs_call'
				else:
					if SNP_matrix_2[i,j]==0:
						tag='call_vs_NA'
					elif SNP_matrix_1[i,j] != SNP_matrix_2[i,j]:
						tag='call_ineq_call'
					elif SNP_matrix_1[i,j] == SNP_matrix_2[i,j]:
						tag='call_eq_call'
				diff_matrix[i,j] = diff_tag_dict[tag]
				diff_tag2counter[tag] += 1
		sys.stderr.write("Done.\n")
		return diff_matrix, diff_tag_dict, diff_tag2counter
	
	def outputDiffType(self, diff_matrix, SNP_matrix_1, SNP_matrix_2, diff_tag_dict, diff_type_to_be_outputted, strain_name_ls, snp_name_ls, output_fname):
		"""
		2007-05-01
		"""
		sys.stderr.write("Outputting diff_matrix where diff code=%s ... "%(diff_type_to_be_outputted))
		for tag, diff_code in diff_tag_dict.iteritems():
			if diff_code==diff_type_to_be_outputted:
				tag_to_be_outputted = tag
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['diff tag = %s'%tag_to_be_outputted])
		writer.writerow(['strain', 'snp', '2010_call', 'justin_call'])
		no_of_rows, no_of_cols = diff_matrix.shape
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				if diff_matrix[i,j] == diff_type_to_be_outputted:
					writer.writerow([strain_name_ls[i], snp_name_ls[j], number2nt[SNP_matrix_1[i,j]], number2nt[SNP_matrix_2[i,j]]])
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-03-29
		2007-04-03
		2007-05-01
			--db_connect()
			--FilterStrainSNPMatrix_instance.read_data()
			if self.comparison_only:
				--FilterStrainSNPMatrix_instance.read_data()
			else:
				--get_SNPpos2index()
				--create_SNP_matrix_2010()
					--get_align_length_from_fname()
						--get_positions_to_be_checked_ls()
					--get_align_matrix_from_fname()
						--get_positions_to_be_checked_ls()
				--get_mapping_info_regarding_strain_acc()
				--shuffle_data_matrix_according_to_strain_acc_ls()
				--FilterStrainSNPMatrix_instance.write_data_matrix()
			
			--extract_sub_data_matrix()
			if self.sub_justin_output_fname:
				--FilterStrainSNPMatrix_instance.write_data_matrix()
			--compare_two_SNP_matrix()
			--outputDiffType()
			
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, src_strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
		if self.comparison_only:
			header, strain_acc_ls, abbr_name_ls_sorted, SNP_matrix_2010_sorted = FilterStrainSNPMatrix_instance.read_data(self.output_fname)
			SNP_matrix_2010_sorted = Numeric.array(SNP_matrix_2010_sorted)
		else:
			#extract data from alignment
			snp_acc_ls = header[2:]
			SNPpos2index = self.get_SNPpos2index(curs, snp_acc_ls, self.snp_locus_table)
			abbr_name_ls, SNP_matrix_2010 = self.create_SNP_matrix_2010(SNPpos2index, self.data_dir_2010)
			strain_acc_ls, strain_acc2abbr_name, strain_acc2index = self.get_mapping_info_regarding_strain_acc(curs, self.strain_info_table, self.strain_info_2010_table, abbr_name_ls)
			SNP_matrix_2010_sorted = self.shuffle_data_matrix_according_to_strain_acc_ls(SNP_matrix_2010, strain_acc_ls, strain_acc2index)
			abbr_name_ls_sorted = []
			for strain_acc in strain_acc_ls:
				abbr_name_ls_sorted.append(strain_acc2abbr_name[strain_acc])
			FilterStrainSNPMatrix_instance.write_data_matrix(SNP_matrix_2010_sorted, self.output_fname, header, strain_acc_ls, abbr_name_ls_sorted)
		
		
		#comparison
		data_matrix = Numeric.array(data_matrix)
		sub_data_matrix = self.extract_sub_data_matrix(src_strain_acc_list, data_matrix, strain_acc_ls)
		if self.sub_justin_output_fname:
			FilterStrainSNPMatrix_instance.write_data_matrix(sub_data_matrix, self.sub_justin_output_fname, header, strain_acc_ls, abbr_name_ls_sorted)
		diff_matrix, diff_tag_dict, diff_tag2counter= self.compare_two_SNP_matrix(SNP_matrix_2010_sorted, sub_data_matrix)
		if self.diff_output_fname:
			self.outputDiffType(diff_matrix, SNP_matrix_2010_sorted, sub_data_matrix, diff_tag_dict, self.diff_type_to_be_outputted, abbr_name_ls_sorted, header[2:], self.diff_output_fname)
		
		summary_result_ls = []
		for tag, counter in diff_tag2counter.iteritems():
			summary_result_ls.append('%s(%s):%s'%(tag, diff_tag_dict[tag], counter))
			print '\t%s(%s)\t%s'%(tag, diff_tag_dict[tag], counter)
		import pylab
		pylab.clf()
		diff_matrix_reverse = list(diff_matrix)
		diff_matrix_reverse.reverse()
		diff_matrix_reverse = Numeric.array(diff_matrix_reverse)
		pylab.imshow(diff_matrix_reverse, interpolation='nearest')
		pylab.title(' '.join(summary_result_ls))
		pylab.colorbar()
		pylab.show()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:a:o:s:t:n:j:e:f:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'dbsnp'
	input_fname = None
	data_dir_2010 = None
	output_fname = None
	strain_info_table = 'strain_info'
	strain_info_2010_table = 'strain_info_2010'
	snp_locus_table = 'snp_locus'
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
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-a",):
			data_dir_2010 = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-t",):
			strain_info_2010_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-j",):
			sub_justin_output_fname = arg
		elif opt in ("-e",):
			diff_type_to_be_outputted = int(arg)
		elif opt in ("-f",):
			diff_output_fname = arg
		elif opt in ("-c",):
			comparison_only = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and data_dir_2010 and output_fname and hostname and dbname and schema:
		instance = Identify2010SNPsGivenSNPset(hostname, dbname, schema, input_fname, data_dir_2010, output_fname,\
			strain_info_table, strain_info_2010_table, snp_locus_table, sub_justin_output_fname, \
			diff_type_to_be_outputted, diff_output_fname, comparison_only, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)