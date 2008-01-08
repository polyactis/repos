#!/usr/bin/env python
"""
Usage: Cmp149SNPVs2010.py [OPTIONS] -i XXX -j XXX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete' or 'at.ecotype2accession'
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-l ...,	calls table, 'stock20071008.calls'(default)
	-i ...,	input_fname1, 149SNP data filename
	-j ...,	input_fname2, 2010 data filename
	-o ...,	latex_output_fname which stores both the summary and detailed difference tables
		specify this latex_output_fname if you want output
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Cmp149SNPVs2010.py -i ~/script/variation/stock20071008/data_y10001101.tsv -j ~/script/variation/data/2010/data_2010_x_149SNP_y1.tsv -p at.ecotype2accession -o ~/script/variation/genotyping/149snp/analysis/stock20071008_data_y10001101Vs2010_ecotype2accession_diff_matrix.tex
	
	Cmp149SNPVs2010.py -i ~/script/variation/stock20071008/data_y10001101.tsv -j ~/script/variation/data/2010/data_2010_x_149SNP_y1.tsv -p at.accession2ecotype_complete -o ~/script/variation/genotyping/149snp/analysis/stock20071008_data_y10001101Vs2010_accession2ecotype_complete_diff_matrix.tex
	
	Cmp149SNPVs2010.py -d stock20071227 -i ~/script/variation/stock20071227/data_y10001101.tsv -j ~/script/variation/data/2010/data_2010_x_149SNP_y1.tsv -p at.ecotype2accession -o ~/script/variation/genotyping/149snp/analysis/stock20071227_data_y10001101Vs2010_ecotype2accession_diff_matrix.tex
	
	Cmp149SNPVs2010.py -d stock20071227 -i ~/script/variation/stock20071227/data_y10001101.tsv -j ~/script/variation/data/2010/data_2010_x_149SNP_y1.tsv -p at.accession2ecotype_complete -o ~/script/variation/genotyping/149snp/analysis/stock20071227_data_y10001101Vs2010_accession2ecotype_complete_diff_matrix.tex
	
Description:
	program to compare the 149snp calls based on the common strains inn 2010 pcr data and Justin's sequenom data.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt
from variation.src.QualityControl import QualityControl

class Cmp149SNPVs2010(QualityControl):
	"""
	2008-01-01
		QC between 149snp and 2010
	"""
	def __init__(self, curs, input_fname1, input_fname2, latex_output_fname, accession2ecotype_table, debug=0, report=0):
		self.curs = curs
		self.input_fname1 = input_fname1
		self.input_fname2 = input_fname2
		self.latex_output_fname = latex_output_fname
		self.accession2ecotype_table = accession2ecotype_table
		self.debug = int(debug)
		self.report = int(report)
	
	def get_row_matching_dstruc(self, strain_acc_list1, category_list1, strain_acc_list2, accession2ecotype_table):
		"""
		2008-01-07
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			ecotype_id = int(strain_acc_list1[i])
			duplicate = int(category_list1[i])
			strain_acc = (ecotype_id, duplicate)
			strain_acc2row_index1[strain_acc] = i
		
		strain_acc2row_index2 = {}
		for i in range(len(strain_acc_list2)):
			strain_acc = strain_acc_list2[i]
			accession_id = int(strain_acc)
			strain_acc2row_index2[accession_id] = i
		
		row_id12row_id2 = {}
		for strain_acc in strain_acc2row_index1:
			ecotype_id = strain_acc[0]
			if accession2ecotype_table=='at.accession2ecotype_complete':
				curs.execute("select accession_id from %s where ecotype_id=%s and is_stockparent_original=1"%(accession2ecotype_table, ecotype_id))
			else:
				curs.execute("select accession_id from %s where ecotype_id=%s"%(accession2ecotype_table, ecotype_id))
			rows = curs.fetchall()
			if rows:
				accession_id = rows[0][0]
				if accession_id in strain_acc2row_index2:	#could be no match
					row_id12row_id2[strain_acc] = accession_id
				else:
					print 'Failure:', strain_acc
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2

	def load_dstruc(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		QualityControl.load_dstruc(self)
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header1, strain_acc_list1, category_list1, self.data_matrix1 = FilterStrainSNPMatrix_instance.read_data(self.input_fname1)
		header2, strain_acc_list2, category_list2, self.data_matrix2 = FilterStrainSNPMatrix_instance.read_data(self.input_fname2)
	 	
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(header1, header2)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(strain_acc_list1, category_list1, strain_acc_list2, self.accession2ecotype_table)
	
	def get_row_id2info(self, row_id_ls, curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype'):
		"""
		2007-12-1
		"""
		row_id2info = {}
		for row_id in row_id_ls:
			ecotypeid, duplicate = row_id
			curs.execute("SELECT e.nativename, e.stockparent FROM %s e where e.id=%s"%(ecotype_table, ecotypeid))
			rows = curs.fetchall()
			if rows:
				nativename, stockparent = rows[0]
				row_id2info[row_id] = nativename
				row_id2info[row_id] = row_id2info[row_id].decode('utf-8', 'ignore')
			else:
				row_id2info[row_id] = '%s'%row_id
		return row_id2info
	
"""
hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

input_fname1 = os.path.expanduser('~/script/variation/stock20071008/data_y10001101.tsv')
input_fname2 = os.path.expanduser('~/script/variation/data/2010/data_2010_x_149SNP.tsv')

latex_output_fname = ''

ins= Cmp149SNPVs2010(curs, input_fname1, input_fname2, latex_output_fname)
ins.load_dstruc()
ins.plot_row_NA_mismatch_rate('149SNP vs 2010 strain-wise')
#ins.cal_row_id2pairwise_dist()
#new_row_id2pairwise_dist = ins.trim_row_id2pairwise_dist(ins.row_id2pairwise_dist, 10)
ins.plot_col_NA_mismatch_rate('149SNP vs 2010 snp-wise')

"""

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:e:p:s:a:n:l:i:j:o:brh", long_options_list)
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
	input_fname1 = None
	input_fname2 = None
	diff_type_to_be_outputted = 5
	latex_output_fname = None
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
		elif opt in ("-i",):
			input_fname1 = arg
		elif opt in ("-j",):
			input_fname2 = arg
		elif opt in ("-o",):
			latex_output_fname = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if ecotype_table and sequence_table and alignment_table and snp_locus_table and calls_table and input_fname1 and input_fname2 and hostname and dbname and schema:
		import MySQLdb
		conn = MySQLdb.connect(db=dbname,host=hostname)
		curs = conn.cursor()
		ins= Cmp149SNPVs2010(curs, input_fname1, input_fname2, latex_output_fname, accession2ecotype_table, debug, report)
		ins.load_dstruc()
		ins.plot_row_NA_mismatch_rate('149SNP vs 2010 strain-wise')
		#ins.cal_row_id2pairwise_dist()
		#new_row_id2pairwise_dist = ins.trim_row_id2pairwise_dist(ins.row_id2pairwise_dist, 10)
		ins.plot_col_NA_mismatch_rate('149SNP vs 2010 snp-wise')
		ins.output_diff_matrix()
	else:
		print __doc__
		sys.exit(2)
