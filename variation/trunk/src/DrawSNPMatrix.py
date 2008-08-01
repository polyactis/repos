#!/usr/bin/env python
"""
Usage: DrawSNPMatrix.py [OPTIONS] -i XXX -o OUTPUT_FNAME_PREFIX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, dbsnp(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input_fname
	-o ...,	output_fname_prefix
	-s ...,	snps_sequenom_info_table, 'dbsnp.snps_sequenom_info'(default)
	-x ...,	row label type, 1 (default)
	-y ...,	drawing_type, 1 (default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	DrawSNPMatrix.py -i ./2010_with_149snps_ecotype2accession.csv -o /tmp/2010_with_149snps_ecotype2accession_y4 -y4

Description:
	Program to draw image for a strain X snp matrix. input_fname's format is like the one outputted by dbSNP2data.py
	two images generated.
		output_fname_prefix.png: the matrix itself
		output_fname_prefix_legend.png: the legend of the color used in drawing
	Drawing Type:
	 1. original, all 12 colors (including '-'/deletion , 'NA', A,C,G,T, all hets)
	 2. into 5 colors, '-', NA, allele1, allele2, heterozygous (allele1/allele2 determined by snps_sequenom_info_table)
	 3. into 4 colors, '-', NA, homozygous, heterozygous
	 4. into 5 colors, '-', NA, allele1, allele2, heterozygous (allele1/allele2 randomly assigned)
	 5. into 6 colors, '-', NA, A,C,G,T (every strain duplicates into two rows (haplotypes))
	Row label type:
	 1. 1st column
	 2. 2nd column.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import numpy
from pymodule.DrawMatrix import drawMatrix, drawLegend
from variation.src.common import nt2number, number2color, number2nt

class DrawSNPMatrix:
	"""
	2007-11-02
	"""
	def __init__(self, hostname='localhost', dbname='dbsnp', schema='dbsnp', input_fname='',\
		output_fname_prefix='', snps_sequenom_info_table='dbsnp.snps_sequenom_info', row_label_type=1, drawing_type=1, debug=0, report=0):
		"""
		2007-11-02
		2007-11-05
			add row_label_type
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_fname_prefix = output_fname_prefix
		self.snps_sequenom_info_table = snps_sequenom_info_table
		self.row_label_type = int(row_label_type)
		self.drawing_type = int(drawing_type)
		self.debug = int(debug)
		self.report = int(report)
	
	def transformMatrixIntoTwoAllelesAndHetero(self, data_matrix, snp_allele2index_ls=None):
		sys.stderr.write("Transforming matrix into Two Alleles and Hetero ...")
		no_of_rows = len(data_matrix)
		no_of_cols = len(data_matrix[0])
		new_data_matrix = numpy.zeros((no_of_rows, no_of_cols), numpy.integer)
		for j in range(no_of_cols):	#column by column
			if not snp_allele2index_ls:	#if not given how to map allele1 and allele2
				snp_allele2index = {}
			for i in range(no_of_rows):
				if data_matrix[i][j] <= 0:	#NA or deletion
					new_data_matrix[i][j] = data_matrix[i][j]
				elif data_matrix[i][j]>4:	#hetero
					new_data_matrix[i][j] = 3
				else:	#homo
					if not snp_allele2index_ls:	#if not given how to map allele1 and allele2
						if data_matrix[i][j] not in snp_allele2index:
							snp_allele2index[data_matrix[i][j]] = len(snp_allele2index)+1
					else:
						snp_allele2index = snp_allele2index_ls[j]
					new_data_matrix[i][j] = snp_allele2index[data_matrix[i][j]]
		sys.stderr.write("Done.\n")
		return new_data_matrix
	
	def transformMatrixIntoHomoAndHetero(self, data_matrix, snp_allele2index_ls=None):
		sys.stderr.write("Transforming matrix into Homo and Hetero ...")
		no_of_rows = len(data_matrix)
		no_of_cols = len(data_matrix[0])
		new_data_matrix = numpy.zeros((no_of_rows, no_of_cols), numpy.integer)
		for j in range(no_of_cols):	#column by column
			for i in range(no_of_rows):
				if data_matrix[i][j] <= 0:	#NA or deletion
					new_data_matrix[i][j] = data_matrix[i][j]
				elif data_matrix[i][j]>4:
					new_data_matrix[i][j] = 2	#hetero
				else:	#homo
					new_data_matrix[i][j] = 1
		sys.stderr.write("Done.\n")
		return new_data_matrix
	
	def transformMatrixIntoFourNucleotides(self, data_matrix):
		sys.stderr.write("Transforming matrix into Four Nucleotides ...")
		no_of_rows = len(data_matrix)
		no_of_cols = len(data_matrix[0])
		new_data_matrix = numpy.zeros((no_of_rows*2, no_of_cols), numpy.integer)
		for j in range(no_of_cols):	#column by column
			for i in range(no_of_rows):
				if data_matrix[i][j] <= 4:	#NA or deletion or homozygous
					new_data_matrix[2*i][j] = data_matrix[i][j]
					new_data_matrix[2*i+1][j] = data_matrix[i][j]
				else:	#heterozygous
					nt = number2nt[data_matrix[i][j]]
					new_data_matrix[2*i][j] = nt2number[nt[0]]
					new_data_matrix[2*i+1][j] = nt2number[nt[1]]
		sys.stderr.write("Done.\n")
		return new_data_matrix
	
	def get_snp_allele2index_ls(self, curs, snp_acc_ls, snps_sequenom_info_table='dbsnp.snps_sequenom_info'):
		sys.stderr.write("Getting snp_allele2index_ls ...")
		snp_allele2index_ls = []
		for snp_acc in snp_acc_ls:
			curs.execute("select ext1_call, ext2_call from %s where snpid='%s'"%(snps_sequenom_info_table, snp_acc))
			rows = curs.fetchall()
			allele1, allele2 = rows[0]
			snp_allele2index = {}
			snp_allele2index[nt2number[allele1]] = 1
			snp_allele2index[nt2number[allele2]] = 2
			snp_allele2index_ls.append(snp_allele2index)
		sys.stderr.write("Done.\n")
		return snp_allele2index_ls
	
	def run(self):
		from pymodule import read_data
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
		snp_acc_ls = header[2:]
		if self.debug:
			import pdb
			pdb.set_trace()
		if self.drawing_type==1:
			matrix_value2label= number2nt
			matrix_value2color = number2color
			data_matrix = numpy.array(data_matrix)
		elif self.drawing_type==2:
			import MySQLdb
			#conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
			conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
			curs = conn.cursor()
			snp_allele2index_ls = self.get_snp_allele2index_ls(curs, snp_acc_ls, self.snps_sequenom_info_table)
			data_matrix = self.transformMatrixIntoTwoAllelesAndHetero(data_matrix, snp_allele2index_ls)
			matrix_value2label= {-1:'deletion', 0:'NA', 1:'allele1', 2:'allele2', 3:'hetero'}
			matrix_value2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(0,255,0), 3:(255,0,0)}
		elif self.drawing_type==3:
			data_matrix = self.transformMatrixIntoHomoAndHetero(data_matrix)
			matrix_value2label= {-1:'deletion', 0:'NA', 1:'homo', 2:'hetero'}
			matrix_value2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(255,0,0)}
		elif self.drawing_type==4:
			data_matrix = self.transformMatrixIntoTwoAllelesAndHetero(data_matrix)
			matrix_value2label= {-1:'deletion', 0:'NA', 1:'allele1', 2:'allele2', 3:'hetero'}
			matrix_value2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(0,255,0), 3:(255,0,0)}
		elif self.drawing_type==5:
			data_matrix = self.transformMatrixIntoFourNucleotides(data_matrix)
			matrix_value2label= {-1: '-', 0: 'NA',	1:'A',	2:'C',	3:'G',	4:'T'}
			matrix_value2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(0,255,0), 3:(255,0,0), 4:(122,0,122)}
			new_strain_acc_list = []
			for strain_acc in strain_acc_list:
				new_strain_acc_list.append(strain_acc)
				new_strain_acc_list.append(strain_acc)
			strain_acc_list = new_strain_acc_list
		else:
			sys.stderr.write("drawing_type %s not supported\n"%drawing_type)
			sys.exit(2)
		row_label_type2label_ls = {1:strain_acc_list,
			2:category_list}
		im = drawLegend(matrix_value2label, matrix_value2color)
		im.save('%s_legend.png'%self.output_fname_prefix)
		im = drawMatrix(data_matrix, matrix_value2color, row_label_type2label_ls[self.row_label_type], snp_acc_ls, with_grid=1)
		im.save('%s.png'%self.output_fname_prefix)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	import getopt
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:x:y:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'dbsnp'
	schema = 'dbsnp'
	input_fname = None
	output_fname_prefix = None
	snps_sequenom_info_table = 'snps_sequenom_info'
	row_label_type = 1
	drawing_type = 1
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
			output_fname_prefix = arg
		elif opt in ("-s",):
			snps_sequenom_info_table = arg
		elif opt in ("-y",):
			drawing_type = int(arg)
		elif opt in ("-x",):
			row_label_type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname_prefix and hostname and dbname and schema:
		instance = DrawSNPMatrix(hostname, dbname, schema, input_fname, output_fname_prefix, \
			snps_sequenom_info_table, row_label_type, drawing_type, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)