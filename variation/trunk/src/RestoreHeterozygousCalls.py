#!/usr/bin/env python
"""
Usage: RestoreHeterozygousCalls.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input file
	-o ...,	output file
	-t ...,	data_table, 'justin_data'(default)
	-s ...,	strain_info table, 'strain_info'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-w,	draw_only, output_fname contains data already
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	RestoreHeterozygousCalls.py -i ../data/justin_data_filtered.csv -o ../data/justin_data_filtered_y.csv

Description:
	Given an input data set, restore those heterozygous calls masked as NA by dbSNP2data.py.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, math
import Numeric, pylab
from sets import Set
from codense.common import db_connect
from common import nt2number

class RestoreHeterozygousCalls:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', \
		input_fname=None, output_fname=None, strain_info_table='strain_info', snp_locus_table='snp_locus',\
		draw_only=0, debug=0, report=0):
		"""
		2007-03-20
		2007-04-03
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.data_table = data_table
		self.strain_info_table = strain_info_table
		self.snp_locus_table = snp_locus_table
		self.draw_only = int(draw_only)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_id2index(self, curs, db_table, acc_ls):
		sys.stderr.write("Getting id2index ... ")
		id2index = {}
		for i in range(len(acc_ls)):
			curs.execute("select id from %s where acc='%s'"%(db_table, acc_ls[i]))
			rows = curs.fetchall()
			id = rows[0][0]
			id2index[id] = i
		sys.stderr.write("Done.\n")
		return id2index
	
	def get_heterozygous_and_coarse_data_matrix(self, data_matrix, heterozygous_number_cutoff=5):
		"""
		2007-03-20
			for heterozygous_data_matrix, heterozygous is 5 to 10, else is 0
			for coarse_data_matrix, heterozygous is 2,
				homozygous is 1
				NA is 0
		"""
		sys.stderr.write("Getting heterozygous data_matrix and coarse_data_matrix...")
		hetero_sign_matrix = data_matrix>=heterozygous_number_cutoff
		homo_sign_matrix = (data_matrix<heterozygous_number_cutoff) * (data_matrix>0)
		heterozygous_data_matrix = hetero_sign_matrix*data_matrix
		coarse_data_matrix = hetero_sign_matrix*2 + homo_sign_matrix
		NA_sign_matrix = data_matrix==0
		no_of_NAs = sum(sum(NA_sign_matrix))
		no_of_fake_NAs = sum(sum(hetero_sign_matrix))
		sys.stderr.write("Done.\n")
		print "no_of_fake_NAs/no_of_NAs: %s"%(float(no_of_fake_NAs)/no_of_NAs)
		return heterozygous_data_matrix, coarse_data_matrix
	
	def displayDataMatrix(self, data_matrix, title):
		"""
		2007-03-20
		"""
		sys.stderr.write("Displaying Data matrix ...")
		pylab.clf()
		pylab.title(title)
		data_matrix_ls = list(data_matrix)	#reverse the data_matrix to match the y index
		data_matrix_ls.reverse()
		data_matrix = Numeric.array(data_matrix_ls)
		pylab.imshow(data_matrix, interpolation='nearest')
		pylab.colorbar()
		pylab.show()
		sys.stderr.write("done.\n")
	
	def run(self):
		"""
		2007-03-20
		2007-04-03
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		if self.draw_only:
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.output_fname)
			data_matrix = Numeric.array(data_matrix)
		else:
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
			
			snp_acc_ls = header[2:]
			strain_id2index = self.get_id2index(curs, self.strain_info_table, strain_acc_list)
			snp_id2index = self.get_id2index(curs, self.snp_locus_table, snp_acc_ls)
			
			from dbSNP2data import dbSNP2data
			dbSNP2data_instance = dbSNP2data(report=self.report)
			data_matrix = dbSNP2data_instance.get_data_matrix(curs, strain_id2index, snp_id2index, nt2number, self.data_table, need_heterozygous_call=1)
			
			FilterStrainSNPMatrix_instance.write_data_matrix(data_matrix, self.output_fname, header, strain_acc_list, category_list)
		
		heterozygous_data_matrix, coarse_data_matrix = self.get_heterozygous_and_coarse_data_matrix(data_matrix)
		self.displayDataMatrix(heterozygous_data_matrix, title='heterozygous_data_matrix, 5-10=hetero, else=0')
		self.displayDataMatrix(coarse_data_matrix, title='coarse_data_matrix, 0=NA, 1=homo, 2=hetero')
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:t:s:n:wbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'dbsnp'
	input_fname = None
	output_fname = None
	data_table = 'justin_data'
	strain_info_table = 'strain_info'
	snp_locus_table = 'snp_locus'
	draw_only = 0
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
			output_fname = arg
		elif opt in ("-t",):
			data_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-w"):
			draw_only = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname and hostname and dbname and schema:
		instance = RestoreHeterozygousCalls(hostname, dbname, schema, input_fname, output_fname,\
			strain_info_table, snp_locus_table, draw_only, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)