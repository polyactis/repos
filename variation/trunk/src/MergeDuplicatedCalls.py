#!/usr/bin/env python
"""
Usage: MergeDuplicatedCalls.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	input file
	-o ...,	output file
	-e ...,	ecotype table 'ecotype'(default)
	-p ...,	ecotype_duplicate2tg_ecotypeid_table, 'ecotype_duplicate2tg_ecotypeid'(default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	MergeDuplicatedCalls.py -i justin_data.csv -o justin_data_filtered.csv
	
Description:
	Merge call data from duplicated strains.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import getopt, math, numpy
from pymodule import process_function_arguments, write_data_matrix
from variation.src.common import get_ecotypeid2nativename
from variation.src.dbSNP2data import dbSNP2data

class MergeDuplicatedCalls(object):
	def __init__(self, **keywords):
		"""
		2008-4-4
		"""
		argument_default_dict = {('hostname',1, ):'localhost',\
								('dbname',1, ):'stock',\
								('schema',0, ):'',\
								('input_fname',1, ):None,\
								('output_fname',1, ):None,\
								('ecotype_table',1, ):'ecotype',\
								('ecotype_duplicate2tg_ecotypeid_table',1,):'ecotype_duplicate2tg_ecotypeid',\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=__doc__, class_to_have_attr=self)
	
	def get_ecotype_duplicate2tg_ecotypeid(cls, curs, ecotype_duplicate2tg_ecotypeid_table):
		"""
		2008-04-04
		"""
		sys.stderr.write("Getting ecotype_duplicate2tg_ecotypeid ... ")
		ecotype_duplicate2tg_ecotypeid = {}
		curs.execute("select ecotypeid, duplicate, tg_ecotypeid from %s"%ecotype_duplicate2tg_ecotypeid_table)
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, duplicate, tg_ecotypeid = row
			key_pair = (ecotypeid, duplicate)
			if key_pair in ecotype_duplicate2tg_ecotypeid:
				sys.stderr.write("Warning: %s already appears in ecotype_duplicate2tg_ecotypeid.\n"%(repr(key_pair)))
			ecotype_duplicate2tg_ecotypeid[key_pair] = tg_ecotypeid
		sys.stderr.write("Done.\n")
		return ecotype_duplicate2tg_ecotypeid
	
	get_ecotype_duplicate2tg_ecotypeid = classmethod(get_ecotype_duplicate2tg_ecotypeid)
	
	def get_tg_ecotypeid2ecotypeid_duplicate_index_ls(cls, strain_acc_list, category_list, ecotype_duplicate2tg_ecotypeid):
		"""
		2008-04-04
		"""
		sys.stderr.write("Getting tg_ecotypeid2ecotypeid_duplicate_index_ls ... ")
		tg_ecotypeid2ecotypeid_duplicate_index_ls = {}
		for i in range(len(strain_acc_list)):
			ecotypeid = int(strain_acc_list[i])
			duplicate = int(category_list[i])
			key_pair = (ecotypeid, duplicate)
			tg_ecotypeid = ecotype_duplicate2tg_ecotypeid.get(key_pair)
			if tg_ecotypeid==None:
				sys.stderr.write("Error: %s not found in ecotype_duplicate2tg_ecotypeid.\n"%(key_pair))
				sys.exit(2)
			if tg_ecotypeid not in tg_ecotypeid2ecotypeid_duplicate_index_ls:
				tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid] = []
			tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid].append(i)
		sys.stderr.write("Done.\n")
		return tg_ecotypeid2ecotypeid_duplicate_index_ls
	
	get_tg_ecotypeid2ecotypeid_duplicate_index_ls = classmethod(get_tg_ecotypeid2ecotypeid_duplicate_index_ls)
	
	def get_merged_matrix(self, tg_ecotypeid2ecotypeid_duplicate_index_ls, data_matrix):
		"""
		"""
		sys.stderr.write("Merging calls from duplicates ... ")
		tg_ecotypeid_ls = tg_ecotypeid2ecotypeid_duplicate_index_ls.keys()
		tg_ecotypeid_ls.sort()
		
		no_of_cols = len(data_matrix[0])
		merge_matrix = numpy.zeros([len(tg_ecotypeid_ls), no_of_cols], numpy.int)
		for i in range(len(tg_ecotypeid_ls)):
			tg_ecotypeid = tg_ecotypeid_ls[i]
			ecotypeid_duplicate_index_ls = tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid]
			if len(ecotypeid_duplicate_index_ls)==1:	#no merging needed. just copy it over
				merge_matrix[i] = data_matrix[ecotypeid_duplicate_index_ls[0]]
			else:
				merge_matrix[i] = self.merge_call_on_one_row(ecotypeid_duplicate_index_ls, data_matrix, no_of_cols)
		sys.stderr.write("Done.\n")
		return tg_ecotypeid_ls, merge_matrix
	
	def merge_call_on_one_row(cls, ecotypeid_duplicate_index_ls, data_matrix, no_of_cols):
		"""
		"""
		one_row = numpy.zeros(no_of_cols)
		for i in range(no_of_cols):
			call_counter_ls = [0]*11
			for index in ecotypeid_duplicate_index_ls:
				call_number = data_matrix[index][i]
				if call_number !=0:	#dont' need NA
					call_counter_ls[call_number] += 1
			one_row[i] = dbSNP2data.get_majority_call_number(call_counter_ls)
		return one_row
	merge_call_on_one_row = classmethod(merge_call_on_one_row)
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		
		ecotype_duplicate2tg_ecotypeid = self.get_ecotype_duplicate2tg_ecotypeid(curs, self.ecotype_duplicate2tg_ecotypeid_table)
		
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_fname)
		
		tg_ecotypeid2ecotypeid_duplicate_index_ls = self.get_tg_ecotypeid2ecotypeid_duplicate_index_ls(strain_acc_list, category_list, ecotype_duplicate2tg_ecotypeid)
		
		tg_ecotypeid_ls, merge_matrix = self.get_merged_matrix(tg_ecotypeid2ecotypeid_duplicate_index_ls, data_matrix)
		
		ecotypeid2nativename = get_ecotypeid2nativename(curs, ecotype_table=self.ecotype_table)
		tg_nativename_ls = []
		for ecotypeid in tg_ecotypeid_ls:
			tg_nativename_ls.append(ecotypeid2nativename[ecotypeid])
		header[1] = 'nativename'
		write_data_matrix(merge_matrix, self.output_fname, header, tg_ecotypeid_ls, tg_nativename_ls)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:e:f:br", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = None
	input_fname = None
	output_fname = None
	ecotype_table = None
	ecotype_duplicate2tg_ecotypeid_table = None
	debug = None
	report = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
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
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-f",):
			ecotype_duplicate2tg_ecotypeid_table = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = MergeDuplicatedCalls(hostname=hostname, dbname=dbname, schema=schema, input_fname=input_fname, \
					output_fname=output_fname, ecotype_table=ecotype_table, ecotype_duplicate2tg_ecotypeid_table=ecotype_duplicate2tg_ecotypeid_table, \
					debug=debug, report=report)
	instance.run()