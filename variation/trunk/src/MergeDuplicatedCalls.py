#!/usr/bin/env python
"""

Examples:
	#use ecotype_duplicate2tg_ecotypeid_table
	./src/MergeDuplicatedCalls.py -i genotyping/149snp/stock_149SNP_y0000110101.tsv -o genotyping/149snp/stock_149SNP_y0000110101_mergedup.csv -t ecotype_duplicate2tg_ecotypeid -z localhost
		
	#use no ecotype_duplicate2tg_ecotypeid_table
	./src/MergeDuplicatedCalls.py -i /stock_149SNP_y0000110101.tsv -o /tmp/stock_149SNP_y0000110101_mergedup.csv
Description:
	Merge call data from duplicated strains.
	duplicates could exist between either different (ecotype_duplicate2tg_ecotypeid_table given) or same ecotypeids.
	
	Input format is strain X snp. 1st two columns are ecotypeid, duplicate.
		delimiter (either tab or comma) is automatically detected. Type of representation of nucleotides is also automatically detected.
	
	Output format is strain X snp. 1st two columns are ecotypeid, nativename.
		same delimiter as input would be used.

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
from pymodule import process_function_arguments, write_data_matrix, read_data, PassingData
from pymodule.SNP import pa_has_characters
from variation.src.common import get_ecotypeid2nativename
from variation.src.dbSNP2data import dbSNP2data
from sets import Set

class MergeDuplicatedCalls(object):
	__doc__ = __doc__
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock', 'd', 1, '', ],\
						('user', 1, ): [None, 'u', 1, 'database username', ],\
						('passwd', 1, ):[None, 'p', 1, 'database password', ],\
						('input_fname', 1, ): [None, 'i', 1, ],\
						('output_fname', 1, ): [None, 'o', 1, 'Output Filename'],\
						('ecotype_table', 1, ): ['stock.ecotype', 'e', 1, 'ecotype Table to get ecotypeid2nativename'],\
						('ecotype_duplicate2tg_ecotypeid_table', 0, ):[None, 't', 1, 'table containing who are duplicates to each other. if not given, use ecotypeid to figure out duplicates'],\
						('processing_bits', 1, ): ['0000111100', 'y', 1, 'processing bits to control which processing step should be turned on.\
							default is 10101101. for what each bit stands, see Description.' ],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-4-4
		"""
		#argument dictionary
		#self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=__doc__, class_to_have_attr=self)
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
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
	
	def get_tg_ecotypeid2ecotypeid_duplicate_index_ls(cls, strain_acc_list, category_list, ecotype_duplicate2tg_ecotypeid=None):
		"""
		2008-08-05
			convert ecotypeid or duplicate into integer if there's no character in it.
		2008-07-11
			stop exitting the program if key_pair not found in ecotype_duplicate2tg_ecotypeid.
			use its immediate ecotypeid as target ecotypeid
		2008-05-12
			if ecotype_duplicate2tg_ecotypeid is None, assume duplicates with same ecotypeid but different duplicate
		2008-04-04
		"""
		sys.stderr.write("Getting tg_ecotypeid2ecotypeid_duplicate_index_ls ... ")
		tg_ecotypeid2ecotypeid_duplicate_index_ls = {}
		for i in range(len(strain_acc_list)):
			if pa_has_characters.search(strain_acc_list[i]):
				ecotypeid = strain_acc_list[i]
			else:
				ecotypeid = int(strain_acc_list[i])
			if pa_has_characters.search(category_list[i]):
				duplicate = category_list[i]
			else:
				duplicate = int(category_list[i])
			key_pair = (ecotypeid, duplicate)
			if isinstance(ecotype_duplicate2tg_ecotypeid, dict):
				tg_ecotypeid = ecotype_duplicate2tg_ecotypeid.get(key_pair)
			else:	#either None or not dictionary
				tg_ecotypeid = ecotypeid
			if tg_ecotypeid==None:
				sys.stderr.write("warning: %s not found in ecotype_duplicate2tg_ecotypeid.\n"%(repr(key_pair)))
				#sys.exit(2)
				tg_ecotypeid = ecotypeid
			if tg_ecotypeid not in tg_ecotypeid2ecotypeid_duplicate_index_ls:
				tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid] = []
			tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid].append(i)
		sys.stderr.write("Done.\n")
		return tg_ecotypeid2ecotypeid_duplicate_index_ls
	
	get_tg_ecotypeid2ecotypeid_duplicate_index_ls = classmethod(get_tg_ecotypeid2ecotypeid_duplicate_index_ls)
	
	def get_merged_matrix(self, tg_ecotypeid2ecotypeid_duplicate_index_ls, data_matrix):
		"""
		2008-07-11
			calculate the inconsistency ratio among duplicates
		"""
		sys.stderr.write("Merging calls from duplicates ... ")
		tg_ecotypeid_ls = tg_ecotypeid2ecotypeid_duplicate_index_ls.keys()
		tg_ecotypeid_ls.sort()
		
		no_of_non_NA_pairs = 0
		no_of_non_NA_inconsistent_pairs = 0
		no_of_duplicated_rows = 0
		
		no_of_cols = len(data_matrix[0])
		merge_matrix = numpy.zeros([len(tg_ecotypeid_ls), no_of_cols], numpy.int)
		for i in range(len(tg_ecotypeid_ls)):
			tg_ecotypeid = tg_ecotypeid_ls[i]
			ecotypeid_duplicate_index_ls = tg_ecotypeid2ecotypeid_duplicate_index_ls[tg_ecotypeid]
			if len(ecotypeid_duplicate_index_ls)==1:	#no merging needed. just copy it over
				merge_matrix[i] = data_matrix[ecotypeid_duplicate_index_ls[0]]
			else:
				passingdata = self.merge_call_on_one_row(ecotypeid_duplicate_index_ls, data_matrix, no_of_cols)
				merge_matrix[i] = passingdata.one_row
				no_of_duplicated_rows += 1
				no_of_non_NA_pairs += passingdata.no_of_non_NA_pairs
				no_of_non_NA_inconsistent_pairs += passingdata.no_of_non_NA_inconsistent_pairs
				
				if passingdata.no_of_non_NA_pairs>0:
					inconsistent_ratio = passingdata.no_of_non_NA_inconsistent_pairs/float(passingdata.no_of_non_NA_pairs)
				else:
					inconsistent_ratio = None
				print '%s\t%s\t%s\t%s'%(tg_ecotypeid, passingdata.no_of_non_NA_inconsistent_pairs, passingdata.no_of_non_NA_pairs, inconsistent_ratio)
		sys.stderr.write("%s/%s=%s inconsistency among %s ecotypes who have duplicates. Done.\n"%\
						(no_of_non_NA_inconsistent_pairs, no_of_non_NA_pairs, \
						no_of_non_NA_inconsistent_pairs/float(no_of_non_NA_pairs), no_of_duplicated_rows))
		return tg_ecotypeid_ls, merge_matrix
	
	def merge_call_on_one_row(cls, ecotypeid_duplicate_index_ls, data_matrix, no_of_cols, NA_set=Set([0, -2])):
		"""
		2008-07-11
			calculate the inconsistency ratio among duplicates
		2008-05-12
			-2 is also ruled out, add NA_set
		"""
		one_row = numpy.zeros(no_of_cols)
		passingdata = PassingData()
		passingdata.no_of_non_NA_pairs = 0
		passingdata.no_of_non_NA_inconsistent_pairs = 0
		for i in range(no_of_cols):
			call_counter_ls = [0]*11
			non_NA_call_number_set = Set()
			for index in ecotypeid_duplicate_index_ls:
				call_number = data_matrix[index][i]
				if call_number not in NA_set:	#dont' need NA and non-touched bit
					call_counter_ls[call_number] += 1
					non_NA_call_number_set.add(call_number)
			if len(non_NA_call_number_set)>0:
				passingdata.no_of_non_NA_pairs += 1
				if len(non_NA_call_number_set)>1:
					passingdata.no_of_non_NA_inconsistent_pairs += 1
			one_row[i] = dbSNP2data.get_majority_call_number(call_counter_ls)
		passingdata.one_row = one_row
		return passingdata
	merge_call_on_one_row = classmethod(merge_call_on_one_row)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user=self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		if self.ecotype_duplicate2tg_ecotypeid_table:
			ecotype_duplicate2tg_ecotypeid = self.get_ecotype_duplicate2tg_ecotypeid(curs, self.ecotype_duplicate2tg_ecotypeid_table)
		else:
			ecotype_duplicate2tg_ecotypeid = None
		from pymodule import figureOutDelimiter
		delimiter = figureOutDelimiter(self.input_fname)
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, delimiter=delimiter)
		
		tg_ecotypeid2ecotypeid_duplicate_index_ls = self.get_tg_ecotypeid2ecotypeid_duplicate_index_ls(strain_acc_list, category_list, ecotype_duplicate2tg_ecotypeid)
		
		tg_ecotypeid_ls, merge_matrix = self.get_merged_matrix(tg_ecotypeid2ecotypeid_duplicate_index_ls, data_matrix)
		
		ecotypeid2nativename = get_ecotypeid2nativename(curs, ecotype_table=self.ecotype_table)
		tg_nativename_ls = []
		for ecotypeid in tg_ecotypeid_ls:
			tg_nativename_ls.append(ecotypeid2nativename[ecotypeid])
		header[1] = 'nativename'
		write_data_matrix(merge_matrix, self.output_fname, header, tg_ecotypeid_ls, tg_nativename_ls, delimiter=delimiter)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MergeDuplicatedCalls
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
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
	"""