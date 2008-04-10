#!/usr/bin/env python
"""
Usage: OutputPhenotype.py [OPTIONS] -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
	-d ..., --dbname=...	the database name, stock_250k(default)
	-u ..., --user=...	the db username, (otherwise it will ask for it).
	-p ..., --passwd=...	the db password, (otherwise it will ask for it).
	-k ..., --schema=...	which schema in the database
	-o ...,	output file
	-e ...,	ecotype table 'stock.ecotype'(default)
	-p ...,	phenotype_avg_table, 'stock_250k.phenotype_avg'(default)
	-m ...,	phenotype_method_table, 'stock_250k.phenotype_method'(default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	OutputPhenotype.py -o /tmp/phenotype.tsv
	
Description:
	program to output phenotype_avg table.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import traceback
from pymodule import process_function_arguments, write_data_matrix

class OutputPhenotype:
	"""
	2008-04-02
	"""
	def __init__(self,  **keywords):
		"""
		2008-4-2
		"""
		argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('output_fname',1, ):None,\
								('ecotype_table',1, ):'stock.ecotype',\
								('phenotype_avg_table',1, ):'stock_250k.phenotype_avg',\
								('phenotype_method_table',1, ):'stock_250k.phenotype_method',\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=__doc__, class_to_have_attr=self)
	
	def get_phenotype_method_id_info(self, curs, phenotype_avg_table, phenotype_method_table ):
		"""
		2008-4-2
		"""
		sys.stderr.write("Getting phenotype_method_id info ... " )
		phenotype_method_id2index = {}	#index of the matrix
		method_id_name_ls = []	#as header for each phenotype
		curs.execute("select distinct method_id from %s order by method_id"%phenotype_avg_table)
		rows = curs.fetchall()
		for row in rows:
			method_id = row[0]
			curs.execute("select short_name from %s where id=%s"%(phenotype_method_table, method_id))
			method_short_name = curs.fetchall()[0][0]
			method_id_name_ls.append('%s_%s'%(method_id, method_short_name))
			phenotype_method_id2index[method_id] = len(phenotype_method_id2index)
		sys.stderr.write("Done\n")
		return phenotype_method_id2index, method_id_name_ls
	
	def get_ecotype_id2info(self, curs, phenotype_avg_table, ecotype_table):
		"""
		2008-4-2
		"""
		sys.stderr.write("Getting ecotype id info ... " )
		ecotype_id2index = {}	#index of the matrix
		ecotype_id_ls = []
		ecotype_name_ls = []
		curs.execute("select distinct ecotype_id from %s order by ecotype_id"%phenotype_avg_table)
		rows = curs.fetchall()
		for row in rows:
			ecotype_id = row[0]
			curs.execute("select nativename from %s where id=%s"%(ecotype_table, ecotype_id))
			nativename = curs.fetchall()[0][0]
			ecotype_name_ls.append(nativename)
			ecotype_id_ls.append(ecotype_id)
			ecotype_id2index[ecotype_id] = len(ecotype_id2index)
		sys.stderr.write("Done\n")
		return ecotype_id2index, ecotype_id_ls, ecotype_name_ls
	
	def get_matrix(self, curs, phenotype_avg_table, ecotype_id2index, phenotype_method_id2index):
		"""
		2008-04-09
			no longer uses numpy matrix. just simple 2-d list.
		2008-4-2
		"""
		sys.stderr.write("Getting matrix ... " )
		#data_matrix = numpy.zeros([len(ecotype_id2index), len(phenotype_method_id2index)], numpy.float)
		data_matrix = [[]]*len(ecotype_id2index)
		for i in range(len(ecotype_id2index)):
			data_matrix[i] = ['NA']*len(phenotype_method_id2index)
		#data_matrix[:] = numpy.nan
		curs.execute("select ecotype_id, method_id, value from %s"%phenotype_avg_table)
		rows = curs.fetchall()
		for row in rows:
			ecotype_id, phenotype_method_id, value = row
			data_matrix[ecotype_id2index[ecotype_id]][phenotype_method_id2index[phenotype_method_id]] = value
		sys.stderr.write("Done\n")
		return data_matrix
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		
		phenotype_method_id2index, method_id_name_ls = self.get_phenotype_method_id_info(curs, self.phenotype_avg_table, self.phenotype_method_table)
		ecotype_id2index, ecotype_id_ls, ecotype_name_ls = self.get_ecotype_id2info(curs, self.phenotype_avg_table, self.ecotype_table)
		data_matrix = self.get_matrix(curs, self.phenotype_avg_table, ecotype_id2index, phenotype_method_id2index)
		
		header = ['ecotype id', 'nativename'] + method_id_name_ls
		write_data_matrix(data_matrix, self.output_fname, header, ecotype_id_ls, ecotype_name_ls)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:o:e:p:m:br", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	output_fname = None
	ecotype_table = None
	phenotype_avg_table = None
	phenotype_method_table = None
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
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-p",):
			phenotype_avg_table = arg
		elif opt in ("-m",):
			phenotype_method_table = arg
		
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = OutputPhenotype(hostname=hostname, dbname=dbname, user=user, passwd=passwd, output_fname=output_fname,
					ecotype_table=ecotype_table, phenotype_avg_table=phenotype_avg_table, \
					phenotype_method_table = phenotype_method_table, debug=debug, report=report)
	instance.run()
