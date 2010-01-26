#!/usr/bin/env python
"""

Examples:
	OutputPhenotype.py -o /tmp/phenotype.tsv
	
	OutputPhenotype.py -o /tmp/phenotype.tsv -e stock.ecotype_usc
	
	#get raw (un-transformed) phenotype
	OutputPhenotype.py -o /tmp/phenotype_g.tsv -g

Description:
	program to output phenotype_avg table.
	The output format is roughly a ecotype_id X phenotype(shortname) matrix.
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
from pymodule import process_function_arguments, write_data_matrix, PassingData, SNPData

class OutputPhenotype(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('output_fname', 1, ): ['', 'o', 1, 'store the pvalue', ],\
							('phenotype_avg_table',1, ):['stock_250k.phenotype_avg', 'q', 1,  ],\
							('phenotype_method_table',1, ):['stock_250k.phenotype_method', 'm', 1, ],\
							('ecotype_table', 1, ): ['stock.ecotype', 'e', 1, 'ecotype table to get name related to each ecotype', ],\
							('get_raw_data', 0, int):[0, 'g', 0, 'whether to output raw phenotype data from db or transform according to column transformation_description in phenotype_method table'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-11-10
			upgrade option handling to ProcessOptions
		2008-4-2
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		#self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=__doc__, class_to_have_attr=self)
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	@classmethod
	def get_phenotype_method_id_info(cls, curs, phenotype_avg_table, phenotype_method_table ):
		"""
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			do two selects in one
			
		2008-4-2
		"""
		sys.stderr.write("Getting phenotype_method_id info ... " )
		phenotype_method_id2index = {}	#index of the matrix
		method_id_name_ls = []	#as header for each phenotype
		phenotype_id_ls = []
		rows = curs.execute("select m.id, m.short_name, m.transformation_description from %s m, (select distinct method_id from %s) p where m.id=p.method_id order by id"%\
					(phenotype_method_table, phenotype_avg_table))
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		phenotype_method_id2transformation_description = {}
		for row in rows:
			if is_elixirdb:
				method_id = row.id
				method_short_name = row.short_name
				transformation_description = row.transformation_description
			else:
				method_id, method_short_name, transformation_description = row[:3]
			"""
			curs.execute("select short_name, transformation_description from %s where id=%s"%(phenotype_method_table, method_id))
			pm_rows = curs.fetchall()
			method_short_name = pm_rows[0][0]
			transformation_description = pm_rows[0][1]
			"""
			phenotype_id_ls.append(method_id)
			method_id_name_ls.append('%s_%s'%(method_id, method_short_name))
			phenotype_method_id2index[method_id] = len(phenotype_method_id2index)
			if transformation_description=='None':
				transformation_description = None
			phenotype_method_id2transformation_description[method_id] = transformation_description
		return_data = PassingData(phenotype_method_id2index=phenotype_method_id2index, method_id_name_ls=method_id_name_ls,\
								phenotype_id_ls=phenotype_id_ls,\
								phenotype_method_id2transformation_description=phenotype_method_id2transformation_description)
		sys.stderr.write("Done\n")
		return return_data
	
	@classmethod
	def get_ecotype_id2info(cls, curs, phenotype_avg_table, ecotype_table):
		"""
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			do two selects in one
		
		2008-4-2
		"""
		sys.stderr.write("Getting ecotype id info ... " )
		ecotype_id2index = {}	#index of the matrix
		ecotype_id_ls = []
		ecotype_name_ls = []
		rows = curs.execute("select e.id, e.nativename from %s e, (select distinct ecotype_id from %s) p where e.id=p.ecotype_id order by id"%\
					(ecotype_table, phenotype_avg_table))
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		for row in rows:
			if is_elixirdb:
				ecotype_id = row.id
				nativename = row.nativename
			else:
				ecotype_id, nativename = row[:2]
			"""
			curs.execute("select nativename from %s where id=%s"%(ecotype_table, ecotype_id))
			nativename = curs.fetchall()[0][0]
			"""
			ecotype_name_ls.append(nativename)
			ecotype_id_ls.append(ecotype_id)
			ecotype_id2index[ecotype_id] = len(ecotype_id2index)
		sys.stderr.write("Done\n")
		return ecotype_id2index, ecotype_id_ls, ecotype_name_ls
	
	@classmethod
	def get_matrix(cls, curs, phenotype_avg_table, ecotype_id2index, phenotype_info, get_raw_data=0, \
				phenotype_method_table='phenotype_method'):
		"""
		2010-1-26
			Comment out the code below, inserted on 2009-9-2, since its purpose is to avoid truncation and 
				the truncation wasn't due to the value being too small.
			It was that the conversion of a 2D list containing character 'NA' into a numpy array renders the whole 2D list
				being converted into a character numpy array. The conversion was default in write_data_matrix() if
				transform_to_numpy=False is not passed on.  
		2009-9-2
			if value>-5e-7 and value<+5e-7:	#beyond float resolution by a python float
				value = 0
			
			without condition above, values like -5.32907e-15 would be taken as -5.32907e, -3.76545e-12 as -3.76545
			
		2009-9-2
			add phenotype_method_table to get stddev, min_value to do certain transformation involving these two variables
		2009-2-2
			curs could be either MySQLdb cursor or elixirdb.metadata.bind.
			average phenotype values among replicates in the same phenotype method
		2008-11-10
			add code to transform phenotype according to phenotype_info.phenotype_method_id2transformation_description
			add option get_raw_data, if True/1, no transformation.
		2008-04-23
			#some db entries (phenotype_avg.value) have nothing there. convert None to 'NA'
		2008-04-09
			no longer uses numpy matrix. just simple 2-d list.
		2008-4-2
		"""
		sys.stderr.write("Getting matrix ... " )
		#data_matrix = numpy.zeros([len(ecotype_id2index), len(phenotype_method_id2index)], numpy.float)
		data_matrix = [[]]*len(ecotype_id2index)
		for i in range(len(ecotype_id2index)):
			data_matrix[i] = ['NA']*len(phenotype_info.phenotype_method_id2index)
		#data_matrix[:] = numpy.nan
		rows = curs.execute("select pa.ecotype_id, pa.method_id, pa.value, pm.min_value, pm.stddev from %s pa, %s pm where pm.id=pa.method_id"%\
						(phenotype_avg_table, phenotype_method_table))
		is_elixirdb = 1
		if hasattr(curs, 'fetchall'):	#2009-2-2 this curs is not elixirdb.metadata.bind
			rows = curs.fetchall()
			is_elixirdb = 0
		
		for row in rows:
			if is_elixirdb:
				ecotype_id = row.ecotype_id
				phenotype_method_id = row.method_id
				value = row.value
				min_value = row.min_value
				stddev = row.stddev
			else:
				ecotype_id, phenotype_method_id, value, min_value, stddev = row
			if value==None:	#some db entries have nothing there. convert None to 'NA'
				value = 'NA'
			elif not get_raw_data:	#2008-11-10
				transformation_description = phenotype_info.phenotype_method_id2transformation_description.get(phenotype_method_id)
				#if value>-5e-7 and value<+5e-7:	#beyond float resolution by a python float
				#	value = 0
				
				if not transformation_description:
					pass
				elif transformation_description.find('Log(x)')!=-1:
					try:
						value = math.log10(value)
					except:
						sys.stderr.write("Ecotype ID %s, phenotype_method_id %s, value %s.\n"%(ecotype_id, phenotype_method_id, value))
						sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()[0]))
						traceback.print_exc()
						print sys.exc_info()
						#raise sys.exc_info()[0]
						sys.exit(2)
				elif transformation_description=="Log(SD/10+x-minVal)":	#2009-9-1 new transformation
					if min_value is not None and stddev is not None:
						value = math.log10(stddev/10. + value-min_value)
					else:
						value = value
				elif transformation_description=='Log(5+x)':
					value = math.log10(5+value)
				elif transformation_description=='Log(0.5+x)':
					value = math.log10(0.5+value)
				elif transformation_description=='(x-3)':
					value = value-3
			col_index = phenotype_info.phenotype_method_id2index[phenotype_method_id]
			data_matrix[ecotype_id2index[ecotype_id]][col_index] = value
		sys.stderr.write("Done\n")
		return data_matrix
	
	@classmethod
	def getPhenotypeData(cls, curs, phenotype_avg_table=None, phenotype_method_table=None, ecotype_table='stock.ecotype', get_raw_data=1):
		"""
		2009-2-2
			wrap up all other 3 methods
		"""
		phenotype_info = cls.get_phenotype_method_id_info(curs, phenotype_avg_table, phenotype_method_table)
		ecotype_id2index, ecotype_id_ls, ecotype_name_ls = cls.get_ecotype_id2info(curs, phenotype_avg_table, ecotype_table)
		data_matrix = cls.get_matrix(curs, phenotype_avg_table, ecotype_id2index, phenotype_info, get_raw_data)
		pheno_data = SNPData(col_id_ls=phenotype_info.phenotype_id_ls, row_id_ls=ecotype_id_ls, data_matrix=data_matrix)
		pheno_data.row_label_ls = ecotype_name_ls
		pheno_data.col_label_ls = phenotype_info.method_id_name_ls
		return pheno_data
	
	def run(self):
		if self.debug==1:
			import pdb
			pdb.set_trace()
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		
		
		pheno_data = self.getPhenotypeData(curs, self.phenotype_avg_table, self.phenotype_method_table, \
										self.ecotype_table, get_raw_data=self.get_raw_data)
		header = ['ecotype id', 'nativename'] + pheno_data.col_label_ls
		write_data_matrix(pheno_data.data_matrix, self.output_fname, header, pheno_data.row_id_ls, pheno_data.row_label_ls, \
						transform_to_numpy=False)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = OutputPhenotype
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:o:e:q:m:br", long_options_list)
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
			sys.exit(2)
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
		elif opt in ("-q",):
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
	"""