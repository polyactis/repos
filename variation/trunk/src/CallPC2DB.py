#!/usr/bin/env python
"""

Examples:
	CallPC2DB.py -l 17 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -e /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.eval -u yh -c
	
Description:
	2009-2-2
		1. Put principal component values of ecotypes/strains calculated by smartpca.perl of EIGENSOFT into database.
			(http://banyan.usc.edu/research/variation/log-2008-11 for how to run smartpca.perl)
		2. Put eigen values into db.
		

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt, numpy
import warnings, traceback
import Stock_250kDB
from Association import Association
from pymodule import figureOutDelimiter

class CallPC2DB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("eigen_vector_fname", 1, ): [None, 'f', 1, 'eigen vector file with PCs outputted by smartpca.perl from EIGENSOFT. first column is the ecotype id.'],\
							('eigen_value_fname', 0, ): [None, 'e', 1, 'eigen value file outputted by smartpca.perl from EIGENSOFT', ],\
							("call_method_id", 1, int): [None, 'l', 1, 'from which call method the ecotypes are derived'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2009-2-2
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def putPCDataIntoDB(self, db, PC_data, call_method_id):
		"""
		2009-2-2
		"""
		sys.stderr.write("Putting PC matrix into db ...")
		no_of_rows, no_of_cols = PC_data.PC_matrix.shape
		CallInfo = Stock_250kDB.CallInfo
		for j in range(no_of_cols):
			ecotype_id2count = {}	#for each PC, record how many one ecotype occurs. if >1, take call info id sequentially. ecotype id supposed to be strictly unique, but some backward cleanup makes some call_method have >1 calls with same ecotype id.
			for i in range(no_of_rows):
				ecotype_id = int(PC_data.id_ls[i])
				call_info_ls = CallInfo.query.filter(CallInfo.array.has(maternal_ecotype_id=ecotype_id)).filter(CallInfo.array.has(paternal_ecotype_id=ecotype_id)).\
					filter_by(method_id=call_method_id).all()
				if ecotype_id not in ecotype_id2count:
					ecotype_id2count[ecotype_id] = 0
				ecotype_id2count[ecotype_id] += 1
				if len(call_info_ls)==0:
					sys.stderr.write("Ignore: No call info found for ecotype_id=%s and call_method_id=%s.\n"%(ecotype_id, call_method_id))
					continue
				elif len(call_info_ls)>1:
					sys.stderr.write("Warning: %s call infos found for ecotype_id=%s and call_method_id=%s. Take the %sth one.\n"%\
									(len(call_info_ls), ecotype_id, call_method_id, ecotype_id2count[ecotype_id]))
				try:
					call_info = call_info_ls[ecotype_id2count[ecotype_id]-1]
				except:
					import pdb
					pdb.set_trace()
				pc_db_entry = Stock_250kDB.CallInfoPCValues(which_pc=j+1, pc_value=PC_data.PC_matrix[i][j])
				pc_db_entry.call_info = call_info
				db.session.save(pc_db_entry)
				db.session.flush()
				
		sys.stderr.write("Done.\n")
	
	def putEigenValuesIntoDB(self, db, eigen_value_ls, explained_var, call_method_id):
		"""
		2009-2-2
		"""
		sys.stderr.write("Putting eigen values into db ...")
		for i in range(len(eigen_value_ls)):
			eigen_value = eigen_value_ls[i]
			variance_perc = explained_var[i]
			db_entry = Stock_250kDB.CallMethodEigenValues(which_eigen=i+1, eigen_value=eigen_value, variance_perc=variance_perc, call_method_id=call_method_id)
			db.session.save(db_entry)
			db.session.flush()
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2009-2-2
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		
		session = db.session
		session.begin()
		
		PC_data = Association.getPCFromFile(self.eigen_vector_fname)
		eigen_value_ls = Association.getEigenValueFromFile(self.eigen_value_fname)
		eigen_value_ls = numpy.array(eigen_value_ls)
		explained_var = eigen_value_ls/numpy.sum(eigen_value_ls)
		
		self.putPCDataIntoDB(db, PC_data, self.call_method_id)
		self.putEigenValuesIntoDB(db, eigen_value_ls, explained_var, self.call_method_id)
		if self.commit:
			session.flush()
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CallPC2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()