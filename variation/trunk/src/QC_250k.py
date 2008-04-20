#!/usr/bin/env python
"""
Usage: QC_250k.py -i DATA_FILE -m QC_method_id

Argument list:
	-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
	-d ..., --dbname=...	the database name, stock_250k(default)
	-u ..., --user=...	the db username, (otherwise it will ask for it).
	-p ..., --passwd=...	the db password, (otherwise it will ask for it).
	-t ...,	call_info_table, 'call_info'(default)
	-i ...,	cmp_data_filename*, the data file to be compared with
	-q ..., call_QC_table, 'call_QC'(default)
	-m ....,	QC_method_id*, check the id in table QC_method
	-c,	commit db transaction
	-b,	toggle debug
	-r, toggle report

Examples:
	#QC between 250k and 2010
	QC_250k.py -i /home/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv -m 1
	
	#QC between 250k and perlegen
	QC_250k.py -i /home/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv -m 2
	
	#QC between 250k and 149SNP
	QC_250k.py -i /home/crocea/script/variation/stock20080403/data_y10001101.tsv -m 3
	
Description:
	QC for 250k call data from call_info_table against 2010, perlegen, 149SNP data.
	It fills the call_QC_table with NA_rate, mismatch_rate and etc.
	
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
import warnings, traceback
from pymodule import process_function_arguments
from variation.src.QualityControl import QualityControl
from variation.src.common import number2nt, nt2number

class SNPData(object):
	def __init__(self, **keywords):
		argument_default_dict = {('header', 1, ): None,\
								('strain_acc_list', 1, ): None,\
								('category_list', 1, ): None,\
								('data_matrix', 1, ): None}
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)

		

class TwoSNPData(QualityControl):
	"""
	QC between SNPData1 and SNPData2. The final NA_rate, mismatch_rate is in terms of row_id in SNPData1.
	"""
	def __init__(self, **keywords):
		argument_default_dict = {('SNPData1', 1, ): None,\
								('SNPData2', 1, ): None}
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)
		self.update_row_col_matching()
	
	def get_row_matching_dstruc(self, strain_acc_list1, category_list1, strain_acc_list2):
		"""
		2008-04-20
			strain_acc_list1 contains call_info_id
			category_list1 contains ecotypeid corresponding to call_info_id
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			call_info_id = int(strain_acc_list1[i])
			ecotypeid = int(category_list1[i])
			strain_acc = (call_info_id, ecotypeid)
			strain_acc2row_index1[strain_acc] = i
		
		strain_acc2row_index2 = {}
		for i in range(len(strain_acc_list2)):
			strain_acc = strain_acc_list2[i]
			ecotypeid = int(strain_acc)
			strain_acc2row_index2[ecotypeid] = i
		
		row_id12row_id2 = {}
		for strain_acc in strain_acc2row_index1:
			call_info_id, ecotypeid = strain_acc
			if ecotypeid in strain_acc2row_index2:
				row_id12row_id2[strain_acc] = ecotypeid
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2
	
	def update_row_col_matching(self):
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.SNPData1.header, self.SNPData2.header)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.SNPData1.strain_acc_list, self.SNPData1.category_list, self.SNPData2.strain_acc_list)
		
	def cmp_row_wise(self):
		return QualityControl.cmp_row_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)



class QC_250k(object):
	__doc__ = __doc__
	def __init__(self,  **keywords):
		"""
		2008-02-28
		"""
		argument_default_dict = {('hostname', 1, ): 'papaya.usc.edu',\
								('dbname', 1, ): 'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('call_info_table',1, ): 'call_info',\
								('cmp_data_filename',1, ): None,\
								('call_QC_table', 1, ): 'call_QC',\
								('QC_method_id', 1, int): None,\
								('update_all_call', 0, int): 0,\
								('experiment_table',1, ):'at.experiment',\
								('phenotype_table',1, ):'stock_250k.phenotype',\
								('phenotype_avg_table',1, ):'stock_250k.phenotype_avg',\
								('method_table',1, ):'stock_250k.method',\
								('commit',0, int):0,\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_call_info_id2fname(cls, curs, call_info_table, call_QC_table, QC_method_id, array_info_table='array_info'):
		sys.stderr.write("Getting call_info_id2fname ... ")
		curs.execute("select distinct c.id, c.filename, a.maternal_ecotype_id, a.paternal_ecotype_id, q.id as qc_id, q.QC_method_id from %s a, %s c left join %s q on c.id=q.call_info_id where c.array_id=a.id and a.maternal_ecotype_id=a.paternal_ecotype_id"%\
					(array_info_table, call_info_table, call_QC_table))
		rows = curs.fetchall()
		call_info_id2fname = {}
		call_info_id_del_ls = []
		for row in rows:
			call_info_id, fname, maternal_ecotype_id, paternal_ecotype_id, qc_id, db_QC_method_id = row
			if db_QC_method_id==QC_method_id:
				call_info_id_del_ls.append(call_info_id)
			call_info_id2fname[call_info_id] = [maternal_ecotype_id, fname]	#take all
		
		for call_info_id in call_info_id_del_ls:	#remove the ones that already have QC from this QC_method_id
			del call_info_id2fname[call_info_id]
		sys.stderr.write("%s call files. Done.\n"%(len(call_info_id2fname)))
		return call_info_id2fname
	
	get_call_info_id2fname = classmethod(get_call_info_id2fname)
	
	def read_call_matrix(cls, call_info_id2fname):
		"""
		2008-04-20
			fake header, strain_acc_list, category_list, data_matrix
		"""
		sys.stderr.write("Creating call matrix ... ")
		header = ['', '']	#1st and 2nd is header for 1st two columns.
		call_matrix = []
		strain_acc_list = []
		category_list = []
		counter = 0
		for call_info_id in call_info_id2fname:
			ecotype_id, fname = call_info_id2fname[call_info_id]
			strain_acc_list.append(call_info_id)
			category_list.append(ecotype_id)
			reader = csv.reader(open(fname), delimiter='\t')
			reader.next()	#throw away the first line
			data_row = []
			
			for row in reader:
				SNP_id, call = row[:2]
				if counter==0:	#first file
					header.append(SNP_id)
				data_row.append(nt2number[call])
			del reader
			call_matrix.append(data_row)
			counter += 1
		sys.stderr.write("Done.\n")
		return header, strain_acc_list, category_list, call_matrix
	
	read_call_matrix = classmethod(read_call_matrix)
	
	def submit_to_call_QC(cls, curs, row_id2NA_mismatch_rate, call_QC_table, QC_method_id, user):
		sys.stderr.write("Submitting row_id2NA_mismatch_rate to %s ..."%(call_QC_table))
		row_id_ls = row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		for row_id in row_id_ls:
			NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
			data_insert_ls = [row_id[0]] + NA_mismatch_ls + [QC_method_id, user]
			NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs = NA_mismatch_ls	#row_id is (call_info_id, ecotypeid)
			curs.execute("insert into " + call_QC_table + " (call_info_id, NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, QC_method_id, created_by)\
				values(%s, %s, %s, %s, %s, %s, %s, %s, %s)", data_insert_ls)
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.cmp_data_filename)
		snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix)
		
		call_info_id2fname = self.get_call_info_id2fname(curs, self.call_info_table, self.call_QC_table, self.QC_method_id)
		
		header, strain_acc_list, category_list, data_matrix = self.read_call_matrix(call_info_id2fname)
		snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix)
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2)
		
		row_id2NA_mismatch_rate = twoSNPData.cmp_row_wise()
		
		self.submit_to_call_QC(curs, row_id2NA_mismatch_rate, self.call_QC_table, self.QC_method_id, self.user)
		if self.commit:
			curs.execute("commit")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:u:p:i:t:q:m:cbrh", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	call_info_table =None
	cmp_data_filename = None
	call_QC_table = None
	QC_method_id = None
	update_all_call = 0
	commit = 0
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
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-i",):
			cmp_data_filename = arg
		elif opt in ("-m",):
			QC_method_id = arg
		elif opt in ("-t",):
			call_info_table = arg
		elif opt in ("-q",):
			call_QC_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = QC_250k(hostname=hostname, dbname=dbname, user=user, passwd=passwd, call_info_table=call_info_table,\
					cmp_data_filename=cmp_data_filename,\
					call_QC_table=call_QC_table, QC_method_id=QC_method_id, update_all_call=update_all_call,\
					commit=commit, debug=debug, report=report)
	instance.run()
