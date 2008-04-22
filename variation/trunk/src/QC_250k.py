#!/usr/bin/env python
"""

Examples:
	#test run (without -c) calculate missing rate for call_info entries. no need for cmp_data_filename. 'nothing' is a place holder.
	QC_250k.py -i nothing -m 0
	
	#test run (without -c) QC between 250k and 2010
	QC_250k.py -i /home/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv -m 1
	
	#test run (without -c) QC between 250k and perlegen
	QC_250k.py -i /home/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv -m 2
	
	#test run (without -c) QC between 250k and 149SNP
	QC_250k.py -i /home/crocea/script/variation/stock20080403/data_y10001101.tsv -m 3
	
Description:
	QC for 250k call data from call_info_table against 2010, perlegen, 149SNP data.
	It will select the call_info entries that haven't been QCed for a particular QC method.
	The output is in the call_QC_table with NA_rate, mismatch_rate and etc.
	
	QC_method_id:
		Except 0, it corresponds to field/column id in table QC_method.
		0 is not a real QC method id. It calculates NA rate for call_info entries whose NA_rates haven't been calculated.
	
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
from pymodule import process_function_arguments, turn_option_default_dict2argument_default_dict
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
								('SNPData2', 1, ): None,\
								('curs', 0, ): None,\
								('snp_locus_table_250k', 0, ): 'stock_250k.snps',\
								('snp_locus_table_149snp', 0, ): 'stock.snps',\
								('QC_method_id', 1, int):1}
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
		if self.QC_method_id==3:	#149SNP data is SNPData2. use database to find out which SNP matches which
			if self.curs==None:
				sys.stderr.write("Error: no database connection but it's required to link SNP ids.\n")
				sys.exit(3)
			from variation.src.Cmp250kVs149SNP import Cmp250kVs149SNP
			self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = Cmp250kVs149SNP.get_col_matching_dstruc(self.SNPData1.header, self.SNPData2.header, self.curs, self.snp_locus_table_250k, self.snp_locus_table_149snp)
		else:	#use the default from QualityControl
			self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.SNPData1.header, self.SNPData2.header)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.SNPData1.strain_acc_list, self.SNPData1.category_list, self.SNPData2.strain_acc_list)
		
	def cmp_row_wise(self):
		return QualityControl.cmp_row_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)



class QC_250k(object):
	__doc__ = __doc__
	option_default_dict = {('z', 'hostname', 1, 'hostname of the db server', 1, ): 'papaya.usc.edu',\
							('d', 'dbname', 1, '', 1, ): 'stock_250k',\
							('u', 'user', 1, '', 1, ):None,\
							('p', 'passwd', 1, '', 1, ):None,\
							('t', 'call_info_table', 1, '', 1, ): 'call_info',\
							('i', 'cmp_data_filename', 1, 'the data file to be compared with.', 1, ): None,\
							('q', 'call_QC_table', 1, '', 1, ): 'call_QC',\
							('m', 'QC_method_id', 1, 'id in table QC_method', 1, int): None,\
							('c', 'commit', 0, 'commit db transaction', 0, int):0,\
							('b', 'debug', 0, 'toggle debug mode', 0, int):0,\
							('r', 'report', 0, 'toggle report, more verbose stdout/stderr.', 0, int):0}
	
	"""
	2008-04-40
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	def __init__(self,  **keywords):
		"""
		2008-04-20
		"""
		argument_default_dict = turn_option_default_dict2argument_default_dict(self.option_default_dict)
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments
			the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_call_info_id2fname(cls, curs, call_info_table, call_QC_table, QC_method_id, array_info_table='array_info'):
		"""
		2008-04-20
			call_info entries that QC haven't been done.
			and its corresponding maternal_ecotype_id=paternal_ecotype_id (no crosses).
			
		"""
		sys.stderr.write("Getting call_info_id2fname ... ")
		curs.execute("select distinct c.id, c.filename, a.maternal_ecotype_id, a.paternal_ecotype_id, q.id as qc_id, q.QC_method_id \
				from %s a, %s c left join %s q on c.id=q.call_info_id where c.array_id=a.id and a.maternal_ecotype_id=a.paternal_ecotype_id"%\
				(array_info_table, call_info_table, call_QC_table))
		rows = curs.fetchall()
		call_info_id2fname = {}
		call_info_id_del_ls = []
		for row in rows:
			call_info_id, fname, maternal_ecotype_id, paternal_ecotype_id, qc_id, db_QC_method_id = row
			if db_QC_method_id==QC_method_id:	#done already, mark it for deletion
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
			data_insert_ls = [row_id[0]] + NA_mismatch_ls + [QC_method_id, user]	#row_id is (call_info_id, ecotypeid)
			curs.execute("insert into " + call_QC_table + " (call_info_id, NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, QC_method_id, created_by)\
				values(%s, %s, %s, %s, %s, %s, %s, %s, %s)", data_insert_ls)
		sys.stderr.write("Done.\n")
	
	submit_to_call_QC = classmethod(submit_to_call_QC)
	
	def cal_independent_NA_rate(cls, curs, call_info_table):
		"""
		2008-04-20
			calculate indepent (no data to be compared) NA rates.
			update it in the db.
		"""
		sys.stderr.write("Calculating indepent NA rate ... \n")
		curs.execute("select c.id, c.filename from %s c where c.NA_rate is null order by id"%\
					(call_info_table))
		rows = curs.fetchall()
		no_of_rows = len(rows)
		sys.stderr.write("\tTotally, %d call_info entries to be processed.\n"%no_of_rows)
		for i in range(no_of_rows):
			sys.stderr.write("%d/%d:\t%s\n"%(i+1, no_of_rows, rows[i][1]))
			call_info_id, fname = rows[i]
			reader = csv.reader(open(fname), delimiter='\t')
			reader.next()	#throw away the first line
			no_of_totals = 0
			no_of_NAs = 0
			for row in reader:
				SNP_id, call = row[:2]
				no_of_totals += 1
				if call=='NA':
					no_of_NAs += 1
			if no_of_totals!=0:
				NA_rate = float(no_of_NAs)/no_of_totals
			else:
				NA_rate = -1
			curs.execute("update " + call_info_table + " set NA_rate=%s where id=%s",\
					(NA_rate, call_info_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		if self.QC_method_id==0:
			self.cal_independent_NA_rate(curs, self.call_info_table)
		else:
			from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.cmp_data_filename)
			snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix)
			
			call_info_id2fname = self.get_call_info_id2fname(curs, self.call_info_table, self.call_QC_table, self.QC_method_id)
			header, strain_acc_list, category_list, data_matrix = self.read_call_matrix(call_info_id2fname)
			snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix)
			twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, curs=curs, QC_method_id=self.QC_method_id)
			
			row_id2NA_mismatch_rate = twoSNPData.cmp_row_wise()
			
			self.submit_to_call_QC(curs, row_id2NA_mismatch_rate, self.call_QC_table, self.QC_method_id, self.user)
		if self.commit:
			curs.execute("commit")


if __name__ == '__main__':
	from pymodule import process_options, generate_program_doc
	opts_dict = process_options(sys.argv, QC_250k.option_default_dict, error_doc=generate_program_doc(sys.argv[0], QC_250k.option_default_dict)+QC_250k.__doc__)
	
	instance = QC_250k(**opts_dict)
	instance.run()
