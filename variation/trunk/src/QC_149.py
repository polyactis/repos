#!/usr/bin/env python
"""

Examples:
	#QC between 149 and perlegen
	QC_149.py -m 2 -o /tmp/149_vs_perlegen.QC.tsv -c
	
	#QC between 149 and 384
	QC_149.py -m 3 -o /tmp/149_vs_384.QC.tsv -c
	
Description:
	QC for 149SNP data from db stock against 2010, perlegen, etc.
	
	The output is in the call_QC_table with NA_rate, mismatch_rate and etc.
	
	QC_method_id:
		All but 0 correspond to field/column id in table QC_method.
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
from pymodule import ProcessOptions, PassingData, SNPData, TwoSNPData, read_data, nt2number, importNumericArray
from variation.src.QualityControl import QualityControl
from variation.src.StockDB import StockDB, QCMethod, CallQC, Calls, SNPs, README, Ecotype, Strain, EcotypeIDStrainID2TGEcotypeID
from pymodule.db import formReadmeObj

num = importNumericArray()

class QC_149(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, 'database name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_dir', 0, ): [None, 'n', 1, 'If given, call data would be read from this directory/file rather than from call info table.\
								Does not work for QC_method_id=0. The data is sorted according to array_id. No call_info_id available. No db submission.', ],\
							('cmp_data_filename', 0,): [None, 'i', 1, 'the data file to be compared with. if not given, it gets figured out by QC_method_id.'],\
							('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputed into it.'],\
							('QC_method_id', 1, int): [None, 'm', 1, 'id in table QC_method'],\
							('run_type', 1, int): [1, 'e', 1, 'QC on 1=accession-wise or 2=snp-wise'],\
							('max_call_info_mismatch_rate', 0, float,):[-1, 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays to calculate snp-wise QC.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-17
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def submit_to_call_QC(self, session, row_id2NA_mismatch_rate, QC_method_id, user, row_id12row_id2, readme):
		"""
		2008-08-18
			to overwrite submit_to_call_QC from QC_250k.py
		"""
		sys.stderr.write("Submitting row_id2NA_mismatch_rate to database ...")
		row_id_ls = row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		for row_id in row_id_ls:
			NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
			ecotypeid, strainid = row_id	#bug here, order changed.
			target_id = row_id12row_id2[row_id]
			NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, relative_NA_rate, relative_no_of_NAs, relative_no_of_totals = NA_mismatch_ls
			#call_QC stores the relative NA rate. call_info already stores the independent NA rate
			NA_rate, no_of_NAs, no_of_totals = relative_NA_rate, relative_no_of_NAs, relative_no_of_totals
			callqc = CallQC(strainid=strainid, ecotypeid=ecotypeid, target_id=target_id,\
						qc_method_id=QC_method_id, NA_rate=NA_rate, mismatch_rate=mismatch_rate,\
						no_of_NAs=no_of_NAs, no_of_totals=no_of_totals, no_of_mismatches=no_of_mismatches, no_of_non_NA_pairs=no_of_non_NA_pairs,\
						created_by=user)
			callqc.readme = readme
			session.save(callqc)
		sys.stderr.write("Done.\n")
	
	def findOutCmpDataFilename(self, cmp_data_filename, QC_method_id):
		"""
		2008-08-18
			copied from QC_250k.py. can't inherit QC_149 from QC_250k because database classes clash
		"""
		# if cmp_data_filename not specified, try to find in the data_description column in table QC_method.
		if not cmp_data_filename and QC_method_id!=0:
			qm = QCMethod.query.get(QC_method_id)
			if qm and qm.data_description:
				data_description_ls = qm.data_description.split('=')
				if len(data_description_ls)>1:
					cmp_data_filename = qm.data_description.split('=')[1].strip()
					return cmp_data_filename
		#after db query, cmp_data_filename is still nothing, exit program.
		if not cmp_data_filename and QC_method_id!=0:
			sys.stderr.write("cmp_data_filename is still nothing even after db query. please specify it on the commandline.\n")
			sys.exit(3)
	
	def output_row_id2NA_mismatch_rate(self, row_id2NA_mismatch_rate, output_fname, file_1st_open=1):
		"""
		2008-08-18
			copied from QC_250k.py
		"""
		sys.stderr.write("Outputting row_id2NA_mismatch_rate to %s ..."%(output_fname))
		if file_1st_open:
			open_flag = 'w'
		else:
			open_flag = 'a'
		writer = csv.writer(open(output_fname, open_flag))
		NA_mismatch_ls_header = ['NA_rate', 'mismatch_rate', 'no_of_NAs', 'no_of_totals', \
				'no_of_mismatches', 'no_of_non_NA_pairs', 'relative_NA_rate', 'relative_no_of_NAs', 'relative_no_of_totals']
		row_id_ls = row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		if len(row_id_ls)>0:
			row_id0 = row_id_ls[0]
			if not isinstance(row_id0, str) and hasattr(row_id0, '__len__'):
				header = ['']*len(row_id0)
			else:
				header = ['']
			header += NA_mismatch_ls_header
			writer.writerow(header)
			for row_id in row_id_ls:
				NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
				if isinstance(row_id, tuple):
					row_id_ls = list(row_id)
				elif isinstance(row_id, list):
					row_id_ls = row_id
				else:
					row_id_ls = [row_id]
				writer.writerow(row_id_ls + NA_mismatch_ls)
		del writer
		sys.stderr.write("Done.\n")
	
	def get_strain_id_info(self, QC_method_id):
		"""
		2008-08-18
			to generate data structure related to strain_id, preparation to get data_matrix
			strainid not QCed yet
			link to tg_ecotypeid
		"""
		sys.stderr.write("Getting strain_id info  ... ")
		strain_id2index = {}
		strain_id_list = []
		strain_id2acc = {}
		strain_id2category = {}
		
		rows = Strain.query.all()
		for row in rows:
			ignore_this = 0
			for call_qc in row.call_qc_ls:
				if call_qc.qc_method_id==QC_method_id:	#QC already done
					ignore_this = 1
					break
			if ignore_this:
				continue
			strain_id = row.id
			strain_index = len(strain_id_list)
			strain_id_list.append(strain_id)
			strain_id2index[strain_id] = strain_index
			strain_id2acc[strain_id] = row.ecotypeid_strainid2tg_ecotypeid.tg_ecotypeid
			strain_id2category[strain_id] = strain_id
		passingdata = PassingData(strain_id2index=strain_id2index, strain_id_list=strain_id_list, strain_id2acc=strain_id2acc,\
								strain_id2category=strain_id2category)
		sys.stderr.write("%s strains. Done.\n"%(len(strain_id_list)))
		return passingdata
	
	def get_data_matrix(self, db, strain_id2index, snp_id2index, call_table_name):
		"""
		2008-08-18
		"""
		sys.stderr.write("Getting data_matrix ...\n")
		data_matrix = num.zeros([len(strain_id2index), len(snp_id2index)], num.int8)
		i = 0
		#block_size = 5000
		#rows = Calls.query.offset(i).limit(block_size)
		rows = db.metadata.bind.execute("select * from %s"%call_table_name)
		for row in rows:
			if row.strainid in strain_id2index:
				data_matrix[strain_id2index[row.strainid], snp_id2index[row.snpid]] = nt2number[row.allele]
			i += 1
			#rows = Calls.query.offset(i).limit(block_size)
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def run(self):
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = StockDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname)
		session = db.session
		session.begin()
		
		self.cmp_data_filename = self.findOutCmpDataFilename(self.cmp_data_filename, self.QC_method_id)
		header, strain_acc_list, category_list, data_matrix = read_data(self.cmp_data_filename)
		strain_acc_list = map(int, strain_acc_list)	#it's ecotypeid, cast it to integer to be compatible to the later ecotype_id_ls from db
		snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix)	#category_list is not used.
		
		readme = formReadmeObj(sys.argv, self.ad, README)
		session.save(readme)
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		from dbSNP2data import dbSNP2data
		snp_id2index, snp_id_list, snp_id2info = dbSNP2data.get_snp_id2index_m(curs, Calls.table.name, SNPs.table.name)
		strain_info_data = self.get_strain_id_info(self.QC_method_id)
		data_matrix = self.get_data_matrix(db, strain_info_data.strain_id2index, snp_id2index, Calls.table.name)
		strain_acc_list = [strain_info_data.strain_id2acc[strain_id] for strain_id in strain_info_data.strain_id_list]
		category_list = [strain_info_data.strain_id2category[strain_id] for strain_id in strain_info_data.strain_id_list]
		header = ['ecotypeid', 'strainid']
		for snp_id in snp_id_list:
			snp_name, chromosome, position = snp_id2info[snp_id]
			header.append(snp_name)
		snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix, \
						snps_table='stock.snps')	#snps_table is set to the stock_250k snps_table
		
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, curs=curs, \
							QC_method_id=self.QC_method_id, user=self.db_user, row_matching_by_which_value=0, debug=self.debug)
		if self.run_type==1:
			row_id2NA_mismatch_rate = twoSNPData.cmp_row_wise()
		elif self.run_type==2:
			#twoSNPData.save_col_wise(session, readme)	#2008-08-18 need to implement a new one for 149SNP
			row_id2NA_mismatch_rate = {}
		else:
			sys.stderr.write("run_type=%s is not supported.\n"%self.run_type)
			sys.exit(5)
		if self.output_fname and self.run_type==1 and row_id2NA_mismatch_rate:
			self.output_row_id2NA_mismatch_rate(row_id2NA_mismatch_rate, self.output_fname)
		
		if self.run_type==1 and self.commit and not self.input_dir and row_id2NA_mismatch_rate:
			#if self.input_dir is given, no db submission. call_info_id2fname here is fake, it's actually keyed by (array_id, ecotypeid)
			#row_id2NA_mismatch_rate might be None if it's method 0.
			self.submit_to_call_QC(session, row_id2NA_mismatch_rate, self.QC_method_id, self.db_user, \
								twoSNPData.row_id12row_id2, readme)
		if self.commit:
			session.commit()
		else:
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = QC_149
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()