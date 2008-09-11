#!/usr/bin/env python
"""

Examples:
	#cross QC between 149 and 149
	QC_149_cross_match.py -m 4 -c
	
	#cross QC between 149 and 384
	QC_149_cross_match.py -m 3 -c
	
Description:
	2008-08-27 A program to do cross match between 149SNP data (from db) and others.
	At this moment QC_method_id=4 (self-cross-match) is supported.

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
from variation.src import StockDB
from pymodule.db import formReadmeObj
from QC_149 import QC_149

class QC_149_cross_match(QC_149):
	__doc__ = __doc__
	option_default_dict = QC_149.option_default_dict.copy()
	option_default_dict.update({('input_fname', 0, ):[None, '', 1, 'Get 149SNP data from this file instead of database.']})
	option_default_dict.pop(('input_dir', 0, ))
	option_default_dict.pop(('max_call_info_mismatch_rate', 0, float,))
	
	def __init__(self,  **keywords):
		"""
		2008-08-26
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def submitToQCCrossMatch(self, session, row_id2pairwise_dist, QC_method_id, readme, commit):
		"""
		2008-08-28
			if QC_method_id!=4
				target_id = row_id2
		2008-8-28
			assign target_id based on QC_method_id
		2008-08-26
		"""
		sys.stderr.write("Submitting row_id2pairwise_dist to database ...")
		counter = 0 
		for row_id, pairwise_dist_ls in row_id2pairwise_dist.iteritems():
			for pairwise_dist in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = pairwise_dist
				#the 2nd position in the both row-id tuples is strain id
				if QC_method_id==4:	#the 2nd position in the row-id2 tuple is strain id
					target_id = row_id2[1]
				else:
					target_id = row_id2
				qc_cross_match = StockDB.QCCrossMatch(strainid=row_id[1], target_id=target_id, mismatch_rate=mismatch_rate, no_of_mismatches=no_of_mismatches,\
									no_of_non_NA_pairs=no_of_non_NA_pairs)
				qc_cross_match.qc_method_id = QC_method_id
				qc_cross_match.readme = readme
				session.save(qc_cross_match)
				if commit:
					session.flush()
				counter += 1
		sys.stderr.write("%s entries. Done.\n"%counter)
	
	def prepareTwoSNPData(self, db):
		"""
		2008-09-10
			if self.input_fname is given, get 149SNP data from it , instead of database
		2008-8-28
			split out of run() so that MpiQC149CrossMatch could call this easily
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		if self.input_fname:
			header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
		else:
			from dbSNP2data import dbSNP2data
			snp_id2index, snp_id_list, snp_id2info = dbSNP2data.get_snp_id2index_m(curs, StockDB.Calls.table.name, StockDB.SNPs.table.name)
			strain_info_data = self.get_strain_id_info(self.QC_method_id, ignore_strains_with_qc=False)
			data_matrix = self.get_data_matrix(db, strain_info_data.strain_id2index, snp_id2index, StockDB.Calls.table.name)
			strain_acc_list = [strain_info_data.strain_id2acc[strain_id] for strain_id in strain_info_data.strain_id_list]	#tg_ecotypeid
			category_list = [strain_info_data.strain_id2category[strain_id] for strain_id in strain_info_data.strain_id_list]	#strainid
			header = ['ecotypeid', 'strainid']
			for snp_id in snp_id_list:
				snp_name, chromosome, position = snp_id2info[snp_id]
				header.append(snp_name)
		snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix, \
						snps_table='stock.snps')	#snps_table is set to the stock_250k snps_table
		if self.QC_method_id==4:
			snpData2 = snpData1
		else:
			self.cmp_data_filename = self.findOutCmpDataFilename(self.cmp_data_filename, self.QC_method_id, StockDB.QCMethod)
			header, strain_acc_list, category_list, data_matrix = read_data(self.cmp_data_filename)
			strain_acc_list = map(int, strain_acc_list)	#it's ecotypeid, cast it to integer to be compatible to the later ecotype_id_ls from db
			snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix)	#category_list is not used to facilitate row-id matching
		
		
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, \
							QC_method_id=self.QC_method_id, user=self.db_user, row_matching_by_which_value=0, debug=self.debug)
		return twoSNPData
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		
		readme = formReadmeObj(sys.argv, self.ad, StockDB.README)
		session.save(readme)
		
		twoSNPData = self.prepareTwoSNPData(db)
		twoSNPData.cal_row_id2pairwise_dist()
		
		self.submitToQCCrossMatch(session, twoSNPData.row_id2pairwise_dist, self.QC_method_id, readme, self.commit)
		
		"""
		if self.commit:
			session.commit()
		else:
			session.rollback()
		"""

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = QC_149_cross_match
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()