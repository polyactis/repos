#!/usr/bin/env python
"""

Examples:
	FigureOut384IlluminaABMapping.py -i genotyping/384-illumina.tsv -j genotyping/250k_l3_w0.15_x0.20_y0.85.tsv -c
	
Description:
	program to figure out mapping between 384-illumina's 'AB' notation and 'ACGT'
	
	2008-05-06

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
from variation.src.dbsnp import DBSNP, SNPs, README, SNPsABAlleleMapping
from pymodule.db import formReadmeObj
from variation.src.QC_250k import QC_250k, SNPData, TwoSNPData
import sqlalchemy

def get_snps_id2mapping(hostname, dbname=None, user=None, passwd=None, readme_id=2):
	"""
	2008-05-07
		only readme_id=2
	2008-05-06
	"""
	sys.stderr.write("Getting snps_id2mapping ... ")
	db = DBSNP(username=user,
			password=passwd, hostname=hostname, database=dbname)
	session = db.session
	snps_id2mapping = {}
	rows = session.query(SNPsABAlleleMapping).filter_by(readme_id=readme_id).order_by(db.tables['snps_ab_allele_mapping'].c.snps_id).list()
	for i in range(0, len(rows), 2):
		a = rows[i]
		b = rows[i+1]
		chosen_obj = a
		if a.mismatch_rate>b.mismatch_rate:
			chosen_obj = b
		mapping = {'N': 'NA', 'A':chosen_obj.allele_A_nt, 'B':chosen_obj.allele_B_nt, 'H':chosen_obj.allele_A_nt+chosen_obj.allele_B_nt}
		snps_id2mapping[a.snps_id] = mapping
	sys.stderr.write("Done.\n")
	return snps_id2mapping

class TwoSNPData384(TwoSNPData):
	"""
	QC between SNPData1 and SNPData2. The final NA_rate, mismatch_rate is in terms of row_id in SNPData1.
	"""
	def __init__(self, **keywords):
		argument_default_dict = {('SNPData1', 1, ): None,\
								('SNPData2', 1, ): None,\
								('curs', 0, ): None,\
								('QC_method_id', 1, int):3,\
								('user', 0,): '',\
								('columns_to_be_selected', 0, ):'s1.name, s2.name'}
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)
		self.columns_to_be_selected = 's1.name, s2.name'
		self.QC_method_id = 3	#2008-05-06 set to 3 to use the Cmp250kVs149SNP.py's col matching function
		self.update_row_col_matching()
	
	def get_row_matching_dstruc(self, strain_acc_list1, category_list1, strain_acc_list2):
		"""
		2008-05-06
		"""
		return QualityControl.get_row_matching_dstruc(strain_acc_list1, category_list1, strain_acc_list2)
	
	def figureOutABMapping(self, session, readme, snps_name2possible_mappings):
		"""
		2008-05-06
		"""
		sys.stderr.write("Comparing col-wise for mismatches ...\n")
		for col_id1 in self.col_id2col_index1:
			col_id2 = self.col_id12col_id2.get(col_id1)
			if col_id2:
				for mapping in snps_name2possible_mappings[col_id1]:
					s = SNPsABAlleleMapping(created_by=self.user, allele_A_nt=number2nt[mapping[1]], allele_B_nt=number2nt[mapping[2]])
					if type(self.SNPData1.col_id2id)==dict:
						s.snps_id = self.SNPData1.col_id2id[col_id1]
					s.tg_snps_name=col_id2
					col_index1 = self.col_id2col_index1[col_id1]
					col_index2 = self.col_id2col_index2[col_id2]
					s.relative_no_of_NAs, s.relative_no_of_totals, s.no_of_mismatches, s.no_of_non_NA_pairs = self.cmp_one_col(self.SNPData1.data_matrix, self.SNPData2.data_matrix, col_index1, col_index2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2, mapping)
					if s.relative_no_of_totals >0:
						s.relative_NA_rate = s.relative_no_of_NAs/float(s.relative_no_of_totals)
					if s.no_of_non_NA_pairs>0:
						s.mismatch_rate = s.no_of_mismatches/float(s.no_of_non_NA_pairs)
					s.no_of_NAs, s.no_of_totals = self.get_NA_rate_for_one_col(self.SNPData1.data_matrix, col_index1)
					if s.no_of_totals >0:
						s.NA_rate = s.no_of_NAs/float(s.no_of_totals)
					
					s.readme = readme
					session.save(s)
		sys.stderr.write("Done.\n")

class FigureOut384IlluminaABMapping(object):
	__doc__ = __doc__
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['dbsnp', 'd', 1, '', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname1', 1, ): [None, 'i', 1, '384-illumina data matrix', ],\
							('input_fname2', 1, ): [None, 'j', 1, '250k or 2010 or 149SNP data matrix', ],\
							('output_fname', 0, ): [None, 'o', 1, 'useless', ],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_snps_name2possible_mappings(self, db):
		"""
		2008-05-06
		"""
		sys.stderr.write("Getting snps_name2possible_mappings ...")
		snps_ls = db.session.query(SNPs).list()
		snps_name2possible_mappings = {}
		snps_name2snps_id = {}
		for snps in snps_ls:
			snps_name2snps_id[snps.name] = snps.id
			#the mapping is between ab_number and nt_number
			snps_name2possible_mappings[snps.name] = [{0:0, 1:nt2number[snps.allele1], 2:nt2number[snps.allele2], 3:nt2number[snps.allele1+snps.allele2]},\
													{0:0, 1:nt2number[snps.allele2], 2:nt2number[snps.allele1], 3:nt2number[snps.allele1+snps.allele2]}]
		sys.stderr.write("Done.\n")
		return snps_name2possible_mappings, snps_name2snps_id
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = DBSNP(username=self.user,
				   password=self.passwd, hostname=self.hostname, database=self.dbname)
		session = db.session
		transaction = session.create_transaction()
		
		snps_name2possible_mappings, snps_name2snps_id = self.get_snps_name2possible_mappings(db)
		
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_fname1)
		snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix, \
							col_id2id=snps_name2snps_id, snps_table='dbsnp.snps')
				
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_fname2)
		snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix,\
						snps_table='stock_250k.snps')
		
		twoSNPData = TwoSNPData384(SNPData1=snpData1, SNPData2=snpData2, curs=curs, user=self.user)
		
		readme = formReadmeObj(sys.argv, self.ad, README)
		session.save(readme)
		twoSNPData.figureOutABMapping(session, readme, snps_name2possible_mappings)
		if self.commit:
			curs.execute("commit")
			transaction.commit()
		else:
			transaction.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = FigureOut384IlluminaABMapping
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()	