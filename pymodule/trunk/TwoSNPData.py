"""
2008-05-18
"""
from SNP import *

class TwoSNPData(QualityControl):
	"""
	2008-05-08
		moved from variation.src.QC_250k
	QC between SNPData1 and SNPData2. The final NA_rate, mismatch_rate is in terms of row_id in SNPData1.
	"""
	argument_default_dict = {('SNPData1', 1, ): None,\
							('SNPData2', 1, ): None,\
							('curs', 0, ): None,\
							('snp_locus_table_250k', 0, ): 'stock_250k.snps',\
							('snp_locus_table_149snp', 0, ): 'stock.snps',\
							('QC_method_id', 1, int):1,\
							('user', 0,): '',\
							('columns_to_be_selected', 0, ):'s1.name, s2.snpid',\
							('row_matching_by_which_value', 0, ): None,\
							('col_matching_by_which_value', 0, ): None,\
							('debug', 0, ): 0}
	def __init__(self, **keywords):
		self.ad = ProcessOptions.process_function_arguments(keywords, self.argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)
		if self.QC_method_id!=3 and self.QC_method_id!=7:	#149SNP uses snpid
			self.columns_to_be_selected = 's1.name, s2.name'
		
		self.update_row_col_matching()
	
	def get_row_matching_dstruc(self, row_id_ls1, row_id_ls2, row_matching_by_which_value=None):
		"""
		2008-05-11
			overhauled
		2008-05-09
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		if row_matching_by_which_value!=None:
			row_matching_by_which_value = row_matching_by_which_value
		else:
			row_matching_by_which_value = self.row_matching_by_which_value
		
		row_id2row_index1= {}
		for i in range(len(row_id_ls1)):
			row_id = row_id_ls1[i]
			row_id2row_index1[row_id] = i
		
		row_id2row_index2 = {}
		for i in range(len(row_id_ls2)):
			row_id = row_id_ls2[i]
			row_id2row_index2[row_id] = i
		
		row_id12row_id2 = {}
		for row_id in row_id2row_index1:
			if isinstance(row_matching_by_which_value, int):
				key = row_id[row_matching_by_which_value]
			else:
				key = row_id
			if key in row_id2row_index2:
				row_id12row_id2[row_id] = key
			else:
				if hasattr(self, 'debug') and getattr(self,'debug'):
					sys.stderr.write('Row Matching Failure: %s.\n'% repr(row_id))
		sys.stderr.write("Done.\n")
		return row_id2row_index1, row_id2row_index2, row_id12row_id2
	
	def get_col_matching_dstruc(self, col_id_ls1, col_id_ls2, col_matching_by_which_value=None):
		"""
		2008-05-11
			copied from QualityControl.py
		2008-01-01
			default version. matching by same names
			copied from get_col_matching_dstruc() of Cmp250kVs2010.py
		"""
		sys.stderr.write("Getting col matching dstruc ...\n")
		if col_matching_by_which_value!=None:
			col_matching_by_which_value = col_matching_by_which_value
		else:
			col_matching_by_which_value = self.col_matching_by_which_value
		
		col_id2col_index1 = {}
		for i in range(len(col_id_ls1)):
			col_id = col_id_ls1[i]
			col_id2col_index1[col_id] = i
		
		col_id2col_index2 = {}
		for i in range(len(col_id_ls2)):
			col_id = col_id_ls2[i]
			col_id2col_index2[col_id] = i
		
		col_id12col_id2 = {}
		for col_id in col_id2col_index1:
			if isinstance(col_matching_by_which_value, int):
				key = col_id[col_matching_by_which_value]
			else:
				key = col_id
			if key in col_id2col_index2:
				col_id12col_id2[col_id] = key
			else:
				if hasattr(self, 'debug') and getattr(self,'debug'):
					sys.stderr.write('Col Matching Failure: %s.\n'% repr(col_id))
		sys.stderr.write("Done.\n")
		return col_id2col_index1, col_id2col_index2, col_id12col_id2
	
	def update_row_col_matching(self):
		"""
		2008-05-11
			fake two headers from col_id_ls
		"""
		if self.QC_method_id==3 or self.QC_method_id==7 or self.QC_method_id==8:	#149SNP data is SNPData2. use database to find out which SNP matches which
			if self.curs==None:
				sys.stderr.write("Error: no database connection but it's required to link SNP ids.\n")
				sys.exit(3)
			from variation.src.Cmp250kVs149SNP import Cmp250kVs149SNP
			header1 = ['', ''] + self.SNPData1.col_id_ls
			header2 = ['', ''] + self.SNPData2.col_id_ls
			self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = Cmp250kVs149SNP.get_col_matching_dstruc(header1, \
					header2, self.curs, self.SNPData1.snps_table, self.SNPData2.snps_table, columns_to_be_selected=self.columns_to_be_selected)
		else:	#use the default from QualityControl
			self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.SNPData1.col_id_ls, self.SNPData2.col_id_ls)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.SNPData1.row_id_ls, self.SNPData2.row_id_ls)
		
	def cmp_row_wise(self):
		return QualityControl.cmp_row_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
	
	def cmp_col_wise(self):
		return QualityControl.cmp_col_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
	
	def save_col_wise(self, session, readme):
		"""
		2008-05-12
			cmp_one_col and get_NA_rate_for_one_col changed interface
		"""
		sys.stderr.write("Comparing col-wise for mismatches ...\n")
		for col_id1 in self.col_id2col_index1:
			col_id2 = self.col_id12col_id2.get(col_id1)
			if col_id2:
				snpsqc = SNPsQC(QC_method_id=self.QC_method_id, min_probability=self.SNPData1.min_probability,\
							call_method_id=self.SNPData1.call_method_id, created_by=self.user, max_call_info_mismatch_rate=self.SNPData1.max_call_info_mismatch_rate)
				if type(self.SNPData1.col_id2id)==dict:
					snpsqc.snps_id = self.SNPData1.col_id2id[col_id1]
				snpsqc.tg_snps_name=col_id2
				col_index1 = self.col_id2col_index1[col_id1]
				col_index2 = self.col_id2col_index2[col_id2]
				snpsqc.relative_NA_rate, snpsqc.mismatch_rate, snpsqc.relative_no_of_NAs, snpsqc.relative_no_of_totals, \
					snpsqc.no_of_mismatches, snpsqc.no_of_non_NA_pairs = self.cmp_one_col(self.SNPData1.data_matrix, \
																						self.SNPData2.data_matrix, col_index1, col_index2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
				if snpsqc.relative_NA_rate==-1:
					snpsqc.relative_NA_rate = None
				if snpsqc.mismatch_rate == -1:
					snpsqc.mismatch_rate = None
				snpsqc.NA_rate, snpsqc.no_of_NAs, snpsqc.no_of_totals = self.get_NA_rate_for_one_col(self.SNPData1.data_matrix, col_index1)
				if snpsqc.NA_rate == -1:
					snpsqc.NA_rate = None
				snpsqc.readme = readme
				session.save(snpsqc)
		sys.stderr.write("Done.\n")
		
