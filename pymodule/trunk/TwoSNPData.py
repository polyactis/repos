#!/usr/bin/env python
"""
Examples:
	TwoSNPData.py -i /Network/Data/250k/db/reference_dataset/data_2010_ecotype_id_y0002_n1c1d110_mergedup.tsv -j /Network/Data/250k/db/reference_dataset/stock_149SNP_y0000110101_mergedup.csv -o /tmp/2010_149.csv -c 1 -w 1
	
Description:
	2009-2-5
		add option row_matching_by_which_value to allow both the 1st and 2nd column from the first file to be retained while matching
			rows from the 2nd file using designated column.
			
	Merge two SNPData.
"""
import sys, os, csv, traceback

from sets import Set
from SNP import get_nt_number2diff_matrix_index, nt2number, number2nt, NA_set, SNPData
from utils import importNumericArray
num = importNumericArray()

class QualityControl(object):
	"""
	2008-05-31 moved from variation.src.QualityControl
	2007-12-19
		an abstract class for more specific classes to inherit
		the functionality is to do quality control between different types of SNP data
	"""
	def __init__(self, **keywords):
		if keywords.has_key('report'):
			self.report = int(keywords['report'])
		else:
			self.report = 0
		if keywords.has_key('debug'):
			self.debug = int(keywords['debug'])
		else:
			self.debug = 0
		if keywords.has_key('latex_output_fname'):
			self.latex_output_fname = int(keywords['latex_output_fname'])
		else:
			self.latex_output_fname = ''
		
		self.diff_details_table = ''
		self.qc_cross_match_table = ''
	
	def calculate_row_NA_rate(self, strain_acc_list, category_list, data_matrix):
		"""
		2007-12-14
		"""
		sys.stderr.write("Calculating row NA_rate ...")
		ecotypeid_duplicate_NA_rate_ls = []
		for i in range(len(strain_acc_list)):
			ecotypeid = int(strain_acc_list[i])
			duplicate = int(category_list[i])
			no_of_NAs = 0
			row = data_matrix[i]
			no_of_SNPs = len(row)
			for call in row:
				if call==0:
					no_of_NAs += 1
			NA_rate = no_of_NAs/float(no_of_SNPs)
			ecotypeid_duplicate_NA_rate_ls.append([ecotypeid, duplicate, NA_rate])
		sys.stderr.write("Done.\n")
		return ecotypeid_duplicate_NA_rate_ls
	
	def calculate_col_NA_rate(self, col_id_ls, data_matrix):
		"""
		2007-12-18
			data_matrix is a list embedded with list
		"""
		sys.stderr.write("Calculating col NA_rate ...")
		col_id_NA_rate_ls = []
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			no_of_NAs = 0
			no_of_rows = len(data_matrix)
			for j in range(no_of_rows):
				call = data_matrix[j][i]
				if call==0:
					no_of_NAs += 1
			NA_rate = no_of_NAs/float(no_of_rows)
			col_id_NA_rate_ls.append([col_id, NA_rate])
		sys.stderr.write("Done.\n")
		return col_id_NA_rate_ls
	
	def get_col_matching_dstruc(self, header1, header2):
		"""
		2008-01-01
			default version. matching by same names
			copied from get_col_matching_dstruc() of Cmp250kVs2010.py
		"""
		sys.stderr.write("Getting col matching dstruc ...\n")
		snpid2col_index1 = {}
		for i in range(2, len(header1)):
			snpid = header1[i]
			snpid2col_index1[snpid] = i-2
		snpid2col_index2 = {}
		for i in range(2, len(header2)):
			snpid = header2[i]
			snpid2col_index2[snpid] = i-2
		
		col_id12col_id2 = {}
		for snpid in snpid2col_index1:
			if snpid in snpid2col_index2:
				col_id12col_id2[snpid] = snpid
		sys.stderr.write("Done.\n")
		return snpid2col_index1, snpid2col_index2, col_id12col_id2
	
	@classmethod
	def get_row_matching_dstruc(cls, strain_acc_list1, category_list1, strain_acc_list2):
		"""
		2008-01-01
			default version.
			copied from get_row_matching_dstruc() of Cmp250kVs2010.py
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			ecotypeid = int(strain_acc_list1[i])
			duplicate = int(category_list1[i])
			strain_acc = (ecotypeid, duplicate)
			strain_acc2row_index1[strain_acc] = i
		
		strain_acc2row_index2 = {}
		for i in range(len(strain_acc_list2)):
			strain_acc = strain_acc_list2[i]
			ecotypeid = int(strain_acc)
			strain_acc2row_index2[ecotypeid] = i
		
		row_id12row_id2 = {}
		for strain_acc in strain_acc2row_index1:
			ecotypeid, duplicate = strain_acc
			if ecotypeid in strain_acc2row_index2:
				row_id12row_id2[strain_acc] = ecotypeid
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2
	
	
	@classmethod
	def get_NA_rate_for_one_row(cls, data_matrix, row_index):
		"""
		2008-05-12
		"""
		return QualityControl.get_NA_rate_for_one_slice(data_matrix, row_index)
	
	
	@classmethod
	def cmp_one_row(cls, data_matrix1, data_matrix2, row_index1, row_index2, col_id2col_index1, col_id2col_index2, col_id12col_id2, \
				mapping_for_data_matrix1=None):
		"""
		2008-05-12
			use NA_set
		2008-05-12
			split from cmp_row_wise(), for one row
		"""
		no_of_mismatches = 0
		no_of_non_NA_pairs = 0
		no_of_NAs = 0
		no_of_totals = 0
		for col_id1, col_index1 in col_id2col_index1.iteritems():
			if col_id1 in col_id12col_id2:
				col_id2 = col_id12col_id2[col_id1]
				col_index2 = col_id2col_index2[col_id2]
				no_of_totals += 1
				if data_matrix1[row_index1][col_index1] in NA_set:
					no_of_NAs += 1
				if data_matrix1[row_index1][col_index1] not in NA_set and data_matrix2[row_index2][col_index2] not in NA_set:	#2008-01-07
					no_of_non_NA_pairs += 1
					if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
						no_of_mismatches += 1
		if no_of_totals >0:
			NA_rate = no_of_NAs/float(no_of_totals)
		else:
			NA_rate = -1
		if no_of_non_NA_pairs>0:
			mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
		else:
			mismatch_rate = -1
		return NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs
	
	@classmethod
	def cmp_row_wise(cls, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, \
					row_id2row_index2, row_id12row_id2):
		"""
		2009-06-10
			calculate some runtime stats
		2008-05-12
			function split up
		2007-12-18
			strain wise
		2007-12-20
			make it more generic
		2007-12-21 even more generic
		2008-01-07
			change the non_NA criteria from !=0 to >0. in 2010 data matrix, there's substantial amount of -2 (not tried).
		"""
		sys.stderr.write("Comparing row-wise for mismatches ...\n")
		row_id2NA_mismatch_rate = {}
		counter = 0
		counter_no_valid_pairs = 0
		counter_no_relative_valid_pairs = 0
		counter_no_valid_non_NA_pairs = 0
		for row_id1, row_id2 in row_id12row_id2.iteritems():
			counter += 1
			row_index1 = row_id2row_index1[row_id1]
			row_index2 = row_id2row_index2[row_id2]
			NA_rate, no_of_NAs, no_of_totals = QualityControl.get_NA_rate_for_one_row(data_matrix1, row_index1)
			if NA_rate==-1:
				counter_no_valid_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid no_of_totals=0 between %s and %s.\n"%(row_id1, row_id2))
				"""
			relative_NA_rate, mismatch_rate, relative_no_of_NAs, relative_no_of_totals, no_of_mismatches, no_of_non_NA_pairs = \
				QualityControl.cmp_one_row(data_matrix1, data_matrix2, row_index1, row_index2, col_id2col_index1, col_id2col_index2, col_id12col_id2)
			if relative_NA_rate==-1:
				counter_no_relative_valid_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid relative_no_of_totals=0 between %s and %s.\n"%(row_id1, row_id2))
				"""
			if mismatch_rate==-1:
				counter_no_valid_non_NA_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
				"""
			row_id2NA_mismatch_rate[row_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs, \
											relative_NA_rate, relative_no_of_NAs, relative_no_of_totals]
		
		sys.stderr.write("%s cmps in total. %s with no valid pairs. %s with no relative valid pairs. %s with no non-NA pairs. Done.\n"%\
						(counter, counter_no_valid_pairs, counter_no_relative_valid_pairs, counter_no_valid_non_NA_pairs))
		return row_id2NA_mismatch_rate
	
	
	def _cal_pairwise_dist(self, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, \
						row_id2row_index2, row_id12row_id2):
		"""
		2009-06-10
			record how many pairs with no_of_non_NA_pairs>0
		2008-08-28
			report no valid pairs only when self.debug=1
		2007-12-21
		2008-01-07
			change the non_NA criteria from !=0 to >0. in 2010 data matrix, there's substantial amount of -2 (not tried).
		"""
		sys.stderr.write("Calculating pairwise distance ...")
		row_id2pairwise_dist = {}
		counter = 0
		counter_with_valid_non_NA_pairs = 0
		for row_id1, row_index1 in row_id2row_index1.iteritems():
			pairwise_dist = []
			for row_id2, row_index2 in row_id2row_index2.iteritems():
				no_of_mismatches = 0
				no_of_non_NA_pairs = 0
				for col_id1, col_id2 in col_id12col_id2.iteritems():
					col_index1 = col_id2col_index1[col_id1]
					col_index2 = col_id2col_index2[col_id2]
					if data_matrix1[row_index1][col_index1]>0 and data_matrix2[row_index2][col_index2]>0:
						no_of_non_NA_pairs += 1
						if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
							no_of_mismatches += 1
				if no_of_non_NA_pairs>0:
					mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
					pairwise_dist.append([mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs])
					counter_with_valid_non_NA_pairs += 1
				#elif self.debug:
				#	sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
				counter += 1
			pairwise_dist.sort()
			row_id2pairwise_dist[row_id1] = pairwise_dist
		sys.stderr.write("%s/%s have no_of_non_NA_pairs>0 Done.\n"%(counter_with_valid_non_NA_pairs, counter))
		return row_id2pairwise_dist
	
	def trim_row_id2pairwise_dist(self, row_id2pairwise_dist, min_no_of_non_NA_pairs=10):
		"""
		2007-12-26
			used to throw away unreliable pairwise comparisons
		"""
		new_row_id2pairwise_dist = {}
		for row_id, pairwise_dist_ls in row_id2pairwise_dist.iteritems():
			new_pairwise_dist_ls = []
			for row in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = row
				if no_of_non_NA_pairs>=min_no_of_non_NA_pairs:
					new_pairwise_dist_ls.append(row)
			new_row_id2pairwise_dist[row_id] = new_pairwise_dist_ls
		return new_row_id2pairwise_dist
	
	@classmethod
	def get_NA_rate_for_one_slice(cls, data_matrix, slice_index, by_col=0):
		"""
		2008-05-12
			for get_NA_rate_for_one_row and get_NA_rate_for_one_col
			use NA_set
		"""
		no_of_NAs = 0
		no_of_totals = 0
		if by_col:
			no_of_cells = len(data_matrix)
		else:
			no_of_cells = len(data_matrix[slice_index])
		for i in range(no_of_cells):
			no_of_totals += 1
			if by_col:
				value = data_matrix[i][slice_index]
			else:
				value = data_matrix[slice_index][i]
			if value in NA_set:
				no_of_NAs += 1
		if no_of_totals >0:
			NA_rate = no_of_NAs/float(no_of_totals)
		else:
			NA_rate = -1
		return NA_rate, no_of_NAs, no_of_totals
	
	
	@classmethod
	def get_NA_rate_for_one_col(cls, data_matrix, col_index):
		"""
		2008-05-06
			calculate independent no_of_NAs, no_of_totals
		"""
		return QualityControl.get_NA_rate_for_one_slice(data_matrix, col_index, by_col=1)
	
	
	@classmethod
	def cmp_one_col(cls, data_matrix1, data_matrix2, col_index1, col_index2, row_id2row_index1, row_id2row_index2, row_id12row_id2, \
				mapping_for_data_matrix1=None):
		"""
		2008-05-12
			use NA_set
		2008-05-06
			add mapping_for_data_matrix1
		2008-05-05
			split from cmp_col_wise(), for one column
		"""
		no_of_mismatches = 0
		no_of_non_NA_pairs = 0
		no_of_NAs = 0
		no_of_totals = 0
		for row_id1, row_index1 in row_id2row_index1.iteritems():
			if row_id1 in row_id12row_id2:
				row_id2 = row_id12row_id2[row_id1]
				row_index2 = row_id2row_index2[row_id2]
				no_of_totals += 1
				if mapping_for_data_matrix1:
					value1 = mapping_for_data_matrix1[data_matrix1[row_index1][col_index1]]
				else:
					value1 = data_matrix1[row_index1][col_index1]
				if value1 in NA_set:
					no_of_NAs += 1
				if value1 not in NA_set and data_matrix2[row_index2][col_index2] not in NA_set:
					no_of_non_NA_pairs += 1
					if value1 != data_matrix2[row_index2][col_index2]:
						no_of_mismatches += 1
		if no_of_totals >0:
			NA_rate = no_of_NAs/float(no_of_totals)
		else:
			NA_rate = -1
		if no_of_non_NA_pairs>0:
			mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
		else:
			mismatch_rate = -1
		return NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs
		
	@classmethod
	def cmp_col_wise(cls, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, \
					row_id2row_index2, row_id12row_id2):
		"""
		2009-06-10
			calculate some runtime stats
		2008-05-12
			return independent/relative NA_rate
		2007-12-18
			SNP wise
		2007-12-20
			make it more generic
		2008-01-07
			change the non_NA criteria from !=0 to >0. in 2010 data matrix, there's substantial amount of -2 (not tried).
		"""
		sys.stderr.write("Comparing col-wise for mismatches ...\n")
		col_id2NA_mismatch_rate = {}
		no_of_rows1 = len(data_matrix1)
		counter = 0
		counter_no_valid_pairs = 0
		counter_no_relative_valid_pairs = 0
		counter_no_valid_non_NA_pairs = 0
		for col_id1, col_id2 in col_id12col_id2.iteritems():
			counter += 0
			col_index1 = col_id2col_index1[col_id1]
			col_index2 = col_id2col_index2[col_id2]
			NA_rate, no_of_NAs, no_of_totals = QualityControl.get_NA_rate_for_one_col(data_matrix1, col_index1)
			if NA_rate==-1:
				counter_no_valid_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid no_of_totals=0 between %s and %s.\n"%(col_id1, col_id2))
				"""
			relative_NA_rate, mismatch_rate, relative_no_of_NAs, relative_no_of_totals, no_of_mismatches, no_of_non_NA_pairs = \
				QualityControl.cmp_one_col(data_matrix1, data_matrix2, col_index1, col_index2, row_id2row_index1, row_id2row_index2, row_id12row_id2)
			if relative_NA_rate==-1:
				counter_no_relative_valid_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid relative_no_of_totals=0 between %s and %s.\n"%(col_id1, col_id2))
				"""
			if mismatch_rate==-1:
				counter_no_valid_non_NA_pairs += 1
				"""
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(col_id1, col_id2))
				"""
			col_id2NA_mismatch_rate[col_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs,  \
											relative_NA_rate, relative_no_of_NAs, relative_no_of_totals]
		sys.stderr.write("%s cmps in total. %s with no valid pairs. %s with no relative valid pairs. %s with no non-NA pairs. Done.\n"%\
						(counter, counter_no_valid_pairs, counter_no_relative_valid_pairs, counter_no_valid_non_NA_pairs))
		return col_id2NA_mismatch_rate
	
	@classmethod
	def get_diff_matrix(cls, data_matrix1, data_matrix2, nt_number2diff_matrix_index, col_id2col_index1, col_id2col_index2, \
					col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2, row_id=-1, col_id=-1, need_diff_code_pair_dict=0):
		"""
		2009-6-12
			become classmethod
		2008-01-24 add flag need_diff_code_pair_dict to output diff_code_pair2diff_details_ls
		2008-01-01 derived from cmp_two_matricies() of CmpAccession2Ecotype.py
		"""
		import numpy
		if (hasattr(cls, 'report') and getattr(cls,'report')) or (hasattr(cls, 'debug') and getattr(cls,'debug')):
			sys.stderr.write("Comparing two matricies ...")
		diff_matrix = numpy.zeros([len(nt_number2diff_matrix_index), len(nt_number2diff_matrix_index)], numpy.integer)
		if row_id!=-1:
			if row_id in row_id12row_id2:
				row_id_ls_to_be_checked = [row_id]
			else:
				return None, None
		else:
			row_id_ls_to_be_checked = row_id12row_id2.keys()
		if col_id!=-1:
			if col_id in col_id12col_id2:
				col_id_ls_to_be_checked = [col_id]
			else:
				return None, None
		else:
			col_id_ls_to_be_checked = col_id12col_id2.keys()
		diff_code_pair2diff_details_ls = {}
		diff_details_ls = []
		for row_id1 in row_id_ls_to_be_checked:
			row_id2 = row_id12row_id2[row_id1]
			row_index1 = row_id2row_index1[row_id1]
			row_index2 = row_id2row_index2[row_id2]
			for col_id1 in col_id_ls_to_be_checked:
				col_id2 = col_id12col_id2[col_id1]
				col_index1 = col_id2col_index1[col_id1]
				col_index2 = col_id2col_index2[col_id2]
				nt1 = data_matrix1[row_index1][col_index1]
				nt2 = data_matrix2[row_index2][col_index2]
				nt_diff_matrix_index1 = nt_number2diff_matrix_index[nt1]
				nt_diff_matrix_index2 = nt_number2diff_matrix_index[nt2]
				diff_matrix[nt_diff_matrix_index1][nt_diff_matrix_index2] += 1
				if need_diff_code_pair_dict:
					diff_code_pair = (nt1, nt2)
					if diff_code_pair not in diff_code_pair2diff_details_ls:
						diff_code_pair2diff_details_ls[diff_code_pair] = []
					diff_code_pair2diff_details_ls[diff_code_pair].append([row_id1, col_id1, nt1, row_id2, col_id2, nt2])
				if nt1!=nt2:
					diff_details_ls.append([row_id1, col_id1, nt1, row_id2, col_id2, nt2])
		if (hasattr(cls, 'report') and getattr(cls,'report')) or (hasattr(cls, 'debug') and getattr(cls,'debug')):
			sys.stderr.write("Done.\n")
		extra_data = [diff_code_pair2diff_details_ls]
		return diff_matrix, diff_details_ls, extra_data
	
	def wrap_diff_matrix_with_row_col_names(self, diff_matrix):
		"""
		2008-01-01
		2008-01-07
			a cleverer way to generate row_name_ls
		"""
		if self.report or self.debug:
			sys.stderr.write("Wrapping diff_matrix with row, column names ...")
		number_nt_ls = []
		for number, nt in number2nt.iteritems():
			number_nt_ls.append([number, nt])
		number_nt_ls.sort()
		row_name_ls = [row[1] for row in number_nt_ls]
		#row_name_ls = ['-', 'NA', 'A', 'C', 'G', 'T', 'AC', 'AG', 'AT', 'CG', 'CT', 'GT']
		wrapped_diff_matrix = [ [''] + row_name_ls]
		for i in range(diff_matrix.shape[0]):
			wrapped_diff_matrix.append([row_name_ls[i]] + diff_matrix[i,:].tolist())
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return wrapped_diff_matrix
	
	def beautify_diff_details_ls(self, diff_details_ls, row_id2info={}):
		"""
		2008-01-01
		2008-01-11
			sort the diff_details_ls and retain the original row_id1 into row_id1_info
		"""
		if self.report or self.debug:
			sys.stderr.write("Beautifying diff_details_ls ...")
		new_diff_details_ls = []
		diff_details_ls.sort()
		for row in diff_details_ls:
			row_id1, col_id1, nt1, row_id2, col_id2, nt2 = row
			if nt1<=0 or nt2<=0:	#skip deletion or missing
				continue
			if row_id1 in row_id2info:
				row_id1_info = '%s/%s'%(repr(row_id1), row_id2info[row_id1])
			else:
				row_id1_info = repr(row_id1)
			new_diff_details_ls.append([row_id1_info, col_id1, number2nt[nt1], row_id2, col_id2, number2nt[nt2]])
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return new_diff_details_ls
	
	def output_diff_matrix(self):
		"""
		2008-01-01
			self.nt_number2diff_matrix_index, self.row_id2info will have to be loaded in ahead.
		2008-01-09 record diff_details_ls into db
		2008-01-11
			add strain-wise and snp-wise diff output
		"""
		self.diff_matrix, self.diff_details_ls, extra_data = self.get_diff_matrix(self.data_matrix1, self.data_matrix2, \
																				self.nt_number2diff_matrix_index, self.col_id2col_index1, \
																				self.col_id2col_index2, self.col_id12col_id2, \
																				self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		print self.diff_matrix
		i = 0
		if self.latex_output_fname:
			outf = open(self.latex_output_fname, 'w')
			outf.write('\\section{Summary} \\label{section_summary}\n')
			from pymodule.latex import outputMatrixInLatexTable, escape_characters
			wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(self.diff_matrix)
			table_label = 'table_dm%s'%i
			outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, '%s vs %s'%\
											(os.path.basename(self.input_fname1), os.path.basename(self.input_fname2)), table_label))
			i += 1
			table_no = i
			
			#output the whole diff_details_ls
			outf.write('\\section{Real Mismatches (deletion/NA excluded)} \\label{section_real_mismatch}\n')
			bea_diff_details_ls = self.beautify_diff_details_ls(self.diff_details_ls, self.row_id2info)
			table_label = 'table_dm%s'%table_no
			caption = 'mismatches (deletion/NA excluded)'
			outf.write(outputMatrixInLatexTable(bea_diff_details_ls, caption, table_label, header_ls=['row id1', 'col id1', 'call1', 'row id2', 'col id2', 'call2']))
			
			#output diff_matrix and diff_details_ls one strain by one
			outf.write('\\section{Diff Matrix for each strain} \\label{section_strain_wise}\n')
			row_id1_ls = self.row_id12row_id2.keys()
			row_id1_ls.sort()
			for row_id1 in row_id1_ls:
				row_id2 = self.row_id12row_id2[row_id1]
				diff_matrix_for_one_row, diff_details_ls_for_one_row, extra_data = self.get_diff_matrix(self.data_matrix1, self.data_matrix2, \
																									self.nt_number2diff_matrix_index, \
																									self.col_id2col_index1, \
																									self.col_id2col_index2, \
																									self.col_id12col_id2, \
																									self.row_id2row_index1, \
																									self.row_id2row_index2, self.row_id12row_id2,\
																									 row_id=row_id1)
				wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(diff_matrix_for_one_row)
				table_no += 1
				table_label = 'table_dm%s'%table_no
				if row_id1 in self.row_id2info:
					row_id1_info = self.row_id2info[row_id1]
				elif type(row_id1)==tuple:
					if row_id1[0] in self.row_id2info:
						row_id1_info = self.row_id2info[row_id1[0]]
					else:
						row_id1_info = ''
				else:
					row_id1_info = ''
				caption = 'row id1=%s (%s) vs row id2=%s'%(row_id1, row_id1_info, row_id2)
				subsection_title = '\\subsection{row id1=%s (%s) vs row id2=%s}\n'%(row_id1, row_id1_info, row_id2)
				subsection_title = escape_characters(subsection_title)
				outf.write(subsection_title)
				outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, caption, table_label))
				if diff_details_ls_for_one_row:
					bea_diff_details_ls = self.beautify_diff_details_ls(diff_details_ls_for_one_row, self.row_id2info)
					table_no += 1
					table_label = 'table_dm%s'%table_no
					caption = 'detailed difference for row id1=%s (%s) vs row id2=%s'%(row_id1, row_id1_info, row_id2)
					if bea_diff_details_ls:
						outf.write(outputMatrixInLatexTable(bea_diff_details_ls, caption, table_label, header_ls=['row id1', 'col id1', 'call1', 'row id2', 'col id2', 'call2']))

			#SNP-wise comparison
			outf.write('\\section{Diff Matrix for each SNP} \\label{section_snp_wise}\n')
			col_id1_ls = self.col_id12col_id2.keys()
			col_id1_ls.sort()
			for col_id1 in col_id1_ls:
				col_id2 = self.col_id12col_id2[col_id1]
				diff_matrix_for_one_col, diff_details_ls_for_one_col, extra_data = self.get_diff_matrix(self.data_matrix1, \
																									self.data_matrix2, self.nt_number2diff_matrix_index, \
																									self.col_id2col_index1, self.col_id2col_index2, \
																									self.col_id12col_id2, self.row_id2row_index1, \
																									self.row_id2row_index2, self.row_id12row_id2, col_id=col_id1)
				wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(diff_matrix_for_one_col)
				table_no += 1
				table_label = 'table_dm%s'%table_no
				caption = 'col id1=%s vs col id2=%s'%(col_id1, col_id2)
				subsection_title = '\\subsection{col id1=%s vs col id2=%s}\n'%(col_id1, col_id2)
				subsection_title = escape_characters(subsection_title)
				outf.write(subsection_title)
				outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, caption, table_label))
				if diff_details_ls_for_one_col:
					bea_diff_details_ls = self.beautify_diff_details_ls(diff_details_ls_for_one_col, self.row_id2info)
					table_no += 1
					table_label = 'table_dm%s'%table_no
					caption = 'detailed difference for col id1=%s vs col id2=%s'%(col_id1, col_id2)
					if bea_diff_details_ls:
						outf.write(outputMatrixInLatexTable(bea_diff_details_ls, caption, table_label, header_ls=['row id1', 'col id1', 'call1', 'row id2', 'col id2', 'call2']))

			del outf
		#2008-01-09 record details into db
		if self.diff_details_table:
			self.create_diff_details_table(self.curs, self.diff_details_table)
			self.submit_diff_details_ls(self.curs, self.diff_details_table, self.diff_details_ls)
	
	def load_dstruc(self):
		"""
		2007-12-21
			let children classes to fill in details
			need to load:
				self.data_matrix1, self.data_matrix2
				self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2
				self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2
		"""
		self.nt_number2diff_matrix_index = get_nt_number2diff_matrix_index(number2nt)

	def get_row_id2info(self, row_id_ls, curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype'):
		"""
		2008-09-18
			calls_250k_duplicate_comment_table is not used anymore.
			only use ecotypeid to find out info of a row
		2007-12-16
		2007-12-20
			make it more generic
		2008-01-01
			'%s'%row_id doesn't work cuz row_id is a tuple. changed to '%s'%repr(row_id)
		"""
		row_id2info = {}
		for row_id in row_id_ls:
			ecotypeid, duplicate = row_id
			try:
				curs.execute("SELECT e.name, e.nativename, e.stockparent FROM %s e where e.id=%s"%(ecotype_table, ecotypeid))
				rows = curs.fetchall()
				if rows:
					name, nativename, stockparent = rows[0]
					#directory = os.path.split(comment)[0]	#take the 1st
					#directory = os.path.split(directory)[-1]	#take the last
					row_id2info[row_id] = '%s,%s,%s'%(nativename,stockparent,name)
					row_id2info[row_id] = row_id2info[row_id].decode('utf-8', 'ignore')
					#ecotypeid_duplicate2info[key_pair] = ecotypeid_duplicate2info[key_pair].decode('latin10')
				else:
					row_id2info[row_id] = '%s'%repr(row_id)
			except:
				import traceback
				traceback.print_exc()
				print sys.exc_info()
				row_id2info[row_id] = '%s'%repr(row_id)
		return row_id2info
	
	def on_click_row(self, event):
		"""
		2007-12-14
		2007-12-16
			need all_ls and ecotypeid_duplicate2info to be global variable
		2007-12-18
			replace all_ls with ecotypeid_duplicate2NA_mismatch_rate
		2007-12-20
			make it more generic
			replace ecotypeid_duplicate2info with row_id2info
			replace ecotypeid_duplicate2NA_mismatch_rate with row_id2NA_mismatch_rate
		2008-01-01 dead
		"""
		# get the x and y coords, flip y from top to bottom
		import pylab
		x, y = event.x, event.y
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for key, value in self.row_id2NA_mismatch_rate.iteritems():
					NA_rate, mismatch_rate = value[:2]
					if abs(NA_rate-event.xdata)<0.005 and abs(mismatch_rate-event.ydata)<0.005:
						pylab.text(event.xdata, event.ydata, self.row_id2info[key], size=8)
						print "row id: %s, NA_mismatch data: %s, info: %s"%(key, value, self.row_id2info[key])
	
	def on_click_col(self, event):
		"""
		2007-12-18
			need ins.snpid2NA_mismatch_rate
		2007-12-20
			make it more generic
			replace snpid2NA_mismatch_rate with col_id2NA_mismatch_rate
		2008-01-01 dead
		"""
		# get the x and y coords, flip y from top to bottom
		import pylab
		x, y = event.x, event.y
		if event.button==1:
			if event.inaxes is not None:
				print 'data coords', event.xdata, event.ydata
				for key, value in self.col_id2NA_mismatch_rate.iteritems():
					NA_rate, mismatch_rate = value[:2]
					if abs(NA_rate-event.xdata)<0.005 and abs(mismatch_rate-event.ydata)<0.005:
						pylab.text(event.xdata, event.ydata, key, size=8)
						print "col id: %s, NA_mismatch data: %s"%(key, value)
	
	def plot_NA_mismatch_rate(self, NA_mismatch_rate_ls, on_click_func, title=''):
		"""
		2007-12-14
		2008-01-01 dead
		"""
		NA_rate_ls = []
		mismatch_rate_ls = []
		for row in NA_mismatch_rate_ls:
			NA_rate, mismatch_rate = row[:2]
			NA_rate_ls.append(NA_rate)
			mismatch_rate_ls.append(mismatch_rate)
		import pylab
		pylab.clf()
		pylab.plot(NA_rate_ls, mismatch_rate_ls, '.')
		#diagonal line give a rough feeling about the notion, more NA, worse calling
		diagonal_start = min(min(NA_rate_ls), min(mismatch_rate_ls))-0.1
		diagonal_end = max(max(NA_rate_ls), max(mismatch_rate_ls))+0.1
		pylab.plot([diagonal_start, diagonal_end],[diagonal_start, diagonal_end])
		if title:
			pylab.title(title)
		pylab.xlabel('NA rate')
		pylab.ylabel('mismatch rate')
		pylab.show()
		pylab.connect('button_press_event', on_click_func)
	
	def plot_row_NA_mismatch_rate(self, title=''):
		"""
		2007-12-20
		"""
		self.row_id2NA_mismatch_rate = self.cmp_row_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, \
														self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 )
		self.row_id2info = self.get_row_id2info(self.row_id2NA_mismatch_rate.keys(), self.curs, \
											calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype')
		from QCVisualize import QCVisualize
		import gtk
		QCVisualize_ins = QCVisualize(self.row_id2NA_mismatch_rate, title, id2info=self.row_id2info, id2index=self.row_id2row_index1, \
									id_is_strain=1, header=self.header1, strain_acc_list=self.strain_acc_list1, category_list=self.category_list1, \
									data_matrix=self.data_matrix1)
		QCVisualize_ins.show_all()
		gtk.main()
	
	def plot_col_NA_mismatch_rate(self, title=''):
		"""
		2007-12-20
		"""
		self.col_id2NA_mismatch_rate = self.cmp_col_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, \
														self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		from QCVisualize import QCVisualize
		import gtk
		QCVisualize_ins = QCVisualize(self.col_id2NA_mismatch_rate, title, id2info={}, id2index=self.col_id2col_index1, id_is_strain=0, \
									header=self.header1, strain_acc_list=self.strain_acc_list1, category_list=self.category_list1, \
									data_matrix=self.data_matrix1)
		QCVisualize_ins.show_all()
		gtk.main()
	
	def cal_row_id2pairwise_dist(self):
		"""
		2008-01-11
			add the part to submit the result to db
		"""
		self.row_id2pairwise_dist = self._cal_pairwise_dist(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, \
														self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		if self.qc_cross_match_table:
			self.create_qc_cross_match_table(self.curs, self.qc_cross_match_table)
			self.submit_row_id2pairwise_dist(self.curs, self.qc_cross_match_table, self.row_id2pairwise_dist)

	def create_diff_details_table(self, curs, diff_details_table):
		"""
		2008-01-09
			store the diff_details_ls into db
		"""
		sys.stderr.write("Creating diff_details_table ...")
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			ecotype_id	integer	not null,\
			duplicate	integer,\
			col_id1	varchar(200) not null,\
			call1	integer,\
			accession_id	integer not null,\
			col_id2	varchar(200) not null,\
			call2	integer)"%diff_details_table)
		sys.stderr.write("Done.\n")
	
	def submit_diff_details_ls(self, curs, diff_details_table, diff_details_ls):
		"""
		2008-01-09
		"""
		sys.stderr.write('Submitting diff_details_ls ...')
		for row in diff_details_ls:
			row_id1, col_id1, number1, row_id2, col_id2, number2 = row
			if type(row_id1)==tuple:
				ecotype_id, duplicate = row_id1
				sql_string = "insert into %s(ecotype_id, duplicate, col_id1, call1, accession_id, col_id2, call2) values(%s, %s, '%s', %s, %s, '%s', %s)"%\
				(diff_details_table, ecotype_id, duplicate, col_id1, number1, row_id2, col_id2, number2)
			else:
				sql_string = "insert into %s(ecotype_id, col_id1, call1, accession_id, col_id2, call2) values(%s, '%s', %s, %s, '%s', %s)"%\
				(diff_details_table, row_id1, col_id1, number1, row_id2, col_id2, number2)
			curs.execute(sql_string)
		sys.stderr.write("Done.\n")
	
	def create_qc_cross_match_table(self, curs, qc_cross_match_table):
		"""
		2009-6-9
			add more columns (qc_method_id, created_by, updated_by, date_created, date_updated)
		2008-07-01
			rename some fields in the cross_match_table
				duplicate => call_info_id
				accession_id => vs_ecotype_id
		2008-01-11
			create the qc_cross_match_table
		"""
		sys.stderr.write("Creating table %s ..."%qc_cross_match_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			ecotype_id	integer	not null,\
			call_info_id	integer,\
			vs_ecotype_id	integer,\
			mismatch_rate	float,\
			no_of_mismatches	integer,\
			no_of_non_NA_pairs	integer,\
			qc_method_id	integer,\
			created_by varchar(200) default NULL,\
			updated_by varchar(200) default NULL,\
			date_created timestamp NOT NULL default CURRENT_TIMESTAMP,\
			date_updated timestamp NOT NULL default '0000-00-00 00:00:00')engine=INNODB"%qc_cross_match_table)
		sys.stderr.write("Done.\n")
	
	def submit_row_id2pairwise_dist(self, curs, qc_cross_match_table, row_id2pairwise_dist, qc_method_id=0, max_mismatch_rate=0.25,\
								min_no_of_non_NA_pairs=5):
		"""
		2009-6-9
			add argument qc_method_id, max_mismatch_rate, min_no_of_non_NA_pairs
			
			max_mismatch_rate & min_no_of_non_NA_pairs are used to filter entries whose mismatch_rate is above that
			
		2009-4-24
			if type(row_id)==tuple:
				check if the entry is already in db before insertion
		2009-3-13
			wrap the insertion sql into "try ... except ..." to allow insertion error (non-unique entries etc)
		2008-07-01
			some fields in qc_cross_match_table renamed
		2008-01-11
			submit row_id2pairwise_dist
		"""
		sys.stderr.write('Submitting row_id2pairwise_dist ...')
		counter = 0
		real_counter = 0
		no_of_entries_in_db = 0
		for row_id, pairwise_dist_ls in row_id2pairwise_dist.iteritems():
			for pairwise_dist in pairwise_dist_ls:
				counter += 1
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = pairwise_dist
				if max_mismatch_rate is not None and max_mismatch_rate>0.0 and mismatch_rate>max_mismatch_rate:
					continue
				if min_no_of_non_NA_pairs is not None and no_of_non_NA_pairs<min_no_of_non_NA_pairs:
					continue
				accession_id = row_id2
				if type(row_id)==tuple:
					ecotype_id, duplicate = row_id
					#2009-4-24 check if it's already in db.
					curs.execute("select id from %s where ecotype_id=%s and call_info_id=%s and vs_ecotype_id=%s and qc_method_id=%s"%\
						(qc_cross_match_table, ecotype_id, duplicate, accession_id, qc_method_id))
					rows = curs.fetchall()
					if len(rows)>0:
						no_of_entries_in_db += 1
						continue
					
					sql_string = "insert into %s(ecotype_id, call_info_id, vs_ecotype_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, \
						qc_method_id) values(%s, %s, %s, %s, %s, %s, %s)"%\
						(qc_cross_match_table, ecotype_id, duplicate, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, qc_method_id)
				else:
					ecotype_id = row_id
					sql_string = "insert into %s(ecotype_id, vs_ecotype_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, qc_method_id) \
							values(%s, %s, %s, %s, %s, %s)"%\
							(qc_cross_match_table, ecotype_id, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, qc_method_id)
				try:
					curs.execute(sql_string)
					real_counter += 1
				except:
					sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
					traceback.print_exc()
		sys.stderr.write("%s/%s already in db. %s/%s inserted into db. Done.\n"%(no_of_entries_in_db, counter, real_counter, counter))


class TwoSNPData(QualityControl):
	"""
	2009-9-24
		add argument (read the pre-CONDITION. table is limited to the shape described in QualityControl.create_qc_cross_match_table().
			max_mismatch_rate: only cross-matching results with mismatch_rate<=max_mismatch_rate are to be stored into db
			min_no_of_non_NA_pairs:  only cross-matching results with no_of_non_NA_pairs above it are to be stored into db
			
			One pre-CONDITION is that TwoSNPData.qc_cross_match_table has to be specified. Otherwise, nothing would go into db
				and these two arguments are USELESS. 
			Toggle TwoSNPData.new_QC_cross_match_table to True to tell the program to construct the table beforehand.
			Table definition is in QualityControl.create_qc_cross_match_table().
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
							('max_mismatch_rate', 0, float): 0.0,\
							('min_no_of_non_NA_pairs', 0, int): 5,\
							('debug', 0, ): 0}
	def __init__(self, **keywords):
		from __init__ import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, \
														howto_deal_with_required_none=2)
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
		2008-08-18
			all reference datasets used in 250k QC use chromsome_position as column header
		2008-05-19
			add row_id2row_index and col_id2col_index to each SNPData
		2008-05-11
			fake two headers from col_id_ls
		"""
		"""
		#2008-08-18
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
		"""
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.SNPData1.col_id_ls, \
																										self.SNPData2.col_id_ls)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.SNPData1.row_id_ls, \
																										self.SNPData2.row_id_ls)
		
		self.SNPData1.row_id2row_index = self.row_id2row_index1
		self.SNPData1.col_id2col_index = self.col_id2col_index1
		self.SNPData2.row_id2row_index = self.row_id2row_index2
		self.SNPData2.col_id2col_index = self.col_id2col_index2
	
	def cmpOneRow(self, row_id1, row_id2):
		"""
		2008-08-28
			for convenience
		"""
		row_index1 = self.SNPData1.row_id2row_index[row_id1]
		row_index2 = self.SNPData2.row_id2row_index[row_id2]
		return QualityControl.cmp_one_row(self.SNPData1.data_matrix, self.SNPData2.data_matrix, row_index1, row_index2, self.col_id2col_index1, \
										self.col_id2col_index2, self.col_id12col_id2)
	
	def cmp_row_wise(self):
		return QualityControl.cmp_row_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, \
										self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
	
	def cmp_col_wise(self, row_id=None, row_id12row_id2=None):
		"""
		2008-10-11
			add option row_id and row_id12row_id2 to do column matching for only one strain and/or within certain pairs
		"""
		if row_id==None:
			row_id2row_index1 = self.row_id2row_index1
		else:
			row_id2row_index1 ={}
			row_index = self.row_id2row_index1.get(row_id)
			if row_index is not None:
				row_id2row_index1[row_id] = row_index
			else:
				return {}	#return Nothing
		
		if row_id12row_id2 is None:
			row_id12row_id2 = self.row_id12row_id2
		return QualityControl.cmp_col_wise(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, self.col_id2col_index2, \
										self.col_id12col_id2, row_id2row_index1, self.row_id2row_index2, row_id12row_id2)
	
	def get_diff_matrix(self, row_id=-1, col_id=-1, need_diff_code_pair_dict=0):
		"""
		2009-6-12
			a wrap up of the namesake verbose method from super class
		"""
		self.nt_number2diff_matrix_index = get_nt_number2diff_matrix_index(number2nt)
		return QualityControl.get_diff_matrix(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.nt_number2diff_matrix_index, \
									self.col_id2col_index1, self.col_id2col_index2, \
									self.col_id12col_id2, self.row_id2row_index1, \
									self.row_id2row_index2, self.row_id12row_id2, \
									row_id=row_id, col_id=col_id, need_diff_code_pair_dict=need_diff_code_pair_dict)
		
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
							call_method_id=self.SNPData1.call_method_id, created_by=self.user, \
							max_call_info_mismatch_rate=self.SNPData1.max_call_info_mismatch_rate)
				if type(self.SNPData1.col_id2id)==dict:
					snpsqc.snps_id = self.SNPData1.col_id2id[col_id1]
				snpsqc.tg_snps_name=col_id2
				col_index1 = self.col_id2col_index1[col_id1]
				col_index2 = self.col_id2col_index2[col_id2]
				snpsqc.relative_NA_rate, snpsqc.mismatch_rate, snpsqc.relative_no_of_NAs, snpsqc.relative_no_of_totals, \
					snpsqc.no_of_mismatches, snpsqc.no_of_non_NA_pairs = self.cmp_one_col(self.SNPData1.data_matrix, \
																						self.SNPData2.data_matrix, col_index1, col_index2, \
																						self.row_id2row_index1, self.row_id2row_index2, \
																						self.row_id12row_id2)
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
	
	def mergeTwoSNPData_no_union(self, priority=1):
		"""
		2008-06-02
			mergeTwoSNPData() renamed to mergeTwoSNPData_no_union() as mergeTwoSNPData2 below becomes functional and is renamed to mergeTwoSNPData().
			in this function, need_row_id_union and need_col_id_union are useless and removed. it's always taking snpData1's row_id_ls and col_id_ls
			priority is useful.
		2008-05-19
			if need_row_id_union==0, take SNPData1's row_id_ls as newSnpData's.
				otherwise take the union of SNPData1 and SNPData2
			same as need_col_id_union
			not implemented the 'otherwise' yet. in progress in mergeTwoSNPData2()
		"""
		sys.stderr.write("Merging two SNPData ...\n")
		
		newSnpData = SNPData(row_id_ls=self.SNPData1.row_id_ls, col_id_ls=self.SNPData1.col_id_ls)
		no_of_rows, no_of_cols = len(newSnpData.row_id_ls), len(newSnpData.col_id_ls)
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		no_of_replaced_snps = 0
		for i in range(no_of_rows):
			row_id = newSnpData.row_id_ls[i]
			if row_id not in self.row_id12row_id2 or priority==1:
				newSnpData.data_matrix[i] = self.SNPData1.data_matrix[i]
			else:
				row_id2 = self.row_id12row_id2[row_id]
				SNPData2_row_index = self.SNPData2.row_id2row_index[row_id2]
				for j in range(no_of_cols):
					col_id = newSnpData.col_id_ls[j]
					if col_id not in self.col_id12col_id2:
						newSnpData.data_matrix[i][j] = self.SNPData1.data_matrix[i][j]
					else:
						col_id2 = self.col_id12col_id2[col_id]
						SNPData2_col_index = self.SNPData2.col_id2col_index[col_id2]
						newSnpData.data_matrix[i][j] = self.SNPData2.data_matrix[SNPData2_row_index][SNPData2_col_index]
						no_of_replaced_snps += 1
		newSnpData.no_of_replaced_snps = no_of_replaced_snps
		no_of_total_snps = no_of_rows * no_of_cols
		sys.stderr.write("%s out of %s replaced. Done.\n"%(no_of_replaced_snps, no_of_total_snps))
		return newSnpData
	
	def mergeUsingDataFromOnlyOneSNPData(self, snpData_ls, newSnpData, new_row_id2old_row_index, new_col_id2old_col_index, which_snpData_index=0):
		"""
		2008-09-17
			called by mergeTwoSNPData()
			
			which_snpData_index==0:
				use SNPData1 no matter what SNPData2 is
			which_snpData_index==1:
				use SNPData2 no matter what SNPData1 is
		"""
		no_of_rows, no_of_cols = len(newSnpData.row_id_ls), len(newSnpData.col_id_ls)
		for i in range(no_of_rows):
			row_id = newSnpData.row_id_ls[i]
			old_row_index_tuple_ls = new_row_id2old_row_index[row_id]
			for j in range(no_of_cols):
				col_id = newSnpData.col_id_ls[j]
				old_col_index_tuple_ls = new_col_id2old_col_index[col_id]
				if len(old_col_index_tuple_ls)==2 and len(old_row_index_tuple_ls)==2:
					old_col_index_tuple = old_col_index_tuple_ls[which_snpData_index]
					old_row_index_tuple = old_row_index_tuple_ls[which_snpData_index]
					which_snpData, old_row_index = old_row_index_tuple
					which_snpData_based_on_col, old_col_index = old_col_index_tuple
				elif len(old_col_index_tuple_ls)==2:	#both SNPData has this column. which snpData is determined by the SNPData that has row.
					old_row_index_tuple = old_row_index_tuple_ls[0]
					if old_row_index_tuple[0]==which_snpData_index:	#it has to match the SNPData we want
						old_col_index_tuple = old_col_index_tuple_ls[old_row_index_tuple[0]]	#which snpData is determined by the one which has rows
					else:
						continue
				elif len(old_row_index_tuple_ls)==2:	#both SNPData has this row. which snpData is determined by the SNPData that has col.
					old_col_index_tuple = old_col_index_tuple_ls[0]
					if old_col_index_tuple[0]==which_snpData_index:	#it has to match the SNPData we want
						old_row_index_tuple = old_row_index_tuple_ls[old_col_index_tuple[0]]
					else:
						continue
				else:
					old_col_index_tuple = old_col_index_tuple_ls[0]
					if old_col_index_tuple[0]!=which_snpData_index:	#skip if it's not the correct SNPData
						continue
					old_row_index_tuple = old_row_index_tuple_ls[0]
					if old_row_index_tuple[0]!=which_snpData_index:	#skip if it's not the correct SNPData
						continue
				
				which_snpData, old_row_index = old_row_index_tuple
				which_snpData_based_on_col, old_col_index = old_col_index_tuple
				if which_snpData!=which_snpData_based_on_col:	#this is region which both SNPData1 and SNPData2 don't have data
					if self.debug:
						sys.stderr.write("Error: which_snpData differs based on row and col structures, %s vs %s.\n"%\
										(which_snpData, which_snpData_based_on_col))
				else:
					newSnpData.data_matrix[i][j] = snpData_ls[which_snpData].data_matrix[old_row_index][old_col_index]
	
	def mergeUsingDataAccordingToPriority(self, snpData_ls, newSnpData, new_row_id2old_row_index, new_col_id2old_col_index, priority=0):
		"""
		2008-09-17
			split off from mergeTwoSNPData so that mergeTwoSNPData can call mergeUsingDataFromOnlyOneSNPData() as an option
		"""
		no_of_rows, no_of_cols = len(newSnpData.row_id_ls), len(newSnpData.col_id_ls)
		for i in range(no_of_rows):
			row_id = newSnpData.row_id_ls[i]
			old_row_index_tuple_ls = new_row_id2old_row_index[row_id]
			for j in range(no_of_cols):
				col_id = newSnpData.col_id_ls[j]
				old_col_index_tuple_ls = new_col_id2old_col_index[col_id]
				if len(old_col_index_tuple_ls)==2 and len(old_row_index_tuple_ls)==2:
					old_col_index_tuple = old_col_index_tuple_ls[priority-1]
					old_row_index_tuple = old_row_index_tuple_ls[priority-1]
					which_snpData, old_row_index = old_row_index_tuple
					which_snpData_based_on_col, old_col_index = old_col_index_tuple
					if snpData_ls[which_snpData].data_matrix[old_row_index][old_col_index] in NA_set:	#it's NA, take the other one then
						old_col_index_tuple = old_col_index_tuple_ls[abs(priority-2)]
						old_row_index_tuple = old_row_index_tuple_ls[abs(priority-2)]
				elif len(old_col_index_tuple_ls)==2:	#both SNPData has this column. which snpData is determined by the SNPData that has rows
					old_row_index_tuple = old_row_index_tuple_ls[0]
					old_col_index_tuple = old_col_index_tuple_ls[old_row_index_tuple[0]]	#which snpData is determined by the one which has rows
				elif len(old_row_index_tuple_ls)==2:	#both SNPData has this row. which snpData is determined by the SNPData that has cols
					old_col_index_tuple = old_col_index_tuple_ls[0]
					old_row_index_tuple = old_row_index_tuple_ls[old_col_index_tuple[0]]
				else:
					old_col_index_tuple = old_col_index_tuple_ls[0]
					old_row_index_tuple = old_row_index_tuple_ls[0]
				
				which_snpData, old_row_index = old_row_index_tuple
				which_snpData_based_on_col, old_col_index = old_col_index_tuple
				if which_snpData!=which_snpData_based_on_col:	#this is region which both SNPData1 and SNPData2 don't have data
					if self.debug:
						sys.stderr.write("Error: which_snpData differs based on row and col structures, %s vs %s.\n"%\
										(which_snpData, which_snpData_based_on_col))
				else:
					newSnpData.data_matrix[i][j] = snpData_ls[which_snpData].data_matrix[old_row_index][old_col_index]
	
	def mergeTwoSNPData(self, row_id_merge_type=0, col_id_merge_type=0, priority=1):
		"""
		2008-09-17
			add two new priorities:
				priority==3:
					use SNPData1 no matter what SNPData2 is
				priority==4:
					use SNPData2 no matter what SNPData1 is
			fix a bug which happens to row_id_merge_type==2 and col_id_merge_type==2
		2008-06-02
			becomes functional. renamed from mergeTwoSNPData2() to mergeTwoSNPData()
			if row_id_merge_type==0:
				new_row_id_ls = SNPData1.row_id_ls
			elif row_id_merge_type==1:
				new_row_id_ls = union of SNPData1.row_id_ls and SNPData2.row_id_ls
			elif row_id_merge_type==2:
				new_row_id_ls = intersection of SNPData1.row_id_ls and SNPData2.row_id_ls
			
			for overlapping row_ids, take the names/ids from 1st SNPData
			
			ditto for col_id_merge_type
			
			if priority==1:
				SNPData1 overwrites SNPData2 if SNPData1 is not NA at that SNP
			elif priority==2:
				SNPData2 overwrites SNPData1 if SNPData2 is not NA at that SNP
			else:
				error
		2008-05-19
			if row_id_merge_type==0, take SNPData1's row_id_ls as newSnpData's.
				otherwise take the union of SNPData1 and SNPData2
			same as col_id_merge_type
		"""
		sys.stderr.write("Merging two SNPData ...\n")
		
		snpData_ls = [self.SNPData1, self.SNPData2]
		no_of_snpDatas = len(snpData_ls)
		
		#construct new row_id_ls
		new_row_id_ls = []
		new_row_id2old_row_index = {}	#2008-09-17 the old_row_index stores a list of 2 tuples, each tuple = (which SNPData, which row index)
		SNPData2_row_id_checked_set = Set()
		for row_id in snpData_ls[0].row_id_ls:
			if row_id in self.row_id12row_id2:	#present in both data, needs in all different row_id_merge_type
				SNPData2_row_id = self.row_id12row_id2[row_id]
				SNPData2_row_id_checked_set.add(SNPData2_row_id)
				new_row_id2old_row_index[row_id] = [(0, snpData_ls[0].row_id2row_index[row_id])]
				new_row_id2old_row_index[row_id].append((1, snpData_ls[1].row_id2row_index[SNPData2_row_id]))
				new_row_id_ls.append(row_id)
			elif row_id_merge_type!=2:	#only in SNPData1, then include only when it's not insection mode
				new_row_id2old_row_index[row_id] = [(0, snpData_ls[0].row_id2row_index[row_id])]
				new_row_id_ls.append(row_id)
		if row_id_merge_type==1:	#to include remaining row_id from SNPData2 back into new_row_id2old_row_index for union
			for row_id in snpData_ls[1].row_id_ls:
				if row_id not in SNPData2_row_id_checked_set:
					new_row_id2old_row_index[row_id] = [(1, snpData_ls[1].row_id2row_index[row_id])]
					new_row_id_ls.append(row_id)
		
		#construct new col_id_ls
		new_col_id_ls = []
		new_col_id2old_col_index = {}	#2008-09-17 the old_col_index stores a list of 2 tuples, each tuple = (which SNPData, which col index)
		SNPData2_col_id_checked_set = Set()
		for col_id in snpData_ls[0].col_id_ls:
			if col_id in self.col_id12col_id2:	#present in both data
				SNPData2_col_id = self.col_id12col_id2[col_id]
				SNPData2_col_id_checked_set.add(SNPData2_col_id)
				new_col_id2old_col_index[col_id] = [(0, snpData_ls[0].col_id2col_index[col_id])]
				new_col_id2old_col_index[col_id].append((1, snpData_ls[1].col_id2col_index[SNPData2_col_id]))
				new_col_id_ls.append(col_id)
			elif col_id_merge_type!=2:	#only in SNPData1, then include only when it's not insection mode
				new_col_id2old_col_index[col_id] = [(0, snpData_ls[0].col_id2col_index[col_id])]
				new_col_id_ls.append(col_id)
		if col_id_merge_type==1:	#to include remaining col_id from SNPData2 back into new_col_id2old_col_index for union
			for col_id in snpData_ls[1].col_id_ls:
				if col_id not in SNPData2_col_id_checked_set:
					new_col_id2old_col_index[col_id] = [(1, snpData_ls[1].col_id2col_index[col_id])]
					new_col_id_ls.append(col_id)
		
		newSnpData = SNPData(row_id_ls=new_row_id_ls, col_id_ls=new_col_id_ls)
		no_of_rows, no_of_cols = len(newSnpData.row_id_ls), len(newSnpData.col_id_ls)
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		if priority==1 or priority==2:
			self.mergeUsingDataAccordingToPriority(snpData_ls, newSnpData, new_row_id2old_row_index, new_col_id2old_col_index, priority=priority)
		else:
			self.mergeUsingDataFromOnlyOneSNPData(snpData_ls, newSnpData, new_row_id2old_row_index, new_col_id2old_col_index, \
												which_snpData_index=priority-3)
		
		sys.stderr.write("Done.\n")
		return newSnpData
	
	def intersectSNPData1_and_SNPData2_row_wise(self):
		"""
		2008-05-19
		"""
		sys.stderr.write("Intersecting row-wise ...")		
		no_of_rows = len(self.row_id12row_id2)
		no_of_cols = len(self.SNPData1.col_id_ls)
		newSnpData = SNPData(col_id_ls=self.SNPData1.col_id_ls, row_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		row_index = 0
		for i in range(len(self.SNPData1.row_id_ls)):
			row_id = self.SNPData1.row_id_ls[i]
			if row_id in self.row_id12row_id2:
				newSnpData.row_id_ls.append(row_id)
				newSnpData.data_matrix[row_index] = self.SNPData1.data_matrix[i]
				row_index += 1
		newSnpData.no_of_rows_removed = len(self.SNPData1.row_id_ls)-no_of_rows
		sys.stderr.write("%s rows removed. Done.\n"%(newSnpData.no_of_rows_removed))
		return newSnpData
	
	def sampleSNPLociFromSNPData1(self, no_of_loci_to_sample_around_on_each_side):
		"""
		2008-06-02
			sample equal amount of loci around ones that match in both SNPData's
		"""
		sys.stderr.write("Sampling %s SNP loci from SNPData1 around each SNPData2 loci ..."%(no_of_loci_to_sample_around_on_each_side))
		no_of_cols = len(self.SNPData1.col_id_ls)
		"""
		for col_id2 in self.SNPData2.col_id_ls:
			col_id2_split_ls = col_id2.split('_')
			chromosome, position = map(int, col_id2_split_ls)
			for i in 
		"""
		col_indices_wanted_set = Set()		
		for col_id1, col_id2 in self.col_id12col_id2.iteritems():
			col_index1 = self.SNPData1.col_id2col_index[col_id1]
			for i in range(max(0, col_index1-no_of_loci_to_sample_around_on_each_side), min(col_index1+no_of_loci_to_sample_around_on_each_side+1, no_of_cols)):
				col_indices_wanted_set.add(i)
		
		cols_to_be_tossed_out = Set(range(no_of_cols)) - col_indices_wanted_set
		sys.stderr.write("%s loci sampled around %s reference loci. Done.\n"%(len(col_indices_wanted_set), len(self.col_id12col_id2)))
		return cols_to_be_tossed_out
	
	def cal_row_id2pairwise_dist(self):
		"""
		2009-6-9
			pass self.QC_method_id & self.max_mismatch_rate, self.min_no_of_non_NA_pairs to submit_row_id2pairwise_dist()
		2008-07-01
			cross-match the two data matricies
		"""
		self.row_id2pairwise_dist = self._cal_pairwise_dist(self.SNPData1.data_matrix, self.SNPData2.data_matrix, self.col_id2col_index1, \
														self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, \
														self.row_id12row_id2)
		if getattr(self, 'qc_cross_match_table', None):
			if getattr(self, 'new_QC_cross_match_table', None):
				self.create_qc_cross_match_table(self.curs, self.qc_cross_match_table)
			self.submit_row_id2pairwise_dist(self.curs, self.qc_cross_match_table, self.row_id2pairwise_dist, self.QC_method_id, \
											self.max_mismatch_rate, self.min_no_of_non_NA_pairs)
	
	@classmethod
	def output_row_id2NA_mismatch_rate(cls, row_id2NA_mismatch_rate, output_fname, file_1st_open=1):
		"""
		2008-10-11
			copied from variation/src/QC_250k.py
		2008-05-12
			smart way to figure out how to deal with row_id and its labels
		2008-05-12
			tsv => csv fromat
		2008-05-11
			add file_1st_open argument
		2008-04-22
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
		
	@classmethod
	def output_col_id2NA_mismatch_rate_InGWRFormat(cls, col_id2NA_mismatch_rate, output_fname, file_1st_open=1):
		"""
		2008-10-11
			output col_id2NA_mismatch_rate in Genome-Wide-Result format (pymodule.SNP.getGenomeWideResultFromFile() could read this output)
			not in strict chromosome/position order
		"""
		sys.stderr.write("Outputting col_id2NA_mismatch_rate to %s ..."%(output_fname))
		if file_1st_open:
			open_flag = 'w'
		else:
			open_flag = 'a'
		writer = csv.writer(open(output_fname, open_flag))
		col_id_ls = col_id2NA_mismatch_rate.keys()
		col_id_ls.sort()	#keep them in at least in chromosome order. position in string/letter order, not in numerical order but 
		if len(col_id_ls)>0:
			for col_id in col_id_ls:
				NA_mismatch_ls = col_id2NA_mismatch_rate[col_id]
				#NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs,  relative_NA_rate, relative_no_of_NAs, relative_no_of_totals
				mismatch_rate = NA_mismatch_ls[1]
				if mismatch_rate>=0 and mismatch_rate<=1:	#out of range means, not enough non-NA pairs
					chr, pos = col_id.split('_')
					writer.writerow([chr, pos, mismatch_rate])
		del writer
		sys.stderr.write("Done.\n")
		
		
class MergeTwoSNPData(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname1',1, ): [None, 'i', 1, 'to form SNPData1'],\
							('input_fname2',1, ): [None, 'j', 1, 'to form SNPData2'],\
							('row_id_merge_type', 1, int): [0, 'w', 1, '0=SNPData1, 1=union, 2=intersection'],\
							('col_id_merge_type', 1, int): [0, 'c', 1, '0=SNPData1, 1=union, 2=intersection'],\
							('priority', 1, int): [1, 'p', 1, '1=SNPData1 with SNPData2 covering where SNPData1 is NA, 2=SNPData2 with SNPData1 covering where SNPData2 is NA, 3=SNPData1 (not care SNPData2), 4=SNPData2'],\
							('output_fname', 1, ): [None, 'o', 1, 'Final output.', ],\
							('row_matching_by_which_value', 0, int):[0, 'm', 1, 'which column in the input_fname1 should be used to establish row-id linking to input_fname2. 0=both inputs discard the 2nd column and use the 1st column, 1=1st column(input_fname1 keeps both columns), 2=2nd column(ditto).'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
							#('input_fname1_format',1,int): [1, 'k', 1, 'Format of input_fname1. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId. 3=snp X strain with arrayId.'],\
							#('input_fname2_format',1,int): [1, 'l', 1, 'Format of input_fname2. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId'],\
	def __init__(self, **keywords):
		"""
		2008-09-17
			allow priority=3 or 4
		2008-06-02
		"""
		from __init__ import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-06-02
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		if self.row_matching_by_which_value==0:
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1, ignore_2nd_column=1)
		else:
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1)
		snpData2 = SNPData(input_fname=self.input_fname2, turn_into_array=1, ignore_2nd_column=1)
		
		if self.row_matching_by_which_value==1 or self.row_matching_by_which_value==2:
			row_matching_by_which_value = self.row_matching_by_which_value-1
		else:
			row_matching_by_which_value = None
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, debug=self.debug, row_matching_by_which_value=row_matching_by_which_value)
		newSnpData= twoSNPData.mergeTwoSNPData(self.row_id_merge_type, self.col_id_merge_type, self.priority)
		newSnpData.tofile(self.output_fname)

if __name__ == '__main__':
	#do simple intersectSNPData1_and_SNPData2_row_wise()
	from __init__ import ProcessOptions
	main_class = MergeTwoSNPData
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()