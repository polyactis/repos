import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


from variation.src.common import get_nt_number2diff_matrix_index, nt2number, number2nt

class QualityControl(object):
	"""
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
	
	get_row_matching_dstruc = classmethod(get_row_matching_dstruc)
	
	def cmp_row_wise(cls, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
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
		for row_id1, row_id2 in row_id12row_id2.iteritems():
			row_index1 = row_id2row_index1[row_id1]
			row_index2 = row_id2row_index2[row_id2]
			no_of_mismatches = 0
			no_of_non_NA_pairs = 0
			no_of_NAs = 0
			no_of_totals = 0
			for col_id1, col_index1 in col_id2col_index1.iteritems():
				if col_id1 in col_id12col_id2:
					col_id2 = col_id12col_id2[col_id1]
					col_index2 = col_id2col_index2[col_id2]
					no_of_totals += 1
					if data_matrix1[row_index1][col_index1] == 0:
						no_of_NAs += 1
					if data_matrix1[row_index1][col_index1] > 0 and data_matrix2[row_index2][col_index2] > 0:	#2008-01-07
						no_of_non_NA_pairs += 1
						if data_matrix1[row_index1][col_index1] != data_matrix2[row_index2][col_index2]:
							no_of_mismatches += 1
			if no_of_totals >0 and no_of_non_NA_pairs>0:
				NA_rate = no_of_NAs/float(no_of_totals)
				mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
				row_id2NA_mismatch_rate[row_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs]
			else:
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
		sys.stderr.write("Done.\n")
		return row_id2NA_mismatch_rate
	
	cmp_row_wise = classmethod(cmp_row_wise)
	
	def _cal_pairwise_dist(self, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
		2007-12-21
		2008-01-07
			change the non_NA criteria from !=0 to >0. in 2010 data matrix, there's substantial amount of -2 (not tried).
		"""
		sys.stderr.write("Calculating pairwise distance ...")
		row_id2pairwise_dist = {}
		counter = 0
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
				else:
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
			pairwise_dist.sort()
			row_id2pairwise_dist[row_id1] = pairwise_dist
		sys.stderr.write("Done.\n")
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
	
	def get_NA_rate_for_one_col(self, data_matrix, col_index):
		"""
		2008-05-06
			calculate independent no_of_NAs, no_of_totals
		"""
		no_of_NAs = 0
		no_of_totals = 0
		for i in range(len(data_matrix)):
			no_of_totals += 1
			if data_matrix[i][col_index] == 0:
				no_of_NAs += 1
		return no_of_NAs, no_of_totals
	
	def cmp_one_col(cls, data_matrix1, data_matrix2, col_index1, col_index2, row_id2row_index1, row_id2row_index2, row_id12row_id2, mapping_for_data_matrix1=None):
		"""
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
				if value1 == 0:
					no_of_NAs += 1
				if value1 > 0 and data_matrix2[row_index2][col_index2] > 0:
					no_of_non_NA_pairs += 1
					if value1 != data_matrix2[row_index2][col_index2]:
						no_of_mismatches += 1
		return no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs
	
	cmp_one_col = classmethod(cmp_one_col)
	
	def cmp_col_wise(cls, data_matrix1, data_matrix2, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2):
		"""
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
		for col_id1, col_id2 in col_id12col_id2.iteritems():
			col_index1 = col_id2col_index1[col_id1]
			col_index2 = col_id2col_index2[col_id2]
			no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs = QualityControl.cmp_one_col(data_matrix1, data_matrix2, col_index1, col_index2, row_id2row_index1, row_id2row_index2, row_id12row_id2)
			if no_of_totals >0 and no_of_non_NA_pairs>0:
				NA_rate = no_of_NAs/float(no_of_totals)
				mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
				col_id2NA_mismatch_rate[col_id1] = [NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs]
			else:
				if hasattr(cls, 'debug') and getattr(cls,'debug'):
					sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(col_id1, col_id2))
		sys.stderr.write("Done.\n")
		return col_id2NA_mismatch_rate
	
	cmp_col_wise = classmethod(cmp_col_wise)
	
	def get_diff_matrix(self, data_matrix1, data_matrix2, nt_number2diff_matrix_index, col_id2col_index1, col_id2col_index2, col_id12col_id2, row_id2row_index1, row_id2row_index2, row_id12row_id2, row_id=-1, col_id=-1, need_diff_code_pair_dict=0):
		"""
		2008-01-24 add flag need_diff_code_pair_dict to output diff_code_pair2diff_details_ls
		2008-01-01 derived from cmp_two_matricies() of CmpAccession2Ecotype.py
		"""
		import numpy
		if self.report or self.debug:
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
		if self.report or self.debug:
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
		self.diff_matrix, self.diff_details_ls, extra_data = self.get_diff_matrix(self.data_matrix1, self.data_matrix2, self.nt_number2diff_matrix_index, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		print self.diff_matrix
		i = 0
		if self.latex_output_fname:
			outf = open(self.latex_output_fname, 'w')
			outf.write('\\section{Summary} \\label{section_summary}\n')
			from pymodule.latex import outputMatrixInLatexTable, escape_characters
			wrapped_diff_matrix = self.wrap_diff_matrix_with_row_col_names(self.diff_matrix)
			table_label = 'table_dm%s'%i
			outf.write(outputMatrixInLatexTable(wrapped_diff_matrix, '%s vs %s'%(os.path.basename(self.input_fname1), os.path.basename(self.input_fname2)), table_label))
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
				diff_matrix_for_one_row, diff_details_ls_for_one_row, extra_data = self.get_diff_matrix(self.data_matrix1, self.data_matrix2, self.nt_number2diff_matrix_index, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2, row_id=row_id1)
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
				diff_matrix_for_one_col, diff_details_ls_for_one_col, extra_data = self.get_diff_matrix(self.data_matrix1, self.data_matrix2, self.nt_number2diff_matrix_index, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2, col_id=col_id1)
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
				curs.execute("SELECT c.comment, e.nativename, e.stockparent FROM %s c, %s e where e.id=c.ecotypeid and e.id=%s and c.duplicate=%s"%(calls_250k_duplicate_comment_table, ecotype_table, ecotypeid, duplicate))
				rows = curs.fetchall()
				if rows:
					comment, nativename, stockparent = rows[0]
					directory = os.path.split(comment)[0]	#take the 1st
					directory = os.path.split(directory)[-1]	#take the last
					row_id2info[row_id] = '%s,%s'%(nativename,directory)
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
		self.row_id2NA_mismatch_rate = self.cmp_row_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 )
		self.row_id2info = self.get_row_id2info(self.row_id2NA_mismatch_rate.keys(), self.curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype')
		from QCVisualize import QCVisualize
		import gtk
		QCVisualize_ins = QCVisualize(self.row_id2NA_mismatch_rate, title, id2info=self.row_id2info, id2index=self.row_id2row_index1, id_is_strain=1, header=self.header1, strain_acc_list=self.strain_acc_list1, category_list=self.category_list1, data_matrix=self.data_matrix1)
		QCVisualize_ins.show_all()
		gtk.main()
	
	def plot_col_NA_mismatch_rate(self, title=''):
		"""
		2007-12-20
		"""
		self.col_id2NA_mismatch_rate = self.cmp_col_wise(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
		from QCVisualize import QCVisualize
		import gtk
		QCVisualize_ins = QCVisualize(self.col_id2NA_mismatch_rate, title, id2info={}, id2index=self.col_id2col_index1, id_is_strain=0, header=self.header1, strain_acc_list=self.strain_acc_list1, category_list=self.category_list1, data_matrix=self.data_matrix1)
		QCVisualize_ins.show_all()
		gtk.main()
	
	def cal_row_id2pairwise_dist(self):
		"""
		2008-01-11
			add the part to submit the result to db
		"""
		self.row_id2pairwise_dist = self._cal_pairwise_dist(self.data_matrix1, self.data_matrix2, self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2, self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2)
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
				sql_string = "insert into %s(ecotype_id, duplicate, col_id1, call1, accession_id, col_id2, call2) values(%s, %s, '%s', %s, %s, '%s', %s)"%(diff_details_table, ecotype_id, duplicate, col_id1, number1, row_id2, col_id2, number2)
			else:
				sql_string = "insert into %s(ecotype_id, col_id1, call1, accession_id, col_id2, call2) values(%s, '%s', %s, %s, '%s', %s)"%(diff_details_table, row_id1, col_id1, number1, row_id2, col_id2, number2)
			curs.execute(sql_string)
		sys.stderr.write("Done.\n")
	
	def create_qc_cross_match_table(self, curs, qc_cross_match_table):
		"""
		2008-01-11
			create the qc_cross_match_table
		"""
		sys.stderr.write("Creating table %s ..."%qc_cross_match_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			ecotype_id	integer	not null,\
			duplicate	integer,\
			accession_id	integer,\
			mismatch_rate	float,\
			no_of_mismatches	integer,\
			no_of_non_NA_pairs	integer)"%qc_cross_match_table)
		sys.stderr.write("Done.\n")
	
	def submit_row_id2pairwise_dist(self, curs, qc_cross_match_table, row_id2pairwise_dist):
		"""
		2008-01-11
			submit row_id2pairwise_dist
		"""
		sys.stderr.write('Submitting row_id2pairwise_dist ...')
		for row_id, pairwise_dist_ls in row_id2pairwise_dist.iteritems():
			for pairwise_dist in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = pairwise_dist
				accession_id = row_id2
				if type(row_id)==tuple:
					ecotype_id, duplicate = row_id
					sql_string = "insert into %s(ecotype_id, duplicate, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs) values(%s, %s, %s, %s, %s, %s)"%(qc_cross_match_table, ecotype_id, duplicate, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs)
				else:
					ecotype_id = row_id
					sql_string = "insert into %s(ecotype_id, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs) values(%s, %s, %s, %s, %s)"%(qc_cross_match_table, ecotype_id, accession_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs)
				curs.execute(sql_string)
		sys.stderr.write("Done.\n")
