#!/usr/bin/env python
"""
2008-05-12
	store SNP data structure related stuff
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray, figureOutDelimiter
import copy

num = importNumericArray()

#2008-05-06 ab2number and number2ab is for 384-illumina data
ab2number = {'N': 0,
	'NA': 0,
	'A': 1,
	'B': 2,
	'H': 3}

number2ab = {0: 'NA',
	1: 'A',
	2: 'B',
	3: 'H'}

#2008-05-12	a common NA set
from sets import Set
NA_set = Set([0, 'NA', 'N', -2, '|'])

nt2number = {'|': -2,	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	None: 0,
	'A': 1,
	'C': 2,
	'G': 3,
	'T': 4,
	'AC':5,
	'CA':5,
	'M':5,
	'AG':6,
	'GA':6,
	'R':6,
	'AT':7,
	'TA':7,
	'W':7,
	'CG':8,
	'GC':8,
	'S':8,
	'CT':9,
	'TC':9,
	'Y':9,
	'GT':10,
	'TG':10,
	'K':10
	}

number2nt = {-2: '|',	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	-1: '-',
	0: 'NA',
	1:'A',
	2:'C',
	3:'G',
	4:'T',
	5:'AC',
	6:'AG',
	7:'AT',
	8:'CG',
	9:'CT',
	10:'GT'
	}

number2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(0,255,0), 3:(255,0,0), 4:(122,0,122), 5:(122,122,0), 6:(122,255,255), 7:(122,122,122), 8:(255,122,0), 9:(255,255,122), 10:(122,122,255) }

#2007-04-16 entry[i,j] means whether nucleotide i and j matches. 0(NA) matches everything. singleton(1-4) matches itself and the doublet containing it. doublet(5-10) matches only itself.
nt_number_matching_matrix = [[1, 1,1,1,1,1, 1,1,1,1,1],
	[1, 1,0,0,0,1, 1,1,0,0,0],
	[1, 0,1,0,0,1, 0,0,1,1,0],
	[1, 0,0,1,0,0, 1,0,1,0,1],
	[1, 0,0,0,1,0, 0,1,0,1,1],
	[1, 1,1,0,0,1, 0,0,0,0,0],
	[1, 1,0,1,0,0, 1,0,0,0,0],
	[1, 1,0,0,1,0, 0,1,0,0,0],
	[1, 0,1,1,0,0, 0,0,1,0,0],
	[1, 0,1,0,1,0, 0,0,0,1,0],
	[1, 0,0,1,1,0, 0,0,0,0,1]]

def get_nt_number2diff_matrix_index(number2nt):
	"""
	2008-05-19
		moved from variation.src.common
	2008-01-01 copied from CmpAccession2Ecotype.py
	2007-10-31/2008-01-07
		nucleotide number ranges from -2 to 10.
		the diff_matrix_index ranges from 0 to 12.
	"""
	sys.stderr.write("Getting nt_number2diff_matrix_index from nt2number ...")
	nt_number2diff_matrix_index = {}
	number_nt_ls = []
	for number, nt in number2nt.iteritems():
		number_nt_ls.append([number,nt])
	number_nt_ls.sort()
	for i in range(len(number_nt_ls)):
		nt_number2diff_matrix_index[number_nt_ls[i][0]] = i
	sys.stderr.write("Done.\n")
	return nt_number2diff_matrix_index

def RawSnpsData_ls2SNPData(rawSnpsData_ls, report=0, use_nt2number=0):
	"""
	2008-05-19
		this returns a SNPData in same orientation as rawSnpsData_ls, SNP (row) X Strain (column).
		apply transposeSNPData after this if another orientation is favored.
	2008-05-19
		swap accession id and array id in col_id. now accession_id in 1st
	2008-05-11
		adapts RawSnpsData(bjarni's SNP data structure) to SNPData
		
		combine all chromsomes together
	"""
	import sys
	if report:
		sys.stderr.write("Converting RawSnpsData_ls to SNPData ...")
	from QC_250k import SNPData
	nt_dict_map = lambda x: nt2number[x]
	snpData = SNPData(row_id_ls = [], col_id_ls=[], data_matrix=[])
	for i in range(len(rawSnpsData_ls)):
		rawSnpsData = rawSnpsData_ls[i]
		chromosome = rawSnpsData.chromosome
		for j in range(len(rawSnpsData.positions)):
			if use_nt2number:
				data_row = map(nt_dict_map, rawSnpsData.snps[j])
			else:
				data_row = rawSnpsData.snps[j]
			snpData.data_matrix.append(data_row)
			snpData.row_id_ls.append((chromosome, rawSnpsData.positions[j]))
		if i==0:	#only need once
			for j in range(len(rawSnpsData.accessions)):
				if rawSnpsData.arrayIds:
					col_id = (rawSnpsData.accessions[j], rawSnpsData.arrayIds[j])
				else:
					col_id = rawSnpsData.accessions[j]
				snpData.col_id_ls.append(col_id)
	if report:
		sys.stderr.write("Done.\n")
	return snpData

def transposeSNPData(snpData, report=0):
	"""
	2008-05-18
		use num.int8 to keep memory small in num.transpose(num.array(snpData.data_matrix, num.int8))
	2008-05-18
		no more copy.deepcopy(snpData), data_matrix takes too long and too much memory
	05/12/08 fix a bug (return snpData)
	2008-05-11
	"""
	import sys
	if report:
		sys.stderr.write("Transposing SNPData ...")
	from pymodule import importNumericArray, SNPData
	num = importNumericArray()
	#copy except data_matrix
	import copy
	newSnpData = SNPData()
	"""
	for option_tuple in SNPData.option_default_dict:
		var_name = option_tuple[0]
		if var_name!='data_matrix':
			setattr(newSnpData, var_name, copy.deepcopy(getattr(snpData, var_name)))
	"""
	newSnpData.row_id_ls = copy.deepcopy(snpData.col_id_ls)
	newSnpData.col_id_ls = copy.deepcopy(snpData.row_id_ls)
	if isinstance(snpData.data_matrix, list):
		newSnpData.data_matrix = num.transpose(num.array(snpData.data_matrix, num.int8))
	else:	#assume it's array type already. Numeric/numarray has ArrayType, but numpy doesn't
		newSnpData.data_matrix = num.transpose(snpData.data_matrix)
	if report:
		sys.stderr.write("Done.\n")
	return newSnpData

def SNPData2RawSnpsData_ls(snpData, use_number2nt=1, need_transposeSNPData=0, report=0, mask_untouched_deleltion_as_NA=1):
	"""
	2008-05-19
		the transformation assumes snpData is in the orientation of SNP(row_id_ls) X Strain (col_id_ls). if not, toggle need_transposeSNPData=1.
	2008-05-18
		add mask_untouched_deleltion_as_NA. turned on by default because bjarni's RawSnpsData structure only recognizes NA, A, T, C, G
		if col_id in newSnpData.col_id_ls is tuple of size >1, the 2nd one  in the tuple is taken as array id.
	2008-05-12
		reverse of RawSnpsData_ls2SNPData
		
		adapts SNPData (Yu's SNP data structure) to RawSnpsData(bjarni's SNP data structure)
		
		split into different chromsomes
	"""
	import sys
	if report:
		sys.stderr.write("Converting SNPData to RawSnpsData_ls ...")
	from snpsdata import RawSnpsData
	
	if need_transposeSNPData:
		newSnpData = transposeSNPData(snpData, report=report)
	else:
		newSnpData = snpData
	
	accessions = []
	arrayIds = []
	accession_id = None	#default
	array_id = None	#default
	for col_id in newSnpData.col_id_ls:
		if isinstance(col_id, tuple):
			if len(col_id)>0:
				accession_id = col_id[0]
			if len(col_id)>1:
				array_id = col_id[1]
		else:
			accession_id = col_id
		accessions.append(accession_id)
		if array_id is not None:
			arrayIds.append(array_id)
	
	if mask_untouched_deleltion_as_NA:
		number2nt[-2] = 'NA'	#mask -2 (untouched) as 'NA'
		number2nt[-1] = 'NA'	#mask -1 (deletion) as 'NA'
	number2nt_dict_map = lambda x: number2nt[x]
	rawSnpsData_ls = []
	
	rawSnpsData = RawSnpsData(accessions=accessions, arrayIds=arrayIds)
	rawSnpsData.snps = []
	rawSnpsData.positions = []
	row_id0 = newSnpData.row_id_ls[0]
	old_chromosome = row_id0.split('_')[:1]
	rawSnpsData.chromosome = old_chromosome
	rawSnpsData_ls.append(rawSnpsData)
	rawSnpsData_ls_index = 0
	for i in range(len(newSnpData.row_id_ls)):
		row_id = newSnpData.row_id_ls[i]
		new_chromosome, position = row_id.split('_')[:2]
		position = int(position)
		if new_chromosome!=old_chromosome:
			rawSnpsData = RawSnpsData(accessions=accessions, arrayIds=arrayIds)
			rawSnpsData.snps = []
			rawSnpsData.positions = []
			rawSnpsData_ls.append(rawSnpsData)
			rawSnpsData_ls_index += 1
			rawSnpsData.chromosome = new_chromosome
			old_chromosome = new_chromosome
		rawSnpsData_ls[rawSnpsData_ls_index].positions.append(position)
		if use_number2nt:
			data_row = map(number2nt_dict_map, newSnpData.data_matrix[i])
		else:
			data_row = newSnpData.data_matrix[i]
		rawSnpsData_ls[rawSnpsData_ls_index].snps.append(data_row)
	if report:
		sys.stderr.write("Done.\n")
	return rawSnpsData_ls

def write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=None, \
					cols_to_be_tossed_out=None, nt_alphabet=0, transform_to_numpy=1,\
					discard_all_NA_rows=0, strain_acc2other_info=None, delimiter='\t', predefined_header_row=['strain', 'duplicate', 'latitude', 'longitude', 'nativename', 'stockparent', 'site', 'country']):
	"""
	strain_acc_list (and category_list) are initial 2 columns in the output.
	
	rows_to_be_tossed_out is a set or dict with row index in it. cols_to_be_tossed_out is similar structure.
	2008-05-19
		add output_fname in report
	2008-05-12
		copied from __init__.py
	2008-05-08
		include more options from dbSNP2data.py's write_data_matrix()
	2008-05-06
		add transform_to_numpy
	2008-04-02
		extracted from variation.src.FilterStrainSNPMatrix to be standalone.
	"""
	from sets import Set
	import sys, csv
	sys.stderr.write("Writing data_matrix to %s ..."%(output_fname))
	if rows_to_be_tossed_out==None:
		rows_to_be_tossed_out = Set()
	if cols_to_be_tossed_out==None:
		cols_to_be_tossed_out = Set()
	
	writer = csv.writer(open(output_fname, 'w'), delimiter=delimiter)
	
	if header:
		new_header = [header[0], header[1]]
		if strain_acc2other_info:
			no_of_fields = len(strain_acc2other_info[strain_acc2other_info.keys()[0]])
			for i in range(no_of_fields):
				new_header.append(predefined_header_row[2+i])
		for i in range(2, len(header)):
			if i-2 not in cols_to_be_tossed_out:
				new_header.append(header[i])
		writer.writerow(new_header)
	
	#figure out no_of_rows, no_of_cols
	if type(data_matrix)==list and transform_to_numpy:	#2008-02-06 transform the 2D list into array
		import numpy
		data_matrix = numpy.array(data_matrix)
		no_of_rows, no_of_cols = data_matrix.shape
	else:
		no_of_rows = len(data_matrix)
		if no_of_rows>0:
			no_of_cols = len(data_matrix[0])
		else:
			no_of_cols = 0
	
	no_of_all_NA_rows = 0
	for i in range(no_of_rows):
		if discard_all_NA_rows and sum(data_matrix[i]==0)==data_matrix.shape[1]:
			no_of_all_NA_rows += 1
			continue
		if i not in rows_to_be_tossed_out:
			new_row = [strain_acc_list[i], category_list[i]]
			if strain_acc2other_info:
				new_row += strain_acc2other_info[strain_acc_list[i]]
			for j in range(no_of_cols):
				if j not in cols_to_be_tossed_out:
					if nt_alphabet:
						new_row.append(number2nt[data_matrix[i][j]])
					else:
						new_row.append(data_matrix[i][j])
			writer.writerow(new_row)
	del writer
	sys.stderr.write("%s NA rows. Done.\n"%no_of_all_NA_rows)

def read_data(input_fname, input_alphabet=0, turn_into_integer=1, double_header=0, delimiter=None, matrix_data_type=int):
	"""
	2008-08-07
		turn_into_integer has to be toggled as well as p_char() detects character before nt2number is used.
	2008-08-03
		if p_char() detects character but dict_map() via nt2number fails to convert every entry in the row, turn data_row back to original un-converted.
	2008-07-11
		use p_char to judge whether there is character in the data, then use nt2number to convert them.
		input_alphabet is deprecated.
	2008-05-21
		add matrix_data_type, default=int. data_row = map(matrix_data_type, data_row)
	2008-05-18
		copied from FilterStrainSNPMatrix.py
		if delimiter not specified, call figureOutDelimiter()
	2008-05-12
		add delimiter
	2008-05-07
		add option double_header
	2007-03-06
		different from the one from SelectStrains.py is map(int, data_row)
	2007-05-14
		add input_alphabet
	2007-10-09
		add turn_into_integer
	"""
	import csv
	sys.stderr.write("Reading data from %s ..."%input_fname)
	if delimiter is None:
		delimiter = figureOutDelimiter(input_fname)
	reader = csv.reader(open(input_fname), delimiter=delimiter)
	header = reader.next()
	if double_header:
		header = [header, reader.next()]
	data_matrix = []
	strain_acc_list = []
	category_list = []
	import re
	p_char = re.compile(r'[a-zA-Z]')
	for row in reader:
		strain_acc_list.append(row[0])
		category_list.append(row[1])
		data_row = row[2:]
		no_of_snps = len(data_row)
		p_char_used = 0	#whether p_char is used to successfully dict_map the data_row
		if p_char.search(data_row[0]) and turn_into_integer:
			data_row = dict_map(nt2number, data_row)
			p_char_used = 1
			if no_of_snps!=len(data_row):
				sys.stderr.write('\n dict_map() via nt2number only maps %s out of %s entries from this row, %s, to integer. Back to original data.\n'%(len(data_row), no_of_snps, repr(row[:5])))
				data_row = row[2:]	#back to original data_row
				p_char_used = 0
		
		if turn_into_integer and not p_char_used:	#if p_char_used ==1, it's already integer.
			data_row = map(matrix_data_type, data_row)
		data_matrix.append(data_row)
	del reader
	sys.stderr.write("Done.\n")
	return header, strain_acc_list, category_list, data_matrix

class SNPData(object):
	"""
	2008-05-19
		either directly specify row_id_ls, col_id_ls, data_matrix
			(if strain_acc_list and category_list is given rather than row_id_ls and col_id_ls, row_id_ls and col_id_ls will be formed from them).
		or give input_fname and associated options to read data from file.
	2008-05-18 moved from variation.src.QC_250k
		add more arguments, input_fname, turn_into_array, ignore_2nd_column
	"""
	option_default_dict = {('row_id_ls', 0, ): None,\
							('col_id_ls', 0, ): None,\
							('header', 0, ): None,\
							('strain_acc_list', 0, ): None,\
							('category_list', 0, ): None,\
							('data_matrix', 0, ): None,\
							
							('input_fname', 0, ): None,\
							('input_alphabet',0, int): 0,\
							('turn_into_integer', 0, int): 1,\
							('double_header', 0, int):0, \
							
							('ignore_2nd_column', 0, int): [0, '', 0, 'Ignore category_list, the 2nd column.'],\
							
							('turn_into_array', 0 ,): [0, '', 0, 'Turn the data_matrix into array'],\
							
							('min_probability', 0, float): -1,\
							('call_method_id', 0, int): -1,\
							('col_id2id', 0, ):None,\
							('max_call_info_mismatch_rate', 0, float): 1,\
							('snps_table', 0, ):None}
	def __init__(self, **keywords):
		from __init__ import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)
		#read it from file
		if not self.data_matrix and isinstance(self.input_fname,str) and os.path.isfile(self.input_fname):
			self.header, self.strain_acc_list, self.category_list, self.data_matrix = read_data(self.input_fname, self.input_alphabet, self.turn_into_integer, self.double_header)
			if self.ignore_2nd_column:
				self.category_list = None
		self.processRowIDColID()
	
	def processRowIDColID(self):
		"""
		2008-06-02
			correct a bug here, opposite judgement of self.data_matrix
		"""
		if self.data_matrix and self.turn_into_array:
			self.data_matrix = num.array(self.data_matrix, num.int8)
		
		if self.row_id_ls is None and self.strain_acc_list is not None:
			self.row_id_ls = []
			for i in range(len(self.strain_acc_list)):
				if self.category_list is not None:
					row_id = (self.strain_acc_list[i], self.category_list[i])
				else:
					row_id = self.strain_acc_list[i]
				self.row_id_ls.append(row_id)
		
		if self.col_id_ls is None and self.header is not None:
			self.col_id_ls = []
			for i in range(2,len(self.header)):
				col_id = self.header[i]
				self.col_id_ls.append(col_id)
	
	def fromfile(self, input_fname, **keywords):
		"""
		2008-05-19
			initialize structure from a file in Strain X SNP format. 1st two columns are about strains.
			for keywords. check read_data()
		"""
		self.header, self.strain_acc_list, self.category_list, self.data_matrix = read_data(input_fname, **keywords)
		if self.ignore_2nd_column:
			self.category_list = None
		self.processRowIDColID()
		
	def tofile(self, output_fname, **keywords):
		"""
		2008-06-02
			either strain_acc_list or category_list could be None
			try to fill them both
		2008-05-18
			keywords is same as write_data_matrix()
		"""
		if self.header is None:
			self.header = ['', '']
			for col_id in self.col_id_ls:
				if isinstance(col_id, tuple):
					if len(col_id)>1:
						col_id = '%s_%s'%(col_id[0], col_id[1])	#chromosome, position
					else:
						col_id = col_id[0]
				self.header.append(col_id)
		if self.strain_acc_list is None or self.category_list is None:
			strain_acc_list = []
			category_list = []
			for row_id in self.row_id_ls:
				strain_acc = None
				category = None
				if isinstance(row_id, tuple):
					if len(row_id)>0:
						strain_acc = row_id[0]
					if len(row_id)>1:
						category = row_id[1]
				else:
					strain_acc = row_id
				if category == None:
					category = strain_acc	#make them same if category is not available
				strain_acc_list.append(strain_acc)
				category_list.append(category)
			if self.strain_acc_list is None:
				self.strain_acc_list = strain_acc_list
			if self.category_list is None:
				self.category_list = category_list
		write_data_matrix(self.data_matrix, output_fname, self.header, self.strain_acc_list, self.category_list, **keywords)
	
	def removeRowsByMismatchRate(cls, snpData, row_id2NA_mismatch_rate, max_mismatch_rate=1):
		"""
		2008-05-19
			more robust handling given max_mismatch_rate
		2008-05-19
		"""
		sys.stderr.write("Removing rows whose mismatch_rate >%s ..."%(max_mismatch_rate))
		if max_mismatch_rate>=0 and max_mismatch_rate<1:
			row_id_wanted_set = Set()	#extra computing time a bit, but to save memory
			for row_id in snpData.row_id_ls:
				if row_id in row_id2NA_mismatch_rate and row_id2NA_mismatch_rate[row_id][1]<=max_mismatch_rate and row_id2NA_mismatch_rate[row_id][1]>=0:
					row_id_wanted_set.add(row_id)

		elif max_mismatch_rate>=1:	#every row is wanted
			row_id_wanted_set = Set(snpData.row_id_ls)
		else:
			row_id_wanted_set = Set()
		no_of_rows = len(row_id_wanted_set)
		
		no_of_cols = len(snpData.col_id_ls)
		newSnpData = SNPData(col_id_ls=copy.deepcopy(snpData.col_id_ls), row_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		row_index = 0
		for i in range(len(snpData.row_id_ls)):
			row_id = snpData.row_id_ls[i]
			if row_id in row_id_wanted_set:
				newSnpData.row_id_ls.append(row_id)
				newSnpData.data_matrix[row_index] = snpData.data_matrix[i]
				row_index += 1
		newSnpData.no_of_rows_filtered_by_mismatch = len(snpData.row_id_ls)-no_of_rows
		sys.stderr.write("%s rows filtered by mismatch. Done.\n"%(newSnpData.no_of_rows_filtered_by_mismatch))
		return newSnpData
	
	removeRowsByMismatchRate = classmethod(removeRowsByMismatchRate)
	
	def removeRowsByNARate(cls, snpData, max_NA_rate=1):
		"""
		2008-05-19
			add is_cutoff_max to FilterStrainSNPMatrix.remove_rows_with_too_many_NAs
		2008-05-19
			if max_NA_rate out of [0,1) range, no calculation
		2008-05-19
		"""
		sys.stderr.write("Removing rows whose NA_rate >%s ..."%(max_NA_rate))
		if max_NA_rate<1 and max_NA_rate>=0:
			from FilterStrainSNPMatrix import FilterStrainSNPMatrix
			remove_rows_data = FilterStrainSNPMatrix.remove_rows_with_too_many_NAs(snpData.data_matrix, max_NA_rate, is_cutoff_max=1)
		
			rows_with_too_many_NAs_set = remove_rows_data.rows_with_too_many_NAs_set
		elif max_NA_rate>=1:	#no row has too many NA sets
			rows_with_too_many_NAs_set = Set()
		else:
			rows_with_too_many_NAs_set = Set(range(len(snpData.row_id_ls)))
		no_of_rows = len(snpData.row_id_ls)-len(rows_with_too_many_NAs_set)
		no_of_cols = len(snpData.col_id_ls)
		newSnpData = SNPData(col_id_ls=copy.deepcopy(snpData.col_id_ls), row_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		row_index = 0
		for i in range(len(snpData.row_id_ls)):
			row_id = snpData.row_id_ls[i]
			if i not in rows_with_too_many_NAs_set:
				newSnpData.row_id_ls.append(row_id)
				newSnpData.data_matrix[row_index] = snpData.data_matrix[i]
				row_index += 1
		newSnpData.no_of_rows_filtered_by_na = len(rows_with_too_many_NAs_set)
		sys.stderr.write("%s rows filtered by NA. Done.\n"%(newSnpData.no_of_rows_filtered_by_na))
		return newSnpData
	removeRowsByNARate = classmethod(removeRowsByNARate)
	
	def removeColsByMismatchRate(cls, snpData, col_id2NA_mismatch_rate, max_mismatch_rate=1):
		"""
		2008-05-19
			more robust handling given max_mismatch_rate
			fix a bug. mismatch rate could be -1 which is no-comparison
		2008-05-19
		"""
		sys.stderr.write("Removing cols whose mismatch_rate >%s ..."%(max_mismatch_rate))
		if max_mismatch_rate>=0 and max_mismatch_rate<1:
			col_id_wanted_set = Set()	#extra computing time a bit, but to save memory
			for col_id in snpData.col_id_ls:
				if col_id in col_id2NA_mismatch_rate and col_id2NA_mismatch_rate[col_id][1]<=max_mismatch_rate and col_id2NA_mismatch_rate[col_id][1]>=0:
					col_id_wanted_set.add(col_id)
		elif max_mismatch_rate>=1:
			col_id_wanted_set = Set(snpData.col_id_ls)
		else:
			col_id_wanted_set = Set()
		
		no_of_cols = len(col_id_wanted_set)
		no_of_rows = len(snpData.row_id_ls)
		newSnpData = SNPData(col_id_ls=[], row_id_ls=snpData.row_id_ls)
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		col_index = 0
		for i in range(len(snpData.col_id_ls)):
			col_id = snpData.col_id_ls[i]
			if col_id in col_id_wanted_set:
				newSnpData.col_id_ls.append(col_id)
				newSnpData.data_matrix[:,col_index] = snpData.data_matrix[:,i]
				col_index += 1
		newSnpData.no_of_cols_filtered_by_mismatch = len(snpData.col_id_ls)-no_of_cols
		sys.stderr.write("%s columns filtered by mismatch. Done.\n"%(newSnpData.no_of_cols_filtered_by_mismatch))
		return newSnpData
	
	removeColsByMismatchRate = classmethod(removeColsByMismatchRate)
	
	def removeColsByNARate(cls, snpData, max_NA_rate=1):
		"""
		2008-05-19
			add is_cutoff_max to FilterStrainSNPMatrix.remove_cols_with_too_many_NAs()
		2008-05-19
			if max_NA_rate out of [0,1) range, no calculation
		2008-05-19
		"""
		sys.stderr.write("Removing cols whose NA_rate >%s ..."%(max_NA_rate))
		if max_NA_rate<1 and max_NA_rate>=0:
			from FilterStrainSNPMatrix import FilterStrainSNPMatrix
			remove_cols_data = FilterStrainSNPMatrix.remove_cols_with_too_many_NAs(snpData.data_matrix, max_NA_rate, is_cutoff_max=1)			
			cols_with_too_many_NAs_set = remove_cols_data.cols_with_too_many_NAs_set
		elif max_NA_rate>=1:
			cols_with_too_many_NAs_set = Set()
		else:
			cols_with_too_many_NAs_set = Set(range(len(snpData.col_id_ls)))
		
		no_of_cols = len(snpData.col_id_ls)-len(cols_with_too_many_NAs_set)
		no_of_rows = len(snpData.row_id_ls)
		newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		col_index = 0
		for i in range(len(snpData.col_id_ls)):
			col_id = snpData.col_id_ls[i]
			if i not in cols_with_too_many_NAs_set:
				newSnpData.col_id_ls.append(col_id)
				newSnpData.data_matrix[:,col_index] = snpData.data_matrix[:,i]
				col_index += 1
		newSnpData.no_of_cols_filtered_by_na = len(cols_with_too_many_NAs_set)
		sys.stderr.write("%s columns filtered by NA. Done.\n"%(newSnpData.no_of_cols_filtered_by_na))
		return newSnpData
	removeColsByNARate = classmethod(removeColsByNARate)

	def removeMonomorphicCols(cls, snpData):
		"""
		2008-05-19
		"""
		sys.stderr.write("Removing monomorphic cols ...")
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		#NA_set = Set([0,-1,-2])
		col_index_wanted_ls = []
		for j in range(no_of_cols):
			non_negative_number_set = Set()
			for i in range(no_of_rows):
				number = snpData.data_matrix[i][j]
				if number >0:
					non_negative_number_set.add(number)
			if len(non_negative_number_set)>1:	#not monomorphic
				col_index_wanted_ls.append(j)
		
		newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, len(col_index_wanted_ls)], num.int8)
		col_index = 0
		for i in col_index_wanted_ls:
			col_id = snpData.col_id_ls[i]
			newSnpData.col_id_ls.append(col_id)
			newSnpData.data_matrix[:,col_index] = snpData.data_matrix[:,i]
			col_index += 1
		newSnpData.no_of_monomorphic_cols = no_of_cols-len(newSnpData.col_id_ls)
		sys.stderr.write("%s monomorphic columns. Done.\n"%(newSnpData.no_of_monomorphic_cols))
		return newSnpData
	removeMonomorphicCols = classmethod(removeMonomorphicCols)
	
	def convertHetero2NA(cls, snpData):
		"""
		2008-06-02
			Convert all heterozygous calls and untouched in the file into NA.
			deletion is not converted
		"""
		sys.stderr.write("Converting Hetero calls to NA ...")
		no_of_hets = 0
		newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=snpData.col_id_ls)
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				if snpData.data_matrix[i][j]<=4 and snpData.data_matrix[i][j]>=1:
					newSnpData.data_matrix[i][j] = snpData.data_matrix[i][j]
				elif snpData.data_matrix[i][j]==-1:	#but not -2 (untouched), -1=deletion
					newSnpData.data_matrix[i][j] = snpData.data_matrix[i][j]
				else:
					no_of_hets += 1
		sys.stderr.write("%s heterozygous calls. Done.\n"%no_of_hets)
		return newSnpData
	convertHetero2NA = classmethod(convertHetero2NA)
	
	def removeSNPsWithMoreThan2Alleles(cls, snpData):
		"""
		2008-08-05
			NA and -2 (not touched) are not considered as an allele
		"""
		sys.stderr.write("Removing SNPs with more than 2 alleles ...")
		no_of_rows, no_of_cols = snpData.data_matrix.shape
		col_index_wanted_ls = []
		for j in range(no_of_cols):
			allele_set = Set(snpData.data_matrix[:,j])
			if 0 in allele_set:	#remove NA if it's there
				allele_set.remove(0)
			if -2 in allele_set:	#remove -2 as well
				allele_set.remove(-2)
			if len(allele_set)==2:	#polymorphic
				col_index_wanted_ls.append(j)
		
		newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, len(col_index_wanted_ls)], num.int8)
		col_index = 0
		for i in col_index_wanted_ls:
			col_id = snpData.col_id_ls[i]
			newSnpData.col_id_ls.append(col_id)
			newSnpData.data_matrix[:,col_index] = snpData.data_matrix[:,i]
			col_index += 1
		newSnpData.no_of_cols_removed = no_of_cols-len(newSnpData.col_id_ls)
		sys.stderr.write("%s columns removed. Done.\n"%(newSnpData.no_of_cols_removed))
		return newSnpData
	removeSNPsWithMoreThan2Alleles = classmethod(removeSNPsWithMoreThan2Alleles)

from db import TableClass
class GenomeWideResults(TableClass):
	genome_wide_result_ls = None
	genome_wide_result_obj_id2index = None
	max_value = 0	#the current top value for all genome wide results
	gap = 1.0	#gap between two genome wide results
	def get_genome_wide_result_by_obj_id(self, obj_id):
		return self.genome_wide_result_ls[self.genome_wide_result_obj_id2index[obj_id]]
	
	def add_genome_wide_result(self, genome_wide_result):
		genome_wide_result_index = len(self.genome_wide_result_ls)
		if genome_wide_result_index==0:	#the first result, no gap necessary
			genome_wide_result.base_value = self.max_value
		else:
			genome_wide_result.base_value = self.max_value + self.gap
		new_max_value = genome_wide_result.base_value + genome_wide_result.max_value - genome_wide_result.min_value
		if self.max_value is None or new_max_value > self.max_value:
			self.max_value = new_max_value
		
		self.genome_wide_result_ls.append(genome_wide_result)
		self.genome_wide_result_obj_id2index[id(genome_wide_result)] = genome_wide_result_index

class GenomeWideResult(TableClass):
	data_obj_ls = None
	data_obj_id2index = None
	name = None
	results_method_id = None
	results_method = None
	min_value = None
	max_value = None
	base_value = 0
	
	def get_data_obj_by_obj_id(self, obj_id):
		return self.data_obj_ls[self.data_obj_id2index[obj_id]]
	
	def get_data_obj_by_obj_index(self, obj_index):
		return self.data_obj_ls[obj_index]
	
	def add_one_data_obj(self, data_obj):
		data_obj_index = len(self.data_obj_ls)
		self.data_obj_ls.append(data_obj)
		
		self.data_obj_id2index[id(data_obj)] = data_obj_index
		if self.min_value is None or data_obj.value<self.min_value:
			self.min_value = data_obj.value
		if self.max_value is None or data_obj.value>self.max_value:
			self.max_value = data_obj.value

class DataObject(TableClass):
	chromosome = None
	position = None
	stop_position = None
	name = None
	value = None
	genome_wide_result_id = None

import re
pa_has_characters = re.compile(r'[a-zA-Z_]')
import math

def getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=False, pdata=None):
	"""
	2008-08-14
		add min_MAF into pdata
	2008-08-03
		add pdata (chromosome, start, stop) to restrain data
	2008-07-17
		moved from GenomeBrowser.py
	2008-05-31
		automatically detect if header exists on the first line.
	2008-05-28
		handle both 3/4-column input file
	"""
	sys.stderr.write("Getting genome wide result from %s ... "%input_fname)
	gwr = GenomeWideResult(name=os.path.basename(input_fname))
	gwr.data_obj_ls = []	#list and dictionary are crazy references.
	gwr.data_obj_id2index = {}
	genome_wide_result_id = id(gwr)
	import csv
	reader = csv.reader(open(input_fname), delimiter='\t')
	no_of_lines = 0
	for row in reader:
		#check if 1st line is header or not
		if no_of_lines ==0 and pa_has_characters.search(row[1]):
			continue
		chr = int(row[0])
		start_pos = int(row[1])
		if len(row)==3:
			stop_pos = None
			score = float(row[2])
			column_4th = None
		elif len(row)>=4:
			score = float(row[2])
			column_4th=row[3]
			stop_pos = None
			#stop_pos = int(row[2])
			#score = float(row[3])
		"""
		else:
			sys.stderr.write("only 3 or 4 columns are allowed in input file.\n")
			return gwr
		"""
		if pdata:	#2008-08-03
			pdata.chromosome = getattr(pdata, 'chromosome', None)
			pdata.start = getattr(pdata, 'start', None)
			pdata.stop = getattr(pdata, 'stop', None)
			pdata.min_MAF = getattr(pdata, 'min_MAF', None)
			if pdata.chromosome!=None and chr!=pdata.chromosome:
				continue
			if pdata.start!=None and start_pos<pdata.start:
				continue
			if pdata.stop!=None and start_pos>pdata.stop:
				continue
			if pdata.min_MAF!=None and column_4th!=None and float(column_4th)<pdata.min_MAF:	#MAF too small
				continue
		if do_log10_transformation:
			score = -math.log10(score)
		if min_value_cutoff is None or score>=min_value_cutoff:
			if stop_pos is not None:
				data_obj = DataObject(chromosome=chr, position=start_pos, stop_position=stop_pos, value =score)
			else:
				data_obj = DataObject(chromosome=chr, position=start_pos, value =score)
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj)
		
		no_of_lines += 1
		
	del reader
	sys.stderr.write("Done.\n")
	return gwr
