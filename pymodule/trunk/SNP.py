#!/usr/bin/env python
"""
2008-05-12
	store SNP data structure related stuff
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray, figureOutDelimiter, PassingData
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

#2008-09-22 A:T, C:G complement group in number
number2complement = {-1:-1, 0:0, 1:4, 4:1, 2:3, 3:2}
#2008-01-06
nt2complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
#2008-05-12	a common NA set
from sets import Set
NA_set = Set([0, 'NA', 'N', -2, '|'])

# 2009-6-5 add lower-case letter to the dictionary
nt2number = {'|': -2,	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	'n': 0,
	'na': 0,
	'X': 0,
	'x': 0,
	None: 0,
	'A': 1,
	'a': 1,
	'C': 2,
	'c': 2,
	'G': 3,
	'g': 3,
	'T': 4,
	't': 4,
	'AC': 5,
	'CA': 5,
	'M': 5,
	'm': 5,
	'AG': 6,
	'GA': 6,
	'R': 6,
	'r': 6,
	'AT': 7,
	'TA': 7,
	'W': 7,
	'w': 7,
	'CG': 8,
	'GC': 8,
	'S': 8,
	's': 8,
	'CT': 9,
	'TC': 9,
	'Y': 9,
	'y': 9,
	'GT': 10,
	'TG': 10,
	'K': 10,
	'k': 10
	}

# 2009-5-19 map nucleotide to number while ignoring heterozygous calls
nt2number_without_het = {'|': -2,	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	None: 0,
	'A': 1,
	'C': 2,
	'G': 3,
	'T': 4,
	'AC':0,
	'CA':0,
	'M':0,
	'AG':0,
	'GA':0,
	'R':0,
	'AT':0,
	'TA':0,
	'W':0,
	'CG':0,
	'GC':0,
	'S':0,
	'CT':0,
	'TC':0,
	'Y':0,
	'GT':0,
	'TG':0,
	'K':0
	}

number2nt = {-2: '|',	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	-1: '-',
	0: 'NA',
	num.nan: 'NA',	#2008-12-03	data_matrix might contain num.nan as NA
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

#2009-6-5 dictionary used in output FASTA format which doesn't recognize two-letter heterozygous representation
number2single_char_nt = {-2: '|',	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	-1: '-',
	0: 'N',
	num.nan: 'N',	#2008-12-03	data_matrix might contain num.nan as NA
	1:'A',
	2:'C',
	3:'G',
	4:'T',
	5:'M',
	6:'R',
	7:'W',
	8:'S',
	9:'Y',
	10:'K'
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
	2009-5-18
		only the integer-type keys in number2nt, which contains num.nan in it now. 
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
		if type(number)==int:	#2009-5-18 number2nt contains num.nan
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
	2008-12-03
		newSnpData.processRowIDColID() after row_id_ls and col_id_ls are given
		also copy strain_acc_list and category_list over to newSnpData
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
	newSnpData.strain_acc_list = copy.deepcopy(snpData.strain_acc_list)
	newSnpData.category_list = copy.deepcopy(snpData.category_list)
	
	newSnpData.processRowIDColID()	#2008-12-02	processRowIDColID() after row_id_ls and col_id_ls are given
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

def read_data(input_fname, input_alphabet=0, turn_into_integer=1, double_header=0, delimiter=None, matrix_data_type=int, ignore_het=0):
	"""
	2009-5-20
		add argument ignore_het, which upon toggled, instructs the function to use nt2number_without_het to map nucleotides to number.
	2009-3-21
		add '-' into vocabulary of p_char and append '$' requiring the entry ends with any characters in that vocabulary
	2009-2-18
		no 'e' or 'E' among the letters to be checked in the 1st column in order to decide whether to cast it into matrix_data_type
		because 'e' or 'E' could be used in scientific number.
		a better version is only to check whether nucleotide letters are in it, because nt2number is used if the check is positive.
	2008-08-29
		put the handling of each row into a "try ... except ..."
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
	if ignore_het:
		nt2number_mapper = nt2number_without_het
	else:
		nt2number_mapper = nt2number
	def ignore_het_func(x):
		if x>4:
			return 0
		else:
			return x
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
	p_char = re.compile(r'[a-df-zA-DF-Z\-]$')	#no 'e' or 'E', used in scientific number, add '-' and append '$'
	i = 0
	for row in reader:
		i += 1
		try:
			strain_acc_list.append(row[0])
			category_list.append(row[1])
			data_row = row[2:]
			no_of_snps = len(data_row)
			p_char_used = 0	#whether p_char is used to successfully dict_map the data_row
			if p_char.search(data_row[0]) and turn_into_integer:
				data_row = dict_map(nt2number_mapper, data_row)
				p_char_used = 1
				if no_of_snps!=len(data_row):
					sys.stderr.write('\n dict_map() via nt2number_mapper only maps %s out of %s entries from this row, %s, to integer. Back to original data.\n'%\
									(len(data_row), no_of_snps, repr(row[:5])))
					data_row = row[2:]	#back to original data_row
					p_char_used = 0
			
			if turn_into_integer and not p_char_used:	#if p_char_used ==1, it's already integer.
				data_row = map(matrix_data_type, data_row)
				if ignore_het:	#2009-5-20 for data that is already in number format, use this function to remove hets
					data_row = map(ignore_het_func, data_row)
			data_matrix.append(data_row)
		except:
			sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			import traceback
			traceback.print_exc()
			sys.stderr.write("Row no: %s. %s.\n"%(i, repr(row)))
			raise
	del reader
	sys.stderr.write("Done.\n")
	return header, strain_acc_list, category_list, data_matrix

class SNPData(object):
	"""
	2009-3-5
		example usages:
		
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1, ignore_2nd_column=1)
			
			snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1)
			
			header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, delimiter=delimiter)
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
	2009-5-20
		add argument ignore_het which is passed to read_data() when data_matrix is not given to ignore heterozygous calls.
		It'll only be functional when data_matrix is not given upon SNPData initialization.
	2008-05-19
		either directly specify row_id_ls, col_id_ls, data_matrix
			(if strain_acc_list and category_list is given rather than row_id_ls and col_id_ls, row_id_ls and col_id_ls will be formed from them).
		or give input_fname and associated options to read data from file.
	2008-05-18 moved from variation.src.QC_250k
		add more arguments, input_fname, turn_into_array, ignore_2nd_column
	"""
	report = 0
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
							('snps_table', 0, ): None,\
							('matrix_data_type', 0, ): int,\
							('ignore_het', 0, int): [0, '', 0, 'Ignore the heterozygous genotypes'],\
							('debug', 0, int): [0, '', 0, 'turn on the debug flag']}
	def __init__(self, **keywords):
		from __init__ import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, \
														class_to_have_attr=self, howto_deal_with_required_none=2)
		#read it from file
		if self.isDataMatrixEmpty(self.data_matrix) and isinstance(self.input_fname,str) and os.path.isfile(self.input_fname):
			self.header, self.strain_acc_list, self.category_list, self.data_matrix = read_data(self.input_fname, self.input_alphabet, \
																							self.turn_into_integer, self.double_header, \
																							matrix_data_type=self.matrix_data_type,\
																							ignore_het=self.ignore_het)
			if self.ignore_2nd_column:
				self.category_list = None
		self.processRowIDColID()
	
	@classmethod
	def isDataMatrixEmpty(cls, data_matrix):
		"""
		2008-08-21
			make it a classmethod
		2008-08-19
			common function to judge whether data_matrix is empty
		"""
		if data_matrix=='':
			return True
		elif data_matrix is None:
			return True
		elif hasattr(data_matrix, '__len__') and len(data_matrix)==0:
			return True
		else:
			return False
	
	def processRowIDColID(self):
		"""
		2008-12-03
			this function is called in __init__(). deal with the case that this class is initialized without any argument,
				which means self.row_id_ls and self.col_id_ls are None.
		2008-09-07
			if turn the data_matrix into array, do not force num.int8 type.
		2008-09-05
			generate id2index for both row and column
		2008-06-02
			correct a bug here, opposite judgement of self.data_matrix
		"""
		if not self.isDataMatrixEmpty(self.data_matrix) and self.turn_into_array:
			self.data_matrix = num.array(self.data_matrix)
		
		if self.row_id_ls is None and self.strain_acc_list is not None:
			self.row_id_ls = []
			for i in range(len(self.strain_acc_list)):
				if self.category_list is not None:
					row_id = (self.strain_acc_list[i], self.category_list[i])
				else:
					row_id = self.strain_acc_list[i]
				self.row_id_ls.append(row_id)
		
		self.row_id2row_index = {}
		if self.row_id_ls:	#2008-12-03
			for i in range(len(self.row_id_ls)):
				row_id = self.row_id_ls[i]
				self.row_id2row_index[row_id] = i
		
		if self.col_id_ls is None and self.header is not None:
			self.col_id_ls = []
			for i in range(2,len(self.header)):
				col_id = self.header[i]
				self.col_id_ls.append(col_id)
		
		self.col_id2col_index = {}
		if self.col_id_ls:	#2008-12-03
			for i in range(len(self.col_id_ls)):
				col_id = self.col_id_ls[i]
				self.col_id2col_index[col_id] = i
	
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
	
	@classmethod
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
	
	@classmethod
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
	
	@classmethod
	def keepRowsByRowID(cls, snpData, row_id_ls):
		"""
		2009-05-19
			keep certain rows in snpData given row_id_ls
		"""
		sys.stderr.write("Keeping rows given row_id_ls ...")
		no_of_rows = len(row_id_ls)
		row_id_wanted_set = set(row_id_ls)
		
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
		sys.stderr.write("%s rows discarded. Done.\n"%(newSnpData.no_of_rows_filtered_by_mismatch))
		return newSnpData
	
	@classmethod
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
	
	
	@classmethod
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
	
	@classmethod
	def removeColsByMAF(cls, snpData, min_MAF=1):
		"""
		2008-05-29
			remove SNPs whose MAF is lower than min_MAF
		"""
		sys.stderr.write("Removing cols whose MAF <%s ..."%(min_MAF))
		
		no_of_rows = len(snpData.data_matrix)
		no_of_cols = len(snpData.data_matrix[0])
		allele2index_ls = []
		allele2count_ls = []
		allele_index2allele_ls = []
		col_id_to_be_kept_ls = []
		no_of_SNPs_with_more_than_2_alleles = 0
		for j in range(no_of_cols):
			col_id = snpData.col_id_ls[j]
			allele2index_ls.append({})
			allele2count_ls.append({})
			allele_index2allele_ls.append({})
			allele_index2allele = allele_index2allele_ls[j]
			for i in range(no_of_rows):
				allele = snpData.data_matrix[i][j]
				if allele==0 or allele==-2:
					allele_index = num.nan	#numpy.nan is better than -2
				elif allele not in allele2index_ls[j]:
					allele_index = len(allele2index_ls[j])
					allele2index_ls[j][allele] = allele_index
					allele2count_ls[j][allele] = 1
					allele_index2allele[allele_index] = allele
				else:
					allele_index = allele2index_ls[j][allele]
					allele2count_ls[j][allele] += 1
				#if cls.report and allele_index>1:
				#	sys.stderr.write("%s (more than 2) alleles at SNP %s (id=%s).\n"%((allele_index+1), j, snpData.col_id_ls[j]))
			if len(allele2count_ls[j])>2:
				no_of_SNPs_with_more_than_2_alleles += 1
				if cls.report:
					sys.stderr.write("more than 2 alleles at SNP %s (id=%s).\n"%(j, snpData.col_id_ls[j]))
			MAF = min(allele2count_ls[j].values())/float(sum(allele2count_ls[j].values()))
			"""
			print MAF
			print j
			print snpData.col_id_ls[j]
			print allele2count_ls[j]
			"""
			if MAF>=min_MAF:
				col_id_to_be_kept_ls.append(col_id)
		
		newSnpData = cls.keepColsByColID(snpData, col_id_to_be_kept_ls)
		newSnpData.no_of_cols_removed = no_of_cols - len(col_id_to_be_kept_ls)
		sys.stderr.write("%s columns filtered by min_MAF and %s SNPs with >2 alleles. Done.\n"%(newSnpData.no_of_cols_removed, no_of_SNPs_with_more_than_2_alleles))
		return newSnpData
	
	@classmethod
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
	
	@classmethod
	def keepColsByColID(cls, snpData, col_id_ls):
		"""
		2009-05-29
			keep certain columns in snpData given col_id_ls
		"""
		sys.stderr.write("Keeping columns given col_id_ls ...")
		no_of_cols = len(col_id_ls)
		col_id_wanted_set = set(col_id_ls)
		
		no_of_rows = len(snpData.row_id_ls)
		newSnpData = SNPData(row_id_ls=copy.deepcopy(snpData.row_id_ls), col_id_ls=[])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		col_index = 0
		for j in range(len(snpData.col_id_ls)):
			col_id = snpData.col_id_ls[j]
			if col_id in col_id_wanted_set:
				newSnpData.col_id_ls.append(col_id)
				newSnpData.data_matrix[:,col_index] = snpData.data_matrix[:, j]
				col_index += 1
		newSnpData.no_of_cols_removed = len(snpData.col_id_ls)-no_of_cols
		sys.stderr.write("%s columns discarded. Done.\n"%(newSnpData.no_of_cols_removed))
		return newSnpData
	
	
	@classmethod
	def convertHetero2NA(cls, snpData):
		"""
		2009-5-29
			bug: previously, no_of_hets includes NA calls. now it doesn't
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
				elif snpData.data_matrix[i][j]>4:	#-2 and 0 all converted to 0 by default
					no_of_hets += 1
		sys.stderr.write("%s heterozygous calls. Done.\n"%no_of_hets)
		return newSnpData
	
	@classmethod
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
	
	def fill_in_snp_allele2index(self, diploid_allele, allele2index):
		"""
		2008-09-05
			used in calLD
		"""
		if diploid_allele>4:
			nt = number2nt[diploid_allele]
			allele1 = nt2number[nt[0]]
			allele2 = nt2number[nt[1]]
		else:
			allele1 = allele2 = diploid_allele
		if allele1 not in allele2index:
			allele2index[allele1] = len(allele2index)
		if allele2 not in allele2index:
			allele2index[allele2] = len(allele2index)
		return allele1, allele2
	
	def calLD(self, col1_id, col2_id):
		"""
		2008-09-05
			adapted from variation.src.misc's LD.calculate_LD class
			only deal with 2-allele loci
			skip if either is NA, or if both are heterozygous (not phased)
		"""
		snp1_index = self.col_id2col_index[col1_id]
		snp2_index = self.col_id2col_index[col2_id]
		counter_matrix = num.zeros([2,2])	#only 2 alleles
		snp1_allele2index = {}
		snp2_allele2index = {}
		for k in range(len(self.row_id_ls)):
			snp1_allele = self.data_matrix[k][snp1_index]
			snp2_allele = self.data_matrix[k][snp2_index]
			if snp1_allele!=0 and snp2_allele!=0 and not (snp1_allele>4 and snp2_allele>4):	#doesn't allow both loci are heterozygous
				snp1_allele1, snp1_allele2 = self.fill_in_snp_allele2index(snp1_allele, snp1_allele2index)
				snp2_allele1, snp2_allele2 = self.fill_in_snp_allele2index(snp2_allele, snp2_allele2index)
				counter_matrix[snp1_allele2index[snp1_allele1],snp2_allele2index[snp2_allele1]] += 1
				counter_matrix[snp1_allele2index[snp1_allele2],snp2_allele2index[snp2_allele2]] += 1
		PA = sum(counter_matrix[0,:])
		Pa = sum(counter_matrix[1,:])
		PB = sum(counter_matrix[:,0])
		Pb = sum(counter_matrix[:,1])
		total_num = float(PA+Pa)
		try:
			PA = PA/total_num
			Pa = Pa/total_num
			PB = PB/total_num
			Pb = Pb/total_num
			PAB = counter_matrix[0,0]/total_num
			D = PAB-PA*PB
			PAPB = PA*PB
			PAPb = PA*Pb
			PaPB = Pa*PB
			PaPb = Pa*Pb
			Dmin = max(-PAPB, -PaPb)
			Dmax = min(PAPb, PaPB)
			if D<0:
				D_prime = D/Dmin
			else:
				D_prime = D/Dmax
			r2 = D*D/(PA*Pa*PB*Pb)
		except:	#2008-01-23 exceptions.ZeroDivisionError, Dmin or Dmax could be 0 if one of(-PAPB, -PaPb)  is >0 or <0
			sys.stderr.write('Unknown except, ignore: %s\n'%repr(sys.exc_info()[0]))
			return None
		allele_freq = (min(PA, Pa),min(PB, Pb))
		return_data = PassingData()
		return_data.D = D
		return_data.D_prime = D_prime
		return_data.r2 = r2
		return_data.allele_freq = allele_freq
		return_data.snp_pair_ls = (col1_id, col2_id)
		return_data.no_of_pairs = total_num
		return return_data
	
	def calRowPairwiseDist(self):
		"""
		2009-4-18
			calculate distance between all rows except itself.
			only calculate half non-redundant pairs.
		"""
		sys.stderr.write("Calculating row-wise pairwise distance ...")
		row_id2pairwise_dist = {}
		counter = 0
		
		for i in range(len(self.row_id_ls)):
			row_id1 = self.row_id_ls[i]
			pairwise_dist = []
			for j in range(i+1, len(self.row_id_ls)):
				row_id2 = self.row_id_ls[j]
				no_of_mismatches = 0
				no_of_non_NA_pairs = 0
				for col_index in range(len(self.col_id_ls)):
					if self.data_matrix[i][col_index] not in NA_set and self.data_matrix[j][col_index] not in NA_set:
						no_of_non_NA_pairs += 1
						if self.data_matrix[i][col_index] != self.data_matrix[j][col_index]:
							no_of_mismatches += 1
				if no_of_non_NA_pairs>0:
					mismatch_rate = no_of_mismatches/float(no_of_non_NA_pairs)
				else:
					mismatch_rate = -1
					if self.debug:
						sys.stderr.write("\t no valid(non-NA) pairs between %s and %s.\n"%(row_id1, row_id2))
				pairwise_dist.append([mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs])
			pairwise_dist.sort()
			row_id2pairwise_dist[row_id1] = pairwise_dist
		sys.stderr.write("Done.\n")
		return row_id2pairwise_dist
	
	def convertSNPAllele2Index(self, report=0):
		"""
		2008-12-03
			set in_major_minor_order to True
		2008-12-02
			code body moved to convert2Binary()
			call convert2Binary with in_major_minor_order=False.
		2008-11-20
			newSnpData also copies self.strain_acc_list and self.category_list over.
		2008-09-07
			Convert SNP matrix into index (0,1,2...) is assigned as first-encounter, first-assign. if only two alleles, it's binary.
			heterozygote is regarded as a different allele.
		"""
		return self.convert2Binary(report=report, in_major_minor_order=True)
	
	def convert2Binary(self, report=0, in_major_minor_order=True):
		"""
		2008-12-02
			based on the old convertSNPAllele2Index()
		"""
		sys.stderr.write("Converting SNP matrix to Binary ...")
		no_of_hets = 0
		newSnpData = SNPData(row_id_ls=self.row_id_ls, col_id_ls=self.col_id_ls)
		newSnpData.strain_acc_list = self.strain_acc_list
		newSnpData.category_list = self.category_list
		no_of_rows = len(self.data_matrix)
		no_of_cols = len(self.data_matrix[0])
		newSnpData.data_matrix = num.zeros([no_of_rows, no_of_cols], num.int8)
		allele2index_ls = []
		allele2count_ls = []
		allele_index2allele_ls = []
		for j in range(no_of_cols):
			allele2index_ls.append({})
			allele2count_ls.append({})
			allele_index2allele_ls.append({})
			allele_index2allele = allele_index2allele_ls[j]
			for i in range(no_of_rows):
				allele = self.data_matrix[i][j]
				if allele==0 or allele==-2:
					allele_index = num.nan	#numpy.nan is better than -2
				elif allele not in allele2index_ls[j]:
					allele_index = len(allele2index_ls[j])
					allele2index_ls[j][allele] = allele_index
					allele2count_ls[j][allele] = 1
					allele_index2allele[allele_index] = allele
				else:
					allele_index = allele2index_ls[j][allele]
					allele2count_ls[j][allele] += 1
				newSnpData.data_matrix[i][j] = allele_index
				if report and allele_index>1:
					sys.stderr.write("%s (more than 2) alleles at SNP %s (id=%s).\n"%((allele_index+1), j, self.col_id_ls[j]))
			#2008-12-02	check if binary-allele codes are in major, minor order
			if in_major_minor_order and len(allele2index_ls[j])==2:	#only two-allele SNPs
				# allele index 0 would be major allele
				# allele index 1 would be minor allele
				allele1 = allele_index2allele[0]
				allele2 = allele_index2allele[1]
				if allele2count_ls[j][allele1]<allele2count_ls[j][allele2]:	#minor allele got assigned the smaller number, reverse it
					newSnpData.data_matrix[:,j] = num.abs(newSnpData.data_matrix[:,j]-1)	#reverse the index. won't affect nan.
					allele2index_ls[j][allele1] = 1
					allele2index_ls[j][allele2] = 0
					allele_index2allele[0] = allele2
					allele_index2allele[1] = allele1
		sys.stderr.write("Done.\n")
		return newSnpData, allele_index2allele_ls
		
from db import TableClass
class GenomeWideResults(TableClass):
	genome_wide_result_ls = None
	genome_wide_result_obj_id2index = None
	max_value = 0	#the current top value for all genome wide results
	#gap = 1.0	#gap between two genome wide results
	def get_genome_wide_result_by_obj_id(self, obj_id):
		return self.genome_wide_result_ls[self.genome_wide_result_obj_id2index[obj_id]]
	
	def add_genome_wide_result(self, genome_wide_result, gap=None):
		"""
		2008-10-12
			specify gap between two genome wide results through option gap or self-guessing, 1/3 of previous (gwr.max_value-gwr.min_value)
		"""
		genome_wide_result_index = len(self.genome_wide_result_ls)
		if genome_wide_result_index==0:	#the first result, no gap necessary
			genome_wide_result.base_value = self.max_value
		else:
			if gap is None:
				prev_gwr = self.genome_wide_result_ls[genome_wide_result_index-1]
				gap = (prev_gwr.max_value - prev_gwr.min_value)/3.	#self-guessed gap
			genome_wide_result.base_value = self.max_value + gap
		new_max_value = genome_wide_result.base_value + genome_wide_result.max_value - genome_wide_result.min_value
		if self.max_value is None or new_max_value > self.max_value:
			self.max_value = new_max_value
		
		self.genome_wide_result_ls.append(genome_wide_result)
		self.genome_wide_result_obj_id2index[id(genome_wide_result)] = genome_wide_result_index

class GenomeWideResult(object):
	"""
	2009-4-24
		add data structure chr2min_max_pos to keep track of min,max chromosomal position all SNPs on that chromosome span
	2008-10-30
		no longer inherit from TableClass
	2008-10-21
		add option construct_data_obj_id2index
	"""
	data_obj_ls = None
	data_obj_id2index = None
	name = None
	results_method_id = None
	results_method = None
	min_value = None
	max_value = None
	chr_pos2index = None
	construct_data_obj_id2index = None	#True if get_data_obj_by_obj_id() is desired.
	construct_chr_pos2index = None	#True if get_data_obj_by_chr_pos() is desired.
	argsort_data_obj_ls = None
	chr2no_of_snps = None
	
	def __init__(self, construct_data_obj_id2index=True, construct_chr_pos2index=False, name=None, base_value = 0):
		"""
		2009-4-24
		2008-10-30
			no longer inherit from TableClass
			add this __init__()
		"""
		self.construct_data_obj_id2index = construct_data_obj_id2index	#True if get_data_obj_by_obj_id() is desired.
		self.construct_chr_pos2index = construct_chr_pos2index	#True if get_data_obj_by_chr_pos() is desired.
		self.name = name
		self.base_value = base_value
		self.chr2min_max_pos = {}
	
	def get_data_obj_by_obj_id(self, obj_id):
		return self.data_obj_ls[self.data_obj_id2index[obj_id]]
	
	def get_data_obj_by_obj_index(self, obj_index):
		return self.data_obj_ls[obj_index]
	
	def get_data_obj_by_chr_pos(self, chromosome, position):
		"""
		2008-09-24
		"""
		if self.chr_pos2index==None:
			return None
		else:
			obj_index = self.chr_pos2index.get((chromosome, position))
			if obj_index is not None:
				return self.data_obj_ls[obj_index]
	
	def add_one_data_obj(self, data_obj, chr_pos2index=None):
		"""
		2009-4-24 update self.chr2min_max_pos
		2008-10-23
			handle chr2no_of_snps
		2008-10-21
			add option chr_pos2index to put data_obj into data_obj_ls with pre-defined order
		2008-09-24
			add snippet to deal with construct_chr_pos2index and chr_pos2index
		"""
		if isinstance(chr_pos2index, dict):
			data_obj_index = chr_pos2index[(data_obj.chromosome, data_obj.position)]
			if not hasattr(self.data_obj_ls, '__len__') or len(self.data_obj_ls)==0:
				self.data_obj_ls = [None]*len(chr_pos2index)
			self.data_obj_ls[data_obj_index] = data_obj
		else:
			data_obj_index = len(self.data_obj_ls)
			self.data_obj_ls.append(data_obj)
		
		if self.construct_data_obj_id2index:	#2008-10-21
			if self.data_obj_id2index == None:
				self.data_obj_id2index = {}
			self.data_obj_id2index[id(data_obj)] = data_obj_index
		if self.construct_chr_pos2index:	#2008-09-24
			if self.chr_pos2index ==None:
				self.chr_pos2index = {}
			chr_pos = (data_obj.chromosome, data_obj.position)
			if chr_pos not in self.chr_pos2index:
				self.chr_pos2index[chr_pos] = data_obj_index
		if self.min_value is None or data_obj.value<self.min_value:
			self.min_value = data_obj.value
		if self.max_value is None or data_obj.value>self.max_value:
			self.max_value = data_obj.value
		
		if self.chr2no_of_snps is None:
			self.chr2no_of_snps = {}
		if data_obj.chromosome not in self.chr2no_of_snps:
			self.chr2no_of_snps[data_obj.chromosome] = 0
		self.chr2no_of_snps[data_obj.chromosome] += 1
		
		#2009-4-24 update self.chr2min_max_pos
		if data_obj.chromosome not in self.chr2min_max_pos:
			self.chr2min_max_pos[data_obj.chromosome] = [data_obj.position, data_obj.position]
		else:	# change the minimum and maximum position if condition is met
			if data_obj.position<self.chr2min_max_pos[data_obj.chromosome][0]:
				self.chr2min_max_pos[data_obj.chromosome][0]=data_obj.position
			if data_obj.position>self.chr2min_max_pos[data_obj.chromosome][1]:
				self.chr2min_max_pos[data_obj.chromosome][1]=data_obj.position
	
	def get_data_obj_at_given_rank(self, rank):
		"""
		2009-2-18
			if rank is beyond reach, return None
		2008-10-02
			rank starts from 1.
		"""
		if self.argsort_data_obj_ls is None:
			self.argsort_data_obj_ls = num.argsort(self.data_obj_ls)	#sort in ascending order
		if rank>len(self.data_obj_ls):
			return None
		else:
			return self.data_obj_ls[self.argsort_data_obj_ls[-rank]]	#value bigger, rank smaller
	
	def get_data_obj_index_given_rank(self, rank):
		"""
		2009-2-18
			if rank is beyond reach, return None
		2008-10-15
			similar to get_data_obj_at_given_rank() but instead of returning data_obj, it returns the index of data_obj in self.data_obj_ls
		"""
		if self.argsort_data_obj_ls is None:
			self.argsort_data_obj_ls = num.argsort(self.data_obj_ls)	#sort in ascending order
		if rank>len(self.data_obj_ls):
			return None
		else:
			return self.argsort_data_obj_ls[-rank]	#value bigger, rank smaller
	
class DataObject(object):
	"""
	2009-1-7
		
	2008-11-20
		add mac
	2008-11-12
		add genotype_var_perc, comment
	"""
	def __init__(self, **keywords):
		self.chromosome = None
		self.position = None
		self.stop_position = None
		self.name = None
		self.value = None
		self.genome_wide_result_id = None
		self.maf = None
		self.mac = None
		self.genotype_var_perc = None
		self.extra_col_ls = []
		self.comment = None
		for key, value in keywords.iteritems():
			setattr(self, key, value)
	
	def __cmp__(self, other):
		"""
		2008-08-20
			define how to compare DataObject
		"""
		return cmp(self.value, other.value)

import re
pa_has_characters = re.compile(r'[a-zA-Z_]')
import math

def getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=False, pdata=None):
	"""
	2009-2-18
		handle filtering by min_MAC (column_5th)
		column 6 and 7 are allowed to be empty placeholder. (previously it causes type cast error if it's empty)
	2008-12-18
		if pdata has attribute 'gwr_name', assign it to GenomeWideResult.name.
		otherwise GenomeWideResult.name = os.path.basename(input_fname)
	2008-11-20
		fix a bug that column_5th was skipped and column_6 was tried when there are only 5 columns
	2008-11-12
		parse lines with column_5th, column_6 and more
	2008-10-28
		handle construct_data_obj_id2index
	2008-10-21
		add score_for_0_pvalue, if pdata doesn't have it, assume 50.
		get chr_pos2index from pdata, which will pre-decide the order of snps in gwr.data_obj_ls
	2008-10-14
		get "is_4th_col_stop_pos" from pdata if it exists to decide how to deal with 4th col
	2008-09-24
		add construct_chr_pos2index
	2008-08-15
		skips non-positive values if need to do_log10_transformation
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
	construct_chr_pos2index = getattr(pdata, 'construct_chr_pos2index', False)	#2008-09-24
	construct_data_obj_id2index = getattr(pdata, 'construct_data_obj_id2index', True)	#2008-10-28 for get_data_obj_by_obj_index()
	is_4th_col_stop_pos = getattr(pdata, 'is_4th_col_stop_pos', False)	#2008-10-14
	chr_pos2index = getattr(pdata, 'chr_pos2index', None)	#2008-10-21
	score_for_0_pvalue = getattr(pdata, 'score_for_0_pvalue', 50)
	gwr_name = getattr(pdata, 'gwr_name', os.path.basename(input_fname))
	gwr = GenomeWideResult(name=gwr_name, construct_chr_pos2index=construct_chr_pos2index, \
						construct_data_obj_id2index=construct_data_obj_id2index)
	gwr.data_obj_ls = []	#list and dictionary are crazy references.
	gwr.data_obj_id2index = {}
	genome_wide_result_id = id(gwr)
	import csv
	delimiter = figureOutDelimiter(input_fname)
	reader = csv.reader(open(input_fname), delimiter=delimiter)
	no_of_lines = 0
	for row in reader:
		#check if 1st line is header or not
		if no_of_lines ==0 and pa_has_characters.search(row[1]):
			continue
		chr = int(row[0])
		start_pos = int(row[1])
		column_4th = None	#it's MAF probably
		column_5th = None	#it's MAC probably
		column_6 = None	#it's genotype_var_perc probably
		rest_of_row = []
		
		stop_pos = None
		if len(row)>=3:
			score = float(row[2])
		
		if len(row)>=4:
			if is_4th_col_stop_pos:
				stop_pos = int(row[3])
			else:
				column_4th=float(row[3])
		if len(row)>=5:
			column_5th = float(row[4])
		if len(row)>=6 and row[5]:
			column_6 = float(row[5])
		if len(row)>=7 and row[6]:
			rest_of_row = row[6:]
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
			pdata.min_MAC = getattr(pdata, 'min_MAC', None)	#2009-1-29
			if pdata.chromosome!=None and chr!=pdata.chromosome:
				continue
			if pdata.start!=None and start_pos<pdata.start:
				continue
			if pdata.stop!=None and start_pos>pdata.stop:
				continue
			if pdata.min_MAF!=None and column_4th!=None and column_4th<pdata.min_MAF:	#MAF too small
				continue
			if pdata.min_MAC!=None and column_5th!=None and column_5th<pdata.min_MAC:	#2009-1-29 MAC too small
				continue
		if do_log10_transformation:
			if score<=0:
				sys.stderr.write("score <=0. can't do log10. row is %s. assign %s to it.\n"%(repr(row), score_for_0_pvalue))
				#continue
				score = score_for_0_pvalue
			else:
				score = -math.log10(score)
		if min_value_cutoff is None or score>=min_value_cutoff:
			if stop_pos is not None:
				data_obj = DataObject(chromosome=chr, position=start_pos, stop_position=stop_pos, value =score)
			else:
				data_obj = DataObject(chromosome=chr, position=start_pos, value =score)
			if column_4th is not None:
				data_obj.maf = column_4th
			if column_5th is not None:
				data_obj.mac = column_5th
			if column_6 is not None:
				data_obj.genotype_var_perc = column_6
			if rest_of_row:
				data_obj.extra_col_ls = rest_of_row
				data_obj.comment = ','.join(rest_of_row)
			data_obj.genome_wide_result_id = genome_wide_result_id
			gwr.add_one_data_obj(data_obj, chr_pos2index)
		
		no_of_lines += 1
		
	del reader
	sys.stderr.write(" %s results. Done.\n"%(len(gwr.data_obj_ls)))
	return gwr


class SNPInfo(object):
	"""
	2009-2-18
		a class to hold chromosome, position, allele, snps_id (db)
		DrawSNPRegion.getSNPInfo(db) does the job of filling it up
	"""
	chr_pos_ls = None
	chr_pos2index = None
	snps_id2index = None
	data_ls = None	#a list of [snps_d, chromosome, position, allele1, allele2]
	
	def __init__(self, **keywords):
		"""
		2009-2-18 allow any type of keywords
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
	
	def getSnpsIDGivenChrPos(self, chromosome, position):
		snp_info_index = self.chr_pos2index.get((chromosome, position))
		if snp_info_index is not None:
			snps_id = self.data_ls[snp_info_index][0]
		else:
			snps_id = None
		return snps_id
	
	def getSnpsAllelesGivenChrPos(self, chromosome, position):
		snp_info_index = self.chr_pos2index.get((chromosome, position))
		if snp_info_index is not None:
			alleles = self.data_ls[snp_info_index][3:5]
		else:
			alleles = None
		return alleles