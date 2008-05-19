#!/usr/bin/env python
"""
2008-05-12
	store SNP data structure related stuff
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray

def write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=None, \
					cols_to_be_tossed_out=None, nt_alphabet=0, transform_to_numpy=1,\
					discard_all_NA_rows=0, strain_acc2other_info=None, delimiter='\t', predefined_header_row=['strain', 'duplicate', 'latitude', 'longitude', 'nativename', 'stockparent', 'site', 'country']):
	"""
	strain_acc_list (and category_list) are initial 2 columns in the output.
	
	rows_to_be_tossed_out is a set or dict with row index in it. cols_to_be_tossed_out is similar structure.
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
	sys.stderr.write("Writing data_matrix ...")
	from variation.src.common import number2nt
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

def read_data(input_fname, input_alphabet=0, turn_into_integer=1, double_header=0, delimiter=None):
	"""
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
	from variation.src.common import nt2number
	#from __init__ import dict_map	#already imported by from __init__ import *
	sys.stderr.write("Reading data ...")
	if delimiter is None:
		delimiter = figureOutDelimiter(input_fname)
	reader = csv.reader(open(input_fname), delimiter=delimiter)
	header = reader.next()
	if double_header:
		header = [header, reader.next()]
	data_matrix = []
	strain_acc_list = []
	category_list = []
	for row in reader:
		strain_acc_list.append(row[0])
		category_list.append(row[1])
		data_row = row[2:]
		no_of_snps = len(data_row)
		if input_alphabet:
			data_row = dict_map(nt2number, data_row)
			if no_of_snps!=len(data_row):
				print row
		else:
			if turn_into_integer:
				data_row = map(int, data_row)
		data_matrix.append(data_row)
	del reader
	sys.stderr.write("Done.\n")
	return header, strain_acc_list, category_list, data_matrix

class SNPData(object):
	"""
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
			if ignore_2nd_column:
				self.category_list = None
				
		if not self.data_matrix and self.turn_into_array:
			num = importNumericArray()
			self.data_matrix = num.array(self.data_matrix)
		
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

try:
	from QualityControl import QualityControl
except:
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
	try:
		from QualityControl import QualityControl
	except:
		QualityControl = object
if QualityControl!=object:
	from TwoSNPData import TwoSNPData