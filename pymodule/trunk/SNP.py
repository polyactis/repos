#!/usr/bin/env python
"""
2008-05-12
	store SNP data structure related stuff
"""

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