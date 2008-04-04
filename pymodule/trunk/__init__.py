"""
2007-10-15. __init__.py for pymodule. pymodule is a concatenation of all common functions/classes.
"""

def process_function_arguments(keywords, argument_default_dict, error_doc='', class_to_have_attr=None):
	"""
	2008-04-02
		add class_to_have_attr to assign the argument and values
	2008-02-28
		process arguments from keywords based on argument_default_dict
		
		keywords is a dictionary of given arguments. 'argument_name':'argument_value'
		
		argument_default_dict is a dictionary of default arguments.
			the key is a tuple, ('argument_name', is_argument_required, argument_type) and argument_type is optional.
	"""
	import sys
	ad = {}
	for argument_key, default_value in argument_default_dict.iteritems():
		argument, is_argument_required = argument_key[0:2]
		if len(argument_key)>2:
			argument_type = argument_key[2]
		else:
			argument_type = None
		if keywords.has_key(argument):
			if keywords[argument]!='' and keywords[argument]!=None:	#only when keywords has this argument and it's not nothing. change default_value
				if argument_type!=None:	#cast to the desired type
					default_value = argument_type(keywords[argument])
				else:
					default_value = keywords[argument]
		if is_argument_required==1 and (default_value=='' or default_value==None):
			if error_doc:
				sys.stderr.write(error_doc)
			sys.stderr.write('Error: %s is required but %s given.\n'%(argument, default_value))
			sys.exit(2)
		else:
			if class_to_have_attr:
				setattr(class_to_have_attr, argument, default_value)
			ad[argument] = default_value
	return ad


def write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=None, cols_to_be_tossed_out=None, nt_alphabet=0):
	"""
	2008-04-02
		extracted from variation.src.FilterStrainSNPMatrix to be standalone.
	"""
	from sets import Set
	import sys, csv
	sys.stderr.write("Writing data_matrix ...")
	if rows_to_be_tossed_out==None:
		rows_to_be_tossed_out = Set()
	if cols_to_be_tossed_out==None:
		cols_to_be_tossed_out = Set()
	from variation.src.common import number2nt
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	
	new_header = [header[0], header[1]]
	for i in range(2, len(header)):
		if i-2 not in cols_to_be_tossed_out:
			new_header.append(header[i])
	writer.writerow(new_header)
	
	if type(data_matrix)==list:	#2008-02-06 transform the 2D list into array
		import numpy
		data_matrix = numpy.array(data_matrix)
	no_of_rows, no_of_cols = data_matrix.shape
	for i in range(no_of_rows):
		if i not in rows_to_be_tossed_out:
			new_row = [strain_acc_list[i], category_list[i]]
			for j in range(no_of_cols):
				if j not in cols_to_be_tossed_out:
					if nt_alphabet:
						new_row.append(number2nt[data_matrix[i][j]])
					else:
						new_row.append(data_matrix[i][j])
			writer.writerow(new_row)
	del writer
	sys.stderr.write("Done.\n")

def dict_map(dict, ls, type=1):
	"""
	2008-04-03 copied from annot.codense.common
	10-13-05
		add type 2 to return item itself if mapping is not available
	2006-12-21
		add type 3 to extract a smaller map
	2007-05-14
		bug, "if value" could miss 0
	"""
	if type==3:
		new_list = {}	#it's a dictionary
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list[item] = value
	else:
		new_list = []
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list.append(value)
			elif type==2:
				new_list.append(item)
	
	return new_list