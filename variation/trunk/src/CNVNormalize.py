#!/usr/bin/env python
"""

Examples:
	CNVNormalize.py -i -o
	CNVNormalize.py -i call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622.tsv -o call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm.tsv
	
Description:
	2008-12-05 program to normalize CNV intensity matrix (outputted by DB_250k2Array.py)
		1. quantile normalize
		2. substract column median
		3. subtract row mean
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter

class CNVNormalize(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): ['', 'i', 1, 'CNV intensity matrix, probe X arrays. 1st column is probe id. 2nd last col is chr. last col is pos.', ],\
							('output_fname', 1, ): ['', 'o', 1, 'self-explanatory', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_input(self, input_fname):
		"""
		
		"""
		sys.stderr.write("Getting input from %s ..."%input_fname)
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		commandline = 'wc -l %s'%input_fname
		command_handler = subprocess.Popen(commandline, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
		stdout_content, stderr_content = command_handler.communicate()
		if stderr_content:
			sys.stderr.write('stderr of %s: %s \n'%(commandline, stderr_content))
		no_of_rows = int(stdout_content.split()[0])-1
		
		header = reader.next()
		no_of_cols = len(header)-3
		data_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.float)
		probe_id_ls = []
		chr_pos_ls = []
		i=0
		for row in reader:
			
			probe_id = row[0]
			probe_id_ls.append(probe_id)
			chr_pos_ls.append(row[-2:])
			for j in range(1, 1+no_of_cols):
				data_matrix[i][j-1] = float(row[j])
			i += 1
		sys.stderr.write("Done.\n")
		return data_matrix, probe_id_ls, chr_pos_ls, header
	
	def quantile_normalize(self, data_matrix):
		sys.stderr.write("Quantile normalizing ...")
		
		sort_data_matrix = numpy.sort(data_matrix, axis=0)	#sort each column
		mean_ar = numpy.mean(sort_data_matrix, axis=1)	#get each row's mean
		del sort_data_matrix
		no_of_rows, no_of_cols = data_matrix.shape
		argsort_data_matrix = numpy.argsort(data_matrix, axis=0)
		
		for j in range(no_of_cols):
			for i in range(no_of_rows):
				which_row_index = argsort_data_matrix[i][j]
				data_matrix[which_row_index][j] = mean_ar[i]	#replace with mean
		del argsort_data_matrix
		sys.stderr.write("Done.\n")
	
	def subtract_col_mean(self, data_matrix):
		sys.stderr.write("Subtracting column mean ...")
		col_mean_ar = numpy.median(data_matrix)	#numpy.median doesn't accept axis. only column-wise, no row-wise
		no_of_rows, no_of_cols = data_matrix.shape
		for j in range(no_of_cols):
			data_matrix[:,j] = data_matrix[:,j]-col_mean_ar[j]
		sys.stderr.write("Done.\n")
	
	def subtract_row_mean(self, data_matrix):
		sys.stderr.write("Subtracting column mean ...")
		row_mean_ar = numpy.mean(data_matrix, axis=1)
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			data_matrix[i,:] = data_matrix[i,:]-row_mean_ar[i]
		sys.stderr.write("Done.\n")
	
	def output(self, data_matrix, probe_id_ls, chr_pos_ls, header, output_fname):
		sys.stderr.write("Outputting ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(header)
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			data_row = [probe_id_ls[i]]
			for j in range(no_of_cols):
				data_row.append(data_matrix[i][j])
			data_row.append(chr_pos_ls[i][0])
			data_row.append(chr_pos_ls[i][1])
			writer.writerow(data_row)
		del writer
		
		sys.stderr.write("Done.\n")
	
	def run(self):
		data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
		self.quantile_normalize(data_matrix)
		self.subtract_col_mean(data_matrix)
		self.subtract_row_mean(data_matrix)
		self.output(data_matrix, probe_id_ls, chr_pos_ls, header, self.output_fname)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVNormalize
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	