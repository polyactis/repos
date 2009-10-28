#!/usr/bin/env python
"""

Examples:	
	#normalize + split data into different chromosomes 
	CNVNormalize.py -i call_method_17_CNV_array_intensity.tsv -o call_method_17_CNV_array_intensity_norm
	
	#2009-9-26 just split data into different chromosomes (output_fname_prefix=input_fname)
	CNVNormalize.py -i call_43_first_311/call_method_43_CNV_intensity.tsv -o call_43_first_311/call_method_43_CNV_intensity -y 00001
	
	
	
Description:
	2009-10-27 program to normalize CNV intensity matrix (outputted by DB_250k2Array.py). Each bit in processing_bits (0=off, 1=on) corresponds to:
		1. quantile normalize (replace the distribution of one array with the mean across all arrays, while preserving the rank of each probe)
		2. subtract reference median
		3. subtract column median
		4. subtract row mean (median if numpy version is not 1.0.x)
		5. split whole-genome data into different chromosomes. data from one chromosome is in one file.
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
from pymodule import figureOutDelimiter, getListOutOfStr

class CNVNormalize(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): ['', 'i', 1, 'CNV intensity matrix, probe X arrays. 1st column is probe id. 2nd last col is chr. last col is pos.', ],\
							('output_fname_prefix', 1, ): ['', 'o', 1, 'Output filename prefix. ".tsv" would be appended if not included.', ],\
							('reference_arrays', 0, ):['1,2,43,139,145', '', 1, 'a string of array ids whose median is used as reference. default is Col arrays.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							('processing_bits', 1, ): ['11001', 'y', 1, 'processing bits to control which processing step should be turned on. check Description for details.' ]}
	
	def __init__(self, **keywords):
		"""
		2008-12-05
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		processing_bits_ls = [1,1,1,1,1]
		
		for i in range(len(self.processing_bits)):
			processing_bits_ls[i] = int(self.processing_bits[i])
		#now pass all values
		self.quantile_normalization,\
		self.need_to_subtract_ref,\
		self.subtract_column_median,\
		self.subtract_row_median,\
		self.split_genome_into_chromosomes = processing_bits_ls
		
		self.reference_arrays = getListOutOfStr(self.reference_arrays, data_type=int)
		
	@classmethod
	def get_input(cls, input_fname, data_type=numpy.float32):
		"""
		2009-10-28
			switch the default data_type to numpy.float32 to save memory on 64bit machines
		2009-9-28
			add argument data_type to specify data type of data_matrix.
			default is numpy.float (numpy.float could be float32, float64, float128 depending on the architecture).
				numpy.double is also fine.
		2009-5-18
			become classmethod
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
		data_matrix = numpy.zeros([no_of_rows, no_of_cols], data_type)
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
	
	def quantile_normalize(self, data_matrix, no_of_probes_per_block=300000):
		"""
		2009-10-28
			use numpy.mean and numpy.median depending on the version
			add argument no_of_probes_per_block and get the median array block by block to save memory.
				median() costs more memory than mean() and program blows up.
		"""
		sys.stderr.write("Quantile normalizing ")
		
		sort_data_matrix = numpy.sort(data_matrix, axis=0)	#sort each column
		no_of_rows, no_of_cols = data_matrix.shape
		if numpy.version.version[:3]=='1.0':	# 2009-10-27, 1.0.x can't get row-wise median 
			sys.stderr.write("by mean ...")
			mean_ar = numpy.mean(sort_data_matrix, axis=1)	#get each row's mean
		else:
			sys.stderr.write("by median ...")
			no_of_blocks = int(math.ceil(no_of_rows/float(no_of_probes_per_block)))
			median_ar_ls = []
			for i in range(no_of_blocks):
				row_start_index = i*no_of_probes_per_block
				row_stop_index = (i+1)*no_of_probes_per_block
				median_ar = numpy.median(sort_data_matrix[row_start_index:row_stop_index,:], axis=1)	#get each row's median
				median_ar_ls.append(median_ar)
			mean_ar = numpy.hstack(median_ar_ls)	#get each row's median
		del sort_data_matrix
		argsort_data_matrix = numpy.argsort(data_matrix, axis=0)
		
		for j in range(no_of_cols):
			for i in range(no_of_rows):
				which_row_index = argsort_data_matrix[i][j]
				data_matrix[which_row_index][j] = mean_ar[i]	#replace with mean
		del argsort_data_matrix
		sys.stderr.write("Done.\n")
	
	def subtract_reference_median(self, data_matrix, array_id_ls, reference_arrays):
		"""
		2009-10-27
			subtract reference median calculated based on reference_arrays from data_matrix
		"""
		sys.stderr.write("Subtracting reference ")
		reference_array_set = set(reference_arrays)
		reference_array_index_ls = []
		for i in range(len(array_id_ls)):
			array_id = array_id_ls[i]
			if array_id in reference_array_set:
				reference_array_index_ls.append(i)
		
		if numpy.version.version[:3]=='1.0':	# numpy version 1.0.x, doesn't accept axis. only column-wise, no row-wise
			sys.stderr.write("mean ...")
			ref_median_ar = numpy.mean(data_matrix[:, reference_array_index_ls], axis=1)
		else:
			sys.stderr.write("median ...")
			ref_median_ar = numpy.median(data_matrix[:, reference_array_index_ls], axis=1)	# version above 1.0. without axis=0, it's getting median of a flattened matrix.
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			data_matrix[i,:] = data_matrix[i,:]-ref_median_ar[i]
		sys.stderr.write("%s reference arrays. Done.\n"%(len(reference_array_index_ls)))
	
	def subtract_col_mean(self, data_matrix):
		sys.stderr.write("Subtracting column ")
		if numpy.version.version[:3]=='1.0':	# numpy version 1.0.x, doesn't accept axis. only column-wise, no row-wise
			sys.stderr.write("mean ...")
			col_mean_ar = numpy.median(data_matrix)
		else:
			sys.stderr.write("median ...")
			col_mean_ar = numpy.median(data_matrix, axis=0)	# version above 1.0. without axis=0, it's getting median of a flattened matrix.
		no_of_rows, no_of_cols = data_matrix.shape
		for j in range(no_of_cols):
			data_matrix[:,j] = data_matrix[:,j]-col_mean_ar[j]
		sys.stderr.write("Done.\n")
	
	def subtract_row_mean(self, data_matrix):
		sys.stderr.write("Subtracting column mean ...")
		if numpy.version.version[:3]=='1.0':
			row_mean_ar = numpy.mean(data_matrix, axis=1)
		else:
			row_mean_ar = numpy.median(data_matrix, axis=1)
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			data_matrix[i,:] = data_matrix[i,:]-row_mean_ar[i]
		sys.stderr.write("Done.\n")
	
	def output(self, data_matrix, probe_id_ls, chr_pos_ls, header, output_fname_prefix, split_genome_into_chromosomes=False):
		"""
		2009-10-11
			add argument split_genome_into_chromosomes
		2009-5-18
			split output into different chromosomes
		"""
		sys.stderr.write("Outputting ...")
		no_of_rows, no_of_cols = data_matrix.shape
		old_chr = None
		old_writer = None
		output_fname_prefix = os.path.splitext(output_fname_prefix)[0]
		for i in range(no_of_rows):
			new_chr = chr_pos_ls[i][0]
			if split_genome_into_chromosomes:
				if old_chr==None or old_chr!=new_chr:
					writer = csv.writer(open('%s_chr%s.tsv'%(output_fname_prefix, new_chr), 'w'), delimiter='\t')
					writer.writerow(header)
					old_chr = new_chr
					del old_writer	#close the old file
					old_writer = writer
			elif old_writer is None:
				old_writer = csv.writer(open('%s.tsv'%(output_fname_prefix), 'w'), delimiter='\t')
				old_writer.writerow(header)
			
			data_row = [probe_id_ls[i]]
			for j in range(no_of_cols):
				data_row.append(data_matrix[i][j])
			data_row.append(chr_pos_ls[i][0])
			data_row.append(chr_pos_ls[i][1])
			old_writer.writerow(data_row)
		del old_writer
		
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
		array_id_ls = header[1:-2]
		array_id_ls = map(int, array_id_ls)
		
		if self.quantile_normalization:
			self.quantile_normalize(data_matrix)
		if self.need_to_subtract_ref:
			self.subtract_reference_median(data_matrix, array_id_ls, self.reference_arrays)
		if self.subtract_column_median:
			self.subtract_col_mean(data_matrix)
		if self.subtract_row_median:
			self.subtract_row_mean(data_matrix)
		
		self.output(data_matrix, probe_id_ls, chr_pos_ls, header, self.output_fname_prefix,\
				split_genome_into_chromosomes=self.split_genome_into_chromosomes)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CNVNormalize
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	