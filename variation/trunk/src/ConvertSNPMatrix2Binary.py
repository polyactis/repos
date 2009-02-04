#!/usr/bin/env python
"""

Examples:
	ConvertSNPMatrix2Binary -i ./panfs/250k/call_method_17.tsv  -o ./tmp/250k/call_method_17_binary.tsv -m ./tmp/250k/call_method_17_binary_map.tsv -r
	
Description:
	2008-09-7 Convert SNP matrix into index (0,1,2...) is assigned as first-encounter, first-assign. if only two alleles, it's binary.
		heterozygote is regarded as a different allele.
			
	Input format is strain X snp. 2nd column ignored by default (change by array_id_2nd_column).
		delimiter (either tab or comma) is automatically detected. 	header is 'chromosome_position'.
	
	Output format is same as input. Another output file holds the mapping between each allele and the index.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
from pymodule import process_function_arguments, write_data_matrix, figureOutDelimiter, read_data, SNPData
import csv

class ConvertSNPMatrix2Binary(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): [None, 'i', 1, 'input file'],\
						('output_fname', 1, ): [None, 'o', 1, 'Output Filename'],\
						('mapping_fname', 1, ): [None, 'm', 1, 'the file to store the mapping between SNP allele and 0/1'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-09-7
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def output_allele2index_ls(self, snpData, allele_index2allele_ls, mapping_fname):
		"""
		2009-2-3
			key/value order is reversed in allele_index2allele_ls. argument allele2index_ls is renamed to allele_index2allele_ls.
		2008-09-07
			output the mapping from allele to index
		"""
		sys.stderr.write("Outputting allele_index2allele_ls ...")
		writer = csv.writer(open(mapping_fname, 'w'), delimiter='\t')
		writer.writerow(['snp_id', 'allele', 'index'])
		for i in range(len(allele_index2allele_ls)):
			snp_id = snpData.col_id_ls[i]
			for allele_index, allele in allele_index2allele_ls[i].iteritems():
				row = [snp_id, allele, allele_index]
				writer.writerow(row)
		del writer
		sys.stderr.write("Done.\n")
		
	
	def run(self):
		"""
		2008-9-7
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		delimiter = figureOutDelimiter(self.input_fname, report=self.report)
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, delimiter=delimiter)
		
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix)
		newSnpData, allele_index2allele_ls = snpData.convert2Binary(self.report)
		
		if self.mapping_fname:	#output allele_index2allele_ls
			self.output_allele2index_ls(snpData, allele_index2allele_ls, self.mapping_fname)
		
		newSnpData.tofile(self.output_fname)

		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConvertSNPMatrix2Binary
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()