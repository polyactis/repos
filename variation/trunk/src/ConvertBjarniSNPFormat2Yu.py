#!/usr/bin/env python
"""

Examples:
	./src/ConvertBjarniSNPFormat2Yu.py -i /mnt/nfs/NPUTE_data/input/2010_149_384_v3.csv -o /tmp/2010_149_384_v3.tsv
	
Description:
	Convert Bjarni's SNP format (input) to Yu's.
	
	Input format is snp X strain, bjarni's format. The extra first row for array id is optional. (change by withArrayIds).
	Output format is strain X snp. tab-delimiter. header is 'chromosome_position'.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import dataParsers
from common import RawSnpsData_ls2SNPData, transposeSNPData

class ConvertBjarniSNPFormat2Yu(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): [None, 'i', 1, ],\
						('output_fname', 1, ): [None, 'o', 1, 'Output Filename'],\
						('withArrayIds', 0, int): [0, 'a', 0, 'does the input_fname come with a line of array ids?'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-05-18
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-5-18
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		snpsd_ls = dataParsers.parseCSVData(self.input_fname, withArrayIds=self.withArrayIds)
		snpData = RawSnpsData_ls2SNPData(snpsd_ls, use_nt2number=1)
		del snpsd_ls
		newSnpData = transposeSNPData(snpData)
		del snpData
		newSnpData.tofile(self.output_fname, transform_to_numpy=0)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConvertBjarniSNPFormat2Yu
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
