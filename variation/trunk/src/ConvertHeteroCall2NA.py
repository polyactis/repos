#!/usr/bin/env python
"""

Examples:
	ConvertHeteroCall2NA.py  -i bin/structure_test/149_popid2ecotypeid_50.csv -o bin/structure_test/149_popid2ecotypeid_50_hets2NA.csv
	
Description:
	Convert all heterozygous calls and untouched in the file into NA.
	
"""
from __init__ import *

class ConvertHeteroCall2NA(object):
	__doc__ = __doc__
	option_default_dict = {	('input_fname',1, ): [None, 'i', 1, ''],\
							('output_fname', 1, ): [None, 'o', 1, '', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-06-02
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		snpData = SNPData(input_fname=self.input_fname, turn_into_array=1)
		newSnpData = SNPData.convertHetero2NA(snpData)
		newSnpData.tofile(self.output_fname)

if __name__ == '__main__':
	main_class = ConvertHeteroCall2NA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()