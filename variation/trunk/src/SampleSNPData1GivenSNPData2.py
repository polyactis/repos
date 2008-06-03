#!/usr/bin/env python
"""

Examples:
	./src/SampleSNPData1GivenSNPData2.py -i /Network/Data/250k/db/reference_dataset/perlegen -j bin/structure_test/149_y0000110101_removeBadSNPs_mergedup_popid2ecotypeid_50_hets2NA.csv -o bin/structure_test/perlegen_50loci_around_149SNP.tsv
	
Description:
	Sample SNP loci from SNPData1 around each SNPData2 loci 
	
"""
from __init__ import *

class SampleSNPData1GivenSNPData2(object):
	__doc__ = __doc__
	option_default_dict = {	('input_fname1',1, ): [None, 'i', 1, 'to form SNPData1'],\
							('input_fname2',1, ): [None, 'j', 1, 'to form SNPData2'],\
							('output_fname', 1, ): [None, 'o', 1, '', ],\
							('no_of_loci_to_sample_around_on_each_side', 1, int ): [50, 'n', 1, 'For each loci in SNPData2, how many around it should be sampled', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-06-02
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		snpData1 = SNPData(input_fname=self.input_fname1, turn_into_array=1, ignore_2nd_column=1)
		snpData2 = SNPData(input_fname=self.input_fname2, turn_into_array=1, ignore_2nd_column=1)
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, debug=self.debug)
		cols_to_be_tossed_out = twoSNPData.sampleSNPLociFromSNPData1(self.no_of_loci_to_sample_around_on_each_side)
		snpData1.tofile(self.output_fname, cols_to_be_tossed_out=cols_to_be_tossed_out)
		

if __name__ == '__main__':
	main_class = SampleSNPData1GivenSNPData2
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()