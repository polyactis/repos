#!/usr/bin/env python
"""

Examples:
	#compare 250k (output by Output250KSNPs.py) and a data merging 2010/149/384
	QC.py -i /mnt/nfs/NPUTE_data/input/250K_m3_70_n1000.csv -k 3 -j /mnt/nfs/NPUTE_data/input/2010_149_384.csv -l 2 -o /tmp/250K_m3_70_n1000_vs_2010_149_384.QC.tsv -b
	
	#compare 250k (output by Output250KSNPs.py) and 384-illumina (output by dbSNP2data.py)
	QC.py -i /mnt/nfs/NPUTE_data/input/250K_m3_70.csv -j /home/crocea/script/variation/genotyping/384-illumina_y0000111111.csv -k 3 -l 2 -o /tmp/250K_m3_70_vs_384-illumina.QC.csv
	
	#compare 250k (output by DB_250k2data.py) with 2010 (output by Output2010InCertainSNPs.py)
	QC.py -i /mnt/nfs/NPUTE_data/250k_l3_y0.85_v0.6_w0.2_x0.2.tsv -j data/2010/data_2010_x_250k.tsv -o /mnt/nfs/NPUTE_data/QC_output/250k_l3_y0.85_v0.6_w0.2_x0.2_vs_2010_x_250k.QC.csv

Description:
	a program to do general QC given two files.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions
from variation.src.QualityControl import QualityControl
import dataParsers, FilterAccessions, FilterSnps, MergeSnpsData
from common import RawSnpsData_ls2SNPData, transposeSNPData
from QC_250k import QC_250k, SNPData, TwoSNPData

class QC(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname1',1, ): [None, 'i', 1, 'File containing the calls'],\
							('input_fname2',1, ): [None, 'j', 1, 'File containing the reference calls to be compared with'],\
							('input_fname1_format',1,int): [1, 'k', 1, 'Format of input_fname1. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId. 3=snp X strain with arrayId.'],\
							('input_fname2_format',1,int): [1, 'l', 1, 'Format of input_fname2. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId'],\
							('output_fname', 1, ): [None, 'o', 1, '', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-05-11
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		
		#to check whether two input file are in different orientation
		file_format2count = {}
		file_format_ls = [self.input_fname1_format, self.input_fname2_format]
		for file_format in file_format_ls:
			if file_format not in file_format2count:
				file_format2count[file_format] = 0
			file_format2count[file_format] += 1
		
		if 1 in file_format2count and file_format2count[1]==1:	#there's one and only one strain x snp format.
			#it needs transpose matrix. only numpy works on character matrix. not sure Numeric or numarray is imported. so transform the input matrix to integer.
			use_nt2number = 1
		else:
			use_nt2number = 1
		
		if self.input_fname1_format==1:
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_fname1)
			snpData1 = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
		elif self.input_fname1_format==2:
			snpsd_ls = dataParsers.parseCSVData(self.input_fname1, withArrayIds=False)
			snpData1 = RawSnpsData_ls2SNPData(snpsd_ls, report=self.report, use_nt2number=use_nt2number)
			del snpsd_ls
		elif self.input_fname1_format==3:
			snpsd_ls = dataParsers.parseCSVData(self.input_fname1, withArrayIds=True)
			snpData1 = RawSnpsData_ls2SNPData(snpsd_ls, report=self.report, use_nt2number=use_nt2number)
			del snpsd_ls
		else:
			sys.stderr.write('Error: unsupported input_fname1 format, %s\n' % self.input_fname1_format)
			sys.exit(2)
		
		if self.input_fname2_format==1:
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(self.input_fname2)
			snpData2 = SNPData(header=header, strain_acc_list=strain_acc_list,\
							data_matrix=data_matrix)
		elif self.input_fname2_format==2:
			snpsd_ls = dataParsers.parseCSVData(self.input_fname2, withArrayIds=False)
			snpData2 = RawSnpsData_ls2SNPData(snpsd_ls, report=self.report, use_nt2number=use_nt2number)
			del snpsd_ls
		else:
			sys.stderr.write('Error: unsupported input_fname2 format, %s\n' % self.input_fname2_format)
			sys.exit(2)
		

		if 1 in file_format2count and file_format2count[1]==1:	#there's one and only one strain x snp format. transpose the 2nd snpData
			snpData2 = transposeSNPData(snpData2, report=self.report)
		
		if self.input_fname1_format == 1:	#row_id for the 1st file = (ecotype_id, duplicate). for 2nd file, row_id=ecotype_id.
			row_matching_by_which_value = 0
			col_matching_by_which_value = None
		elif self.input_fname1_format == 2:	#col_id for the 1st file = accession. for 2nd file, col_id=accession.
			row_matching_by_which_value = None
			col_matching_by_which_value = None
		elif self.input_fname1_format == 3:	#col_id for the 1st file = (array_id, accession). for 2nd file, col_id=accession.
			row_matching_by_which_value = None
			col_matching_by_which_value = 1
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, row_matching_by_which_value=row_matching_by_which_value,\
							col_matching_by_which_value=col_matching_by_which_value, debug=self.debug)
		row_id2NA_mismatch_rate = twoSNPData.cmp_row_wise()
		col_id2NA_mismatch_rate = twoSNPData.cmp_col_wise()
		if row_id2NA_mismatch_rate:
			QC_250k.output_row_id2NA_mismatch_rate(row_id2NA_mismatch_rate, self.output_fname, file_1st_open=1)
		if col_id2NA_mismatch_rate:
			QC_250k.output_row_id2NA_mismatch_rate(col_id2NA_mismatch_rate, self.output_fname, file_1st_open=0)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = QC
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()