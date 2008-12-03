#!/usr/bin/env python
"""

Examples:
	ConvertYuSNPFormat2EigenStrat.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -o /tmp/call_method_17_binary_eigenstrat

Description:
	Convert Yu's SNP format (input) to EigenStrat's (Price2006, http://genepath.med.harvard.edu/~reich/Software.htm).
	
	EigenStrat is one program in the Price2006's package, EigenSoft. smartpca.perl in EigenSoft uses the same format as EigenStrat.

	Input format is strain X snp. 2nd column ignored by default (change by array_id_2nd_column).
		delimiter (either tab or comma) is automatically detected. 	header is 'chromosome_position'.
	
	There are 3 output files.
		.geno
			The genotype file contains 1 line per SNP.
				Each line contains 1 character per individual:
				0 means zero copies of reference allele.
				1 means one copy of reference allele.
				2 means two copies of reference allele.
				9 means missing data.
		
		.ind
			The indiv file contains 1 line per individual.  There are 3 columns:
				1st column is sample ID
				2nd column is gender (M or F).  If unknown, ok to set to U for Unknown.
				3rd column is a label which might refer to Case or Control status, or
					might be a population group label.  If this entry is set to "Ignore",
					then that individual and all genotype data from that individual will be
					removed from the data set in all convertf output.
		.snp
			The snp file contains 1 line per SNP.  There are 6 columns (last 2 optional):
				1st column is SNP name
				2nd column is chromosome.  X chromosome is encoded as 23.
					Also, Y is encoded as 24, mtDNA is encoded as 90, and XY is encoded as 91.
					Note: SNPs with illegal chromosome values, such as 0, will be removed
				3rd column is genetic position (in Morgans).  If unknown, ok to set to 0.0.
				4th column is physical position (in bases)
					Optional 5th and 6th columns are reference and variant alleles.
					For monomorphic SNPs, the variant allele can be encoded as X.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
from pymodule import process_function_arguments, write_data_matrix, figureOutDelimiter, read_data
from pymodule.SNP import transposeSNPData
from QC_250k import SNPData
from common import SNPData2RawSnpsData_ls
import snpsdata, csv

class ConvertYuSNPFormat2EigenStrat(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): [None, 'i', 1, ],\
						('output_fname_prefix', 1, ): [None, 'o', 1, 'Output Filename prefix'],\
						('array_id_2nd_column', 0, int): [0, 'a', 0, 'whether 2nd column in input_fname is array id or not'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-12-02
			modelled after ConvertYuSNPFormat2Bjarni.py
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-12-02
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		delimiter = figureOutDelimiter(self.input_fname, report=self.report)
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, delimiter=delimiter)
		
		if self.array_id_2nd_column:
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
		else:
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list,\
							data_matrix=data_matrix)	#ignore category_list
		
		newSnpData, allele_index2allele_ls = snpData.convert2Binary(self.report)
		
		genotype_f = open('%s.geno'%self.output_fname_prefix, 'w')
		ind_writer = csv.writer(open('%s.ind'%self.output_fname_prefix, 'w'), delimiter='\t')
		snp_writer = csv.writer(open('%s.snp'%self.output_fname_prefix, 'w'), delimiter='\t')
		
		#transpose it
		newSnpData = transposeSNPData(newSnpData)
		
		no_of_rows = len(newSnpData.data_matrix)
		no_of_cols = len(newSnpData.data_matrix[0])
		for i in range(no_of_rows):
			snp_id = newSnpData.row_id_ls[i]
			chr, pos = snp_id.split('_')
			allele1 = allele_index2allele_ls[i][0]	#major allele
			allele2 = allele_index2allele_ls[i][1]	#minor allele
			snp_writer.writerow([snp_id, chr, 0.0, pos, allele1, allele2])
			geno_line = ''
			for j in range(no_of_cols):
				if i==0:	#write out the accessions
					ind_writer.writerow([newSnpData.col_id_ls[j], 'U', 'Case'])
				allele = newSnpData.data_matrix[i][j]
				if allele==0:
					geno_line += '0'
				elif allele==1:
					geno_line += '2'
				else:
					geno_line += '9'
			geno_line += '\n'
			genotype_f.write(geno_line)
		
		del genotype_f, ind_writer, snp_writer

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConvertYuSNPFormat2EigenStrat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
