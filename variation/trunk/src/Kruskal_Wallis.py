#!/usr/bin/env python
"""

Examples:
	./src/Kruskal_Wallis.py -i /tmp/250K_method_5_after_imputation_noRedundant_051908.tsv -p /Network/Data/250k/finalData_051808/phenotypes.tsv -e -r -o /tmp/250K_method_5_after_imputation_noRedundant_051908.LD.pvalue

Description:
	class to do kruskal wallis test on SNP data.
	
	Input genotype file format is Strain X SNP format (Yu's format, Output by DB_250k2data.py Or Output250KSNPs.py + ConvertBjarniSNPFormat2Yu.py).
	Input phenotype file format is Strain X phenotype format (Output by OutputPhenotype.py). 
	
	It requires a minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis test.
	
	It will automatically match strains in two files. NO worry for missing/extra data in either input file.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy
from pymodule import read_data, ProcessOptions, PassingData

class Kruskal_Wallis:
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): ['', 'i', 1, 'input genotype matrix. Strain X SNP format.', ],\
							('output_fname', 1, ): ['', 'o', 1, 'store the pvalue', ],\
							('phenotype_fname', 1, ): [None, 'p', 1, 'phenotype file', ],\
							('minus_log_pvalue', 0, ): [0, 'e', 0, 'toggle -log(pvalue)', ],\
							('which_phenotype', 1, int): [0, 'w', 1, 'which phenotype, 0=first phenotype (3rd column in phenotype_fname) and so on.',],\
							('min_data_point', 1, int): [3, 'm', 1, 'minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis test'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-02-14
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		"""
		if not input_fname or not phenotype_fname or not output_fname:
			print self.__doc__
			sys.exit(2)
		self.input_fname = input_fname
		self.phenotype_fname = phenotype_fname
		self.output_fname = output_fname
		self.log_pvalue = int(log_pvalue)
		self.debug = int(debug)
		self.report = int(report)
		"""
	
	def get_phenotype_ls_in_data_matrix_order(self, strain_acc_list, strain_acc_list_phen, category_list_phen, data_type=float):
		"""
		2008-05-21
			strain_acc could be missing in strain_acc_phen2index
			phenotype_value could be 'NA'
		2008-02-14
		"""
		sys.stderr.write("Getting phenotype_ls in order with data_matrix ...")
		strain_acc_phen2index = dict(zip(strain_acc_list_phen, range(len(strain_acc_list_phen)) ) )
		phenotype_ls = []
		for strain_acc in strain_acc_list:
			if strain_acc in strain_acc_phen2index:
				phen_index = strain_acc_phen2index[strain_acc]
				phenotype_value = category_list_phen[phen_index]
				if phenotype_value !='NA':
					phenotype_ls.append(data_type(phenotype_value))
				else:
					phenotype_ls.append(None)
			else:
				phenotype_ls.append(None)
		sys.stderr.write("Done.\n")
		return phenotype_ls
	
	def _kruskal_wallis(self, data_matrix, phenotype_ls, min_data_point=3):
		"""
		2008-05-21
			phenotype_ls could have None value, skip them
			each kw result is wrapped in PassingData
		2008-02-14
		"""
		sys.stderr.write("Doing kruskal wallis test ...\n")
		import rpy
		rpy.r.as_factor.local_mode(rpy.NO_CONVERSION)
		
		if type(data_matrix)==list:
			import numpy
			data_matrix = numpy.array(data_matrix)
		no_of_rows, no_of_cols = data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		for j in range(no_of_cols):
			genotype_ls = data_matrix[:,j]
			non_NA_genotype_ls = []
			non_NA_phenotype_ls = []
			non_NA_genotype2count = {}
			for i in range(no_of_rows):
				if genotype_ls[i]!=0 and phenotype_ls[i]!=None:
					non_NA_genotype = genotype_ls[i]
					non_NA_genotype_ls.append(non_NA_genotype)
					non_NA_phenotype_ls.append(phenotype_ls[i])
					if non_NA_genotype not in non_NA_genotype2count:
						non_NA_genotype2count[non_NA_genotype] = 0
					non_NA_genotype2count[non_NA_genotype] += 1
			count_ls = non_NA_genotype2count.values()
			if len(count_ls)>=2 and count_ls[0]>=min_data_point and count_ls[1]>=min_data_point:
				pvalue = rpy.r.kruskal_test(x=non_NA_phenotype_ls, g=rpy.r.as_factor(non_NA_genotype_ls))['p.value']
				pdata = PassingData(snp_index=j, pvalue=pvalue, count_ls=count_ls)
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		sys.stderr.write("Done.\n")
		return results
	
	def output_kw_results(self, kw_results, SNP_header, output_fname, log_pvalue=0):
		"""
		2008-05-27
			log10
		2008-05-21
			more stuff in kw_results
			each kw result is wrapped in PassingData
		2008-02-14
		"""
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for i in range(len(kw_results)):
			pdata = kw_results[i]
			snp_index = pdata.snp_index
			pvalue = pdata.pvalue
			SNP_name = SNP_header[snp_index]
			SNP_pos_ls = SNP_name.split('_')
			if log_pvalue:
				if pvalue>0:
					pvalue = -math.log10(pvalue)
				else:
					pvalue = 'NA'
			writer.writerow([SNP_pos_ls[0], SNP_pos_ls[1], pvalue] + pdata.count_ls)
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
		data_matrix_phen = numpy.array(data_matrix_phen)
		#2008-05-21 take the 1st phenotype
		phenotype_ls = self.get_phenotype_ls_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, data_matrix_phen[:, self.which_phenotype], data_type=float)
		kw_results = self._kruskal_wallis(data_matrix, phenotype_ls, self.min_data_point)
		self.output_kw_results(kw_results, header[2:], self.output_fname, self.minus_log_pvalue)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Kruskal_Wallis
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()