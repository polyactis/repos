#!/usr/bin/env python
"""

Examples:
	#try kruskal wallis
	Association.py -i /tmp/250K_method_5_after_imputation_noRedundant_051908.tsv -p /Network/Data/250k/finalData_051808/phenotypes.tsv -e -r -o /tmp/250K_method_5_after_imputation_noRedundant_051908.LD.pvalue
	
	#try linear model
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -p ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_lm.tsv -y2
	
	#try emma (linear mixture model) on 1st 7 phenotypes
	Association.py -i ./mnt2/panfs/250k/call_method_17.tsv -p ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_y3.tsv  -y3 -w 0-6

Description:
	class to do association test on SNP data. option 'test_type' decides which test to run.
	
	Input genotype file format is Strain X SNP format (Yu's format, Output by DB_250k2data.py Or Output250KSNPs.py + ConvertBjarniSNPFormat2Yu.py).
	Input phenotype file format is Strain X phenotype format (Output by OutputPhenotype.py). 
	
	It requires a minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis or linear model test.
	
	For kruskal wallis & linear model, it will automatically match strains in two files.
		NO worry for missing/extra data in either input file.
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
from pymodule import read_data, ProcessOptions, PassingData, SNPData, getListOutOfStr
from numpy import linalg

from Kruskal_Wallis import Kruskal_Wallis
import rpy

class Association(Kruskal_Wallis):
	__doc__ = __doc__
	option_default_dict = Kruskal_Wallis.option_default_dict.copy()
	option_default_dict.pop(("which_phenotype", 1, int))
	option_default_dict.update({('which_phenotype_ls', 1, ): ['0', 'w', 1, 'list of index indicating which phenotype, 0-1,3 = 1st,2nd, and 4th phenotype (starting from 3rd column in phenotype_fname) and so on.',]})
	option_default_dict.update({('test_type', 1, int): [1, 'y', 1, 'Which type of test to do. 1:Kruskal_Wallis, 2:linear model(y=xb+e), 3:Emma']})
	def __init__(self, **keywords):
		"""
		2008-11-10
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.which_phenotype_ls = getListOutOfStr(self.which_phenotype_ls, data_type=int)
		self.run_whole_matrix = {1:self._kruskal_wallis_whole_matrix,
								2:self.LM_whole_matrix,
								3:self.Emma_whole_matrix}
		
		self.output_results = {1:self.output_kw_results,
							2:self.output_lm_results,
							3:self.output_lm_results}
	
	def pure_linear_model(cls, non_NA_genotype_ls, non_NA_phenotype_ls):
		"""
		2008-11-10
			split out of linear_model()
		"""
		#rpy.r.as_factor.local_mode(rpy.NO_CONVERSION)
		genotype_matrix = numpy.array(non_NA_genotype_ls, numpy.int)
		genotype_var = numpy.var(genotype_matrix[:,1]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		"""
		#2008-11-10 do linear regression by numpy
		(p, residuals, rank, s) = linalg.lstsq(genotype_matrix, non_NA_phenotype_ls)	#all 4 returned elements are arrays
		F = p[1]*p[1]/((residuals[0]/(n-2))/(genotype_var*(n-1)))	#page 108 of Seber2003LRA. F=r^2(n-2)/(1-r^2). T=sqrt(F) is t-stat with n-2 df. it's used in r.glm(). F-test gives slightly bigger p-value than t-test.
		pvalue = rpy.r.pf(F,1, n-2, lower_tail=rpy.r.FALSE)
		geno_effect_var = genotype_var*p[1]*p[1]*(n-1)
		var_perc = geno_effect_var/(residuals[0]+geno_effect_var)
		#var_perc2 = geno_effect_var/numpy.var(non_NA_phenotype_ls)	#this is also same as var_perc.
		"""
		
		#2008-11-10 do linear regression by R
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		#data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":rpy.r.as_factor(genotype_matrix[:,1])})
		data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":genotype_matrix[:,1]})
		if len(non_NA_phenotype2count)==2:	#binary phenotype, use logistic regression
			lm_result = rpy.r.glm(rpy.r("phenotype~genotype"), data=data_frame, family=rpy.r("binomial"))
		else:
			lm_result = rpy.r.glm(rpy.r("phenotype~genotype"), data=data_frame)	#06-30-05 use formula_list
		rpy.set_default_mode(rpy.BASIC_CONVERSION) #04-07-05
		#04-07-05 r.summary() requires lm_result in NO_CONVERSION state
		summary_stat = rpy.r.summary(lm_result)
		
		#06-30-05	index 0 in summary_stat['coefficients'] is intercept
		coeff_list = []
		coeff_p_value_list = []
		for i in range(len(summary_stat['coefficients'])):
			coeff_list.append(summary_stat['coefficients'][i][0])	#0 is the coefficient
			coeff_p_value_list.append(summary_stat['coefficients'][i][-1])	#-1 is the corresponding p-value
		#06-30-05	fill in other efficients based on bit_string, NOTE i+1
		pvalue = coeff_p_value_list[1]
		residuals = summary_stat['deviance']
		geno_effect_var = genotype_var*coeff_list[1]*coeff_list[1]*(n-1)
		var_perc = geno_effect_var/(residuals+geno_effect_var)
		
		#pvalue = rpy.r.kruskal_test(x=non_NA_phenotype_ls, g=rpy.r.as_factor(non_NA_genotype_ls))['p.value']
		#2008-08-06 try wilcox
		#pvalue = rpy.r.wilcox_test(non_NA_genotype2phenotype_ls[top_2_allele_ls[0]], non_NA_genotype2phenotype_ls[top_2_allele_ls[1]], conf_int=rpy.r.TRUE)['p.value']
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	pure_linear_model = classmethod(pure_linear_model)
	def linear_model(cls, genotype_ls, phenotype_ls, min_data_point=3, snp_index=None, kinship_matrix=None, eig_L=None, run_type=1):
		"""
		
		2008-11-13
			run_type
				1: pure_linear_model
				2: emma
			genotype_ls is already in binary (or integer starting from 0)
		2008-11-10
			similar to _kruskal_wallis() of Kruskal_Wallis
		"""
		non_NA_genotype_ls = []
		non_NA_phenotype_ls = []
		non_NA_genotype2count = {}
		#non_NA_genotype2allele = {}
		non_NA_phenotype2count = {}
		#non_NA_genotype2phenotype_ls = {}	#2008-08-06 try wilcox
		non_NA_index_ls = []
		for i in range(len(genotype_ls)):
			if genotype_ls[i]!=-2 and not numpy.isnan(phenotype_ls[i]):
				non_NA_index_ls.append(i)
				non_NA_genotype = genotype_ls[i]
				non_NA_phenotype = phenotype_ls[i]
				if non_NA_phenotype not in non_NA_phenotype2count:
					non_NA_phenotype2count[non_NA_phenotype] = 0
				non_NA_phenotype2count[non_NA_phenotype] += 1
				non_NA_phenotype_ls.append(non_NA_phenotype)
				
				if non_NA_genotype not in non_NA_genotype2count:
					non_NA_genotype2count[non_NA_genotype] = 0
					#non_NA_genotype2allele[non_NA_genotype] = len(non_NA_genotype2allele)
					#non_NA_genotype2phenotype_ls[non_NA_genotype] = []	#2008-08-06 try wilcox
				#allele = non_NA_genotype2allele[non_NA_genotype]
				allele = genotype_ls[i]
				non_NA_genotype_ls.append([1, allele])
				non_NA_genotype2count[non_NA_genotype] += 1
				#non_NA_genotype2phenotype_ls[non_NA_genotype].append(phenotype_ls[i])	#2008-08-06 try wilcox
		"""
		#2008-08-06 try wilcox
		new_snp_allele2index = returnTop2Allele(non_NA_genotype2count)
		top_2_allele_ls = new_snp_allele2index.keys()
		non_NA_genotype2count = {top_2_allele_ls[0]: non_NA_genotype2count[top_2_allele_ls[0]],
								top_2_allele_ls[1]: non_NA_genotype2count[top_2_allele_ls[1]]}
		"""
		count_ls = non_NA_genotype2count.values()
		n = len(non_NA_phenotype_ls)
		if len(count_ls)>=2 and min(count_ls)>=min_data_point:	#require all alleles meet the min data point requirement
			if run_type==1:
				pdata = cls.pure_linear_model(non_NA_genotype_ls, non_NA_phenotype_ls)
			elif run_type==2:
				if kinship_matrix.shape[0]!=n:	#there is NA and need slicing
					new_kinship_matrix = kinship_matrix[non_NA_index_ls, non_NA_index_ls]
				else:
					new_kinship_matrix = kinship_matrix
				pdata = cls.emma(non_NA_genotype_ls, non_NA_phenotype_ls, new_kinship_matrix, eig_L)
			else:
				sys.stderr.write("run_type=%s not supported.\n"%run_type)
				return None
			pdata.snp_index = snp_index
			pdata.count_ls = count_ls
		else:
			pdata = None
		return pdata
	
	linear_model = classmethod(linear_model)
	
	def LM_whole_matrix(self, data_matrix, phenotype_ls, min_data_point=3):
		"""
		2008-11-10
			adapted from _kruskal_wallis_whole_matrix() of Kruskal_Wallis
		"""
		sys.stderr.write("Association by pure linear model ...\n")
		no_of_rows, no_of_cols = data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		for j in range(no_of_cols):
			genotype_ls = data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, phenotype_ls, min_data_point, snp_index=j)
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		sys.stderr.write("Done.\n")
		return results
	
	def get_kinship_matrix(self, data_matrix):
		"""
		2008-11-11
		"""
		no_of_rows, no_of_cols = data_matrix.shape
		kinship_matrix = numpy.identity(no_of_rows, numpy.float)
		for i in range(no_of_rows):
			for j in range(i+1, no_of_rows):
				identity_vector = data_matrix[i,:]*data_matrix[j,:]+ (1-data_matrix[i,:])*(1-data_matrix[j,:])	#only for binary data_matrix. identity_vector[k]=1 only when data_matrix[i,k]=data_matrix[j,k]
				kinship_matrix[i,j] = sum(identity_vector)/float(no_of_cols)
				kinship_matrix[j,i] = kinship_matrix[i,j]
		return kinship_matrix
	
	def emma(cls, non_NA_genotype_ls, non_NA_phenotype_ls, kinship_matrix, eig_L = None):
		"""
		2008-11-13
			call emma.REMLE()
		"""
		n = len(non_NA_genotype_ls)
		genotype_matrix =  numpy.array(non_NA_genotype_ls, numpy.int)	#non_NA_genotype_ls is a list of (1, X)-tuples
		#numpy.hstack(numpy.ones([n], numpy.int), ...)
		one_marker_rs = rpy.r.emma_REMLE(non_NA_phenotype_ls, genotype_matrix, kinship_matrix, eig_L=eig_L, cal_pvalue=rpy.r.TRUE)
		coeff_list=[beta[0] for beta in one_marker_rs['beta']]
		coeff_p_value_list=[None]*len(coeff_list)
		pdata = PassingData(pvalue=one_marker_rs['pvalue'], var_perc=one_marker_rs['genotype_var_perc'][0][0], \
						coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	emma = classmethod(emma)
	def Emma_whole_matrix(self, data_matrix, phenotype_ls, min_data_point=3):
		"""
		2008-11-13
			iteratively call rpy.r.emma_REMLE() through self.linear_model() in order to get MAF, MAC etc
			stop calling rpy.r.emma_REML_t
		2008-11-11
			assume:
				1. no NA in data_matrix (imputed) but allow NA in phenotype_ls
				2. data_matrix is binary, used in get_kinship_matrix()
			procedure:
			
			1. remove rows that have NA in phenotype_ls
			2. calculate kinship
			3. run emma
		"""
		sys.stderr.write("Association by emma (linear mixture model) ...\n")
		
		#remove non-NA phenotype
		non_phenotype_NA_row_index_ls = []
		non_NA_phenotype_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_phenotype_NA_row_index_ls.append(i)
				non_NA_phenotype_ls.append(phenotype_ls[i])
		
		non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
		new_data_matrix = data_matrix[non_phenotype_NA_row_index_ls,:]
		kinship_matrix = self.get_kinship_matrix(new_data_matrix)
		
		no_of_rows, no_of_cols = data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		
		eig_L = rpy.r.emma_eigen_L(None, kinship_matrix)	#to avoid repetitive computing of eig_L inside emma.REMLE
		for j in range(no_of_cols):
			genotype_ls = data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, phenotype_ls, min_data_point, snp_index=j, \
									kinship_matrix=kinship_matrix, eig_L=eig_L, run_type=2)
			if pdata is not None:
				results.append(pdata)
				real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		"""
		#
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		#rpy.set_default_mode(rpy.NO_CONVERSION)
		non_NA_phenotype_ar.resize([len(non_NA_phenotype_ar),1])	#transform it into 2-D for emma
		results = rpy.r.emma_REML_t(non_NA_phenotype_ar.transpose(), new_data_matrix.transpose(), kinship_matrix)	#in emma, row is marker. col is strain.
		"""
		
		sys.stderr.write("Done.\n")
		return results
	
	def output_emma_results(self, results, SNP_header, output_fname, log_pvalue=0):
		"""
		2008-11-12
			this is written when rpy.r.emma_REML_t() is used in Emma_whole_matrix().
		"""
		#collect results
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		
		no_of_snps = len(results['ps'])
		counter = 0
		real_counter = 0
		for i in range(no_of_snps):
			SNP_name = SNP_header[i]
			SNP_pos_ls = SNP_name.split('_')
			
			pvalue = results['ps'][i][0]
			var_perc = results['genotype_var_perc'][i][0]
			coeff_list = [results['beta0_est'][i][0], results['beta1_est'][i][0]]
			
			#pdata = PassingData(snp_index=i, pvalue=pvalue, count_ls=[0.5,20], var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
			MAF = 0.5
			MAC = 20
			writer.writerow([SNP_pos_ls[0], SNP_pos_ls[1], pvalue, MAF, MAC, var_perc] + coeff_list)
			
			real_counter += 1
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		del writer
		sys.stderr.write("Done.\n")
		
	def output_lm_results(self, results, SNP_header, output_fname, log_pvalue=0):
		"""
		2008-11-11
			adapted from output_kw_results() of Kruskal_Wallis.py
		2008-05-27
			log10
		2008-05-21
			more stuff in results
			each kw result is wrapped in PassingData
		2008-02-14
		"""
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for i in range(len(results)):
			pdata = results[i]
			snp_index = pdata.snp_index
			pvalue = pdata.pvalue
			SNP_name = SNP_header[snp_index]
			SNP_pos_ls = SNP_name.split('_')
			if log_pvalue:
				if pvalue>0:
					pvalue = -math.log10(pvalue)
				else:
					pvalue = 'NA'
			MAC = min(pdata.count_ls)
			MAF = float(MAC)/sum(pdata.count_ls)
			coeff_and_pvalue_ls = []
			for j in range(len(pdata.coeff_list)):
				coeff = pdata.coeff_list[j]
				coeff_pvalue = pdata.coeff_p_value_list[j]
				coeff_and_pvalue = '%s'%coeff
				if coeff_pvalue:
					coeff_and_pvalue += ':%s'%coeff_pvalue
				coeff_and_pvalue_ls.append(coeff_and_pvalue)
			writer.writerow([SNP_pos_ls[0], SNP_pos_ls[1], pvalue, MAF, MAC, pdata.var_perc] + coeff_and_pvalue_ls)
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
		#if type(data_matrix)==list:
		#	data_matrix = numpy.array(data_matrix)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix)
		
		#convert SNP matrix into binary
		if self.test_type!=1:	#kruskal wallis doesn't need binary data matrix
			newSnpData, allele2index_ls = snpData.convertSNPAllele2Index(self.report)	#0 (NA) or -2 (untouched) is all converted to -2 as 0 is used to denote allele
		else:
			newSnpData = snpData
		
		data_matrix_phen = self.get_phenotype_matrix_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, data_matrix_phen)
		for which_phenotype in self.which_phenotype_ls:
			phenotype_name = header_phen[2+which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			output_fname='%s_pheno_%s.tsv'%(os.path.splitext(self.output_fname)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			results = self.run_whole_matrix[self.test_type](newSnpData.data_matrix, data_matrix_phen[:, which_phenotype], self.min_data_point)
			self.output_results[self.test_type](results, header[2:], output_fname, self.minus_log_pvalue)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Association
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
