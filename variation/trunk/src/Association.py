#!/usr/bin/env python
"""

Examples:
	#try kruskal wallis
	Association.py -i /tmp/250K_method_5_after_imputation_noRedundant_051908.tsv -p /Network/Data/250k/finalData_051808/phenotypes.tsv -e -r -o /tmp/250K_method_5_after_imputation_noRedundant_051908.LD.pvalue
	
	#try linear model
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -p ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_lm.tsv -y2
	
	#try emma (linear mixture model) on 1st 7 phenotypes
	Association.py -i ./mnt2/panfs/250k/call_method_17.tsv -p ./banyan_fs/tmp/phenotype.tsv -o /tmp/call_method_17_y3.tsv  -y3 -w 1-7

	#linear model with principal components 0 to 9, phenotype from 1 to 7
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -y4 -o /Network/Data/250k/tmp-yh/eigenstrat//call_method_17_lm_with_pc0_9 -W 0-9 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -r -w 1-7
	
	#linear model with PCs 0 to 1, phenotype from 1 to 5. the PCs are calculated on the fly according to the snp input file.
	Association.py -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_chr4_100000_700000.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_chr4_100000_700000_y4_pc0_1 -W 0-1 -y 4 -w 1-5 -r
	
	#y = SNP + environment + noise, for phenotype 1 & 2
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y5.tsv -y5 -w 1,2 -r
	
	#y = SNP + environment + PC1 + PC2 + noise, for phenotype 1 & 2
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm_with_PC12/call_method_17_y5.tsv -W 0-1 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -y5 -w 1,2 -r
	
	#y = SNP + environment + SNP X environ + noise, for phenotype 1 & 2
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6.tsv -y6 -w 1,2 -r
	
	#y = SNP + environment + SNP X environ + PC1 + PC2 + noise, for phenotype 1 & 2
	Association.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -p /Network/Data/250k/tmp-yh/phenotype.tsv -o /Network/Data/250k/tmp-yh/association_results/lm_with_PC12/call_method_17_y6.tsv -W 0-1 -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -y6 -w 1,2 -r
	
	# 2010-2-1 EMMAX
	~/script/variation/src/Association.py -i /Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv -p /Network/Data/250k/tmp-yh//phenotype.tsv -o /tmp/call_method_17_y8.tsv  -y8 -w 1
	
Description:
	class to do association test on SNP data. option 'test_type' decides which test to run.
	
	Input genotype file format is Strain X SNP format (Yu's format, Output by DB_250k2data.py Or Output250KSNPs.py + ConvertBjarniSNPFormat2Yu.py).
	Input phenotype file format is Strain X phenotype format (Output by OutputPhenotype.py). 
	
	It requires a minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis or linear model test.
	
	For all methods, it will automatically match strains in two files.
		NO worry for missing/extra data in either input file.
	
	All methods iterate through phenotypes given by '-w' except that Method "5: LM two phenotypes with PCs" takes two phenotypes from '-w'.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback
from pymodule import read_data, ProcessOptions, PassingData, SNPData, getListOutOfStr
from numpy import linalg

from Kruskal_Wallis import Kruskal_Wallis
import rpy
from PlotGroupOfSNPs import PlotGroupOfSNPs
from sets import Set
#from DrawEcotypeOnMap import DrawEcotypeOnMap

class Association(Kruskal_Wallis):
	debug = 0
	report = 0
	__doc__ = __doc__
	option_default_dict = Kruskal_Wallis.option_default_dict.copy()
	option_default_dict.pop(("which_phenotype", 1, int))
	option_default_dict.update({('phenotype_method_id_ls', 0, ): ['1', 'w', 1, 'which phenotypes to work on. a comma-dash-separated list phenotype_method ids in the phenotype file. Check db Table phenotype_method. \
		if not available, take all phenotypes in the phenotype_fname.',]})
	option_default_dict.update({('test_type', 1, int): [1, 'y', 1, 'Which type of test to do. 1:Kruskal_Wallis, 2:linear model(y=xb+e), 3:Emma, 4:LM with PCs, 5: LM two phenotypes with PCs, \
		6: LM two phenotypes with PCs, GeneXEnvironment Interaction, 7: Emma for genotype matrix without NA (no MAC and MAF output), 8: EMMAX (variance matrix is estimated once)']})
	option_default_dict.update({('eigen_vector_fname', 0, ): [None, 'f', 1, 'eigen vector file with PCs outputted by smartpca.perl from EIGENSOFT', ]})
	option_default_dict.update({('which_PC_index_ls', 0, ): [None, 'W', 1, 'list of indices indicating which PC(s) from eigen_vector_fname should be used. format: 0,1-3', ]})
	def __init__(self, **keywords):
		"""
		2008-11-10
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.phenotype_method_id_ls:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
		
		self.which_PC_index_ls = getListOutOfStr(self.which_PC_index_ls, data_type=int)
		
		self.run_whole_matrix = {1:self._kruskal_wallis_whole_matrix,
								2:self.LM_whole_matrix,
								3:self.Emma_whole_matrix,
								4:self.LM_with_PCs_whole_matrix,
								5:self.LM_with_PCs_whole_matrix,
								6:self.LM_with_PCs_whole_matrix,
								7:self.Emma_whole_matrixForNoNAGenotypeMatrix,
								8:self.EMMAX}
		
		self.output_results = {1:self.output_kw_results,
							2:self.output_lm_results,
							3:self.output_lm_results,
							7:self.output_emma_results}
	
	@classmethod
	def getEigenValueFromFile(cls, eigen_value_fname):
		"""
		2008-12-03
		"""
		sys.stderr.write("Getting eigen values from %s ..."%eigen_value_fname)
		inf = open(eigen_value_fname)
		eigen_value_ls = []
		for line in inf:
			eigen_value = float(line.strip())
			eigen_value_ls.append(eigen_value)
		sys.stderr.write("Done.\n")
		return eigen_value_ls
	
	
	@classmethod
	def getPCFromFile(cls, eigen_vector_fname):
		"""
		2009-2-2
			add id_ls in the returning data
			wrap the data under PC_data
		2008-12-03
			
		"""
		sys.stderr.write("Getting principal components from %s ..."%eigen_vector_fname)
		inf = open(eigen_vector_fname)
		inf.next()	#skip first row (corresponding eigen values)
		PC_matrix = []
		id_ls = []
		for line in inf:
			row = line.split()
			PCs_for_one_entity = map(float, row[1:-1])	#1st entry in row is individual label. last entry is Case/Control.
			PC_matrix.append(PCs_for_one_entity)
			id_ls.append(row[0])
		PC_matrix = numpy.array(PC_matrix)
		PC_data = PassingData(PC_matrix=PC_matrix, id_ls=id_ls)
		sys.stderr.write("Done.\n")
		return PC_data
	
	@classmethod
	def multi_linear_model(cls, genotype_matrix, phenotype_ls, min_data_point=3, snp_index=None, kinship_matrix=None, eig_L=None, run_type=1):
		"""
		
		2008-11-13
			run_type
				1: pure_linear_model
				2: emma
			genotype_ls is already in binary (or integer starting from 0)
		2008-11-10
			similar to _kruskal_wallis() of Kruskal_Wallis
		"""
		non_NA_genotype_matrix = []
		non_NA_phenotype_ls = []
		non_NA_genotype2count = {}
		#non_NA_genotype2phenotype_ls = {}	#2008-08-06 try wilcox
		non_NA_index_ls = []
		for i in range(len(genotype_matrix)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_index_ls.append(i)
				non_NA_phenotype = phenotype_ls[i]
				non_NA_phenotype_ls.append(non_NA_phenotype)
				
				non_NA_genotype_matrix.append(genotype_matrix[i,:])
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
	
	@classmethod
	def pure_linear_model(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None):
		"""
		2009-12-18
			non_NA_phenotype2count is not used.
		2008-11-10
			split out of linear_model()
		"""
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls, add_intercept=True)
		no_of_rows, no_of_cols = genotype_matrix.shape
		
		genotype_matrix2 = numpy.transpose(genotype_matrix)
		D = linalg.inv(numpy.inner(genotype_matrix2, genotype_matrix2))	#2nd matrix's shape is opposite to real matrix algebra.

		(p, residuals, rank, s) = linalg.lstsq(genotype_matrix, non_NA_phenotype_ls)	#all 4 returned elements are arrays
		#coeff_list = [p[0]]	#put the intercept beta there first
		#coeff_p_value_list = [-1]
		coeff_list = []
		coeff_p_value_list = []
		for i in range(len(p)):
			coeff_list.append(p[i])
			F = p[i]*p[i]/((residuals[0]/(no_of_rows-len(p)))*D[i,i])	#page 106 of Seber2003LRA
				#page 108 of Seber2003LRA. F=r^2(n-2)/(1-r^2). 
				#T=sqrt(F) is t-stat with n-2 df. it's used in r.glm(). F-test gives slightly bigger p-value than t-test.
			pvalue = rpy.r.pf(F,1, no_of_rows-len(p), lower_tail=rpy.r.FALSE)
			coeff_p_value_list.append(pvalue)
		pvalue = coeff_p_value_list[1]	#this is the pvalue to return, corresponding to the genotype vector
		genotype_var = numpy.var(genotype_matrix[:,1]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		geno_effect_var = genotype_var*p[1]*p[1]*(no_of_rows-1)
		var_perc = geno_effect_var/(residuals[0]+geno_effect_var)
		#var_perc2 = geno_effect_var/numpy.var(non_NA_phenotype_ls)	#this is also same as var_perc.
		
		
		#pvalue = rpy.r.kruskal_test(x=non_NA_phenotype_ls, g=rpy.r.as_factor(non_NA_genotype_ls))['p.value']
		#2008-08-06 try wilcox
		#pvalue = rpy.r.wilcox_test(non_NA_genotype2phenotype_ls[top_2_allele_ls[0]], non_NA_genotype2phenotype_ls[top_2_allele_ls[1]], conf_int=rpy.r.TRUE)['p.value']
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	@classmethod
	def pure_linear_model_via_R(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None):
		"""
		2010-2-25
			use createDesignMatrix() to generate a design matrix
		2009-8-28
			split out of pure_linear_model(). same functionality as pure_linear_model(), but invoke R to run regression.
		"""
		
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls)
		#2008-11-10 do linear regression by R
		genotype_var = numpy.var(genotype_matrix[:,0]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		#data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":rpy.r.as_factor(genotype_matrix[:,1])})
		formula_list = []
		data_frame_dict = {"phenotype":non_NA_phenotype_ls}
		for i in range(genotype_matrix.shape[1]):
			var_name = 'genotype%s'%i
			formula_list.append(var_name)
			data_frame_dict.update({var_name: genotype_matrix[:,i]})
		data_frame = rpy.r.as_data_frame(data_frame_dict)
		formula = 'phenotype~%s'%'+'.join(formula_list)
		
		if non_NA_phenotype2count and len(non_NA_phenotype2count)==2:	#binary phenotype, use logistic regression
			lm_result = rpy.r.glm(rpy.r(formula), data=data_frame, family=rpy.r("binomial"))
		else:
			lm_result = rpy.r.glm(rpy.r(formula), data=data_frame)
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
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
		geno_effect_var = genotype_var*coeff_list[1]*coeff_list[1]*(no_of_rows-1)
		var_perc = geno_effect_var/(residuals+geno_effect_var)
		
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	@classmethod
	def createDesignMatrix(cls, genotype_ls, add_intercept=False,n=None):
		"""
		2009-12-23
			argument add_intercept determines whether a constant vector would be added to the design matrix or not
		"""
		if genotype_ls is None:
			if add_intercept:
				if n is None:
					sys.stderr.write("n (number of accessions) must be specified to create a design matrix with genotype_ls=None.\n")
					sys.exit(3)
				design_matrix = (numpy.ones([n ,1], numpy.int))
		else:
			#rpy.r.as_factor.local_mode(rpy.NO_CONVERSION)
			if not isinstance(genotype_ls, numpy.ndarray):
				design_matrix = numpy.array(genotype_ls)
			else:
				design_matrix = genotype_ls
			if len(design_matrix.shape)==1:	#transform into 2D array
				design_matrix = numpy.resize(design_matrix, [len(design_matrix), 1])				
			no_of_rows, no_of_cols = design_matrix.shape
			
			
			if add_intercept:
				#2008-11-10 do linear regression by numpy
				design_matrix = numpy.hstack((numpy.ones([len(design_matrix),1], numpy.int), design_matrix))	#the design variable for intercept has to be included.
			
		return design_matrix
	
	@classmethod
	def gls_via_R(cls, non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count=None, variance_matrix=None):
		"""
		2009-12-23
			general least square model via calling equivalent function in R.
			
		"""
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls)  # no need to add a constant vector.
		
		if hasattr(cls, 'corStruct'):
			corStruct = cls.corStruct
		else:
			if variance_matrix is not None:
				corStruct = cls.generateCorStructForGLSFromVarianceMatrix(variance_matrix)
				setattr(cls, "corStruct", corStruct)
			else:
				corStruct = None
		#2008-11-10 do linear regression by R
		genotype_var = numpy.var(genotype_matrix[:,0]) 	#2008-11-10 var=\sum(x_i-\bar{x})^2/(n-1)
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		rpy.r.library("nlme")
		
		#data_frame = rpy.r.as_data_frame({"phenotype":non_NA_phenotype_ls, "genotype":rpy.r.as_factor(genotype_matrix[:,1])})
		formula_list = []
		data_frame_dict = {"phenotype":non_NA_phenotype_ls}
		for i in range(genotype_matrix.shape[1]):
			var_name = 'genotype%s'%i
			formula_list.append(var_name)
			data_frame_dict.update({var_name: genotype_matrix[:,i]})
		data_frame = rpy.r.as_data_frame(data_frame_dict)
		formula = 'phenotype~%s'%'+'.join(formula_list)
		
		lm_result = rpy.r.gls(rpy.r(formula), data=data_frame, correlation=corStruct)
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
		#04-07-05 r.summary() requires lm_result in NO_CONVERSION state
		summary_stat = rpy.r.summary(lm_result)
		
		rpy.set_default_mode(rpy.NO_CONVERSION)
		summary_stat1 = rpy.r.summary(lm_result)
		
		rpy.set_default_mode(rpy.VECTOR_CONVERSION)
		summary_stat2 = rpy.r.summary(lm_result)
		
		rpy.set_default_mode(rpy.TOP_CONVERSION)
		summary_stat3 = rpy.r.summary(lm_result)
		
		#06-30-05	index 0 in summary_stat['coefficients'] is intercept
		coeff_list = []
		coeff_p_value_list = []
		for i in range(len(summary_stat['coefficients'])):
			coeff_list.append(summary_stat['coefficients'][i][0])	#0 is the coefficient
			coeff_p_value_list.append(summary_stat['coefficients'][i][-1])	#-1 is the corresponding p-value
		#06-30-05	fill in other efficients based on bit_string, NOTE i+1
		pvalue = coeff_p_value_list[1]
		residuals = summary_stat['deviance']
		geno_effect_var = genotype_var*coeff_list[1]*coeff_list[1]*(no_of_rows-1)
		var_perc = geno_effect_var/(residuals+geno_effect_var)
		
		pdata = PassingData(pvalue=pvalue, var_perc=var_perc, coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list)
		return pdata
	
	
	@classmethod
	def linear_model(cls, genotype_ls, phenotype_ls, min_data_point=3, snp_index=None, kinship_matrix=None, eig_L=None, \
					run_type=1, counting_and_NA_checking=True, variance_matrix=None):
		"""
		2009-12-23
			add argument variance_matrix for gls()
		2009-12-18
			add argument counting_and_NA_checking, which if enabled, this function counts the number of accessions with either allele
				(, which allows calculation of MAF and makes sure that MAF wasn't too small)  
				and exclude the NA genotypes or phenotypes.
			If it's disabled, pure_linear_model_via_R() would always call glm() with no family specification. Others are same.
		2009-12-23
			run_type
				1: pure_linear_model
				2: emma
				3: pure_linear_model via R
				4: gls with specified variance matrix
		2009-5-15
			report full error information if cls.debug is set.
		2009-2-11
			simply report the snp_index when numpy.linalg.linalg.LinAlgError is caught
		2008-12-04
			genotype_ls could be either 1D or 2D matrix
			if it's 2D, only the 1st column is genotype. 2nd or further are PCs (so far).
		2008-11-13
			run_type
				1: pure_linear_model
				2: emma
			genotype_ls is already in binary (or integer starting from 0)
		2008-11-10
			similar to _kruskal_wallis() of Kruskal_Wallis
		"""
		if counting_and_NA_checking:
			non_NA_genotype_ls = []
			non_NA_phenotype_ls = []
			non_NA_genotype2count = {}
			#non_NA_genotype2allele = {}
			non_NA_phenotype2count = {}
			#non_NA_genotype2phenotype_ls = {}	#2008-08-06 try wilcox
			non_NA_index_ls = []
			for i in range(len(genotype_ls)):
				if isinstance(genotype_ls[i], numpy.ndarray):
					genotype = genotype_ls[i][0]
				else:
					genotype = genotype_ls[i]
				if genotype!=-2 and not numpy.isnan(genotype) and not numpy.isnan(phenotype_ls[i]):
					non_NA_index_ls.append(i)
					non_NA_genotype = genotype
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
					if not isinstance(allele, numpy.ndarray):	#make non_NA_genotype_ls 2D
						allele = [allele]
					non_NA_genotype_ls.append(allele)
					non_NA_genotype2count[non_NA_genotype] += 1
					#non_NA_genotype2phenotype_ls[non_NA_genotype].append(phenotype_ls[i])	#2008-08-06 try wilcox
		else:
			non_NA_genotype_ls = genotype_ls
			non_NA_phenotype_ls = phenotype_ls
			non_NA_genotype2count = {0:min_data_point, 1:min_data_point}	# fake one and allow it to pass the condition below
			non_NA_phenotype2count = {}	# based on its length, pure_linear_model_via_R would decide whether to call glm() with 'binomial'
										# if phenotypes are binary or not. So if this is empty, envoke glm() with no family specification.
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
			#try:
			if run_type==1:
				pdata = cls.pure_linear_model(non_NA_genotype_ls, non_NA_phenotype_ls)
			elif run_type==3:
				pdata = cls.pure_linear_model_via_R(non_NA_genotype_ls, non_NA_phenotype_ls, non_NA_phenotype2count)
			elif run_type==2:
				if kinship_matrix.shape[0]!=n:	#there is NA and need slicing
					new_kinship_matrix = kinship_matrix[non_NA_index_ls, non_NA_index_ls]
					eig_L = rpy.r.emma_eigen_L(None, new_kinship_matrix)	#2008-1-5 generate new eig_L
				else:
					new_kinship_matrix = kinship_matrix
				pdata = cls.emma(non_NA_genotype_ls, non_NA_phenotype_ls, new_kinship_matrix, eig_L)
			elif run_type==4:
				pdata = cls.gls_via_R(non_NA_genotype_ls, non_NA_phenotype_ls, variance_matrix=variance_matrix)
			else:
				sys.stderr.write("run_type=%s not supported.\n"%run_type)
				return None
			pdata.snp_index = snp_index
			pdata.count_ls = count_ls
			"""
			except numpy.linalg.linalg.LinAlgError:
				sys.stderr.write("Except while running pure_linear_model on snp_index=%s:\n"%(snp_index))
				sys.stderr.write('\t%s.\n'%repr(sys.exc_info()))
				pdata = None
			except:
				if cls.debug:	# 2009-5-15
					sys.stderr.write("Except while running pure_linear_model on snp_index=%s with non_NA_genotype_ls=%s, \
					non_NA_phenotype_ls=%s, non_NA_phenotype2count=%s.\n"%(snp_index, repr(non_NA_genotype_ls), \
																			repr(non_NA_phenotype_ls), repr(non_NA_phenotype2count)))
					traceback.print_exc()
					sys.stderr.write('%s.\n'%repr(sys.exc_info()))
				else:	#2009-4-5 simpilify output if not debug
					sys.stderr.write("Except (%s) while running pure_linear_model on snp_index=%s.\n"%(repr(sys.exc_info()), snp_index))
				pdata = None
			"""
		else:
			pdata = None
		return pdata
	
	
	def LM_whole_matrix(self, data_matrix, phenotype_ls, min_data_point=3, **keywords):
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
	
	def LM_with_PCs_whole_matrix(self, data_matrix, phenotype_ls, min_data_point=3, **keywords):
		"""
		2008-01-04
			add code to deal with environment_matrix
		2008-12-04
		"""
		sys.stderr.write("Association by linear model with principle components ...\n")
		no_of_rows, no_of_cols = data_matrix.shape
		which_PC_index_ls = keywords.get('which_PC_index_ls') or [0]
		PC_matrix = keywords.get('PC_matrix')
		environment_matrix = keywords.get('environment_matrix')
		gene_environ_interaction = keywords.get('gene_environ_interaction')
		
		results = []
		counter = 0
		real_counter = 0
		if PC_matrix is not None:
			sub_PC_matrix = PC_matrix[:, which_PC_index_ls]
		else:
			sub_PC_matrix = None
		for j in range(no_of_cols):
			genotype_ls = data_matrix[:,j]
			genotype_ls = numpy.resize(genotype_ls, [len(genotype_ls),1])	#make it 2D , so able to hstack with sub_PC_matrix
			
			if gene_environ_interaction:
				gene_environ_matrix = genotype_ls*environment_matrix
			
			if environment_matrix is not None:
				genotype_ls = numpy.hstack((genotype_ls, environment_matrix))
			
			if gene_environ_interaction:
				genotype_ls = numpy.hstack((genotype_ls, gene_environ_matrix))
			if sub_PC_matrix is not None:
				genotype_ls = numpy.hstack((genotype_ls, sub_PC_matrix))
			
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
	
	@classmethod
	def get_kinship_matrix(cls, data_matrix):
		"""
		2009-2-9
			make it classmethod
		2008-11-11
			only for binary data_matrix. identity_vector[k]=1 only when data_matrix[i,k]=data_matrix[j,k]
		"""
		no_of_rows, no_of_cols = data_matrix.shape
		kinship_matrix = numpy.identity(no_of_rows, numpy.float)
		for i in range(no_of_rows):
			for j in range(i+1, no_of_rows):
				identity_vector = data_matrix[i,:]*data_matrix[j,:]+ (1-data_matrix[i,:])*(1-data_matrix[j,:])	#only for binary data_matrix. identity_vector[k]=1 only when data_matrix[i,k]=data_matrix[j,k]
				kinship_matrix[i,j] = sum(identity_vector)/float(no_of_cols)
				kinship_matrix[j,i] = kinship_matrix[i,j]
		return kinship_matrix
	
	@classmethod
	def emma(cls, non_NA_genotype_ls=None, non_NA_phenotype_ls=None, kinship_matrix=None, eig_L = None):
		"""
		2009-12-19
			non_NA_genotype_ls could be None
		2009-8-26
			update coeff_p_value_list with pvalues returned from EMMA
		2008-11-13
			call emma.REMLE()
		"""
		if non_NA_genotype_ls is None:
			if non_NA_phenotype_ls is not None:
				n = len(non_NA_phenotype_ls)
			elif kinship_matrix is not None:
				n = kinship_matrix.shape[0]
			else:
				n = None
		else:
			n = len(non_NA_genotype_ls)
		genotype_matrix = cls.createDesignMatrix(non_NA_genotype_ls, add_intercept=True, n=n)
		one_marker_rs = rpy.r.emma_REMLE(non_NA_phenotype_ls, genotype_matrix, kinship_matrix, eig_L=eig_L, cal_pvalue=rpy.r.TRUE)
		coeff_list = [beta[0] for beta in one_marker_rs['beta']]
		if non_NA_genotype_ls is None:    # no covariate used, so no pvalues
			coeff_p_value_list=[0.]
		else:
			coeff_p_value_list = [pvalue for pvalue in one_marker_rs['pvalues']]
		
		"""# 2009-3-24 temporary, get random_effect+residual
		coeff_ar = numpy.array(coeff_list)
		random_effect_and_residual_ls = list(numpy.array(non_NA_phenotype_ls)-numpy.inner(genotype_matrix, coeff_ar))
		coeff_list.extend(random_effect_and_residual_ls)
		"""
		pdata = PassingData(pvalue=one_marker_rs['pvalue'], var_perc=one_marker_rs['genotype_var_perc'][0][0], \
						coeff_list=coeff_list, coeff_p_value_list=coeff_p_value_list, ve=one_marker_rs['ve'], \
						vg=one_marker_rs['vg'])	# vg is the variance for the random effect (multi-small-effect-gene effect), ve is the residual variance. 
		return pdata
	
	def Emma_whole_matrix(self, data_matrix, phenotype_ls, min_data_point=3, **keywords):
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
				
		#remove NA phenotype and genotypes from the corresponding accessions
		non_NA_phenotype_row_index_ls = []
		non_NA_phenotype_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_phenotype_row_index_ls.append(i)
				non_NA_phenotype_ls.append(phenotype_ls[i])
		
		non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
		new_data_matrix = data_matrix[non_NA_phenotype_row_index_ls,:]
		kinship_matrix = self.get_kinship_matrix(new_data_matrix)
		
		no_of_rows, no_of_cols = new_data_matrix.shape
		results = []
		counter = 0
		real_counter = 0
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		
		eig_L = rpy.r.emma_eigen_L(None, kinship_matrix)	#to avoid repeating the computation of eig_L inside emma.REMLE
		for j in range(no_of_cols):
			genotype_ls = new_data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, non_NA_phenotype_ls, min_data_point, snp_index=j, \
									kinship_matrix=kinship_matrix, eig_L=eig_L, run_type=2)
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
	
	def Emma_whole_matrixForNoNAGenotypeMatrix(self, data_matrix, phenotype_ls, min_data_point=3, **keywords):
		"""
		2009-3-18
			carved out of old Emma_whole_matrix() for data_matrix that has no NA in it.
			no MAC and MAF information.
		"""
		sys.stderr.write("Association by monolithic-emma (linear mixture model) ...\n")
		
		#remove NA phenotype
		non_NA_phenotype_row_index_ls = []
		non_NA_phenotype_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_phenotype_row_index_ls.append(i)
				non_NA_phenotype_ls.append(phenotype_ls[i])
		
		non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
		new_data_matrix = data_matrix[non_NA_phenotype_row_index_ls,:]
		kinship_matrix = self.get_kinship_matrix(new_data_matrix)
		
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		#rpy.set_default_mode(rpy.NO_CONVERSION)
		non_NA_phenotype_ar.resize([len(non_NA_phenotype_ar),1])	#transform it into 2-D for emma
		results = rpy.r.emma_REML_t(non_NA_phenotype_ar.transpose(), new_data_matrix.transpose(), kinship_matrix)	#in emma, row is marker. col is strain.
		
		sys.stderr.write("Done.\n")
		return results
	
	@classmethod
	def generateCorStructForGLSFromVarianceMatrix(cls, variance_matrix):
		"""
		2009-12-23
			generate the corStruct for gls()
		"""
		sys.stderr.write("Generating corStruct for gls() from variance_matrix ...")
		
		rpy.set_default_mode(rpy.NO_CONVERSION) #04-07-05
		rpy.r.library("nlme")
		
		# bring the lower-triangle of variance_matrix into a list, row by row
		no_of_rows, no_of_cols = variance_matrix.shape
		lower_triangle_cor_vector = []
		for i in range(1, no_of_rows):
			for j in range(i):
				lower_triangle_cor_vector.append(variance_matrix[i][j]/math.sqrt(variance_matrix[i][i]*variance_matrix[j][j]))
		
		csSymm = rpy.r.corSymm(value=lower_triangle_cor_vector)
		
		data_frame = rpy.r.as_data_frame({"fakedata":[1]*no_of_rows})
		csSymm = rpy.r.Initialize(csSymm, data=data_frame)
		rpy.set_default_mode(rpy.BASIC_CONVERSION)
		sys.stderr.write("Done.\n")
		return csSymm
	
	def EMMAX(self, data_matrix, phenotype_ls, min_data_point=3, **keywords):
		"""
		2009-12-18
			fast-version EMMA
				estimate variance first by calling EMMA with only intercept. then call regular glm() with the new variance matrix.
		"""
		sys.stderr.write("Association by EMMAX (fast-version EMMA) ...\n")
		
		#remove NA phenotype and genotypes from the corresponding accessions
		non_NA_phenotype_row_index_ls = []
		non_NA_phenotype_ls = []
		for i in range(len(phenotype_ls)):
			if not numpy.isnan(phenotype_ls[i]):
				non_NA_phenotype_row_index_ls.append(i)
				non_NA_phenotype_ls.append(phenotype_ls[i])
		
		non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
		new_data_matrix = data_matrix[non_NA_phenotype_row_index_ls,:]
		kinship_matrix = self.get_kinship_matrix(new_data_matrix)
		
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		# run EMMA here to get vg & ve, scalars for the two variance matrices (random effect, residual) 
		one_emma_rs = self.emma(non_NA_genotype_ls=None, non_NA_phenotype_ls=non_NA_phenotype_ls, kinship_matrix=kinship_matrix, eig_L=None)
		vg = one_emma_rs.vg
		ve = one_emma_rs.ve
		
		no_of_rows, no_of_cols = new_data_matrix.shape
		
		variance_matrix = vg*kinship_matrix + ve*numpy.identity(no_of_rows)
		
		results = []
		counter = 0
		real_counter = 0
		
		for j in range(no_of_cols):
			genotype_ls = new_data_matrix[:,j]
			pdata = self.linear_model(genotype_ls, non_NA_phenotype_ls, min_data_point, snp_index=j, \
									kinship_matrix=None, eig_L=None, run_type=4, counting_and_NA_checking=False,\
									variance_matrix=variance_matrix)
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
	
	@classmethod
	def output_emma_results(cls, results, SNP_header, output_fname, log_pvalue=0):
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
			if cls.report and counter%2000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		if cls.report:
			sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
		
		del writer
		sys.stderr.write("Done.\n")
	
	@classmethod
	def output_lm_results(cls, results, SNP_header, output_fname, log_pvalue=0):
		"""
		2009-9-19 fix a bug
			if coeff_pvalue is not None:	#2009-9-19 different from "if coeff_pvalue:", which would skip 0.0
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
				if coeff_pvalue is not None:	#2009-9-19 different from "if coeff_pvalue:", which would skip 0.0
					coeff_and_pvalue += ':%s'%coeff_pvalue
				coeff_and_pvalue_ls.append(coeff_and_pvalue)
			writer.writerow([SNP_pos_ls[0], SNP_pos_ls[1], pvalue, MAF, MAC, pdata.var_perc] + coeff_and_pvalue_ls)
		del writer
		sys.stderr.write("Done.\n")
	
	@classmethod
	def removeUnPhenotypedSNPData(clf, snpData, header_phen, strain_acc_list_phen, data_matrix_phen, phenotype_method_id_ls):
		"""
		2010-2-25
			remove un-phenotyped ecotypes from the SNP data in order to keep the snp dataset small 
		"""
		sys.stderr.write("Removing un-phenotyped ecotypes from the SNP data ...")
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, data_matrix=data_matrix_phen)
		if phenotype_method_id_ls:
			which_phenotype_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(phenData, Set(phenotype_method_id_ls))
		else:	#if not available, take all phenotypes
			which_phenotype_ls = range(len(phenData.col_id_ls))
		
		phenotyped_ecotype_id_set = set()
		for i in range(len(phenData.row_id_ls)):
			ecotype_id = phenData.row_id_ls[i]
			keep_this_ecotype = False
			for col_index in which_phenotype_ls:
				if phenData.data_matrix[i][col_index]!='NA':	# 2010-2-25 phenotype values are in raw string.
					keep_this_ecotype = True
					break
			if keep_this_ecotype:
				phenotyped_ecotype_id_set.add(ecotype_id)
		
		row_ids_to_be_kept = set()	# 2010-2-21
		no_of_ecotypes_in_total = len(snpData.row_id_ls)
		for row_id in snpData.row_id_ls:
			ecotype_id = row_id[0]	#1st column is ecotype_id, 2nd is array id
			if ecotype_id in phenotyped_ecotype_id_set:
				row_ids_to_be_kept.add(row_id)
		snpData = SNPData.keepRowsByRowID(snpData, row_ids_to_be_kept)
		no_of_removed = no_of_ecotypes_in_total - len(row_ids_to_be_kept)
		sys.stderr.write("%s removed. Done.\n"%(no_of_removed))
		return snpData
	
	@classmethod
	def readInData(cls, phenotype_fname, input_fname, eigen_vector_fname, phenotype_method_id_ls, test_type=1, report=0):
		"""
		2010-2-25
			call removeUnPhenotypedSNPData() to shrink the snp dataset by removing un-phenotyped ecotypes
		2009-3-20
			refactored out of run(), easy for MpiAssociation.py to call
		"""
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix)
		
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
		snpData = cls.removeUnPhenotypedSNPData(snpData, header_phen, strain_acc_list_phen, data_matrix_phen, phenotype_method_id_ls)
		
		newSnpData, allele2index_ls = snpData.convertSNPAllele2Index(report)	#0 (NA) or -2 (untouched) is all converted to -2 as 0 is used to denote allele
		newSnpData.header = snpData.header
		
		data_matrix_phen = cls.get_phenotype_matrix_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, data_matrix_phen)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)
		
		if eigen_vector_fname:
			PC_data = cls.getPCFromFile(eigen_vector_fname)
			PC_matrix = PC_data.PC_matrix
		else:
			if test_type==4:	#eigen_vector_fname not given for this test_type. calcualte PCs.
				import pca_module
				T, P, explained_var = pca_module.PCA_svd(newSnpData.data_matrix, standardize=False)
				PC_matrix = T
			else:
				PC_matrix = None
		
		del snpData
		if phenotype_method_id_ls:
			which_phenotype_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(phenData, Set(phenotype_method_id_ls))
		else:	#if not available, take all phenotypes
			which_phenotype_ls = range(len(phenData.col_id_ls))
		pdata = PassingData(snpData=newSnpData, phenData=phenData, PC_matrix=PC_matrix, which_phenotype_ls=which_phenotype_ls, \
						phenotype_method_id_ls=phenotype_method_id_ls)
		return pdata
	
	def preprocessForTwoPhenotypeAsso(self, snpData, which_phenotype_ls, data_matrix_phen):
		"""
		2009-3-20
			refactored out of run(), easy for MpiAssociation.py to call
		"""
		if len(which_phenotype_ls)<2:
			sys.stderr.write("Error: Require to specify 2 phenotypes in order to carry out this test type.\n")
			sys.exit(3)
		snpData.data_matrix = numpy.vstack((snpData.data_matrix, snpData.data_matrix))
		#stack the two phenotypes, fake a data_matrix_phen, which_phenotype_ls
		no_of_strains = len(strain_acc_list)
		which_phenotype1, which_phenotype2 = which_phenotype_ls[:2]
		phenotype1 = data_matrix_phen[:, which_phenotype1]
		phenotype1 = numpy.resize(phenotype1, [no_of_strains, 1])	#phenotype1.resize([10,1]) doesn't work. ValueError: 'resize only works on single-segment arrays'
		phenotype2 = data_matrix_phen[:, which_phenotype2]
		phenotype2 = numpy.resize(phenotype2, [no_of_strains, 1])
		data_matrix_phen = numpy.vstack((phenotype1, phenotype2))
		which_phenotype_index_ls = [0]
		#stack the PC_matrix as well
		if PC_matrix is not None:
			PC_matrix = numpy.vstack((PC_matrix, PC_matrix))
		
		#create an environment variable
		a = numpy.zeros(no_of_strains)
		a.resize([no_of_strains, 1])
		b = numpy.ones(no_of_strains)
		b.resize([no_of_strains, 1])
		environment_matrix = numpy.vstack((a, b))
		return which_phenotype_index_ls, environment_matrix, data_matrix_phen
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		initData = self.readInData(self.phenotype_fname, self.input_fname, self.eigen_vector_fname, self.phenotype_method_id_ls, self.test_type, self.report)
		
		if self.test_type==5 or self.test_type==6:
			which_phenotype_index_ls, environment_matrix, initData.phenData.data_matrix = self.preprocessForTwoPhenotypeAsso(initData.snpData, \
																							initData.which_phenotype_ls,\
																							initData.phenData.data_matrix)
		else:
			which_phenotype_index_ls = initData.which_phenotype_ls
			environment_matrix = None
		
		if self.test_type==6:
			gene_environ_interaction = True
		else:
			gene_environ_interaction = False
		
		for which_phenotype in which_phenotype_index_ls:
			if self.test_type==5 or self.test_type==6:
				which_phenotype1, which_phenotype2 = initData.which_phenotype_ls[:2]
				phenotype1_name = initData.phenData.col_id_ls[which_phenotype1]
				phenotype2_name = initData.phenData.col_id_ls[which_phenotype2]
				phenotype_name = phenotype1_name+'_'+phenotype2_name
			else:
				phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			output_fname='%s_pheno_%s.tsv'%(os.path.splitext(self.output_fname)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			results = self.run_whole_matrix[self.test_type](initData.snpData.data_matrix, initData.phenData.data_matrix[:, which_phenotype], \
														self.min_data_point, PC_matrix=initData.PC_matrix, \
														which_PC_index_ls=self.which_PC_index_ls, environment_matrix=environment_matrix,\
														gene_environ_interaction=gene_environ_interaction)
			output_results_func = self.output_results.get(self.test_type)
			if output_results_func is None:
				output_results_func = self.output_lm_results
			output_results_func(results, initData.snpData.col_id_ls, output_fname, self.minus_log_pvalue)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Association
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
