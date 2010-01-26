#!/usr/bin/env python
"""
2008-11-18
A module to do PCA.
"""

import sys, os, math
import numpy

class PCA(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1,):['', 'i', 1, 'input file', ],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-11-18
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def normalize(cls, data_matrix, divide_variance=True):
		"""
		2008-11-18
			this is based on Patterson2006a.
			
			assume no missing data in data_matrix. SNP allele in data_matrix is 0 or 1. two chromosomes are homozygous.
			so each strain has either allele 0 or 1 on both chromosomes. (selfing plants, like arabidopsis thaliana)
			
			normalizing steps:
				1. subtract mean from each column
				2. make variance equal among columns by dividing each column by the estimated stdev
					if variance is 0, skip this step.
		"""
		sys.stderr.write("Normalizing ...")
		no_of_rows, no_of_cols = data_matrix.shape
		new_data_matrix = data_matrix.astype(numpy.float)	#change the data type
		
		for j in range(no_of_cols):
			genotype_ls = data_matrix[:,j]
			col_mean = numpy.mean(genotype_ls)
			col_var = numpy.var(genotype_ls)
			#col_var = col_mean*(1-col_mean)	#2009-9-3 only good for binary matrix
			if col_mean!=0 and divide_variance and col_var!=0:
				new_data_matrix[:,j] = (new_data_matrix[:,j]-col_mean)/numpy.sqrt(col_var)
			else:
				new_data_matrix[:,j] = new_data_matrix[:,j]-col_mean
		sys.stderr.write("Done.\n")
		return new_data_matrix
	normalize = classmethod(normalize)
	
	def eig(cls, data_matrix, normalize=True):
		"""
		2008-11-18
			numpy.inner(matrix1, matrix2) is weird. matrix1's 2nd dimension is equal to matrix2's 2nd dimension.
				Not as traditional, matrix1's 2nd dimension is equal to matrix2's 1st dimension.
		"""
		if normalize:
			new_data_matrix = cls.normalize(data_matrix)
			new_data_matrix2 = numpy.transpose(new_data_matrix)	#2008-11-20 transpose only to get cov_matrix. a later mulitplication between eigen_vector and new_data_matrix
			cov_matrix = 1.0/new_data_matrix2.shape[1]*numpy.inner(new_data_matrix2, new_data_matrix2)	#2008-11-19 mysteriously, numpy.tranpose() is not required on the 2nd new_data_matrix, it'll cause ValueError('matrices are not aligned',)
		else:
			new_data_matrix = cls.normalize(data_matrix, divide_variance=False)
			#new_data_matrix = data_matrix
			cov_matrix = numpy.cov(new_data_matrix, rowvar=0)	#2009-9-3 numpy.cov or numpy.corrcoef
		#eigen_values, eigen_vectors = numpy.linalg.eig(cov_matrix)	#get complex values out of this
		import rpy
		eigen_result = rpy.r.eigen(cov_matrix)
		eigen_values = numpy.array(eigen_result['values'])
		eigen_vectors = eigen_result['vectors']
		explained_var = eigen_values/numpy.sum(eigen_values)
		pc_matrix = numpy.inner(numpy.transpose(eigen_vectors), new_data_matrix)	#eigen_vectors has to be transposed in row-vector form
		pc_matrix = numpy.transpose(pc_matrix)	#transpose again. maybe pc_matrix = numpy.inner(new_data_matrix, eigen_vectors)
		return pc_matrix, eigen_vectors, explained_var
	
	eig = classmethod(eig)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawFTPathway
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()