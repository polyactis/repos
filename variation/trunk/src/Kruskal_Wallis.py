import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv


class Kruskal_Wallis:
	"""
	2008-02-14
		class to do kruskal wallis test on SNP data
	Argument list:
		-i ..., input_fname*
		-o ..., output_fname*
		-a ..., phenotype_fname(argument1)*
		-e ..., minus log the pvalue(argument2), 0 or 1
		-b,	toggle debug
		-r, toggle report
	"""
	def __init__(self, input_fname, phenotype_fname, output_fname, log_pvalue=0, debug=0, report=0):
		"""
		2008-02-14
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
	
	def get_phenotype_ls_in_data_matrix_order(self, strain_acc_list, strain_acc_list_phen, category_list_phen, data_type=float):
		"""
		2008-02-14
		"""
		sys.stderr.write("Getting phenotype_ls in order with data_matrix ...")
		strain_acc_phen2index = dict(zip(strain_acc_list_phen, range(len(strain_acc_list_phen)) ) )
		phenotype_ls = []
		for strain_acc in strain_acc_list:
			phen_index = strain_acc_phen2index[strain_acc]
			phenotype_ls.append(float(category_list_phen[phen_index]))
		sys.stderr.write("Done.\n")
		return phenotype_ls
	
	def _kruskal_wallis(self, data_matrix, phenotype_ls, min_data_point=3):
		"""
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
				if genotype_ls[i]!=0:
					non_NA_genotype = genotype_ls[i]
					non_NA_genotype_ls.append(non_NA_genotype)
					non_NA_phenotype_ls.append(phenotype_ls[i])
					if non_NA_genotype not in non_NA_genotype2count:
						non_NA_genotype2count[non_NA_genotype] = 0
					non_NA_genotype2count[non_NA_genotype] += 1
			count_ls = non_NA_genotype2count.values()
			if len(count_ls)>=2 and count_ls[0]>=min_data_point and count_ls[1]>=min_data_point:
				pvalue = rpy.r.kruskal_test(x=non_NA_phenotype_ls, g=rpy.r.as_factor(non_NA_genotype_ls))['p.value']
				results.append([j, pvalue])
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
		2008-02-14
		"""
		sys.stderr.write("Outputting pvalue results ...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for column_index, pvalue in kw_results:
			SNP_name = SNP_header[column_index]
			SNP_pos_ls = SNP_name.split('_')
			if log_pvalue:
				if pvalue>0:
					pvalue = -math.log(pvalue)
				else:
					pvalue = 'NA'
			writer.writerow([SNP_pos_ls[0], SNP_pos_ls[1], pvalue])
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = FilterStrainSNPMatrix_instance.read_data(self.phenotype_fname)
		phenotype_ls = self.get_phenotype_ls_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, category_list_phen, data_type=float)
		kw_results = self._kruskal_wallis(data_matrix, phenotype_ls)
		self.output_kw_results(kw_results, header[2:], self.output_fname, self.log_pvalue)
		