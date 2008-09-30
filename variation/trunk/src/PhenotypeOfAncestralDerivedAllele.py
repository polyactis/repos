#!/usr/bin/env python
"""

Examples:
	#
	PhenotypeOfAncestralDerivedAllele.py -a 1 -n /Network/Data/A_lyrata/250K.ancestral -e /Network/Data/250k/dataFreeze_080608/phenotypes_091408_raw.tsv -g /tmp/call_method_17_1_1000.tsv -l 17 -t 1 -o /tmp/ancestral_derived
	
	#run it on hpc-cmb cluster, specify input_dir
	PhenotypeOfAncestralDerivedAllele.py -a 1 -n /Network/Data/A_lyrata/250K.ancestral -e /Network/Data/250k/dataFreeze_080608/phenotypes_091408_raw.tsv -g /tmp/call_method_17_1_1000.tsv -l 17 -t 1 -o /tmp/ancestral_derived -i /Network/Data/250k/db/results/type_1
	
Description:
	Program to check phenotype of accessions harboring ancestral/derived alleles.

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
from pymodule import SNPData, PassingData
from pymodule.SNP import number2complement
from variation.src.common import number2nt, nt2number
import Stock_250kDB
from GeneListRankTest import GeneListRankTest

"""
2008-09-19
	check phenotype histogram of ancestral/derived alleles among Top SNPs (ranked by pvalues)
"""
class PhenotypeOfAncestralDerivedAllele(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_dir', 0, ): [None, 'i', 1, 'If given is directory, results_method.filename is assumed to be in this directory.'],\
							('call_method_id', 1, int): [None, 'l', 1, 'id in table call_method', ],\
							('analysis_method_id', 1, int): [None, 'a', 1, 'If given is directory, call_info.filename is assumed to be in this directory.'],\
							('phenotype_method_id', 1, int): [None, 't', 1, 'id in table call_method', ],\
							('genotype_fname', 1, ): [None, 'g', 1, 'File containing the genotypes in StrainXSNP format', ],\
							('phenotype_fname', 1, ): [None, 'e', 1, '', ],\
							('ancestral_allele_fname', 1, ): [None, 'n', 1, 'File containing the ancestral alleles for all SNPs', ],\
							('output_fname_prefix', 0, ): [None, 'o', 1, 'figure filename prefix (no extension, .svg will be added). supply if you want the figure to be saved in svg format.', ],\
							('no_of_top_snps', 1, int): [1000, 'f', 1, 'how many number of top snps ranked by score or -log(pvalue).'],\
							('min_MAF', 1, float): [0.1, '', 1, 'minimum Minor Allele Frequency.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	"""
	2008-09-19
		add option take_unique_ecotype
	"""
	def __init__(self, **keywords):
		"""
		2008-09-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_chr_pos2ancestral_allele(self, ancestral_allele_fname):
		"""
		2008-09-19
			/Network/Data/A_lyrata/250K.ancestral, file generated by Tina
		"""
		sys.stderr.write("Getting chr_pos2ancestral_allele ... ")
		reader = csv.reader(open(ancestral_allele_fname), delimiter='\t')
		reader.next()
		chr_pos2ancestral_allele = {}
		for row in reader:
			chr_pos = row[0]
			ancestral_allele = row[1]
			chromosome = int(chr_pos[0])
			position = int(chr_pos[1:])
			chr_pos = '%s_%s'%(chromosome, position)
			chr_pos2ancestral_allele[chr_pos] = nt2number[ancestral_allele]
		sys.stderr.write("Done.\n")
		return chr_pos2ancestral_allele
	
	def process_phenotype_data(self, pheno_data):
		"""
		2008-09-19
			the header of each column in phenotype_fname is 'phenotype_method_id_phenotype_method_short_name'
			so by parsing the column header, get the phenotype_method_id out
		"""
		sys.stderr.write("Processing phenotype_header to get phenotype_method_id2col_index ...")
		phenotype_method_id2col_index = {}
		for i in range(len(pheno_data.col_id_ls)):
			phenotype_label = pheno_data.col_id_ls[i]
			phenotype_method_id = int(phenotype_label.split('_')[0])
			phenotype_method_id2col_index[phenotype_method_id] = i
		pheno_data.phenotype_method_id2col_index = phenotype_method_id2col_index
		
		
		sys.stderr.write("Done.\n")
		return pheno_data
		
	def get_phenotype_ls(self, rm, no_of_top_snps, chr_pos2ancestral_allele, pheno_data, geno_data, min_MAF, results_directory=None):
		"""
		2008-09-26
			consider the complement of the ancestral allele as ancestral as well.
		2008-09-19
			differentiate phenotype between ancestral/derived alleles
		"""
		sys.stderr.write("Getting phenotype_ls for ancestral/derived alleles ... ")
		genome_wide_result = GeneListRankTest.getResultMethodContent(rm, results_directory, min_MAF)
		genome_wide_result.data_obj_ls.sort()	#in value descending order. each SNP object has a defined method for comparison based on its value
		genome_wide_result.data_obj_ls.reverse()
		ancestral_allele_phenotype_ls = []
		derived_allele_phenotype_ls = []
		no_of_genotype_accessions = len(geno_data.row_id_ls)
		pheno_data_col_index = pheno_data.phenotype_method_id2col_index[rm.phenotype_method_id]
		no_of_polarized_snps = 0
		for i in range(no_of_top_snps):
			data_obj = genome_wide_result.data_obj_ls[i]
			chr_pos = '%s_%s'%(data_obj.chromosome, data_obj.position)
			geno_data_col_index = geno_data.col_id2col_index.get(chr_pos)
			if chr_pos in chr_pos2ancestral_allele and geno_data_col_index is not None:
				no_of_polarized_snps += 1
				for i in range(no_of_genotype_accessions):
					allele = geno_data.data_matrix[i][geno_data_col_index]
					pheno_row_index = pheno_data.row_id2row_index.get(geno_data.row_id_ls[i])	#find corresponding accession row index in phenotype matrix
					if pheno_row_index is not None and allele>=1 and allele<=4:	#no heterozygote or NA or deletion
						phenotype = pheno_data.data_matrix[pheno_row_index][pheno_data_col_index]
						if phenotype!='NA':
							phenotype = float(phenotype)
							if allele==chr_pos2ancestral_allele[chr_pos] or allele==number2complement[chr_pos2ancestral_allele[chr_pos]]:
								ancestral_allele_phenotype_ls.append(phenotype)
							else:
								derived_allele_phenotype_ls.append(phenotype)
		phenotype_ls_data = PassingData()
		phenotype_ls_data.ancestral_allele_phenotype_ls = ancestral_allele_phenotype_ls
		phenotype_ls_data.derived_allele_phenotype_ls = derived_allele_phenotype_ls
		sys.stderr.write("no_of_polarized_snps/no_of_top_snps=%s/%s=%.3f. Done.\n"%(no_of_polarized_snps, no_of_top_snps, no_of_polarized_snps/float(no_of_top_snps)))
		return phenotype_ls_data
		
	
	def run(self):
		"""
		
		"""
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.user,
				   password=self.passwd, hostname=self.hostname, database=self.dbname)
		db.setup(create_tables=False)
		session = db.session
		if self.debug:
			import pdb
			pdb.set_trace()
		chr_pos2ancestral_allele = self.get_chr_pos2ancestral_allele(self.ancestral_allele_fname)
		pheno_data = SNPData(input_fname=self.phenotype_fname, turn_into_integer=0, ignore_2nd_column=1)
		pheno_data = self.process_phenotype_data(pheno_data)
		
		geno_data = SNPData(input_fname=self.genotype_fname, turn_into_array=1, matrix_data_type=int, ignore_2nd_column=1)
		
		query = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=self.call_method_id).filter_by(analysis_method_id=self.analysis_method_id).filter_by(phenotype_method_id=self.phenotype_method_id)
		if query.count()==1:
			rm = query.first()
		elif query.count()>1:
			sys.stderr.write("Warning: more than 1 results_method for call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%(self.call_method_id, self.analysis_method_id, self.phenotype_method_id))
			rm = query.first()
		else:
			sys.stderr.write("Error: no results_method for call_method_id=%s, analysis_method_id=%s, phenotype_method_id=%s.\n"%(self.call_method_id, self.analysis_method_id, self.phenotype_method_id))
			sys.exit(3)
		
		phenotype_ls_data = self.get_phenotype_ls(rm, self.no_of_top_snps, chr_pos2ancestral_allele, pheno_data, geno_data, \
												self.min_MAF, results_directory=self.input_dir)
		
		import pylab
		pylab.clf()
		hist_patch_ls = []
		legend_ls = []
		if len(phenotype_ls_data.ancestral_allele_phenotype_ls)>2:
			n1 = pylab.hist(phenotype_ls_data.ancestral_allele_phenotype_ls, 100, alpha=0.4, normed=1)
			hist_patch_ls.append(n1[2][0])	#first patch in all patches of a histogram
			legend_ls.append('ancestral allele')
		if len(phenotype_ls_data.derived_allele_phenotype_ls)>2:
			n2 = pylab.hist(phenotype_ls_data.derived_allele_phenotype_ls, 100, alpha=0.4, normed=1, facecolor='r')
			hist_patch_ls.append(n2[2][0])
			legend_ls.append('derived allele')
		pylab.legend(hist_patch_ls, legend_ls)
		if self.output_fname_prefix:
			pylab.savefig('%s.svg'%self.output_fname_prefix, dpi=300)
		#pylab.show()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PhenotypeOfAncestralDerivedAllele
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
