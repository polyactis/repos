#!/usr/bin/env python
"""
Examples:
	#plots of comparison between analysis method 1 and 7 on phenotype 1-7 with gene list 28
	PlotCmpTwoAnalysisMethods.py -e 1-7 -l 28 -u yh -p yh324 -o /Network/Data/250k/tmp-yh/PlotCmpTwoAnalysisMethods/ -m 500 -s ./mnt2/panfs/250k/snps_context_g0_m500 -a 1,7
	
Description:
	Program that makes plots comparing pvalues from two analysis methods. Highlight the SNPs that are close to genes in candidate gene list.
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")
import time, csv, cPickle, numpy, random
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr
import Stock_250kDB
from Stock_250kDB import ResultsByGene, ResultsMethod
from sets import Set
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
from DrawSNPRegion import DrawSNPRegion
#from sqlalchemy.orm import join
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 4
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 4
rcParams['axes.titlesize'] = 6
rcParams['xtick.labelsize'] = 4
rcParams['ytick.labelsize'] = 4
import pylab
import StringIO
from common import get_total_gene_ls
from CheckCandidateGeneRank import CheckCandidateGeneRank

class PlotCmpTwoAnalysisMethods(CheckCandidateGeneRank):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("phenotype_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7'],\
							("min_distance", 0, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 0, float): [0, 'n', 1, 'minimum Minor Allele Frequency. deprecated.'],\
							('min_sample_size', 0, int): [5, 'i', 1, 'minimum size for candidate gene sets to draw histogram'],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							("snps_context_picklef", 0, ): [None, 's', 1, 'given the option, if the file does not exist yet, to store a pickled snps_context_wrapper into it, min_distance and flag get_closest will be attached to the filename. If the file exists, load snps_context_wrapper out of it.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("results_type", 1, int): [1, 'w', 1, 'which type of results. 1; ResultsMethod, 2: ResultsByGene'],\
							("output_dir", 0, ): [None, 'o', 1, 'directory to store output'],\
							('call_method_id', 0, int):[17, 'j', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id_ls', 0, ):['1,7', 'a', 1, 'two analysis_method_ids separated by comma. Compare results from these two analysis methods. If more than two are given, only 1st two is taken.'],\
							("allow_two_sample_overlapping", 1, int): [0, 'x', 0, 'whether to allow one SNP to be assigned to both candidate and non-candidate gene group'],\
							('null_distribution_type_id', 0, int):[1, 'y', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. in db table null_distribution_type'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-11-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.phenotype_id_ls = getListOutOfStr(self.phenotype_id_ls, data_type=int)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		if self.output_dir and not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
	
	def matchPvaluesFromTwoResults(self, genome_wide_result_ls, snps_context_wrapper, candidate_gene_set):
		"""
		2008-11-19
		"""
		sys.stderr.write("Matching pvalues from two results ...")
		pvalue_ls1_in_candidate = []
		pvalue_ls2_in_candidate = []
		pvalue_ls1_in_non_candidate = []
		pvalue_ls2_in_non_candidate = []
		gwr1, gwr2 = genome_wide_result_ls[:2]
		chr_pos_ls1 = gwr1.chr_pos2index.keys()
		chr_pos_ls2 = gwr2.chr_pos2index.keys()
		chr_pos_set = Set(chr_pos_ls1)&Set(chr_pos_ls2)	#intersection of all SNPs
		for chr, pos in chr_pos_set:
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				
				if gene_id in candidate_gene_set:
					assign_snp_candidate_gene = 1
					break
			data_obj1 = gwr1.get_data_obj_by_chr_pos(chr, pos)
			data_obj2 = gwr2.get_data_obj_by_chr_pos(chr, pos)
			if assign_snp_candidate_gene:
				pvalue_ls1_in_candidate.append(data_obj1.value)
				pvalue_ls2_in_candidate.append(data_obj2.value)
			else:
				pvalue_ls1_in_non_candidate.append(data_obj1.value)
				pvalue_ls2_in_non_candidate.append(data_obj2.value)
		
		return_data = PassingData(pvalue_ls1_in_candidate=pvalue_ls1_in_candidate, pvalue_ls2_in_candidate=pvalue_ls2_in_candidate,\
								pvalue_ls1_in_non_candidate=pvalue_ls1_in_non_candidate, pvalue_ls2_in_non_candidate=pvalue_ls2_in_non_candidate)
		sys.stderr.write("Done.\n")
		return return_data
	
	def plot(self, phenotype_method, list_type, rm_ls, pvalue_matching_data, output_dir=None, commit=0):
		"""
		2008-11-19
		"""
		sys.stderr.write("Making plots ...")
		pylab.clf()
		pylab.subplots_adjust(left=0.08, right=0.92,bottom = 0.05)
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		no_of_rows = 1
		rm1, rm2 = rm_ls
		pylab.subplot(no_of_rows,2, 1, frameon=False)
		#pylab.title()
		pylab.xlabel(rm1.analysis_method.short_name)
		pylab.ylabel(rm2.analysis_method.short_name)
		pylab.grid(True, alpha=0.3)
		legend_ls = []
		legend_patch_ls = []
		s2 = pylab.scatter(pvalue_matching_data.pvalue_ls1_in_non_candidate, pvalue_matching_data.pvalue_ls2_in_non_candidate, c='b', alpha=0.3, linewidth=0)
		legend_patch_ls.append(s2)
		legend_ls.append('non-candidate (%s)'%(len(pvalue_matching_data.pvalue_ls1_in_non_candidate)))
		
		s1 = pylab.scatter(pvalue_matching_data.pvalue_ls1_in_candidate, pvalue_matching_data.pvalue_ls2_in_candidate, c='r', alpha=0.3, linewidth=0)
		legend_ls.append('candidate (%s)'%(len(pvalue_matching_data.pvalue_ls1_in_candidate)))
		legend_patch_ls.append(s1)
		
		pylab.legend(legend_patch_ls, legend_ls, loc='upper left', handlelen=0.02)
						
		ax = pylab.axes([0.1, 0.1, 0.8,0.8], frameon=False)
		ax.set_xticks([])
		ax.set_yticks([])
		title = 'Phenotype %s %s by %s %s Analysis %s vs %s'%(phenotype_method.id, phenotype_method.short_name, list_type.id, \
													list_type.short_name, rm1.analysis_method.short_name, rm2.analysis_method.short_name)
		ax.set_title(title)
		
		png_data = None
		svg_data = None
		if commit:
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=300)
			#pylab.savefig(svg_data, format='svg', dpi=300)
		elif output_dir:
			output_fname_prefix = os.path.join(output_dir, title.replace('/', '_'))
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
			#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		
		"""
		#2008-11-22 this part is commented out right now. which looks at distribution of (pvalue1, pvalue2) tuples by a 2D histogram
		import rpy
		rpy.r.library('gplots')
		pvalue_ls1 = pvalue_matching_data.pvalue_ls1_in_non_candidate + pvalue_matching_data.pvalue_ls1_in_candidate
		pvalue_ls1 = [10**(-value) for value in pvalue_ls1]
		pvalue_ls2 = pvalue_matching_data.pvalue_ls2_in_non_candidate + pvalue_matching_data.pvalue_ls2_in_candidate
		pvalue_ls2 = [10**(-value) for value in pvalue_ls2]
		rpy.r.postscript('%s.ps'%output_fname_prefix)
		h2d = rpy.r.hist2d(pvalue_ls1, pvalue_ls2, same_scale=rpy.r.TRUE, nbins=[50,50])	#show=rpy.r.FALSE, 
		#2008-11-21 3D plot
		#persp( h2d$x, h2d$y, h2d$counts, ticktype="detailed", theta=30, phi=30, expand=0.5, shade=0.5, col="cyan", ltheta=-30)
		rpy.r.dev_off()
		rpy.r.postscript('%s_3D.ps'%output_fname_prefix)
		rpy.r.persp(h2d['x'], h2d['y'], h2d['counts'], xlab=rm1.analysis_method.short_name, ylab=rm2.analysis_method.short_name, \
				zlab='counts', ticktype="detailed", theta=30, phi=30, expand=0.5, shade=0.5, col="cyan", ltheta=-30)
		rpy.r.dev_off()
		"""
		sys.stderr.write("Done.\n")
		return png_data, svg_data
	
	def run(self):
		"""
		2008-11-19
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		
		if self.results_type==1:
			ResultsClass = Stock_250kDB.ResultsMethod
			snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		elif self.results_type==2:
			ResultsClass = Stock_250kDB.ResultsByGene
		else:
			sys.stderr.write("Invalid results type : %s.\n"%self.results_type)
			return None
		
		hist_type = self.getHistType(self.call_method_id, self.min_distance, self.get_closest, self.min_MAF, \
									self.allow_two_sample_overlapping, self.results_type, self.null_distribution_type_id)
		
		candidate_gene_list = self.dealWithCandidateGeneList(self.list_type_id)
		candidate_gene_set = self.dealWithCandidateGeneList(self.list_type_id, return_set=True)
		list_type = Stock_250kDB.GeneListType.get(self.list_type_id)
		
		analysis_method_id_set = Set(self.analysis_method_id_ls[:2])	#only first two
		phenotype_id2results_id_ls = self.getResultsIDLs(db, ResultsClass, self.results_type, self.phenotype_id_ls, \
														self.min_distance, self.get_closest, self.min_MAF, self.call_method_id,\
														analysis_method_id_set=analysis_method_id_set)
		
		param_data = PassingData(results_directory=self.results_directory, candidate_gene_list=candidate_gene_list, \
			min_MAF=self.min_MAF, allow_two_sample_overlapping=self.allow_two_sample_overlapping, need_the_value=1, \
			construct_chr_pos2index=True)
			#need_the_value means to get the pvalue/score
		
		for phenotype_id, results_id_ls in phenotype_id2results_id_ls.iteritems():
			"""
			if hist_type.id:	#hist_type already in database
				rows = Stock_250kDB.MatchTwoAnalysisMethodPlot.query.filter_by(phenotype_method_id=phenotype_id).\
					filter_by(list_type_id=self.list_type_id).filter_by(hist_type_id=hist_type.id)
				if rows.count()>0:
					row = rows.first()
					sys.stderr.write("Histogram already in database. id=%s, phenotype_id=%s, list_type_id=%s, hist_type_id=%s.\n"%\
									(row.id, row.phenotype_method_id, row.list_type_id, row.hist_type_id))
					continue
			"""
			phenotype_method = Stock_250kDB.PhenotypeMethod.get(phenotype_id)
			if not phenotype_method:
				continue
			
			sys.stderr.write("Checking phenotype %s (%s) on list_type %s (%s) ...\n"%\
							(phenotype_method.id, phenotype_method.short_name, list_type.id, list_type.short_name))
			rm_ls = []
			genome_wide_result_ls = []
			results_id_ls.reverse()	#in descending analysis method to make KW behind and on the y-axis
			for results_id in results_id_ls:
				rm = ResultsClass.get(results_id)
				genome_wide_result = self.getResultMethodContent(rm, self.results_directory, self.min_MAF, pdata=param_data)
				if not genome_wide_result:
					continue
				rm_ls.append(rm)
				genome_wide_result_ls.append(genome_wide_result)
			
			if len(genome_wide_result_ls)!=2:	#not two results, skip
				continue
			
			pvalue_matching_data = self.matchPvaluesFromTwoResults(genome_wide_result_ls, snps_context_wrapper, candidate_gene_set)
			self.plot(phenotype_method, list_type, rm_ls, pvalue_matching_data, output_dir=self.output_dir, commit=self.commit)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotCmpTwoAnalysisMethods
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()