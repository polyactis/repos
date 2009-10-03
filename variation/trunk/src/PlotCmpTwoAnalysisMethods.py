#!/usr/bin/env python
"""
Examples:
	#plots of comparison between analysis method 1 and 7 on phenotype 1-7 with gene list 28
	PlotCmpTwoAnalysisMethods.py -e 1-7 -l 28 -u yh -p passw** -o /Network/Data/250k/tmp-yh/PlotCmpTwoAnalysisMethods/ -m 500 -s ./mnt2/panfs/250k/snps_context_g0_m500 -a 1,7
	
	#save the enrichment ratio by pvalue-based thresholds (-A, -B) into db
	~/script/variation/src/PlotCmpTwoAnalysisMethods.py -e 1-7,39-48,57-59,80-82 -l 145 -o /Network/Data/250k/tmp-yh/MatchedPvalueData -m 20000 
	-s /Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m20000 -a 1,7 -q 1 -A 3.5 -B 5.5 -j 32 -c
	
	#output the association from both results SNP by SNP
	~/script/variation/src/PlotCmpTwoAnalysisMethods.py -e 1-7,39-48,57-59,80-82 -l 145
	-o /Network/Data/250k/tmp-yh/MatchedPvalueData -m 20000 -s /Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m20000
	-a 1,7 -q 2 -A 2.5 -B 3 -j 32

	#save the enrichment ratio by rank-based thresholds (-A, -B) into db
	~/script/variation/src/PlotCmpTwoAnalysisMethods.py -e 1-7,39-48,57-59,80-82 -l 145
	-o /Network/Data/250k/tmp-yh/MatchedPvalueData -m 20000 -s /Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m20000
	-a 1,7 -q 3 -A 700 -B 700 -j 32 -c
	
	#null distribution type 2, to get enrichment ratio & permutation pvalues for published flowering phenotypes 
	PlotCmpTwoAnalysisMethods.py -e 1-7,39-48,57-59,80-82 -l 145 -o /Network/Data/250k/tmp-yh/MatchedPvalueData -m 20000
	-s /Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m20000 -a 1,7 -q 1 -A 3 -B 5 -j 32 -c -y2
	
Description:
	Program that makes 2 figures in one plot.
	Figure 1 compares pvalues from two analysis methods. Highlight the SNPs that are close to genes in candidate gene list.
	Figure 2 (on the right) has same x,y axis with space paritioned into sub-spaces by pvalue cutoffs. Each sub-space is painted
	in a color representing the enrichment ratio (= candidate-ratio/non-candidate-ratio).
	
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
import pylab
import StringIO
from common import get_total_gene_ls
from CheckCandidateGeneRank import CheckCandidateGeneRank
from TopSNPTest import TopSNPTest
from matplotlib.patches import Polygon, CirclePolygon, Ellipse, Wedge
from matplotlib import rcParams
rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 8
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 8
rcParams['axes.titlesize'] = 10
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8

class PlotCmpTwoAnalysisMethods(CheckCandidateGeneRank, TopSNPTest):
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
							('min_MAF', 0, float): [0, 'n', 1, 'minimum Minor Allele Frequency. If above 1, it is minor allele count. only applied to Emma-based methods'],\
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
							('run_type', 0, int):[1, 'q', 1, 'type of running. 1: just get all pvalues into different bins, \
	2: return data above threshold from both results in matched fashion, 3: r1_pvalue_cutoff and r2_pvalue_cutoff become rank-based thresholds.'],\
							('pvalue_int_gap', 0, int):[1., 'f', 1, 'the size of bin to partition the pvalues. It is used only when both r1_pvalue_cutoff and r2_pvalue_cutoff are zero.'],\
							('r1_pvalue_cutoff', 0, float):[3., 'A', 1, '-log pvalue cutoff for the analysis method with bigger id (emma if -a 1,7)'],\
							('r2_pvalue_cutoff', 0, float):[5., 'B', 1, '-log pvalue cutoff for the analysis method with smaller id (kw if -a 1,7)'],\
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
	
	@classmethod
	def matchPvaluesFromTwoResults(cls, genome_wide_result_ls, snps_context_wrapper, candidate_gene_set, pvalue_int_gap=1.0, \
								r1_pvalue_cutoff=3, r2_pvalue_cutoff=5, run_type=1, need_snp_index=False, need_chr_pos_ls=True):
		"""
		2009-10-2
			add argument need_snp_index. toggled to True to return some data, necessary for permutation pvalue of enrichment_ratio 
		2009-7-8
			become a classmethod
		2009-5-1
			add argument run_type
				1: just get all pvalues into different bins
				2: return data above threshold from both results in matched fashion
		2009-4-16
			genome_wide_result_ls is sorted descendingly by analysis_method_id (check how it's got in run())
			r1_pvalue_cutoff is cutoff for 1st gwa
			r2_pvalue_cutoff is cutoff for 2nd gwa
		2009-1-7
			get both raw pvalue and -log pvalue
		2008-11-19
		"""
		sys.stderr.write("Matching pvalues from two results ...")
		
		
		pvalue_ls1_in_candidate = []
		pvalue_ls2_in_candidate = []
		logpvalue_ls1_in_candidate = []
		logpvalue_ls2_in_candidate = []
		pvalue_ls1_in_non_candidate = []
		pvalue_ls2_in_non_candidate = []
		logpvalue_ls1_in_non_candidate = []
		logpvalue_ls2_in_non_candidate = []
		
		data_matrix = []	# to store reu
		
		gwa1, gwa2 = genome_wide_result_ls[:2]
		chr_pos_ls1 = gwa1.chr_pos2index.keys()
		chr_pos_ls2 = gwa2.chr_pos2index.keys()
		chr_pos_set = Set(chr_pos_ls1)&Set(chr_pos_ls2)	#intersection of all SNPs
		
		total_chr_pos_ls = list(chr_pos_set)
		total_chr_pos_ls.sort()
		no_of_total_snps = len(chr_pos_set)
		
		pvalue_int_pair2count_candidate = {}
		pvalue_int_pair2count_non_candidate = {}
		pvalue_int_pair2cand_snp_index_ls = {}
		pvalue_int_pair2non_cand_snp_index_ls = {}
		
		for i in range(int(no_of_total_snps)):
			#for chr, pos in chr_pos_set:
			chr, pos = total_chr_pos_ls[i]
			
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				
				if gene_id in candidate_gene_set:
					assign_snp_candidate_gene = 1
					break
			data_obj1 = gwa1.get_data_obj_by_chr_pos(chr, pos)
			data_obj2 = gwa2.get_data_obj_by_chr_pos(chr, pos)
			
			if r1_pvalue_cutoff==0.0 and r2_pvalue_cutoff==0.0:
				int(data_obj1.value/pvalue_int_gap)
				int(data_obj2.value/pvalue_int_gap)
			else:
				pvalue_int1 = int(data_obj1.value>=r1_pvalue_cutoff)
				pvalue_int2 = int(data_obj2.value>=r2_pvalue_cutoff)
			pvalue_int_pair = (pvalue_int1, pvalue_int2)
			
			if run_type==2:
				#if pvalue_int_pair!=(0,0):
				data_row = [chr, pos, data_obj1.value, data_obj2.value, assign_snp_candidate_gene]
				data_matrix.append(data_row)
				continue	#ignore the rest
			
			if assign_snp_candidate_gene:
				logpvalue_ls1_in_candidate.append(data_obj1.value)
				logpvalue_ls2_in_candidate.append(data_obj2.value)
				pvalue_ls1_in_candidate.append(10**(-data_obj1.value))
				pvalue_ls2_in_candidate.append(10**(-data_obj2.value))
				if pvalue_int_pair not in pvalue_int_pair2count_candidate:
					pvalue_int_pair2count_candidate[pvalue_int_pair] = 0
					pvalue_int_pair2cand_snp_index_ls[pvalue_int_pair] = []
				pvalue_int_pair2count_candidate[pvalue_int_pair] += 1
				if need_snp_index:
					pvalue_int_pair2cand_snp_index_ls[pvalue_int_pair].append(i)
			else:
				logpvalue_ls1_in_non_candidate.append(data_obj1.value)
				logpvalue_ls2_in_non_candidate.append(data_obj2.value)
				pvalue_ls1_in_non_candidate.append(10**(-data_obj1.value))
				pvalue_ls2_in_non_candidate.append(10**(-data_obj2.value))
				if pvalue_int_pair not in pvalue_int_pair2count_non_candidate:
					pvalue_int_pair2count_non_candidate[pvalue_int_pair] = 0
					pvalue_int_pair2non_cand_snp_index_ls[pvalue_int_pair] = []
				pvalue_int_pair2count_non_candidate[pvalue_int_pair] += 1
				if need_snp_index:
					pvalue_int_pair2non_cand_snp_index_ls[pvalue_int_pair].append(i)
		
		if need_snp_index:
			#turn them into numpy arrays
			for pvalue_int_pair in pvalue_int_pair2cand_snp_index_ls.keys():
				snp_index_ls = pvalue_int_pair2cand_snp_index_ls[pvalue_int_pair]
				pvalue_int_pair2cand_snp_index_ls[pvalue_int_pair] = numpy.array(snp_index_ls, numpy.int)
			for pvalue_int_pair in pvalue_int_pair2non_cand_snp_index_ls.keys():
				snp_index_ls = pvalue_int_pair2non_cand_snp_index_ls[pvalue_int_pair]
				pvalue_int_pair2non_cand_snp_index_ls[pvalue_int_pair] = numpy.array(snp_index_ls, numpy.int)
		
		total_chr_pos_ar = None
		if need_chr_pos_ls:
			total_chr_pos_ar = numpy.array(total_chr_pos_ls)
		
		
		return_data = PassingData(pvalue_ls1_in_candidate=pvalue_ls1_in_candidate, \
								pvalue_ls2_in_candidate=pvalue_ls2_in_candidate,\
								logpvalue_ls1_in_candidate=logpvalue_ls1_in_candidate, \
								logpvalue_ls2_in_candidate=logpvalue_ls2_in_candidate,\
								pvalue_ls1_in_non_candidate=pvalue_ls1_in_non_candidate, \
								pvalue_ls2_in_non_candidate=pvalue_ls2_in_non_candidate,\
								logpvalue_ls1_in_non_candidate=logpvalue_ls1_in_non_candidate, \
								logpvalue_ls2_in_non_candidate=logpvalue_ls2_in_non_candidate,\
								pvalue_int_pair2count_candidate=pvalue_int_pair2count_candidate, \
								pvalue_int_pair2count_non_candidate=pvalue_int_pair2count_non_candidate,\
								
								pvalue_int_pair2cand_snp_index_ls = pvalue_int_pair2cand_snp_index_ls,\
								pvalue_int_pair2non_cand_snp_index_ls = pvalue_int_pair2non_cand_snp_index_ls,\
								no_of_total_snps = no_of_total_snps,\
								total_chr_pos_ar = total_chr_pos_ar,\
								
								data_matrix=data_matrix)
		sys.stderr.write("Done.\n")
		return return_data
	
	def plotEnrichmentRatio(self, rm_ls, pvalue_matching_data, pvalue_int_gap=1.0, output_fname_prefix=None, \
						no_of_rows=1, no_of_cols=3, which_figure=3,\
						r1_pvalue_cutoff=3, r2_pvalue_cutoff=5):
		"""
		2008-11-25
			plot the enrichment ratio for each pvalue window, represented by colors
		"""
		sys.stderr.write("Plotting enrichment ratio ...")
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		rm1, rm2 = rm_ls
		
		pvalue_int_pair2enrichment_ratio = {}
		for pvalue_int_pair in pvalue_matching_data.pvalue_int_pair2count_candidate:
			if pvalue_int_pair in pvalue_matching_data.pvalue_int_pair2count_non_candidate:
				candidate_ratio = pvalue_matching_data.pvalue_int_pair2count_candidate[pvalue_int_pair]/float(len(pvalue_matching_data.pvalue_ls1_in_candidate))
				non_candidate_ratio = pvalue_matching_data.pvalue_int_pair2count_non_candidate[pvalue_int_pair]/float(len(pvalue_matching_data.pvalue_ls1_in_non_candidate))
				enrichment_ratio = candidate_ratio/non_candidate_ratio
				pvalue_int_pair2enrichment_ratio[pvalue_int_pair] = enrichment_ratio
		
		phenotype_cmap = mpl.cm.jet
		max_phenotype = max(pvalue_int_pair2enrichment_ratio.values())
		min_phenotype = min(pvalue_int_pair2enrichment_ratio.values())
		phenotype_gap = max_phenotype - min_phenotype
		phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype, vmax=max_phenotype)
		
		if output_fname_prefix:
			pylab.savefig('%s_%s.png'%(output_fname_prefix,1), dpi=300)
		
		ax = pylab.subplot(no_of_rows, no_of_cols, which_figure, frameon=False)
		#pylab.title()
		ax.set_xlabel(rm1.analysis_method.short_name)
		ax.set_ylabel(rm2.analysis_method.short_name)
		pylab.grid(True, alpha=0.3)
		
		for pvalue_int_pair, enrichment_ratio in pvalue_int_pair2enrichment_ratio.iteritems():
			pvalue_int1, pvalue_int2 = pvalue_int_pair
			if r1_pvalue_cutoff==0.0 and r2_pvalue_cutoff==0.0:			
				xs = [pvalue_int1*pvalue_int_gap, (pvalue_int1+1)*pvalue_int_gap, (pvalue_int1+1)*pvalue_int_gap, pvalue_int1*pvalue_int_gap]	#anti-clockwise starting from lower left
				ys = [pvalue_int2*pvalue_int_gap, pvalue_int2*pvalue_int_gap, (pvalue_int2+1)*pvalue_int_gap, (pvalue_int2+1)*pvalue_int_gap]	#anti-clockwise starting from lower left
			else:
				x0 = pvalue_int1*r1_pvalue_cutoff
				x1 = (pvalue_int1+1)*r1_pvalue_cutoff
				xs = [x0, x1, x1, x0]	#anti-clockwise starting from lower left
				y0 = pvalue_int2*r2_pvalue_cutoff
				y1 = (pvalue_int2+1)*r2_pvalue_cutoff
				ys = [y0, y0, y1, y1]	#anti-clockwise starting from lower left
			facecolor = phenotype_cmap(phenotype_norm(enrichment_ratio))
			patch = Polygon(zip(xs,ys), facecolor=facecolor, linewidth=0.5)
			ax.add_patch(patch)
		
		axe_map_phenotype_legend = pylab.axes([0.92, 0.08, 0.05, 0.3], frameon=False)
		
		#draw after ax, to avoid overwriting
		cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
									norm=phenotype_norm,
									orientation='vertical')
		cb.set_label('Enrichment Ratio')
		axe_map_phenotype_legend.set_title('Enrichment Ratio')
		sys.stderr.write("Done.\n")
		return ax
	
	@classmethod
	def plot_scatter(cls, ax, pvalue_matching_data, rm_ls, data_type=1):
		"""
		2008-12-17
			refactored out of plot()
		"""
		rm1, rm2 = rm_ls
		ax.set_xlabel(rm1.analysis_method.short_name)
		ax.set_ylabel(rm2.analysis_method.short_name)
		ax.grid(True, alpha=0.3)
		legend_ls = []
		legend_patch_ls = []
		if data_type==1:
			x_ls = pvalue_matching_data.pvalue_ls1_in_non_candidate
			y_ls = pvalue_matching_data.pvalue_ls2_in_non_candidate
		else:
			x_ls = pvalue_matching_data.logpvalue_ls1_in_non_candidate
			y_ls = pvalue_matching_data.logpvalue_ls2_in_non_candidate
		s2 = ax.scatter(x_ls, y_ls, c='b', alpha=0.3, linewidth=0)
		legend_patch_ls.append(s2)
		legend_ls.append('non-candidate (%s)'%(len(pvalue_matching_data.pvalue_ls1_in_non_candidate)))
		
		if data_type==1:
			x_ls = pvalue_matching_data.pvalue_ls1_in_candidate
			y_ls = pvalue_matching_data.pvalue_ls2_in_candidate
		else:
			x_ls = pvalue_matching_data.logpvalue_ls1_in_candidate
			y_ls = pvalue_matching_data.logpvalue_ls2_in_candidate
		s1 = ax.scatter(x_ls, y_ls, c='r', alpha=0.3, linewidth=0)
		legend_ls.append('candidate (%s)'%(len(pvalue_matching_data.pvalue_ls1_in_candidate)))
		legend_patch_ls.append(s1)
		
		ax.legend(legend_patch_ls, legend_ls, loc='upper left', handlelen=0.02)
	
	def plot(self, phenotype_method, list_type, rm_ls, pvalue_matching_data, output_dir=None, commit=0, pvalue_int_gap=1.0,\
			r1_pvalue_cutoff=3, r2_pvalue_cutoff=5):
		"""
		2009-5-1
			add argument r1_pvalue_cutoff, r2_pvalue_cutoff
		2009-1-7
			add a scatter plot of raw pvalues 
		2008-11-19
		"""
		sys.stderr.write("Making plots ...")
		pylab.clf()
		pylab.subplots_adjust(left=0.08, right=0.92,bottom = 0.05)
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		no_of_rows = 1
		no_of_cols = 3
		rm1, rm2 = rm_ls
		ax_scatter_pvalue = pylab.subplot(no_of_rows, no_of_cols, 1, frameon=False)
		self.plot_scatter(ax_scatter_pvalue, pvalue_matching_data, rm_ls, data_type=1)
		
		ax_scatter_logpvalue = pylab.subplot(no_of_rows, no_of_cols, 2, frameon=False)
		self.plot_scatter(ax_scatter_logpvalue, pvalue_matching_data, rm_ls, data_type=2)
		
		ax_enrichment_ratio = self.plotEnrichmentRatio(rm_ls, pvalue_matching_data, pvalue_int_gap,\
													no_of_rows=no_of_rows, no_of_cols=no_of_cols, which_figure=3)
		ax_enrichment_ratio.set_xlim(ax_scatter_logpvalue.get_xlim())
		ax_enrichment_ratio.set_ylim(ax_scatter_logpvalue.get_ylim())
		
		ax = pylab.axes([0.1, 0.1, 0.8,0.8], frameon=False)
		ax.set_xticks([])
		ax.set_yticks([])
		title = 'Call Method %s Phenotype %s %s by %s %s Analysis %s vs %s Pvalue int gap %s r1 r2 pvalue cutoff %s %s'%\
				(rm1.call_method_id, phenotype_method.id, phenotype_method.short_name, list_type.id, \
				list_type.short_name, rm1.analysis_method.short_name, rm2.analysis_method.short_name, pvalue_int_gap,\
				r1_pvalue_cutoff, r2_pvalue_cutoff)
		ax.set_title(title)
		
		png_data = None
		svg_data = None
		"""
		if commit:	#don't need this. never saved into db.
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=300)
			#pylab.savefig(svg_data, format='svg', dpi=300)
		elif output_dir:
		"""
		output_fname_prefix = os.path.join(output_dir, title.replace('/', '_').replace(' ', '_'))
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
	
	def saveDataIntoDB(self, session, genome_wide_result_ls, hist_type, threshold_type, pvalue_matching_data, list_type_id, \
					r1_pvalue_cutoff=3, r2_pvalue_cutoff=5, null_distribution_type_id=1, candidate_gene_set=set(), snps_context_wrapper=None):
		"""
		2009-10-3
			If null_distribution_type_id=2, calculate the permutation pvalue before saving into db.
			If null_distribution_type_id=1, pvalue = None. maybe 2X2 table test (Fisher test) 
		2009-4-16
		"""
		sys.stderr.write("Saving enrichment data into db ...\n")
		results_id1 = genome_wide_result_ls[0].results_id
		results_id2 = genome_wide_result_ls[1].results_id
		
		pvalue_int_pair_set = set(pvalue_matching_data.pvalue_int_pair2count_non_candidate.keys())
		pvalue_int_pair_set.update(set(pvalue_matching_data.pvalue_int_pair2count_candidate.keys()))
		for pvalue_int_pair in pvalue_int_pair_set:
			if pvalue_int_pair[0]==0:
				r1_min_score = 0
				r1_max_score = r1_pvalue_cutoff
			else:
				r1_min_score = r1_pvalue_cutoff
				r1_max_score = None
			if pvalue_int_pair[1]==0:
				r2_min_score = 0
				r2_max_score = r2_pvalue_cutoff
			else:
				r2_min_score = r2_pvalue_cutoff
				r2_max_score = None
			candidate_sample_size = pvalue_matching_data.pvalue_int_pair2count_candidate.get(pvalue_int_pair)
			non_candidate_sample_size = pvalue_matching_data.pvalue_int_pair2count_non_candidate.get(pvalue_int_pair)
			candidate_gw_size = len(pvalue_matching_data.pvalue_ls1_in_candidate)
			non_candidate_gw_size = len(pvalue_matching_data.pvalue_ls1_in_non_candidate)
			if candidate_sample_size is not None and non_candidate_sample_size is not None and candidate_gw_size>0 and non_candidate_sample_size>0:
				enrichment_ratio = (candidate_sample_size*non_candidate_gw_size)/float(non_candidate_sample_size*candidate_gw_size)
			else:
				enrichment_ratio = None
			
			### 2009-10-2
			if null_distribution_type_id==1:
				pvalue = None	# need to figure out a way to calculate the pvalue , maybe 2X2 table test (Fisher test) 
			elif null_distribution_type_id==2:
				cand_snp_index_ls = pvalue_matching_data.pvalue_int_pair2cand_snp_index_ls.get(pvalue_int_pair, numpy.array([], numpy.int))	# without numpy.array around [], [] would be treated as numpy.float64 by default.
				non_cand_snp_index_ls = pvalue_matching_data.pvalue_int_pair2non_cand_snp_index_ls.get(pvalue_int_pair, numpy.array([], numpy.int))
				top_snp_index_ls = numpy.hstack((cand_snp_index_ls, non_cand_snp_index_ls))
				if len(top_snp_index_ls)==0:
					pvalue = 1.
				else:
					return_data = self.get_enrichment_pvalue_by_gw_looping(candidate_sample_size, top_snp_index_ls, candidate_gene_set, \
											snps_context_wrapper, \
											pvalue_matching_data.no_of_total_snps, total_chr_pos_ar=pvalue_matching_data.total_chr_pos_ar, \
											no_of_permutations=20000, no_of_min_breaks=30)
					pvalue = return_data.pvalue
					no_of_tests = return_data.no_of_tests
					no_of_tests_passed = return_data.no_of_tests_passed
			
			entry = Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods(results_id1=results_id1, results_id2=results_id2, list_type_id=list_type_id,\
													type=hist_type, r1_min_score=r1_min_score, r1_max_score=r1_max_score,\
													r2_min_score=r2_min_score, r2_max_score=r2_max_score,\
													candidate_sample_size=candidate_sample_size, non_candidate_sample_size=non_candidate_sample_size,\
													candidate_gw_size=candidate_gw_size, non_candidate_gw_size=non_candidate_gw_size,\
													enrichment_ratio=enrichment_ratio, pvalue=pvalue)
			entry.threshold_type=threshold_type
			session.save(entry)
			#session.flush()
		sys.stderr.write("Done.\n")
	
	def outputMatchedPvalueData(self, phenotype_method, rm_ls, list_type, pvalue_matching_data, output_dir=None):
		"""
		2009-5-1
			output the matched pvalue data (above threshold) under run_type=2
		"""
		sys.stderr.write("Outputting matched pvalue data ...")
		rm1, rm2 = rm_ls
		title = 'Call Method %s Phenotype %s %s by %s %s Analysis %s vs %s'%(rm1.call_method_id, phenotype_method.id, \
													phenotype_method.short_name, list_type.id, \
													list_type.short_name, rm1.analysis_method.short_name, \
													rm2.analysis_method.short_name)
		output_fname = os.path.join(output_dir, title.replace('/', '_').replace(' ', '_')+'.tsv')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['chr', 'pos', 'Emma-pvalue', 'KW-pvalue', 'close_to_candidate_gene']
		writer.writerow(header)
		for data_row in pvalue_matching_data.data_matrix:
			writer.writerow(data_row)
		
		del writer
		
		sys.stderr.write("Done.\n")
	
	def getThresholdType(self, r1_pvalue_cutoff, r2_pvalue_cutoff):
		"""
		2009-5-2
			
		"""
		threshold_type = Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethodsType.query.filter_by(r1_threshold=r1_pvalue_cutoff).\
			filter_by(r2_threshold=r2_pvalue_cutoff).first()
		if not threshold_type:
			threshold_type = Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethodsType(r1_threshold=r1_pvalue_cutoff, r2_threshold=r2_pvalue_cutoff)
		return threshold_type
	
	def checkWhetherTwoResultsBeenCompared(self, hist_type, threshold_type, results_id1, results_id2, list_type_id):
		"""
		2009-5-2
			
		"""
		if hist_type.id and threshold_type.id:	# check if it's already in db
			db_query = Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods.query.filter_by(results_id1=results_id1).\
						filter_by(results_id2=results_id2).\
						filter_by(list_type_id=list_type_id).filter_by(type_id=hist_type.id).\
						filter_by(threshold_type_id=threshold_type.id).first()
			if db_query:
				return True
		return False
	
	def findPvalueCutoffThruRankCutoff(self, genome_wide_result_ls, r1_rank_cutoff, r2_rank_cutoff):
		"""
		2009-5-2
		"""
		sys.stderr.write("Finding pvalue cutoff based on rank cutoff ... \n")
		gwa1, gwa2 = genome_wide_result_ls[:2]
		
		data_obj1 = gwa1.get_data_obj_at_given_rank(int(r1_rank_cutoff))
		if data_obj1 is not None:
			r1_pvalue_cutoff = data_obj1.value
		else:
			sys.stderr.write("\tFail to find pvalue cutoff at rank %s for r1.\n"%(r1_rank_cutoff))
			r1_pvalue_cutoff = 0
		data_obj2 = gwa2.get_data_obj_at_given_rank(int(r2_rank_cutoff))
		if data_obj2 is not None:
			r2_pvalue_cutoff = data_obj2.value
		else:
			sys.stderr.write("\tFail to find pvalue cutoff at rank %s for r2.\n"%(r2_rank_cutoff))
			r2_pvalue_cutoff = 0
		sys.stderr.write("Done.\n")
		return (r1_pvalue_cutoff, r2_pvalue_cutoff)
	
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
		session.begin()
		
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
		threshold_type = self.getThresholdType(self.r1_pvalue_cutoff, self.r2_pvalue_cutoff)
		
		candidate_gene_list = self.dealWithCandidateGeneList(self.list_type_id)
		candidate_gene_set = self.dealWithCandidateGeneList(self.list_type_id, return_set=True)
		list_type = Stock_250kDB.GeneListType.get(self.list_type_id)
		
		analysis_method_id_set = Set(self.analysis_method_id_ls[:2])	#only first two
		phenotype_id2results_id_ls = self.getResultsIDLs(db, ResultsClass, self.results_type, self.phenotype_id_ls, \
														self.min_distance, self.get_closest, self.min_MAF, self.call_method_id,\
														analysis_method_id_set=analysis_method_id_set)
		
		param_data = PassingData(results_directory=self.results_directory, candidate_gene_list=candidate_gene_list, \
			allow_two_sample_overlapping=self.allow_two_sample_overlapping, need_the_value=1, \
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
			
			if self.run_type!=2:	# check whether this has been done before only when run_type is not 2 (matched pvalues output)
				if self.checkWhetherTwoResultsBeenCompared(hist_type, threshold_type, results_id_ls[0], results_id_ls[1], self.list_type_id):
					sys.stderr.write('Phenotype %s (%s) with list_type %s under type %s and threshold_type %s already in db. Skip.\n'%\
								(phenotype_method.id, phenotype_method.short_name, list_type.id, hist_type.id, threshold_type.id))
					continue
			
			
			for results_id in results_id_ls:
				rm = ResultsClass.get(results_id)
				
				if self.min_MAF <1:
					if rm.analysis_method_id==7 or rm.analysis_method_id==4:	#Emma gets the min_MAF set by the user. others get 0 (no cutoff). 
						min_MAF = self.min_MAF
					else:
						min_MAF = 0
					min_MAC = None
				else:	# 2009-5-3 if <1, self.min_MAF is minor allele frequency; otherwise, it's minor allele count.
					min_MAF = None
					if rm.analysis_method_id==7 or rm.analysis_method_id==4:	#Emma gets the cutoff set by the user. others get None (no cutoff). 
						min_MAC = self.min_MAF
					else:
						min_MAC = None
				#param_data.min_MAF = min_MAF	#min_MAF is passed to param_data in getResultMethodContent()
				param_data.min_MAC = min_MAC
				genome_wide_result = self.getResultMethodContent(rm, self.results_directory, min_MAF, pdata=param_data)
				if not genome_wide_result:
					continue
				rm_ls.append(rm)
				genome_wide_result.results_id = results_id
				genome_wide_result_ls.append(genome_wide_result)
			
			if len(genome_wide_result_ls)!=2:	#not two results, skip
				continue
			
			if self.run_type==3:	#self.r1_pvalue_cutoff and self.r2_pvalue_cutoff is rank cutoff, so find corresponding pvalue cutoff
				r1_pvalue_cutoff, r2_pvalue_cutoff = self.findPvalueCutoffThruRankCutoff(genome_wide_result_ls, \
																					self.r1_pvalue_cutoff, self.r2_pvalue_cutoff)
			else:
				r1_pvalue_cutoff, r2_pvalue_cutoff = self.r1_pvalue_cutoff, self.r2_pvalue_cutoff
			
			if self.null_distribution_type_id!=1:
				need_snp_index = True
			else:
				need_snp_index = False
			pvalue_matching_data = self.matchPvaluesFromTwoResults(genome_wide_result_ls, snps_context_wrapper, candidate_gene_set,\
																pvalue_int_gap=self.pvalue_int_gap, r1_pvalue_cutoff=r1_pvalue_cutoff,\
																r2_pvalue_cutoff=r2_pvalue_cutoff, run_type=self.run_type, \
																need_snp_index=need_snp_index)
			
			if self.run_type==1 or self.run_type==3:
				if r1_pvalue_cutoff!=0 or r2_pvalue_cutoff!=0:	#2009-5-2 run only when necessary
					self.saveDataIntoDB(session, genome_wide_result_ls, hist_type, threshold_type, pvalue_matching_data, self.list_type_id, \
									r1_pvalue_cutoff, \
									r2_pvalue_cutoff, null_distribution_type_id=self.null_distribution_type_id, \
									candidate_gene_set=candidate_gene_set, snps_context_wrapper=snps_context_wrapper)
				#2009-5-20 comment out temporarily
				self.plot(phenotype_method, list_type, rm_ls, pvalue_matching_data, output_dir=self.output_dir, commit=self.commit,\
						pvalue_int_gap=self.pvalue_int_gap, r1_pvalue_cutoff=r1_pvalue_cutoff,\
						r2_pvalue_cutoff=r2_pvalue_cutoff)
			elif self.run_type==2:
				self.outputMatchedPvalueData(phenotype_method, rm_ls, list_type, pvalue_matching_data, output_dir=self.output_dir)
		if self.commit:
			session.flush()
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotCmpTwoAnalysisMethods
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()