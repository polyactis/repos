#!/usr/bin/env python
"""
Examples:
	#SNP pairs from two different genes
	PlotGenePairAssoResult.py -i /Network/Data/250k/tmp-yh/FTGenePair/ft_gene_pairs.tsv -s ./mnt2/panfs/250k/snps_context_g0_m5000 -I ./mnt2/panfs/250k/boolean_snp_pair_ft_gene_phenotype1_4/SNPpair_1_LD.tsv -u yh -p pass*** -q 1  -o /Network/Data/250k/tmp-yh/FTGenePair/ -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf
	
	#SNP pairs from the same gene.
	PlotGenePairAssoResult.py -i /Network/Data/250k/tmp-yh/FTGenePair/ft_gene_pair_itself.tsv -s ./mnt2/panfs/250k/snps_context_g0_m5000 -I /Network/Data/250k/db/results/type_1/2559_results.tsv -u yh -p pass*** -q 1  -o /Network/Data/250k/tmp-yh/FTGenePair/ -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf

Description:
	2008-11-30 program to plot association significance results of SNP pairs (BooleanSNPPair method) from two genes.
	
	The input file contains at least 2 columns, gene1_id, gene2_id. tab or comma-delimited (auto detected).
		The data starts from the 1st row (no header).
		
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, GeneModel, read_data, SNPData
import Stock_250kDB
from sets import Set
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment
import matplotlib as mpl
from GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()
from PlotGroupOfSNPs import PlotGroupOfSNPs

from matplotlib.patches import Polygon, CirclePolygon
import pylab
from pymodule.yh_matplotlib_artists import ExonIntronCollection
import ImageColor
import numpy, StringIO
from pymodule.latex import outputMatrixInLatexTable
from pymodule.latex import outputFigureInLatex
from common import getEcotypeInfo
from Kruskal_Wallis import Kruskal_Wallis
from DrawSNPRegion import DrawSNPRegion
from GeneListRankTest import SnpsContextWrapper


from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6

class PlotGenePairAssoResult(DrawSNPRegion):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 0, float): [0, 'n', 1, 'minimum Minor Allele Frequency.'],\
							("list_type_id", 0, int): [0, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('rbg_results_directory', 0, ):[os.path.expanduser('~/panfs/db/results_by_gene/'), '', 1, 'The rbg results directory. Default is None. use the one given by db.'],\
							("input_fname", 1, ): [None, 'i', 1, 'Filename which contains at least 2 columns: gene_id1, gene_id2'],\
							('boolean_pair_fname', 0, ): ['', 'I', 1, 'file containing boolean pair result', ],\
							('snps_context_picklef',1, ): [None, 's', 1, 'a file containing a pickled snps_context_wrapper. outputted by GeneListRankTest.constructDataStruc()'],\
							('phenotype_fname', 0, ): [None, 'N', 1, 'phenotype file, if snp_matrix_fname is given, this is needed as well.', ],\
							("output_dir", 1, ): [None, 'o', 1, 'directory to store all images'],\
							('call_method_id', 0, int):[17, '', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, ):['1', 'a', 1, 'Restrict results based on this list of analysis_method. Default is no such restriction.'],\
							('phenotype_method_id', 1, ):[None, 'q', 1, 'Restrict results based on this list of phenotype_method. Default is no such restriction.'],\
							('no_of_top_hits', 1, int): [1000, 'f', 1, 'how many number of top hits based on score or -log(pvalue).'],\
							("LD_info_picklef", 0, ): [None, 'D', 1, 'given the option, If the file does not exist yet, store a pickled LD_info into it (min_MAF and min_distance will be attached to the filename). If the file exists, load LD_info out of it.'],\
							("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
							("draw_LD_relative_to_center_SNP", 0, int): [0, '', 0, 'toggle this to draw LD between other SNPs and the center SNP'],\
							("plot_type_short_name", 0, ): [None, 'y', 1, 'A short name (<256 chars) to characterize the type of all the plots. it could have existed in table SNPRegionPlotType'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-11-25
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	
	def get_gene_id2snps_id_ls(self, snps_context_wrapper):
		"""
		2008-11-25 from MpiIntraGeneSNPPairAsso.py
		"""
		sys.stderr.write("Getting gene_id2snps_id_ls  ...")
		gene_id2snps_id_ls = {}
		for chrpos_key, snps_id in snps_context_wrapper.chrpos2snps_id.iteritems():
			chr_pos_label = (chrpos_key[0], chrpos_key[1])	#2008-11-25 tuple is better
			if 'touch' in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id]:	#'touch' is closest, forget about all others if there's 'touch'
				for disp_pos, gene_id in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id]['touch']:
					if gene_id not in gene_id2snps_id_ls:
						gene_id2snps_id_ls[gene_id] = []
					gene_id2snps_id_ls[gene_id].append(chr_pos_label)
			else:
				for left_or_right, gene_id_ls in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id].iteritems():
					for disp_pos, gene_id in gene_id_ls:
						if gene_id not in gene_id2snps_id_ls:
							gene_id2snps_id_ls[gene_id] = []
						gene_id2snps_id_ls[gene_id].append(chr_pos_label)
		sys.stderr.write("Done.\n")
		return gene_id2snps_id_ls
	
	def read_input_fname(self, input_fname):
		sys.stderr.write("Getting gene pairs from %s ..."%input_fname)
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		gene_id_pair_ls = []
		for row in reader:
			gene1_id = int(row[0])
			gene2_id = int(row[1])
			gene_id_pair_ls.append((gene1_id, gene2_id))
		
		sys.stderr.write("Done.\n")
		return gene_id_pair_ls
	
	def loadDataStructure(self, gene_annotation_picklef, min_MAF=0, min_distance=20000, \
						list_type_id=None):
		"""
		2008-11-25
		2008-10-01
			wrap a few functions up, convenient for both run() and drawSNPRegion()
		"""
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		self.db = db
		snp_info = self.getSNPInfo(db)
		gene_annotation = self.dealWithGeneAnnotation(gene_annotation_picklef)
		if list_type_id:
			candidate_gene_list = self.getGeneList(list_type_id)
			candidate_gene_set = Set(candidate_gene_list)
		else:
			candidate_gene_set = Set()
		
		return_data = PassingData(gene_annotation=gene_annotation, snp_info=snp_info, \
								candidate_gene_set=candidate_gene_set)
		return return_data
	
	def get_snp_pair2value_type(self, boolean_pair_fname):
		"""
		2008-11-25
		"""
		sys.stderr.write("Getting snp_pair2value_type ...")
		snp_pair2value_type = {}
		reader = csv.reader(open(boolean_pair_fname), delimiter=figureOutDelimiter(boolean_pair_fname))
		reader.next()
		min_value = None
		max_value = None
		for row in reader:
			snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row
			if not snp2_id:
				snp2_id = snp1_id
				continue	#2008-11-26 skip a row if it's pvalue from single SNP.
			
			if not gene2_id:
				gene2_id = gene1_id
			
			snp1_id = snp1_id.split('_')
			snp1_id = map(int, snp1_id)
			
			snp2_id = snp2_id.split('_')
			snp2_id = map(int, snp2_id)
			
			snp_pair = [tuple(snp1_id), tuple(snp2_id)]
			snp_pair.sort()
			snp_pair = tuple(snp_pair)
			pvalue = -math.log10(float(pvalue))
			
			value = pvalue
			if min_value is None:
				min_value =value
			elif value<min_value:
				min_value = value
				
			if max_value is None:
				max_value = value
			elif value>max_value:
				max_value = value
			
			bool_type = int(bool_type)
			if snp_pair not in snp_pair2value_type:
				snp_pair2value_type[snp_pair] = (pvalue, bool_type)
			else:
				if pvalue>snp_pair2value_type[snp_pair][0]:	#only take maximum
					snp_pair2value_type[snp_pair] = (pvalue, bool_type)
		del reader
		sys.stderr.write("Done.\n")
		return snp_pair2value_type, min_value, max_value
	
	bool_type2marker_and_name = {1: ('v', 'AND'),
					2: ('+', 'Inhibition'),
					4: ('o', 'Inhibition'),
					6: ('^', 'XOR'),
					7: ('D', 'OR')}
		
	def drawSNPPairPvalue(self, axe_pvalue, snps_region1, snps_region2, snp_pair2value_type, phenotype_cmap, phenotype_norm):
		"""
		"""
		sys.stderr.write("Drawing snp pair pvalue.\n")
		pscatter_ls = []
		legend_ls = []
		bool_type2counter = {}
		for bool_type, marker in self.bool_type2marker_and_name.iteritems():
			for chr_pos1 in snps_region1.chr_pos2adjacent_window:
				for chr_pos2 in snps_region2.chr_pos2adjacent_window:
					snp_pair = (min(chr_pos1, chr_pos2),max(chr_pos1, chr_pos2))
					if snp_pair in snp_pair2value_type:
						value, bool_type = snp_pair2value_type[snp_pair]
						if bool_type not in bool_type2counter:
							bool_type2counter[bool_type] = 0
						bool_type2counter[bool_type] += 1
						
						facecolor = phenotype_cmap(phenotype_norm(value))
						marker, bool_type_name = self.bool_type2marker_and_name[bool_type]
						pscatter = axe_pvalue.plot([chr_pos1[1]], [chr_pos2[1]], marker=marker, markersize=3, linewidth=0.6, markeredgecolor=facecolor, markerfacecolor=facecolor)
						#markersize=8
						if bool_type2counter[bool_type] ==1:	#first time to draw this bool type
							legend_ls.append(bool_type_name)
							pscatter_ls.append(pscatter[0])
		sys.stderr.write("Done.\n")
		return legend_ls, pscatter_ls
	
	def drawSNPLinesCrossGenesAndPvalue(self, axe_gene_model, snps_within_this_region, gene_model_min_y,\
				gene_model_max_y, gene_width, gwr, phenotype_cmap, phenotype_norm, rotate_xy=False):
		"""
		2008-11-25
		"""
		sys.stderr.write("\t Drawing LD info  ...")
		no_of_snps = len(snps_within_this_region.chr_pos_ls)
		left_chr, left_pos = snps_within_this_region.chr_pos_ls[0]
		right_chr, right_pos = snps_within_this_region.chr_pos_ls[-1]
		#ax1.hlines(y_value, left_pos, right_pos, linewidth=0.3)
		for i in range(no_of_snps):
			#draw the SNP line first
			chr_pos1 = snps_within_this_region.chr_pos_ls[i]
			if rotate_xy:
				axe_gene_model.hlines(chr_pos1[1], gene_model_min_y, gene_model_max_y, linestyle='dashed', alpha=0.3, linewidth=0.3)
			else:
				axe_gene_model.vlines(chr_pos1[1], gene_model_min_y, gene_model_max_y, linestyle='dashed', alpha=0.3, linewidth=0.3)
			
			#put the pvalue dot there
			data_obj = gwr.get_data_obj_by_chr_pos(chr_pos1[0], chr_pos1[1])
			if data_obj is not None:
				edgecolor = phenotype_cmap(phenotype_norm(data_obj.value))
				if rotate_xy:
					xpos_ls = [gene_model_max_y-gene_width/4.]	#-gene_width/4. makes it a bit down
					ypos_ls = [chr_pos1[1]]
				else:
					ypos_ls = [gene_model_max_y-gene_width/4.]
					xpos_ls = [chr_pos1[1]]
				axe_gene_model.scatter(xpos_ls, ypos_ls, s=10, linewidth=0.6, edgecolor=edgecolor, facecolor='w')
			
		sys.stderr.write("Done.\n")
	
	def drawOneGenePair(self, snp_pair2value_type, gene_id2snps_id_ls, gene_id_pair, \
						phenotype_method_id, candidate_gene_set, gene_annotation, snp_info, gwr, \
						output_dir, min_value, max_value, snp_region=None, min_distance=40000, list_type_id=None, label_gene=0,
						draw_LD_relative_to_center_SNP=0, commit=0):
		"""
		
		"""
		sys.stderr.write("drawing one pair ...")
		phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
		gene1_id, gene2_id = gene_id_pair
		gene1_model = gene_annotation.gene_id2model.get(gene1_id)
		gene2_model = gene_annotation.gene_id2model.get(gene2_id)
		gene1_symbol = getattr(gene1_model, 'gene_symbol', '')
		gene2_symbol = getattr(gene2_model, 'gene_symbol', '')
		
		#list_type = Stock_250kDB.GeneListType.get(list_type_id)
		fname_basename = 'gene_pair_%s_%s_vs_%s_%s_phenotype_%s_%s'%\
									(gene1_id, gene1_symbol, gene2_id, gene2_symbol, phenotype.short_name, phenotype.id)
		fname_basename = fname_basename.replace('/', '_')
		output_fname_prefix = os.path.join(output_dir, fname_basename)
		
		snps_id_ls1 = gene_id2snps_id_ls[gene1_id]
		snps_id_ls1.sort()
		snps_id_ls2 = gene_id2snps_id_ls[gene2_id]
		snps_id_ls2.sort()
		
		snps_region1 = self.findSNPsInRegion(snp_info, snps_id_ls1[0][0], snps_id_ls1[0][1], snps_id_ls1[-1][1])
		snps_region2 = self.findSNPsInRegion(snp_info, snps_id_ls2[0][0], snps_id_ls2[0][1], snps_id_ls2[-1][1])
		
		pylab.clf()
		
		gap = 0.02
		
		axe_y_offset1 = 0.02	#y_offset for
		axe_height1 = 0.15	#height of axe_gene_model1
		axe_y_offset2 = axe_y_offset1+axe_height1+gap
		axe_height2 = 0.7	#height of axe_gene_model2 and axe_pvalue
		axe_y_offset3 = axe_y_offset2+axe_height2
		
		axe_x_offset1 = 0.02	#
		axe_width1 = 0.15	#width of axe_gene_model2
		axe_x_offset2 = axe_x_offset1 + axe_width1
		axe_width2 = 0.7	#width of axe_pvalue, axe_gene_model1
		axe_x_offset3 = axe_x_offset2 + axe_width2+gap
		axe_width3 = 0.06	#width of axe_score
		axe_x_offset4 = axe_x_offset3 + axe_width3
		no_of_axes_drawn = 0
		
		axe_pvalue = pylab.axes([axe_x_offset2, axe_y_offset2, axe_width2, axe_height2], frameon=False)	#left gap, bottom gap, width, height, axes for pvalue, gene models
		axe_pvalue.grid(True, alpha=0.3)
		axe_pvalue.set_xticklabels([])	#remove xtick labels on axe_pvalue because axe_gene_model1 covers this.
		axe_gene_model1 = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, axe_height1], frameon=False)
		#do not share x-axis with axe_pvalue because no xticks on axe_pvalue, while it's on axe_gene_model1
		#manually set their xlim same
		axe_gene_model1.set_yticks([])
		
		axe_gene_model2 = pylab.axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2], frameon=False, sharey=axe_pvalue)
		axe_gene_model2.set_xticks([])
				
		phenotype_cmap = mpl.cm.jet
		max_phenotype = max_value
		min_phenotype = min_value
		phenotype_gap = max_phenotype - min_phenotype
		phenotype_jitter = phenotype_gap/10.
		phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
		axe_score = pylab.axes([axe_x_offset3, axe_y_offset2, axe_width3, axe_height2/3.], frameon=False)
		cb = mpl.colorbar.ColorbarBase(axe_score, cmap=phenotype_cmap,
									norm=phenotype_norm,
									orientation='vertical')
		cb.set_label('Score')
		axe_score.title.set_text('Score')
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		fig_title = 'Gene %s(%s) vs %s(%s) '%(gene1_symbol, gene1_id, gene2_symbol, gene2_id)
		fig_title += "Phenotype %s (id=%s)."%(phenotype.short_name, phenotype.id)
		axe_pvalue.title.set_text(fig_title)	#main title using this snp.
		
		legend_ls, pscatter_ls = self.drawSNPPairPvalue(axe_pvalue, snps_region1, snps_region2, snp_pair2value_type, phenotype_cmap, phenotype_norm)
		
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		gene_position_cycle = 5
		base_y_value = 1
		gene_width=0.8
		gene_box_text_gap = min_distance*2*0.005
		matrix_of_gene_descriptions1 = self.drawGeneModel(axe_gene_model1, snps_region1, gene_annotation, \
														candidate_gene_set, gene_width=gene_width, \
						gene_position_cycle=gene_position_cycle, base_y_value=base_y_value, gene_box_text_gap=gene_box_text_gap,\
						label_gene=label_gene)
		
		gene_model_min_y = base_y_value-gene_width
		gene_model_max_y = gene_position_cycle + base_y_value -1 + gene_width	#"-1" because genes never sit on y=gene_position_cycle + base_y_value
		self.drawSNPLinesCrossGenesAndPvalue(axe_gene_model1, snps_region1, gene_model_min_y,\
				gene_model_max_y, gene_width, gwr, phenotype_cmap, phenotype_norm)
		
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		matrix_of_gene_descriptions2 = self.drawGeneModel(axe_gene_model2, snps_region2, gene_annotation, \
														candidate_gene_set, gene_width=gene_width, \
						gene_position_cycle=gene_position_cycle, base_y_value=base_y_value, gene_box_text_gap=gene_box_text_gap,\
						label_gene=label_gene, rotate_xy=True)
		
		gene_model_min_y = base_y_value-gene_width
		gene_model_max_y = gene_position_cycle + base_y_value -1 + gene_width	#"-1" because genes never sit on y=gene_position_cycle + base_y_value
		self.drawSNPLinesCrossGenesAndPvalue(axe_gene_model2, snps_region2, gene_model_min_y,\
				gene_model_max_y, gene_width, gwr, phenotype_cmap, phenotype_norm, rotate_xy=True)
		
		axe_gene_model1.set_xlim(axe_pvalue.get_xlim())
		
		axe_cover = pylab.axes([0.03, 0.03, 0.9,0.9], frameon=False)
		axe_cover.set_xticks([])
		axe_cover.set_yticks([])
		axe_cover.legend(pscatter_ls, legend_ls, loc='lower left', handlelen=0.02)
		
		
		png_data = None
		svg_data = None
		png_output_fname = None
		if commit:	#2008-10-24
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=400)
			pylab.savefig(svg_data, format='svg', dpi=300)
		else:
			png_output_fname = '%s.png'%output_fname_prefix
			pylab.savefig(png_output_fname, dpi=500)
			
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		if self.debug:
			pylab.show()
		
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
			
		grand_dataStructure = self.loadDataStructure(self.gene_annotation_picklef, \
							self.min_MAF, self.min_distance, self.list_type_id)
		param_data = grand_dataStructure
		
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		
		gene_id2snps_id_ls = self.get_gene_id2snps_id_ls(snps_context_wrapper)
		del snps_context_wrapper
		
		gene_id_pair_ls = self.read_input_fname(self.input_fname)
		snp_pair2value_type, min_value, max_value = self.get_snp_pair2value_type(self.boolean_pair_fname)
		
		analysis_method_id2gwr = self.getSimilarGWResultsGivenResultsByGene(self.phenotype_method_id, self.call_method_id, \
																		self.results_directory, analysis_method_id_ls=[1])
		gwr = analysis_method_id2gwr.values()[0]
		if gwr.min_value<min_value:
			min_value = gwr.min_value
		if gwr.max_value>max_value:
			max_value = gwr.max_value
		
		for gene_id_pair in gene_id_pair_ls:
			#self.drawOneGenePair(snp_pair2value_type, gene_id2snps_id_ls, gene_id_pair, param_data)
			self.drawOneGenePair(snp_pair2value_type, gene_id2snps_id_ls, gene_id_pair, \
						self.phenotype_method_id, param_data.candidate_gene_set, param_data.gene_annotation, param_data.snp_info, gwr, \
						self.output_dir, min_value, max_value, min_distance=self.min_distance, list_type_id=None, label_gene=1,
						draw_LD_relative_to_center_SNP=0, commit=self.commit)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotGenePairAssoResult
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()