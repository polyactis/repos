#!/usr/bin/env python
"""
Examples:
	PlotGroupOfSNPs.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -n /tmp/phenotype.tsv -o /tmp/phenotype_1_top_200_snps_a$j -s /Network/Data/A_lyrata/250K.ancestral -y 1 -a $j

Description:
	Program to plot groups of SNPs for one phenotype/one analysis method. both Strains and SNPs go under PCA clustering.
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, SNPData, read_data,\
	assignMatPlotlibHueColorToLs, drawName2FCLegend
import Stock_250kDB, StockDB
from sets import Set
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment
from GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6
from matplotlib.patches import Polygon, CirclePolygon, Ellipse, Wedge
import pylab
import ImageColor
import numpy
from Kruskal_Wallis import Kruskal_Wallis
from PhenotypeOfAncestralDerivedAllele import PhenotypeOfAncestralDerivedAllele
from common import get_chr_id2size, get_ecotypeid2pos, get_ecotypeid2nativename, get_chr_id2size, get_chr_id2cumu_size
from pymodule.SNP import number2complement

class PlotGroupOfSNPs(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'input genotype matrix. Strain X SNP format.', ],\
							('output_fname_prefix', 1, ): ['', 'o', 1, 'store the pvalue', ],\
							('phenotype_fname', 1, ): [None, 'n', 1, 'phenotype file', ],\
							('ancestral_allele_fname', 0, ): [None, 's', 1, 'File containing the ancestral alleles to polarize SNPs', ],\
							('minus_log_pvalue', 0, ): [0, 'e', 0, 'toggle -log(pvalue)', ],\
							('phenotype_method_id', 1, int): [None, 'y', 1, 'which phenotype',],\
							('min_MAF', 1, float): [0.1, 'm', 1, 'minimum Minor Allele Frequency.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('call_method_id', 0, int):[17, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							('no_of_top_hits', 1, int): [200, 'f', 1, 'how many number of top hits based on score or -log(pvalue).'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-09-24
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def drawSNPMtrix(self, axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
					ecotype_info, strain_id_label_x_offset=0.9, snp_id_label_y_offset=0.9):
		"""
		2008-10-07
			1. allele=1 no color (white), 2(red)
			2. snpData.snp_id2index_type. index_type=1(polarized), shape=square. index_type=2(MAF), circle. index_type=3(weirdo), triangle
			3. alpha of each cell determined by the score in top_snp_data.score_ls
			4. label strain, snp plainlys
		"""
		sys.stderr.write("Drawing SNP matrix ...")
		no_of_snps = len(top_snp_data.snp_id_ls)
		max_score = max(top_snp_data.score_ls)
		min_score = min(top_snp_data.score_ls)
		cell_x_len = SNPID2PCAPosInfo.step
		cell_y_len = StrainID2PCAPosInfo.step
		radius = min(cell_x_len, cell_y_len)/2.	#half of the smaller one
		for i in range(no_of_snps):
			snp_id = top_snp_data.snp_id_ls[i]
			#draw snp label
			snp_img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			axe_snp_matrix.text(snp_img_x_pos, snp_id_label_y_offset, snp_id, rotation='vertical', \
							horizontalalignment ='left', verticalalignment='bottom', size=2)
			col_index = subSNPData.col_id2col_index[snp_id]
			index_type = subSNPData.snp_id2index_type[snp_id]
			score = top_snp_data.score_ls[i]
			alpha = (score-min_score)/(max_score-min_score)*(1-0.2)+0.2	#alpha from 0.2 to 1. can't be low.
			for strain_id in StrainID2PCAPosInfo.strain_id_ls:
				strain_img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
				row_index = subSNPData.row_id2row_index[strain_id]
				allele = subSNPData.data_matrix[row_index][col_index]
				if allele==1:
					facecolor = 'w'
				elif allele==2:
					facecolor = 'r'
				else:
					facecolor = 'b'
				if i ==0:	#draw strain label, first SNP
					ecotype_id = int(strain_id)
					nativename = ecotype_info.ecotypeid2nativename[ecotype_id]
					strain_label = '%s_%s'%(strain_id, nativename)
					axe_snp_matrix.text(strain_id_label_x_offset, strain_img_y_pos, strain_label, \
							horizontalalignment ='left', verticalalignment='bottom', size=2)
				center = (snp_img_x_pos+cell_x_len/2., strain_img_y_pos+cell_y_len/2.)
				if index_type==1:
					xs = [snp_img_x_pos, snp_img_x_pos+cell_x_len, snp_img_x_pos+cell_x_len, snp_img_x_pos]
					ys = [strain_img_y_pos, strain_img_y_pos, strain_img_y_pos+cell_y_len, strain_img_y_pos+cell_y_len]
					patch = Polygon(zip(xs,ys), facecolor=facecolor, linewidth=0, alpha=alpha)	#
				elif index_type==2:
					patch = Ellipse(center, cell_x_len, cell_y_len, facecolor=facecolor, linewidth=0, alpha=alpha)
				else:
					patch = Wedge(center, radius, -90, 90, facecolor=facecolor, linewidth=0, alpha=alpha)
				axe_snp_matrix.add_patch(patch)
		
		sys.stderr.write("Done.\n")
			
	def drawStrainPCA(self, axe_strain_pca, axe_strain_pca_legend, StrainID2PCAPosInfo, ecotype_info, rightmost_x_value=1.05):
		"""
		2008-10-07
			1. find all distinct countries
			2. assign colors to the countries and draw a country legend
			3. plot each strain on the axe by pca x/y value, color by its country
			4. draw line linking each strain from axe_snp_matrix to the point in axe_strain_pca
		"""
		sys.stderr.write("Drawing strain PCA ...")
		country_abbr_set = Set()
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			ecotype_id = int(strain_id)
			if ecotype_id in ecotype_info.ecotypeid2country:
				country = ecotype_info.ecotypeid2country[ecotype_id]
			else:
				country = 'Unknown'
			country_abbr_set.add(country)
		country_ls = list(country_abbr_set)
		country_ls.sort()
		country2fc = assignMatPlotlibHueColorToLs(country_ls)
		drawName2FCLegend(axe_strain_pca_legend, country_ls, country2fc, no_face_color=False, font_size=2)
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			ecotype_id = int(strain_id)
			if ecotype_id in ecotype_info.ecotypeid2country:
				country = ecotype_info.ecotypeid2country[ecotype_id]
			else:
				country = 'Unknown'
			x_value = StrainID2PCAPosInfo.strain_id2pca_x[strain_id]
			y_value = StrainID2PCAPosInfo.strain_id2pca_y[strain_id]
			country_fc = country2fc[country]
			
			axe_strain_pca.scatter([x_value],[y_value], s=10, linewidth=0.6, edgecolor=country_fc, facecolor='w', alpha=0.6)
			
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			axe_strain_pca.plot([rightmost_x_value, x_value], [img_y_pos, y_value], '--', alpha=0.2, linewidth=0.3)
		axe_strain_pca.title.set_text('Strain by PC2(x-axis) vs PC1(y-axis)')
		sys.stderr.write("Done.\n")
	
	def drawSNPPCA(self, axe_snp_pca, axe_snp_pca_legend, SNPID2PCAPosInfo, lowest_y_value=-0.05):
		"""
		2008-10-07
		"""
		sys.stderr.write("Drawing SNP PCA ...")
		chr_set = Set()
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr_set.add(chr)
		chr_ls = list(chr_set)
		chr_ls.sort()
		chr2fc = assignMatPlotlibHueColorToLs(chr_ls)
		drawName2FCLegend(axe_snp_pca_legend, chr_ls, chr2fc, no_face_color=True, font_size=4)
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr_fc = chr2fc[chr]
			x_value = SNPID2PCAPosInfo.snp_id2pca_x[snp_id]
			y_value = SNPID2PCAPosInfo.snp_id2pca_y[snp_id]
			axe_snp_pca.scatter([x_value],[y_value], s=10, linewidth=0.6, edgecolor=chr_fc, facecolor='w', alpha=0.6)
			
			img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			axe_snp_pca.plot([img_x_pos, x_value], [lowest_y_value, y_value], linestyle='--', alpha=0.2, linewidth=0.3)
		axe_snp_pca.title.set_text('SNP by PC1(x-axis) vs PC2(y-axis)')
		sys.stderr.write("Done.\n")
	
	def drawMap(self, axe_map, StrainID2PCAPosInfo, phenData, phenotype_col_index, ecotype_info):
		"""
		2008-10-07
			1. draw each phenotype with a bar (temporary solution)
			2. draw map, locate ecotype and color according to phentoype (later)
		"""
		sys.stderr.write("Drawing phenotype and map ...")
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			axe_map.hlines(img_y_pos, 0, phenotype, linewidth=0.3, alpha=0.6)
		axe_map.title.set_text('Phenotype')
		sys.stderr.write("Done.\n")
	
	def scale_chr_size(self, chr_id2size, total_cumu_size, scale_range=[0,1]):
		new_chr_id2size = {}
		for chr_id, size in chr_id2size.iteritems():
			scale_min, scale_max = scale_range
			new_size = size/float(total_cumu_size)*(scale_max-scale_min) + scale_min
			new_chr_id2size[chr_id] = new_size
		return new_chr_id2size
	
	def drawChromosome(self, axe_chromosome, SNPID2PCAPosInfo, chr_id2size, chr_y_value = 0, scale_range=[0,1]):
		"""
		2008-10-07
			lay down five chromosomes
		"""
		sys.stderr.write("Drawing chromosomes ...")
		"""
		chr_size_ls = chr_id2size.values()
		min_chr_size = min(chr_size_ls)
		chr_gap = min_chr_size/3.
		"""
		chr_id2cumu_size, chr_gap, chr_id_ls = get_chr_id2cumu_size(chr_id2size)
		last_chr_id = chr_id_ls[-1]
		
		#rescale everything within [0,1]
		total_cumu_size = float(chr_id2cumu_size[last_chr_id])
		scale_min, scale_max = scale_range
		scale_chr_id2size = self.scale_chr_size(chr_id2size, total_cumu_size, scale_range)
		scale_chr_id2cumu_size = self.scale_chr_size(chr_id2cumu_size, total_cumu_size, scale_range)
		scale_chr_gap = chr_gap/total_cumu_size*(scale_max-scale_min) + scale_min
		
		
		start_pos = 0
		chr_index = 0
		
		axe_chromosome.hlines(chr_y_value, start_pos, scale_chr_id2cumu_size[chr_id_ls[chr_index]])
		
		for chr_index in range(1, len(chr_id_ls)):
			prev_chr_id = chr_id_ls[chr_index-1]
			start_pos = scale_chr_id2cumu_size[prev_chr_id] + scale_chr_gap
			curr_chr_id = chr_id_ls[chr_index]
			stop_pos = scale_chr_id2cumu_size[curr_chr_id]
			axe_chromosome.hlines(chr_y_value, start_pos, stop_pos)
			axe_chromosome.text((start_pos+stop_pos)/2., chr_y_value, 'chr %s'%curr_chr_id, \
							horizontalalignment ='center', verticalalignment='top', size=4)
		
		for snp_id in SNPID2PCAPosInfo.snp_id_ls:
			chr, pos = snp_id.split('_')
			chr = int(chr)
			pos = int(pos)
			img_x_pos = SNPID2PCAPosInfo.snp_id2img_x_pos[snp_id]
			chr_size = float(chr_id2size[chr])
			
			cumu_pos = scale_chr_id2cumu_size[chr]-scale_chr_id2size[chr]+pos/total_cumu_size
			
			axe_chromosome.plot([img_x_pos, cumu_pos], [1, chr_y_value], '--', alpha=0.2, linewidth=0.3)
		sys.stderr.write("Done.\n")
	
	def findOutWhichPhenotypeColumn(self, phenData, phenotype_method_id):
		"""
		2008-10-07
		"""
		col_index_to_return = None
		for col_id, col_index in phenData.col_id2col_index.iteritems():
			col_id_ls = col_id.split('_')
			if int(col_id_ls[0])==phenotype_method_id:
				col_index_to_return = col_index
		return col_index_to_return
		
	def assignAncestralAlleleSmallerNumber(self, snpData, chr_pos2ancestral_allele):
		"""
		2008-10-07
			index_type = 1: Ancestral-allele=1 if it's in chr_pos2ancestral_allele and snpData has that allele.
			index_type = 2: no ancestral in chr_pos2ancestral_allele
			index_type = 3: ancestral exists in chr_pos2ancestral_allele but doesn't appear in snpData
		"""
		snp_id2index_type = {}
		ancestral_allele_used = 1
		for snp_id, snp_allele2index in snpData.snp_id2snp_allele2index.iteritems():
			col_index = snpData.col_id2col_index[snp_id]
			snp_index2allele = snpData.snp_id2snp_allele2index[snp_id]
			if snp_id in chr_pos2ancestral_allele:
				ancestral_allele = chr_pos2ancestral_allele[snp_id]
				ancestral_allele_complement=number2complement[ancestral_allele]	#get the complement
				if ancestral_allele in snp_allele2index:
					if snp_allele2index[ancestral_allele]!=1:	#ancestral allele is not assigned the smallest number. swap. assume only 2 alleles
						#max_index_in_this_col = numpy.max(snpData.data_matrix[:,col_index])
						snpData.data_matrix[:,col_index] = numpy.abs(snpData.data_matrix[:,col_index]-3)	#reverse the index
					index_type = 1
				elif ancestral_allele_complement in snp_allele2index:
					if snp_allele2index[ancestral_allele_complement]!=1:
						snpData.data_matrix[:,col_index] = numpy.abs(snpData.data_matrix[:,col_index]-3)
					index_type = 1
				else:
					index_type = 3
				
			else:
				index_type = 2
			snp_id2index_type[snp_id] = index_type
		snpData.snp_id2index_type = snp_id2index_type
		return snpData
	
	def getSubStrainSNPMatrix(self, snpData, phenData, phenotype_method_id, phenotype_col_index, snp_id_ls, chr_pos2ancestral_allele=None):
		"""
		2008-10-07
			Ancestral-allele=1 if it exists
			Major allele=1 if no ancestral allele info
			remove all NA to do PCA (assume snpData is imputed and has no NA, so use phenData to kick out strain with NA phenotype.)
		"""
		sys.stderr.write("Getting a sub strain-snp matrix ...")
		if phenotype_col_index is None:
			sys.stderr.write("Phenotype %s not in phenotype file.\n"%phenotype_method_id)
			return None
		sub_row_id_ls = []
		
		for row_id, row_index in phenData.row_id2row_index.iteritems():
			if phenData.data_matrix[row_index][phenotype_col_index]!=numpy.nan:	#skip rows with NA
				sub_row_id_ls.append(row_id)
		no_of_sub_rows = len(sub_row_id_ls)
		no_of_sub_snps = len(snp_id_ls)
		sub_matrix = numpy.zeros([no_of_sub_rows, no_of_sub_snps], numpy.int8)
		snp_id2snp_index2allele = {}
		snp_id2snp_allele2index = {}
		for j in range(no_of_sub_snps):
			snp_id = snp_id_ls[j]
			col_index = snpData.col_id2col_index[snp_id]
			snp_allele2index = {}
			snp_allele2count = {}
			snp_index2allele = {}
			for i in range(no_of_sub_rows):
				row_id = sub_row_id_ls[i]
				row_index = snpData.row_id2row_index[row_id]
				allele = snpData.data_matrix[row_index][col_index]
				if allele not in snp_allele2index:	#if not given how to map allele1 and allele2
					snp_allele2index[allele] = len(snp_allele2index)+1
					snp_index = snp_allele2index[allele]
					snp_allele2count[allele] = 0
					snp_index2allele[snp_index] = allele
				#increase the counter
				snp_index = snp_allele2index[allele]
				snp_allele2count[allele] += 1
				sub_matrix[i][j] = snp_index
			if len(snp_allele2count)>2:
				sys.stderr.write("Warnings: This SNP %s has %s (>2) alleles.\n"%(snp_id, len(snp_allele2count)))
			
			#make sure minor allele gets bigger number (assume newly-derived allele)
			if len(snp_allele2count)==2:
				allele1 = snp_index2allele[1]
				allele2 = snp_index2allele[2]
				if snp_allele2count[allele1]<snp_allele2count[allele2]:	#minor allele got assigned the smaller number, reverse it
					sub_matrix[:,j] = numpy.abs(sub_matrix[:,j]-3)	#reverse the index
					snp_allele2index[allele1] = 2
					snp_allele2index[allele2] = 1
					snp_index2allele[2] = allele1
					snp_index2allele[1] = allele2
			
			snp_id2snp_index2allele[snp_id] = snp_index2allele
			snp_id2snp_allele2index[snp_id] = snp_allele2index
		subSNPData = SNPData(col_id_ls=snp_id_ls, row_id_ls=sub_row_id_ls, data_matrix=sub_matrix)
		subSNPData.snp_id2snp_index2allele = snp_id2snp_index2allele
		subSNPData.snp_id2snp_allele2index = snp_id2snp_allele2index
		subSNPData = self.assignAncestralAlleleSmallerNumber(subSNPData, chr_pos2ancestral_allele)
		
		sys.stderr.write("Done.\n")
		return subSNPData
	
	def getTopSNPData(self, genome_wide_result, no_of_top_hits):
		"""
		2008-10-07
		"""
		sys.stderr.write("Getting Top %s SNPs ..."%(no_of_top_hits))
		top_snp_data = PassingData()
		snp_id_ls = []
		score_ls = []
		for i in range(no_of_top_hits):
			data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)	#rank starts from 1
			snp_id_ls.append('%s_%s'%(data_obj.chromosome, data_obj.position))
			score_ls.append(data_obj.value)
		top_snp_data.snp_id_ls = snp_id_ls
		top_snp_data.score_ls = score_ls
		sys.stderr.write("Done.\n")
		return top_snp_data
	
	def sortObjByPCA_value_ls(self, obj_id2index, pca_value_ls, pca_range):
		"""
		2008-10-07
			position on axe_strain_pca or axe_snp_pca, the constraint axis (adjacent to axe_snp_matrix)
		"""
		tuple_to_sort_ls = []
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			tuple_to_sort_ls.append((pca_value, obj_id))
		tuple_to_sort_ls.sort()
		obj_id2pos = {}
		obj_id_ls = []
		pca_value_min = tuple_to_sort_ls[0][0]
		pca_value_max = tuple_to_sort_ls[-1][0]
		for pca_value, obj_id in tuple_to_sort_ls:
			if pca_range:
				pca_range_min, pca_range_max = pca_range
				pos = ((pca_value-pca_value_min)/(pca_value_max-pca_value_min))*(pca_range_max-pca_range_min)+pca_range_min
			else:
				pos = pca_value
			obj_id_ls.append(obj_id)
			obj_id2pos[obj_id] = pos
		return obj_id_ls, obj_id2pos
	
	def getObj2posFromPCA_value_ls(self, obj_id2index, pca_value_ls, pca_range):
		"""
		2008-10-07
			position on axe_strain_pca or axe_snp_pca, the non-constraint axis (not adjacent to axe_snp_matrix)
		"""
		obj_id2pos = {}
		tuple_to_sort_ls = []
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			tuple_to_sort_ls.append((pca_value, obj_id))
		tuple_to_sort_ls.sort()
		pca_value_min = tuple_to_sort_ls[0][0]
		pca_value_max = tuple_to_sort_ls[-1][0]
		
		for obj_id, index in obj_id2index.iteritems():
			pca_value = pca_value_ls[index]
			if pca_range:
				pca_range_min, pca_range_max = pca_range
				pos = ((pca_value-pca_value_min)/(pca_value_max-pca_value_min))*(pca_range_max-pca_range_min)+pca_range_min
			else:
				pos = pca_value
			obj_id2pos[obj_id] = pos
		return 	obj_id2pos
	
	def getObj2ImgPos(self, obj_id_ls, img_pos_range=[0., 0.9]):
		"""
		2008-10-07
			position on axe_snp_matrix (either x or y axis)
		"""
		obj_id2img_pos = {}
		no_of_objs = len(obj_id_ls)
		step = (img_pos_range[1]-img_pos_range[0])/no_of_objs
		for i in range(no_of_objs):
			obj_id = obj_id_ls[i]
			obj_id2img_pos[obj_id] = i*step
		return obj_id2img_pos, step
			
	
	def getStrainID2PCAPosInfo(self, subSNPData, pca_range=[0,1]):
		"""
		2008-10-07
		"""
		sys.stderr.write("Getting  StrainID2PCAPosInfo ... ")
		StrainID2PCAPosInfo = PassingData()
		import pca_module
		T, P, explained_var = pca_module.PCA_svd(subSNPData.data_matrix, standardize=True)
		#T[:,1] and T[:,2] are new positions for strain
		strain_id_ls, strain_id2pca_y = self.sortObjByPCA_value_ls(subSNPData.row_id2row_index, T[:,1], pca_range)
		strain_id2pca_x = self.getObj2posFromPCA_value_ls(subSNPData.row_id2row_index, T[:,2], pca_range)
		StrainID2PCAPosInfo.strain_id_ls = strain_id_ls
		StrainID2PCAPosInfo.strain_id2pca_y = strain_id2pca_y
		StrainID2PCAPosInfo.strain_id2pca_x = strain_id2pca_x
		obj_id2img_pos, step = self.getObj2ImgPos(strain_id_ls, img_pos_range=[0., 0.9])
		StrainID2PCAPosInfo.strain_id2img_y_pos = obj_id2img_pos
		StrainID2PCAPosInfo.step = step
		sys.stderr.write("Done.\n")
		return StrainID2PCAPosInfo
	
	def getSNPID2PCAPosInfo(self, subSNPData, pca_range=[0,1]):
		sys.stderr.write("Getting  SNPID2PCAPosInfo ... ")
		SNPID2PCAPosInfo = PassingData()
		import pca_module
		T, P, explained_var = pca_module.PCA_svd(numpy.transpose(subSNPData.data_matrix), standardize=True)
		#T[:,1] and T[:,2] are new positions for strain
		snp_id_ls, snp_id2pca_x = self.sortObjByPCA_value_ls(subSNPData.col_id2col_index, T[:,1], pca_range)
		snp_id2pca_y = self.getObj2posFromPCA_value_ls(subSNPData.col_id2col_index, T[:,2], pca_range)
		SNPID2PCAPosInfo.snp_id_ls = snp_id_ls
		SNPID2PCAPosInfo.snp_id2pca_x = snp_id2pca_x
		SNPID2PCAPosInfo.snp_id2pca_y = snp_id2pca_y
		obj_id2img_pos, step = self.getObj2ImgPos(snp_id_ls, img_pos_range=[0., 0.9])
		SNPID2PCAPosInfo.snp_id2img_x_pos = obj_id2img_pos
		SNPID2PCAPosInfo.step = step
		sys.stderr.write("Done.\n")
		return SNPID2PCAPosInfo
	
	def getEcotypeInfo(self, db):
		"""
		2008-10-07
		"""
		sys.stderr.write("Getting  Ecotype info ... ")
		ecotype_info = PassingData()
		rows = db.metadata.bind.execute("select e.id, e.nativename, e.latitude, e.longitude, c.abbr from stock.%s e, stock.%s s, stock.%s a, stock.%s c \
				where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id"%(StockDB.Ecotype.table.name, \
				StockDB.Site.table.name, StockDB.Address.table.name, StockDB.Country.table.name))
		ecotypeid2pos = {}
		ecotypeid2nativename = {}
		ecotypeid2country = {}
		for row in rows:
			ecotypeid2pos[row.id] = (row.latitude, row.longitude)
			ecotypeid2nativename[row.id] = row.nativename
			ecotypeid2country[row.id] = row.abbr
		ecotype_info.ecotypeid2pos = ecotypeid2pos
		ecotype_info.ecotypeid2nativename = ecotypeid2nativename
		ecotype_info.ecotypeid2country = ecotypeid2country
		sys.stderr.write("Done.\n")
		return ecotype_info
	
	def run(self):
		"""
		2008-10-07
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		rm = Stock_250kDB.ResultsMethod.query.filter_by(phenotype_method_id=self.phenotype_method_id, \
													call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id).first()
		
		if rm:
			genome_wide_result = self.getResultMethodContent(rm, self.results_directory, min_MAF=0, construct_chr_pos2index=True)
		else:
			sys.stderr.write("No genome association results for phenotype_method_id=%s, call_method_id=%s, analysis_method_id=%s.\n"%\
							(self.phenotype_method_id, self.call_method_id, self.analysis_method_id))
			sys.exit(3)
		
		snpData = SNPData(input_fname=self.input_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1)
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)	#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
		phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, phenData.row_id_ls, phenData.data_matrix)
		if self.ancestral_allele_fname:
			chr_pos2ancestral_allele = PhenotypeOfAncestralDerivedAllele.get_chr_pos2ancestral_allele(self.ancestral_allele_fname)
		else:
			chr_pos2ancestral_allele = {}
		
		phenotype_col_index = self.findOutWhichPhenotypeColumn(phenData, self.phenotype_method_id)
		
		top_snp_data = self.getTopSNPData(genome_wide_result, self.no_of_top_hits)
		subSNPData = self.getSubStrainSNPMatrix(snpData, phenData, self.phenotype_method_id, phenotype_col_index, top_snp_data.snp_id_ls, chr_pos2ancestral_allele)
		
		"""
		ecotypeid2pos = 
		ecotypeid2nativename =
		ecotypeid2country = 
		"""
		ecotype_info = self.getEcotypeInfo(db)
		chr_id2size = get_chr_id2size(db.metadata.bind)
		strain_id_label_x_offset=0.9	#not used yet
		snp_id_label_y_offset=0.9
		
		StrainID2PCAPosInfo = self.getStrainID2PCAPosInfo(subSNPData, pca_range=[0,1])
		SNPID2PCAPosInfo = self.getSNPID2PCAPosInfo(subSNPData, pca_range=[0,1])
		
		no_of_axes_drawn = 0
		
		axe_chromosome = pylab.axes([0.3, 0.05, 0.4, 0.2], frameon=False)
		axe_chromosome.set_xticks([])
		axe_chromosome.set_yticks([])
		self.drawChromosome(axe_chromosome, SNPID2PCAPosInfo, chr_id2size, chr_y_value = 0, scale_range=[0,1])
		#axe_chromosome.set_xlim([0,1])
		#axe_chromosome.set_ylim([0,1])
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		axe_strain_pca = pylab.axes([0.0, 0.25, 0.3, 0.4], frameon=False)
		axe_strain_pca.grid(True, alpha=0.3)
		axe_strain_pca.set_xticks([])
		axe_strain_pca.set_yticks([])
		axe_strain_pca_legend = pylab.axes([0.1, 0.7, 0.1, 0.3], frameon=False)
		axe_strain_pca_legend.set_xticks([])
		axe_strain_pca_legend.set_yticks([])
		axe_strain_pca_xlim = [-0.05,1.05]
		self.drawStrainPCA(axe_strain_pca, axe_strain_pca_legend, StrainID2PCAPosInfo, ecotype_info, rightmost_x_value=axe_strain_pca_xlim[1])
		axe_strain_pca.set_xlim(axe_strain_pca_xlim)
		axe_strain_pca.set_ylim([0,1])
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		axe_snp_matrix = pylab.axes([0.3, 0.25, 0.4, 0.4], frameon=False, sharex=axe_chromosome, sharey=axe_strain_pca)
		axe_snp_matrix.set_xticks([])
		axe_snp_matrix.set_yticks([])
		self.drawSNPMtrix(axe_snp_matrix, subSNPData, top_snp_data, StrainID2PCAPosInfo, SNPID2PCAPosInfo, \
						ecotype_info, strain_id_label_x_offset, snp_id_label_y_offset)
		axe_snp_matrix.set_xlim([0,1])
		#axe_snp_matrix.set_ylim([0,1])
		no_of_axes_drawn += 1
		#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		axe_map = pylab.axes([0.7, 0.25, 0.28, 0.4], frameon=False, sharey=axe_snp_matrix)
		axe_map.set_yticks([])
		self.drawMap(axe_map, StrainID2PCAPosInfo, phenData, phenotype_col_index, ecotype_info)
		no_of_axes_drawn += 1
		axe_map.set_ylim([0,1])	#without this, ylim of all 3 axes are set to [0,0.9] because axe_map automatically adjust to 0-0.9
		#pylab.savefig('%s_%s.png'%(self.output_fname_prefix, no_of_axes_drawn), dpi=400)
		
		axe_snp_pca = pylab.axes([0.3, 0.65, 0.4, 0.3], frameon=False, sharex=axe_snp_matrix)
		axe_snp_pca.grid(True, alpha=0.3)
		axe_snp_pca.set_xticks([])
		axe_snp_pca.set_yticks([])
		axe_snp_pca_legend = pylab.axes([0.8, 0.8, 0.1, 0.15], frameon=False)
		axe_snp_pca_legend.set_xticks([])
		axe_snp_pca_legend.set_yticks([])
		axe_snp_pca_ylim = [-0.05,1.05]
		self.drawSNPPCA(axe_snp_pca, axe_snp_pca_legend, SNPID2PCAPosInfo, lowest_y_value=axe_snp_pca_ylim[0])
		axe_snp_pca.set_ylim(axe_snp_pca_ylim)
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=600)
		#pylab.savefig('%s.svg'%self.output_fname_prefix)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PlotGroupOfSNPs
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()