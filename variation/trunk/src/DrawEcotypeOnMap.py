#!/usr/bin/env python
"""
Examples:
	DrawEcotypeOnMap.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -n /tmp/phenotype.tsv -o /tmp/ecotype_map_with_phenotype_1 -y 1

	DrawEcotypeOnMap.py  -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -n /tmp/phenotype.tsv -o /tmp/ecotype_map_with_phenotype_1 -y 1 -u yh
	
	#use PCA results (eigen values, principal components from smartpca.perl of EIGENSOFT)
	#the input SNP matrix is good as far as the order of strains match the PCA result files.
	DrawEcotypeOnMap.py -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -n /Network/Data/250k/tmp-yh/phenotype.tsv  -o ./banyan_fs/tmp/ecotype_map_with_phenotype_1_smartpca -y 1 -e /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.eval -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec
	
	#ditto, but covering whole world
	DrawEcotypeOnMap.py -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -n /Network/Data/250k/tmp-yh/phenotype.tsv  -o ./banyan_fs/tmp/ecotype_map_with_phenotype_1_smartpca_world -y 1 -e /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.eval -f /Network/Data/250k/tmp-yh/eigenstrat/call_method_17_eigenstrat.pca.evec -u yh -p yh324 -j -125,-38,180,66
	
	#ditto but with 149SNP data
	DrawEcotypeOnMap.py -i /Network/Data/250k/db/reference_dataset/stock_149SNP_y0000110101_mergedup.csv -n /Network/Data/250k/tmp-yh/phenotype.tsv  -o /Network/Data/250k/tmp-yh/eigenstrat/stock149_ecotype_map_with_phenotype_1_smartpca_world -y 1  -u yh -p yh324 -j -125,-38,180,66
	
Description:
	Program to plot ecotypes on a map. Color according to phenotype. also linking the map to a plot by PC1 and PC2 from PCA clustering.
	

	2008-12-08 add another output file (name like '%s_%s'%(output_fname_prefix, 'LatLonPhenVsPC') ) which
		make all kinds of two-variable (longitude vs PC1 or latitude vs PC2 or phenotype vs longitude ...) plots with linear model specifics
		every plot comprises one subplot in a giant plot.
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import matplotlib as mpl; mpl.use("Agg")	#to avoid popup and collapse in X11-disabled environment
from matplotlib import rcParams
from PlotGroupOfSNPs import PlotGroupOfSNPs
rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 8
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 6
rcParams['xtick.labelsize'] = 4
rcParams['ytick.labelsize'] = 4

import time, csv, cPickle, random
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, SNPData, read_data,\
	assignMatPlotlibHueColorToLs, drawName2FCLegend
from pymodule.SNP import number2complement
import Stock_250kDB, StockDB
from sets import Set
from GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()

import pylab
import ImageColor
import numpy
from Kruskal_Wallis import Kruskal_Wallis
from PhenotypeOfAncestralDerivedAllele import PhenotypeOfAncestralDerivedAllele
from common import get_chr_id2size, get_chr_id2cumu_size, getEcotypeInfo
from Association import Association

class DrawEcotypeOnMap(PlotGroupOfSNPs, Association):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'input genotype matrix. Strain X SNP format.', ],\
							('output_fname_prefix', 1, ): ['', 'o', 1, 'output filename prefix, no suffix needs to be attached.', ],\
							('phenotype_fname', 1, ): [None, 'n', 1, 'phenotype file', ],\
							('phenotype_method_id', 1, int): [None, 'y', 1, 'which phenotype',],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('call_method_id', 0, int):[17, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							('country_order_type', 1, int): [1, 'g', 1, 'How to order countries from where strains are from. 1: order by latitude. 2: by longitude.'],\
							('eigen_value_fname', 0, ): [None, 'e', 1, 'eigen value file outputted by smartpca.perl from EIGENSOFT', ],\
							('eigen_vector_fname', 0, ): [None, 'f', 1, 'eigen vector file outputted by smartpca.perl from EIGENSOFT', ],\
							('pic_area', 1, ): ['-15,30,38,66', 'j', 1, 'left longitude, bottom latitude, right longitude, upper latitude. #pic_area=[-15,30,38,66]	#[-125,-38,180,66] covers everything, [-125,10,90,66] covers US and europe.',],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-11-14
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.pic_area =  self.pic_area.split(',')
		self.pic_area = map(int, self.pic_area)
	
	
	def subplotLatLonPhenVsPC(self, ecotype_info, StrainID2PCAPosInfo, phenData, phenotype_col_index, \
							phenotype_cmap, phenotype_norm, which_figure=0, sub_title='',\
							ydata='latitude', which_PC_index=0, no_of_rows=2, dot_size=10, alpha=0.7):
		"""
		2008-12-08
			one single subplot in plotLatLonPhenVsPC()
		"""
		ax = pylab.subplot(no_of_rows,2,which_figure, frameon=False)
		pylab.ylabel(ydata)
		pylab.grid(True, alpha=0.3)
		
		x_ls = []
		y_ls = []
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			ecotype_id = int(strain_id)
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj:
				lat, lon = ecotype_obj.latitude, ecotype_obj.longitude
				y_value = getattr(ecotype_obj, ydata, None)
			else:
				sys.stderr.write("Warning: Ecotype %s not in ecotype_info (fetched from stock db).\n"%ecotype_id)
				continue
			
			#x_value = StrainID2PCAPosInfo.PC_matrix[row_index][which_PC_index]
			if which_PC_index==1:
				x_value = StrainID2PCAPosInfo.strain_id2pca_x[strain_id]
				xlabel = 'PC%s'%(which_PC_index+1)
			elif which_PC_index==0:
				x_value = StrainID2PCAPosInfo.strain_id2pca_y[strain_id]
				xlabel = 'PC%s'%(which_PC_index+1)
			elif which_PC_index=='longitude':
				x_value = lon
				xlabel = which_PC_index
			elif which_PC_index =='latitude':
				x_value = lat
				xlabel = which_PC_index
			
			#img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]

			#strain color according to phenotype
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			strain_fc = phenotype_cmap(phenotype_norm(phenotype))
			if numpy.isnan(phenotype):
				linewidth=0.5
				strain_fc = 'w'
				edgecolor = 'k'
				_alpha = 0	#facecolor gets very transparent
			else:
				linewidth=0
				strain_fc = strain_fc
				edgecolor = 'k'
				_alpha = alpha
				
			if ydata=='phenotype':
				y_value = phenotype
				if numpy.isnan(phenotype):	#can't do regression or plot
					continue
			if y_value is None or numpy.isnan(y_value):
				continue
			pylab.scatter([x_value],[y_value], s=dot_size, linewidth=linewidth, facecolor=strain_fc, alpha=_alpha, zorder=10)
			x_ls.append(x_value)
			y_ls.append(y_value)
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()
		
		lm_result = Association.pure_linear_model(x_ls, y_ls)
		try:
			pvalue = '%.2f'%-math.log10(lm_result.pvalue)
		except:
			pvalue = '%s'%lm_result.pvalue
		beta0 = lm_result.coeff_list[0]
		beta = lm_result.coeff_list[1]
		#draw a line showing the trend
		x_lm_ls = [min(x_ls), max(x_ls)]
		lm_func = lambda x: lm_result.coeff_list[0]+lm_result.coeff_list[1]*x
		y_lm_ls = map(lm_func, x_lm_ls)
		ax.plot(x_lm_ls, y_lm_ls, alpha=0.5)
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		
		pylab.xlabel('%s, pvalue=%s, beta=%.3f'%(xlabel, pvalue, beta))
		#pylab.title('%s, pvalue=%.2f, beta=%.3f'%(sub_title, pvalue, beta))
	
	def plotLatLonPhenVsPC(self, ecotype_info, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_cmap, phenotype_norm, 
						output_fname_prefix=None, commit=0):
		"""
		2008-12-08
			make all kinds of two-variable (longitude vs PC1 or latitude vs PC2 or phenotype vs longitude ...) plots with linear model specifics
			
			all comprises one subplot in a giant plot
		"""
		sys.stderr.write("Drawing Latitude/Longitude/Phenotype vs PCs ...")
		pylab.clf()
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		pylab.subplots_adjust(left=0.05, right=0.95, bottom = 0.05, top=0.95, wspace=0.2, hspace=0.25)	#wspace is the amount of width reserved for blank space between subplots
		
		ydata_which_PC_index_ls = [('longitude', 0), ('longitude', 1), ('latitude', 0), ('latitude', 1), ('phenotype', 0), ('phenotype',1),
								('phenotype', 'longitude'), ('phenotype','latitude')]
		no_of_rows = len(ydata_which_PC_index_ls)/2.
		if no_of_rows%1>0:
			no_of_rows = int(no_of_rows)+1
		else:
			no_of_rows = int(no_of_rows)
		
		for i in range(len(ydata_which_PC_index_ls)):
			ydata, which_PC_index = ydata_which_PC_index_ls[i]
			self.subplotLatLonPhenVsPC(ecotype_info, StrainID2PCAPosInfo, phenData, phenotype_col_index, \
							phenotype_cmap, phenotype_norm, which_figure=i+1, sub_title='',\
							ydata=ydata, which_PC_index=which_PC_index, no_of_rows=no_of_rows)
		
		ax = pylab.axes([0.1, 0.1, 0.8,0.8], frameon=False)
		ax.set_xticks([])
		ax.set_yticks([])
		#title = 'Phenotype %s %s by %s %s'%(phenotype_method.id, phenotype_method.short_name, list_type.id, \
		#											list_type.short_name)
		#ax.set_title(title)
		png_data = None
		svg_data = None
		if commit:
			png_data = StringIO.StringIO()
			svg_data = StringIO.StringIO()
			pylab.savefig(png_data, format='png', dpi=300)
			pylab.savefig(svg_data, format='svg', dpi=300)
		elif output_fname_prefix:
			#output_fname_prefix = os.path.join(output_dir, title.replace('/', '_'))
			output_fname_prefix = '%s_%s'%(output_fname_prefix, 'LatLonPhenVsPC')
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		return png_data, svg_data
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		snpData = SNPData(input_fname=self.input_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1)
		
		
		if self.eigen_vector_fname and self.eigen_value_fname:
			eigen_value_ls = self.getEigenValueFromFile(self.eigen_value_fname)
			eigen_value_ls = numpy.array(eigen_value_ls)
			explained_var = eigen_value_ls/numpy.sum(eigen_value_ls)
			PC_data = self.getPCFromFile(self.eigen_vector_fname)
			PC_matrix = PC_data.PC_matrix
		else:
			max_no_of_snps = 10000
			if len(snpData.col_id_ls)>max_no_of_snps:	#2008-12-01 randomly pick max_no_of_snps SNPs
				picked_col_index_ls = random.sample(range(len(snpData.col_id_ls)), max_no_of_snps)
				new_col_id_ls = [snpData.col_id_ls[i] for i in picked_col_index_ls]
				newSnpData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=new_col_id_ls, strain_acc_list=snpData.strain_acc_list,\
								category_list=snpData.category_list)
				newSnpData.data_matrix = snpData.data_matrix[:, picked_col_index_ls]
				snpData = newSnpData
		
			snpData, allele_index2allele_ls = snpData.convertSNPAllele2Index()
			explained_var = None
			PC_matrix = None
		
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)	#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
		phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, strain_acc_list_phen, phenData.data_matrix)	#tricky, using strain_acc_list_phen
		
		phenotype_col_index = self.findOutWhichPhenotypeColumn(phenData, Set([self.phenotype_method_id]))[0]
		
		
		ecotype_info = getEcotypeInfo(db, self.country_order_type)
		
		#the offset below decides where the label of strains/snps should start in axe_snp_matrix
		#2008-11-14 only for PlotGroupOfSNPs.py. you can set it to 1 cuz we dont' draw axe_snp_matrix here.
		snp_id_label_y_offset = 0.95
		StrainID2PCAPosInfo = self.getStrainID2PCAPosInfo(snpData, pca_range=[0,1], snp_id_label_y_offset=snp_id_label_y_offset, explained_var=explained_var, T=PC_matrix)
		
		axe_y_offset1 = 0.03
		axe_height1 = 0.45	#height of axe_chromosome, twice height of axe_map_phenotype_legend
		axe_y_offset2 = axe_y_offset1+axe_height1
		axe_height2 = 0.5	#height of axe_strain_pca, axe_snp_matrix, axe_map
		axe_y_offset3 = axe_y_offset2+axe_height2
		
		axe_x_offset1 = 0.05
		axe_width1 = 0.8	#width of axe_strain_pca
		axe_x_offset2 = axe_x_offset1 + 0.02 + axe_width1
		axe_width2 = 0.05	#width of axe_chromosome, axe_snp_matrix, axe_snp_pca
		axe_x_offset3 = axe_x_offset2 + axe_width2
		axe_width3 = 0.02	#width of axe_phenotype
		
		phenotype_method = Stock_250kDB.PhenotypeMethod.get(self.phenotype_method_id)
		
		phenotype_cmap = mpl.cm.jet
		max_phenotype = numpy.nanmax(phenData.data_matrix[:,phenotype_col_index])	#nanmax ignores the nan elements
		min_phenotype = numpy.nanmin(phenData.data_matrix[:,phenotype_col_index])	#nanmin ignores the nan elements
		phenotype_gap = max_phenotype - min_phenotype
		phenotype_jitter = phenotype_gap/10.
		phenotype_norm = mpl.colors.Normalize(vmin=min_phenotype-phenotype_jitter, vmax=max_phenotype+phenotype_jitter)
		axe_map_phenotype_legend = pylab.axes([axe_x_offset2, axe_y_offset1, axe_width2, 0.3], frameon=False)
		cb = mpl.colorbar.ColorbarBase(axe_map_phenotype_legend, cmap=phenotype_cmap,
									norm=phenotype_norm,
									orientation='vertical')
		cb.set_label('Legend Of Phenotype %s %s'%(phenotype_method.id, phenotype_method.short_name))
		
		axe_strain_map = pylab.axes([axe_x_offset1, axe_y_offset2, axe_width1, axe_height2], frameon=False)
		axe_strain_pca = pylab.axes([axe_x_offset1, axe_y_offset1, axe_width1, axe_height1], frameon=False)
		axe_strain_map_pca_cover = pylab.axes([axe_x_offset1, axe_y_offset1, axe_width1, axe_height1+axe_height2], frameon=False, \
											sharex=axe_strain_pca)	#cover both axe_strain_map and axe_strain_pca
		axe_strain_map_pca_cover.set_yticks([])
		axe_strain_pca_xlim = [-0.05,1.05]
		axe_strain_pca_ylim = [0, 1.05]
		axe_strain_pca.set_xlim(axe_strain_pca_xlim)
		axe_strain_pca.set_ylim(axe_strain_pca_ylim)
		axe_strain_map_pca_cover_ylim = [0, (axe_height1+axe_height2)/axe_height1]	#set it accordingly
		axe_strain_map_pca_cover.set_ylim(axe_strain_map_pca_cover_ylim)
				
		axe_strain_pca.grid(True, alpha=0.3)
		axe_strain_pca.set_xticks([])
		axe_strain_pca.set_yticks([])
		axe_strain_pca_legend = None	#no pca legend
		self.drawStrainPCA(axe_strain_pca, axe_strain_map, axe_strain_map_pca_cover, axe_strain_pca_legend, StrainID2PCAPosInfo, \
						ecotype_info, phenData, \
					phenotype_col_index, phenotype_cmap, phenotype_norm, rightmost_x_value=axe_strain_pca_xlim[1],\
					strain_color_type=2, pca2map_line_color=None, ecotype_width_on_map=10,\
					draw_lines_to_axe_snp_matrix = False, strain_size_on_axe_strain_pca=14, pic_area=self.pic_area,\
					map_pca_line_alpha=0.2, map_pca_linewidth=0.2)	#customize a couple of things
		
		axe_strain_pca.set_xlim(axe_strain_pca_xlim)
		axe_strain_pca.set_ylim(axe_strain_pca_ylim)
		axe_strain_map_pca_cover.set_ylim(axe_strain_map_pca_cover_ylim)
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=400)
		pylab.savefig('%s.svg'%self.output_fname_prefix)
		
		self.plotLatLonPhenVsPC(ecotype_info, StrainID2PCAPosInfo, phenData, phenotype_col_index, phenotype_cmap, phenotype_norm, 
						self.output_fname_prefix, commit=self.commit)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawEcotypeOnMap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()