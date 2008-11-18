#!/usr/bin/env python
"""
Examples:
	DrawEcotypeOnMap.py -i /Network/Data/250k/tmp-yh/call_method_17.tsv -n /tmp/phenotype.tsv -o /tmp/ecotype_map_with_phenotype_1 -y 1

	DrawEcotypeOnMap.py  -i /Network/Data/250k/tmp-yh/call_method_17_test.tsv -n /tmp/phenotype.tsv -o /tmp/ecotype_map_with_phenotype_1 -y 1 -u yh
Description:
	Program to plot ecotypes on a map. Color according to phenotype. also linking the map to a plot by PC1 and PC2 from PCA clustering.
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
rcParams['axes.labelsize'] = 8
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6

import time, csv, cPickle
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

class DrawEcotypeOnMap(PlotGroupOfSNPs):
	__doc__ = __doc__
	
	def __init__(self,  **keywords):
		"""
		2008-11-14
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		snpData = SNPData(input_fname=self.input_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1)
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)	#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
		phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, strain_acc_list_phen, phenData.data_matrix)	#tricky, using strain_acc_list_phen
		
		phenotype_col_index = self.findOutWhichPhenotypeColumn(phenData, self.phenotype_method_id)
		
		
		ecotype_info = getEcotypeInfo(db, self.country_order_type)
		
		#the offset below decides where the label of strains/snps should start in axe_snp_matrix
		#2008-11-14 only for PlotGroupOfSNPs.py. you can set it to 1 cuz we dont' draw axe_snp_matrix here.
		snp_id_label_y_offset = 0.95
		StrainID2PCAPosInfo = self.getStrainID2PCAPosInfo(snpData, pca_range=[0,1], snp_id_label_y_offset=snp_id_label_y_offset)
		
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
		max_phenotype = max(phenData.data_matrix[:,phenotype_col_index])
		min_phenotype = min(phenData.data_matrix[:,phenotype_col_index])
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
					draw_lines_to_axe_snp_matrix = False, strain_size_on_axe_strain_pca=14)	#customize a couple of things
		
		axe_strain_pca.set_xlim(axe_strain_pca_xlim)
		axe_strain_pca.set_ylim(axe_strain_pca_ylim)
		axe_strain_map_pca_cover.set_ylim(axe_strain_map_pca_cover_ylim)
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=400)
		pylab.savefig('%s.svg'%self.output_fname_prefix)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawEcotypeOnMap
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()