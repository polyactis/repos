#!/usr/bin/env python
"""
Examples:
	ContrastTwoPhenotypes.py -n /tmp/phenotype.tsv -o /tmp/phenotype_2_vs_4 -x 2 -y4

Description:
	Program to contrast phenotype 1 vs phenotype 2. To see how correlated they are and which accessions are not.
	
	todo: add a map showing the geographic source of accessions.
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
rcParams['axes.labelsize'] = 10
rcParams['axes.titlesize'] = 12
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8

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

class ContrastTwoPhenotypes(PlotGroupOfSNPs):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('output_fname_prefix', 1, ): ['', 'o', 1, 'store the pvalue', ],\
							('phenotype_fname', 1, ): [None, 'n', 1, 'phenotype file', ],\
							('phenotype_method_id1', 1, int): [None, 'x', 1, 'phenotype on the X-axis',],\
							('phenotype_method_id2', 1, int): [None, 'y', 1, 'phenotype on the Y-axis',],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-11-30
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
		
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, data_matrix=data_matrix_phen)	#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
		phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(phenData.row_id_ls, strain_acc_list_phen, phenData.data_matrix)	#tricky, using strain_acc_list_phen
		
		phenotype_col_index1 = self.findOutWhichPhenotypeColumn(phenData, Set([self.phenotype_method_id1]))[0]
		phenotype_col_index2 = self.findOutWhichPhenotypeColumn(phenData, Set([self.phenotype_method_id2]))[0]
		
		x_ls = []
		y_ls = []
		for i in range(phenData.data_matrix.shape[0]):
			if not numpy.isnan(phenData.data_matrix[i][phenotype_col_index1]) and not numpy.isnan(phenData.data_matrix[i][phenotype_col_index2]):
				x_ls.append(phenData.data_matrix[i][phenotype_col_index1])
				y_ls.append(phenData.data_matrix[i][phenotype_col_index2])
		
		pylab.clf()
		pylab.title('Phenotype Contrast')
		pylab.plot(x_ls, y_ls, '.', alpha=0.6)
		pylab.grid(alpha=0.3)
		phenotype_method1 = Stock_250kDB.PhenotypeMethod.get(self.phenotype_method_id1)
		phenotype_method2 = Stock_250kDB.PhenotypeMethod.get(self.phenotype_method_id2)
		pylab.xlabel(phenotype_method1.short_name)
		pylab.ylabel(phenotype_method2.short_name)
		
		#draw diagonal line to show perfect correlation
		max_min_value = max(min(x_ls), min(y_ls))
		min_max_value = min(max(x_ls), max(y_ls))
		pylab.plot([max_min_value, min_max_value], [max_min_value, min_max_value], c='g', alpha=0.7)
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=400)
		pylab.savefig('%s.svg'%self.output_fname_prefix)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ContrastTwoPhenotypes
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()