#!/usr/bin/env python
"""

Examples:
	OutputTestResultInMatrix.py  -l 6 -o /tmp/rank_test_call_method_6_m40000.tsv -m 40000
	
	OutputTestResultInMatrix.py  -l 6,7 -o /tmp/rank_test_call_method_6_7_a.tsv -a
	
	#check TopSNP test results
	OutputTestResultInMatrix.py  -l 17 -x /tmp/top_snp_test_call_method_17_g_f200.png -g  -y 2 -f 200
Description:
	Output Test result in data matrix. either CandidateGeneRankSumTestResult or CandidateGeneTopSNPTest
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, write_data_matrix, SNPData
from TopSNPTest import TopSNPTest
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, \
	CandidateGeneTopSNPTest, CandidateGeneRankSumTestResult, AnalysisMethod, PhenotypeMethod
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages

num = importNumericArray()

class OutputTestResultInMatrix(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
							('max_pvalue_per_gene', 0, int): [0, 'a', 0, 'take the most significant among all SNPs associated with one gene'],\
							('min_sample_size', 0, int): [5, 'i', 1, 'minimum size for both candidate and non-candidate sets to do wilcox.test'],\
							('no_of_top_snps', 1, int): [100, 'f', 1, 'how many number of top snps based on score or -log10(pvalue).'],\
							('call_method_id_ls', 1, ):[None, 'l', 1, 'Restrict results based on these call_methods, coma-separated list of ids.'],\
							('analysis_method_id_ls', 0, ):[None, 'j', 1, 'Restrict results based on these analysis, coma-separated. Default is no such restriction.'],\
							("result_type", 1, int): [1, 'y', 1, 'Which test result to output. 1: CandidateGeneRankSumTestResult, 2: CandidateGeneTopSNPTest'],\
							("super_type_id", 0, int): [1, 'q', 1, 'Super Type ID for candidate gene lists'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [5, 't', 1, 'Number of ticks on the legend'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-24
			split results_method_id_ls if it exists, to accomodate MpiGeneListRankTest which removed this option
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getListTypeInfo(self, db, result_class,  where_condition):
		sys.stderr.write("Getting list type info ...")
		rows = db.metadata.bind.execute("select distinct c.list_type_id, g.biology_category_id from %s order by g.biology_category_id, list_type_id"\
								%(where_condition))
		list_type_id_ls = []
		
		for row in rows:
			list_type_id_ls.append(row.list_type_id)
		sys.stderr.write("Done.\n")
		return list_type_id_ls
			
	
	def getAnalysisMethodInfo(self, db, result_class, where_condition):
		sys.stderr.write("Getting analysis method info ...")
		rows = db.metadata.bind.execute("select distinct r.analysis_method_id from %s order by analysis_method_id"\
								%(where_condition))
		analysis_method_id_ls = []
		for row in rows:
			analysis_method_id_ls.append(row.analysis_method_id)
		sys.stderr.write("Done.\n")
		return analysis_method_id_ls
	
	def getPhenotypeInfo(self, db, result_class,  where_condition):
		sys.stderr.write("Getting phenotype method info ...")
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id, p.biology_category_id from %s p, %s and p.id=r.phenotype_method_id order by p.biology_category_id, r.phenotype_method_id"\
								%(PhenotypeMethod.table.name, where_condition))
		phenotype_method_id_ls = []
		phenotype_method_id2index = {}
		phenotype_method_label_ls = []
		for row in rows:
			phenotype_method_id2index[row.phenotype_method_id] = len(phenotype_method_id_ls)
			phenotype_method_id_ls.append(row.phenotype_method_id)
			pm = PhenotypeMethod.get(row.phenotype_method_id)
			phenotype_method_label_ls.append('%s_%s'%(pm.id, pm.short_name))
		phenotype_info = PassingData()
		phenotype_info.phenotype_method_id2index = phenotype_method_id2index
		phenotype_info.phenotype_method_id_ls = phenotype_method_id_ls
		phenotype_info.phenotype_method_label_ls = phenotype_method_label_ls
		sys.stderr.write("Done.\n")
		return phenotype_info
	
	def orderListTypeAnalysisMethodID(self, list_type_id_ls, analysis_method_id_ls):
		sys.stderr.write("Orderinig list type id and analysis_method id ... ")
		list_type_id_analysis_method_id_ls = []
		list_type_id_analysis_method_id2index = {}
		list_type_analysis_method_label_ls = []
		for list_type_id in list_type_id_ls:
			list_type_short_name = GeneListType.get(list_type_id).short_name
			for analysis_method_id in analysis_method_id_ls:
				analysis_method_short_name = AnalysisMethod.get(analysis_method_id).short_name
				tup = (list_type_id, analysis_method_id)
				list_type_id_analysis_method_id2index[tup] = len(list_type_id_analysis_method_id_ls)
				list_type_id_analysis_method_id_ls.append(tup)
				list_type_analysis_method_label_ls.append('%s_%s_%s'%(analysis_method_short_name, list_type_short_name, list_type_id))
		return_data = PassingData()
		return_data.list_type_id_analysis_method_id_ls = list_type_id_analysis_method_id_ls
		return_data.list_type_id_analysis_method_id2index = list_type_id_analysis_method_id2index
		return_data.list_type_analysis_method_label_ls = list_type_analysis_method_label_ls
		sys.stderr.write("Done.\n")
		return return_data
	
	def get_data_matrix(self, db, phenotype_info, list_type_analysis_method_info, result_class,  where_condition):
		sys.stderr.write("Getting data matrix ...")
		data_matrix = num.zeros([len(list_type_analysis_method_info.list_type_id_analysis_method_id2index), len(phenotype_info.phenotype_method_id2index)], num.float)
		data_matrix[:] = -2
		i = 0
		rows = db.metadata.bind.execute("select r.analysis_method_id, r.phenotype_method_id, c.* from %s order by analysis_method_id"\
								%(where_condition))
		min_value = None
		max_value = None
		for row in rows:
			tup = (row.list_type_id, row.analysis_method_id)
			row_index = list_type_analysis_method_info.list_type_id_analysis_method_id2index[tup]
			col_index = phenotype_info.phenotype_method_id2index[row.phenotype_method_id]
			if row.pvalue>0:
				data_value = -math.log10(row.pvalue)
				if min_value==None:
					min_value = data_value
				elif data_value<min_value:
					min_value = data_value
				
				if max_value==None:
					max_value=data_value
				elif data_value>max_value:
					max_value =data_value
			else:
				data_value = -1	#0 pvalue
			data_matrix[row_index, col_index] = data_value
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		return return_data
	
	max_hue_value = 255	#In Inkscape, the maximum possible hue value, 255, looks almost same as hue=0. cut off before reaching 255.
	#but it's not the case in PIL.
	def value2RGBcolor(cls, value, min_value=0., max_value=255.):
		"""
		2008-08-28
			use Hue-Saturation-Lightness (HSL) color to replace the simple gray gradient represented by (R,G,B)
		2008-08-21
			color span is (0,0,0) to (255,255,255).
		"""
		if value==-1:	#pvalue=0
			return "red"	#(255,0,0)
		elif value==-2:	#NA
			return (255,255,255)	#white #(0,255,0)
		else:
			Y = (value-min_value)/(max_value-min_value)*(cls.max_hue_value-0)
			hue_value = cls.max_hue_value-int(Y)	#the smaller the value is, the higher hue_value is.
			#in (R,G,B) mode, the bigger R/G/B is, the darker the color is
			#R_value = int(Y/math.pow(2,8))
			#G_value = int(Y- R_value*math.pow(2,8))
			return "hsl(%s"%(hue_value)+",100%,50%)"
		
	
	def run(self):	
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		if self.result_type==1:
			result_class = CandidateGeneRankSumTestResult
		elif self.result_type==2:
			result_class = CandidateGeneTopSNPTest
		else:
			sys.stderr.write(" result_type %s not supported.\n"%(self.result_type))
			sys.exit(2)
		
		where_condition = "%s r, %s c, %s g where g.id=c.list_type_id and r.analysis_method_id is not null and r.id=c.results_method_id and c.get_closest=%s and c.min_distance=%s and abs(c.min_MAF-%s)<0.00001"\
				%(ResultsMethod.table.name, result_class.table.name, GeneListType.table.name, self.get_closest, self.min_distance, self.min_MAF)	#the condition for min_MAF is tricky because of the floating precision.
		
		if self.call_method_id_ls:
			where_condition += " and r.call_method_id in (%s)"%self.call_method_id_ls
		
		if self.analysis_method_id_ls:
			where_condition += " and r.analysis_method_id in (%s)"%self.analysis_method_id_ls
		if self.super_type_id:
			where_condition += " and g.super_type_id=%s"%self.super_type_id
		
		if self.result_type==1:
			result_class = CandidateGeneRankSumTestResult
			where_condition += " and c.max_pvalue_per_gene=%s"%(self.max_pvalue_per_gene)
		elif self.result_type==2:
			result_class = CandidateGeneTopSNPTest
			where_condition += " and c.no_of_top_snps=%s"%(self.no_of_top_snps)
		else:
			sys.stderr.write(" result_type %s not supported.\n"%(self.result_type))
			sys.exit(2)
		
		list_type_id_ls = self.getListTypeInfo(db, result_class, where_condition)
		analysis_method_id_ls = self.getAnalysisMethodInfo(db, result_class, where_condition)
		list_type_analysis_method_info = self.orderListTypeAnalysisMethodID(list_type_id_ls, analysis_method_id_ls)
		phenotype_info = self.getPhenotypeInfo(db, result_class, where_condition)
		rdata = self.get_data_matrix(db, phenotype_info, list_type_analysis_method_info, result_class, where_condition)
		
		header = ['list_type_analysis_method', ''] + phenotype_info.phenotype_method_label_ls
		strain_acc_list = list_type_analysis_method_info.list_type_analysis_method_label_ls
		category_list = list_type_analysis_method_info.list_type_id_analysis_method_id_ls
		if SNPData.isDataMatrixEmpty(rdata.data_matrix):
			sys.stderr.write("Nothing fetched from database.\n")
			sys.exit(3)
		if self.output_fname:
			write_data_matrix(rdata.data_matrix, self.output_fname, header, strain_acc_list, category_list)
		
		if self.fig_fname:
			font = get_font(self.font_path, font_size=self.font_size)	#2008-08-01
			value2color_func = lambda x: self.value2RGBcolor(x, rdata.min_value, rdata.max_value)
			im_legend = drawContinousLegend(rdata.min_value, rdata.max_value, self.no_of_ticks, value2color_func, font)
			#im.save('%s_legend.png'%self.fig_fname_prefix)
			im = drawMatrix(rdata.data_matrix, value2color_func, list_type_analysis_method_info.list_type_analysis_method_label_ls,\
						phenotype_info.phenotype_method_label_ls, with_grid=1, font=font)
			im = combineTwoImages(im, im_legend, font=font)
			im.save(self.fig_fname)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = OutputTestResultInMatrix
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()