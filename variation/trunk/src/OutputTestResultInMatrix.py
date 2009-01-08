#!/usr/bin/env python
"""

Examples:
	OutputTestResultInMatrix.py  -l 6 -o /tmp/rank_test_call_method_6_m40000.tsv -m 40000
	
	OutputTestResultInMatrix.py  -l 6,7 -o /tmp/rank_test_call_method_6_7_a.tsv -a
	
	#check TopSNP test results
	OutputTestResultInMatrix.py  -l 17 -x /tmp/top_snp_test_call_method_17_g_f200.png -g  -y 2 -f 200
	
Description:
	Output Test result in data matrix. either CandidateGeneRankSumTestResult or CandidateGeneTopSNPTest
	In the outputted matrix, -1 = NA. -2 = (pvalue=0).
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, write_data_matrix, SNPData
import Stock_250kDB
from Stock_250kDB import Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, \
	CandidateGeneTopSNPTest, CandidateGeneRankSumTestResult, AnalysisMethod, PhenotypeMethod, ResultsByGene
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color

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
							('no_of_top_snps', 1, int): [100, 'f', 1, 'For test_result_type=2. how many number of top snps based on score or -log10(pvalue).'],\
							('call_method_id_ls', 1, ):[None, 'l', 1, 'Restrict results based on these call_methods, coma-separated list of ids.'],\
							('analysis_method_id_ls', 0, ):[None, 'j', 1, 'Restrict results based on these analysis, coma-separated. Default is no such restriction.'],\
							("test_result_type", 1, int): [1, 'y', 1, 'Which test result to output. 1: CandidateGeneRankSumTestResult, 2: CandidateGeneTopSNPTest, 3: CandidateGeneRankSumTestResultMethod'],\
							("super_type_id", 0, int): [1, 'q', 1, 'Super Type ID for candidate gene lists'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [15, 't', 1, 'Number of ticks on the legend'],\
							("test_type", 1, int): [0, 'w', 1, 'which type of test in each result_type. column test_type. Check GeneListRankTest.py or TopSNPTest.py for info.'],\
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
	
	def getListTypeInfo(self, db,  where_condition):
		"""
		2008-08-29
			add -1 as a separator into list_type_id_ls
		"""
		sys.stderr.write("Getting list type info ...")
		rows = db.metadata.bind.execute("select distinct c.list_type_id, g.biology_category_id from %s order by g.biology_category_id, list_type_id"\
								%(where_condition))
		list_type_id_ls = []
		prev_biology_category_id = None
		for row in rows:
			if prev_biology_category_id == None:
				prev_biology_category_id = row.biology_category_id
			elif row.biology_category_id!=prev_biology_category_id:
				prev_biology_category_id = row.biology_category_id
				list_type_id_ls.append(-1)
			list_type_id_ls.append(row.list_type_id)
		sys.stderr.write("Done.\n")
		return list_type_id_ls
			
	
	def getAnalysisMethodInfo(self, db, where_condition):
		sys.stderr.write("Getting analysis method info ...")
		rows = db.metadata.bind.execute("select distinct r.analysis_method_id from %s order by analysis_method_id"\
								%(where_condition))
		analysis_method_id_ls = []
		for row in rows:
			analysis_method_id_ls.append(row.analysis_method_id)
		sys.stderr.write("Done.\n")
		return analysis_method_id_ls
	
	def getPhenotypeInfo(self, db,  where_condition):
		"""
		2008-08-29
			add -1 as a separator into phenotype_method_id_ls and others
		"""
		sys.stderr.write("Getting phenotype method info ...")
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id, p.biology_category_id from %s p, %s and p.id=r.phenotype_method_id order by p.biology_category_id, r.phenotype_method_id"\
								%(PhenotypeMethod.table.name, where_condition))
		phenotype_method_id_ls = []
		phenotype_method_id2index = {}
		phenotype_method_label_ls = []
		prev_biology_category_id = None
		no_of_separators = 0
		for row in rows:
			if prev_biology_category_id == None:
				prev_biology_category_id = row.biology_category_id
			elif row.biology_category_id!=prev_biology_category_id:
				prev_biology_category_id = row.biology_category_id
				#add a blank phenotype id as separator
				no_of_separators += 1
				phenotype_method_id2index[-no_of_separators] = len(phenotype_method_id_ls)
				phenotype_method_id_ls.append(-no_of_separators)
				phenotype_method_label_ls.append('')
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
		"""
		2008-08-29
			deal with separator (list_type_id=-1) in list_type_id_ls
		"""
		sys.stderr.write("Orderinig list type id and analysis_method id ... ")
		list_type_id_analysis_method_id_ls = []
		list_type_id_analysis_method_id2index = {}
		list_type_analysis_method_label_ls = []
		no_of_separators = 0
		for list_type_id in list_type_id_ls:
			if list_type_id==-1:	#separator
				no_of_separators += 1
				tup = (-no_of_separators,-1)
				list_type_id_analysis_method_id2index[tup] = len(list_type_id_analysis_method_id_ls)
				list_type_id_analysis_method_id_ls.append(tup)
				list_type_analysis_method_label_ls.append('')
				continue
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
	
	def get_data_matrix(self, db, phenotype_info, list_type_analysis_method_info, where_condition):
		sys.stderr.write("Getting data matrix ...")
		data_matrix = num.zeros([len(list_type_analysis_method_info.list_type_id_analysis_method_id2index), len(phenotype_info.phenotype_method_id2index)], num.float)
		data_matrix[:] = -1
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
				data_value = -2	#0 pvalue
			data_matrix[row_index, col_index] = data_value
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		return return_data
	
	def markDataMatrixBoundary(self, data_matrix, phenotype_info, list_type_analysis_method_info):
		"""
		2008-09-01
			mark boundary between different groups of rows and columns
			all those separators have id <0 and mark them with a value different from NA_value, -3.
		"""
		sys.stderr.write("Marking data matrix boundaries ...")
		for row_id in list_type_analysis_method_info.list_type_id_analysis_method_id2index:
			if row_id[0]<0:	#row_id is a tuple
				row_index = list_type_analysis_method_info.list_type_id_analysis_method_id2index[row_id]
				data_matrix[row_index,:] = -3
		for col_id in phenotype_info.phenotype_method_id2index:
			if col_id<0:
				col_index = phenotype_info.phenotype_method_id2index[col_id]
				data_matrix[:, col_index] = -3
		sys.stderr.write("Done.\n")
		return data_matrix
	
	def run(self):	
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		if self.test_result_type==1:
			test_result_class_table = CandidateGeneRankSumTestResult.table.name
			test_result_class_table = 'candidate_gene_rank_sum_test_result_2008_09_15'
		elif self.test_result_type==2:
			test_result_class_table = CandidateGeneTopSNPTest.table.name
		elif self.test_result_type==3:
			test_result_class_table = Stock_250kDB.CandidateGeneRankSumTestResultMethod.table.name
		else:
			sys.stderr.write(" test_result_type %s not supported.\n"%(self.test_result_type))
			sys.exit(2)

		#the condition for min_MAF is tricky because of the floating precision.
		if self.test_result_type==1:
			where_condition = "%s r, %s c, %s g where g.id=c.list_type_id and r.analysis_method_id is not null \
				and c.results_id=r.id and c.get_closest=%s and c.min_distance=%s and abs(c.min_MAF-%s)<0.00001"\
				%(ResultsMethod.table.name, test_result_class_table, GeneListType.table.name, self.get_closest, self.min_distance, self.min_MAF)
		elif self.test_result_type==2:
			where_condition = "%s r, %s rg, %s c, %s g where g.id=c.list_type_id and r.analysis_method_id is not null and r.id=rg.results_method_id \
				and c.results_id=rg.id and c.get_closest=%s and c.min_distance=%s and abs(c.min_MAF-%s)<0.00001"\
				%(ResultsMethod.table.name, ResultsByGene.table.name, test_result_class_table, GeneListType.table.name, self.get_closest, self.min_distance, self.min_MAF)
		elif self.test_result_type==3:
			where_condition = "%s r, %s c, %s g where g.id=c.list_type_id and r.analysis_method_id is not null \
				and c.results_id=r.id and c.get_closest=%s and c.min_distance=%s and abs(c.min_MAF-%s)<0.00001"\
				%(ResultsMethod.table.name, test_result_class_table, GeneListType.table.name, self.get_closest, self.min_distance, self.min_MAF)
		if self.call_method_id_ls:
			where_condition += " and r.call_method_id in (%s)"%self.call_method_id_ls
		
		if self.analysis_method_id_ls:
			where_condition += " and r.analysis_method_id in (%s)"%self.analysis_method_id_ls
		if self.super_type_id:
			where_condition += " and g.super_type_id=%s"%self.super_type_id
		
		if self.test_type:
			where_condition += " and c.test_type=%s"%self.test_type
		
		if self.test_result_type==1:
			pass
			where_condition += " and c.max_pvalue_per_gene=%s"%(self.max_pvalue_per_gene)
		elif self.test_result_type==2:
			where_condition += " and c.no_of_top_snps=%s"%(self.no_of_top_snps)		
		
		list_type_id_ls = self.getListTypeInfo(db, where_condition)
		analysis_method_id_ls = self.getAnalysisMethodInfo(db, where_condition)
		list_type_analysis_method_info = self.orderListTypeAnalysisMethodID(list_type_id_ls, analysis_method_id_ls)
		phenotype_info = self.getPhenotypeInfo(db, where_condition)
		rdata = self.get_data_matrix(db, phenotype_info, list_type_analysis_method_info, where_condition)
		
		rdata.data_matrix = self.markDataMatrixBoundary(rdata.data_matrix, phenotype_info, list_type_analysis_method_info)
		
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
			value2color_func = lambda x: Value2Color.value2HSLcolor(x, rdata.min_value, rdata.max_value)
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