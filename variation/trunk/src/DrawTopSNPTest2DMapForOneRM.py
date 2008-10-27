#!/usr/bin/env python
"""

Examples:
	#draw the 2D significance map of a result (results_id=2316 emma on LD) on candidate_gene_list 28(FT_short).
	DrawTopSNPTest2DMapForOneRM.py -i 2316 -l 28 -x /tmp/rm_2316_l1_C2.png -g  -C 2

Description:
	2008-10-23
		draw 2D map based on top snp test results (-log pvalue). row is no_of_top_snps(score cutoff). col is minimum distance to associate SNP to genes.
	
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
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color

num = importNumericArray()


class DrawTopSNPTest2DMapForOneRM(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene. Deprecated.'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0, 'n', 1, 'minimum Minor Allele Frequency.'],\
							('max_pvalue_per_gene', 0, int): [0, 'a', 0, 'take the most significant among all SNPs associated with one gene'],\
							('no_of_top_snps', 1, int): [100, 'f', 1, 'For test_result_type=2. how many number of top snps based on score or -log10(pvalue). Deprecated.'],\
							('results_id', 1, ):[None, 'i', 1, 'data from which ResultsMethod/ResultsByGene to draw'],\
							('list_type_id', 1, ):[None, 'l', 1, 'data from which list_type to look at'],\
							("test_type_id", 1, int): [15, 'y', 1, 'which type of tests. check db table analysis_method. likely be 14,15 etc.'],\
							("results_type", 1, int): [1, 'w', 1, 'data from Which result to output. 1: ResultsMethod, 2: ResultsByGene'],\
							('null_distribution_type_id', 0, int):[1, 'C', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. in db table null_distribution_type'],\
							("allow_two_sample_overlapping", 1, int): [0, '', 0, 'whether to allow one SNP to be assigned to both candidate and non-candidate gene group'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [15, 't', 1, 'Number of ticks on the legend'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-23
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getTopSNPTestType_id_ls(self, get_closest, min_MAF, allow_two_sample_overlapping, results_type, \
							test_type_id, null_distribution_type_id):
		"""
		"""
		sys.stderr.write("Getting TopSNPTestType_id_ls ...")
		TopSNPTestType_id_ls = []
		rows = Stock_250kDB.CandidateGeneTopSNPTestRMType.query.filter_by(get_closest=get_closest).\
			filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF>=min_MAF-0.0001).\
			filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF<=min_MAF+0.0001).\
			filter_by(allow_two_sample_overlapping = allow_two_sample_overlapping).\
			filter_by(results_type=results_type).\
			filter_by(test_type_id=test_type_id).\
			filter_by(null_distribution_type_id=null_distribution_type_id)
		for row in rows:
			TopSNPTestType_id_ls.append(row.id)
		sys.stderr.write("Done.\n")
		return TopSNPTestType_id_ls
	
	def get_no_of_top_snps_info(cls, db, from_where_clause):
		"""
		2008-10-23
		"""
		sys.stderr.write("Getting no_of_top_snps_info ...")
		rows = db.metadata.bind.execute("select distinct t.no_of_top_snps, t.min_score %s order by no_of_top_snps"%from_where_clause)
		id_ls = []
		id2index = {}
		label_ls = []
		no_of_separators = 0
		for row in rows:
			id2index[row.no_of_top_snps] = len(id_ls)
			id_ls.append(row.no_of_top_snps)
			label_ls.append('%s %s'%(row.no_of_top_snps, row.min_score))
		list_info = PassingData()
		list_info.id2index = id2index
		list_info.id_ls = id_ls
		list_info.label_ls = label_ls
		sys.stderr.write("Done.\n")
		return list_info
	get_no_of_top_snps_info = classmethod(get_no_of_top_snps_info)
	
	def get_min_distance_info(cls, db, from_where_clause):
		"""
		2008-10-23
		"""
		sys.stderr.write("Getting min_distance_info ...")
		rows = db.metadata.bind.execute("select distinct t.min_distance %s order by min_distance"%from_where_clause)
		id_ls = []
		id2index = {}
		label_ls = []
		no_of_separators = 0
		for row in rows:
			id2index[row.min_distance] = len(id_ls)
			id_ls.append(row.min_distance)
			label_ls.append('%s'%(row.min_distance))
		list_info = PassingData()
		list_info.id2index = id2index
		list_info.id_ls = id_ls
		list_info.label_ls = label_ls
		sys.stderr.write("Done.\n")
		return list_info
	
	get_min_distance_info = classmethod(get_min_distance_info)
	
	def get_data_matrix(cls, db, row_info, col_info, from_where_clause, need_other_values=False):
		sys.stderr.write("Getting data matrix ...")
		data_matrix = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
		data_matrix[:] = -1
		if need_other_values:
			data_matrix_candidate_sample_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_candidate_sample_size[:] = -1
			data_matrix_non_candidate_sample_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_non_candidate_sample_size[:] = -1
			data_matrix_candidate_gw_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_candidate_gw_size[:] = -1
			data_matrix_non_candidate_gw_size = num.zeros([len(row_info.id2index), len(col_info.id2index)], num.float)
			data_matrix_non_candidate_gw_size[:] = -1
		else:
			data_matrix_candidate_sample_size = None
			data_matrix_non_candidate_sample_size = None
			data_matrix_candidate_gw_size = None
			data_matrix_non_candidate_gw_size = None
		rows = db.metadata.bind.execute("select t.no_of_top_snps, t.min_distance, t.pvalue, t.candidate_sample_size, \
			t.non_candidate_sample_size, t.candidate_gw_size, t.non_candidate_gw_size %s"%from_where_clause)
		min_value = None
		max_value = None
		for row in rows:
			row_index = row_info.id2index[row.no_of_top_snps]
			col_index = col_info.id2index[row.min_distance]
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
			if need_other_values:
				data_matrix_candidate_sample_size[row_index, col_index] = row.candidate_sample_size
				data_matrix_non_candidate_sample_size[row_index, col_index] = row.non_candidate_sample_size
				data_matrix_candidate_gw_size[row_index, col_index] = row.candidate_gw_size
				data_matrix_non_candidate_gw_size[row_index, col_index] = row.non_candidate_gw_size
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.data_matrix_candidate_sample_size = data_matrix_candidate_sample_size
		return_data.data_matrix_non_candidate_sample_size = data_matrix_non_candidate_sample_size
		return_data.data_matrix_candidate_gw_size = data_matrix_candidate_gw_size
		return_data.data_matrix_non_candidate_gw_size = data_matrix_non_candidate_gw_size
		return_data.min_value = min_value
		return_data.max_value = max_value
		return return_data
	get_data_matrix = classmethod(get_data_matrix)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		TopSNPTestType_id_ls = self.getTopSNPTestType_id_ls(self.get_closest, self.min_MAF, self.allow_two_sample_overlapping, self.results_type, \
							self.test_type_id, self.null_distribution_type_id)
		if not TopSNPTestType_id_ls:
			sys.stderr.write("No TopSNPTestType matches the input requirements. Exit.\n")
			sys.exit(3)
		TopSNPTestType_id_ls_str = map(str, TopSNPTestType_id_ls)
		from_where_clause = "from %s t, %s y where t.type_id=y.id and t.results_id=%s and t.list_type_id=%s and y.id in (%s)"%\
			(Stock_250kDB.CandidateGeneTopSNPTestRM.table.name, Stock_250kDB.CandidateGeneTopSNPTestRMType.table.name,\
			self.results_id, self.list_type_id, ','.join(TopSNPTestType_id_ls_str))
		
		no_of_top_snps_info = self.get_no_of_top_snps_info(db, from_where_clause)
		min_distance_info = self.get_min_distance_info(db, from_where_clause)
		rdata = self.get_data_matrix(db, no_of_top_snps_info, min_distance_info, from_where_clause)
		
		header = ['no_of_top_snps', ''] + min_distance_info.label_ls
		strain_acc_list = no_of_top_snps_info.label_ls
		category_list = no_of_top_snps_info.label_ls
		
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
			im = drawMatrix(rdata.data_matrix, value2color_func, no_of_top_snps_info.label_ls,\
						min_distance_info.label_ls, with_grid=1, font=font)
			im = combineTwoImages(im, im_legend, font=font)
			im.save(self.fig_fname)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawTopSNPTest2DMapForOneRM
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()