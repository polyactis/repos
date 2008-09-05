#!/usr/bin/env python
"""

Examples:
	OutputOverlappingResultsByGene.py -l 17 -o overlapping_results.tsv -x overlapping_results.png -r  -s 10 -e ./FreeSerif.ttf  -t 40
	
Description:
	Output a matrix with each cell showing how many genes are overlapping between two results.
	It has the option to draw the matrix or defer it to pymodule/DrawMatrix.py.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, write_data_matrix, SNPData
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, \
	CandidateGeneTopSNPTest, CandidateGeneRankSumTestResult, AnalysisMethod, PhenotypeMethod, ResultsByGene
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color

num = importNumericArray()

class OutputOverlappingResultsByGene(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('call_method_id_ls', 1, ):[None, 'l', 1, 'Restrict results based on these call_methods, coma-separated list of ids.'],\
							('analysis_method_id_ls', 0, ):[None, 'j', 1, 'Restrict results based on these analysis, coma-separated. Default is no such restriction.'],\
							("result_type", 1, int): [1, 'y', 1, 'Which test result to output. 1: CandidateGeneRankSumTestResult, 2: CandidateGeneTopSNPTest'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 1, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [5, 't', 1, 'Number of ticks on the legend'],\
							("max_rank", 1, int): [1000, 'm', 1, 'only take genes above this rank'],\
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
	
	def getResultsMethodIDInfo(self, db):
		sys.stderr.write("Gettiing ResultsMethodIDInfo ...")
		results_method_id_info = PassingData()
		results_method_id_ls = []
		results_method_id2index = {}
		results_method_id_label_ls = []
		rows = db.metadata.bind.execute("select distinct rg.results_method_id, r.analysis_method_id, r.phenotype_method_id, p.biology_category_id from %s rg, %s r, %s p \
			where rg.results_method_id=r.id and p.id=r.phenotype_method_id \
			order by p.biology_category_id, r.phenotype_method_id, r.analysis_method_id"%(ResultsByGene.table.name, ResultsMethod.table.name, PhenotypeMethod.table.name))
		prev_phenotype_method_id = None
		prev_biology_category_id = None
		no_of_separators = 0
		for row in rows:
			if prev_biology_category_id==None:
				prev_biology_category_id = row.biology_category_id
			elif row.biology_category_id!=prev_biology_category_id:
				prev_biology_category_id = row.biology_category_id
				no_of_separators += 1
				results_method_id2index[-no_of_separators] = len(results_method_id_ls)
				results_method_id_ls.append(-no_of_separators)
				results_method_id_label_ls.append('')
			
			if prev_phenotype_method_id == None:
				prev_phenotype_method_id = row.phenotype_method_id
			elif row.phenotype_method_id!=prev_phenotype_method_id:
				prev_phenotype_method_id = row.phenotype_method_id
				#add a blank phenotype id as separator
				no_of_separators += 1
				results_method_id2index[-no_of_separators] = len(results_method_id_ls)
				results_method_id_ls.append(-no_of_separators)
				results_method_id_label_ls.append('')
			results_method_id2index[row.results_method_id] = len(results_method_id_ls)
			results_method_id_ls.append(row.results_method_id)
			am = AnalysisMethod.get(row.analysis_method_id)
			pm = PhenotypeMethod.get(row.phenotype_method_id)
			results_method_id_label_ls.append('%s_%s_%s'%(am.short_name, pm.short_name, pm.id))
		results_method_id_info.results_method_id_ls = results_method_id_ls
		results_method_id_info.results_method_id2index = results_method_id2index
		results_method_id_info.results_method_id_label_ls = results_method_id_label_ls
		sys.stderr.write("Done.\n")
		return results_method_id_info
	
	def getResultsMethodID2GeneSet(self, db, max_rank=1000):
		sys.stderr.write("Gettiing results_method_id2gene_set ... \n")
		results_method_id2gene_set = {}
		i = 0
		block_size = 10000
		rows = ResultsByGene.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				if row.rank<=max_rank:
					if row.results_method_id not in results_method_id2gene_set:
						results_method_id2gene_set[row.results_method_id] = Set()
					results_method_id2gene_set[row.results_method_id].add(row.gene_id)
				i += 1
			
			if self.debug and i>50000:
				break
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*40, i))
			rows = ResultsByGene.query.offset(i).limit(block_size)
		sys.stderr.write("Done.\n")
		return results_method_id2gene_set
	
	def getDataMatrix(self, results_method_id2gene_set, results_method_id_info):
		sys.stderr.write("Gettiing data matrix ...")
		data_matrix = num.zeros([len(results_method_id_info.results_method_id_ls), len(results_method_id_info.results_method_id_ls)], num.float)
		data_matrix[:] = -1
		min_value = None
		max_value = None
		no_of_results = len(results_method_id_info.results_method_id_ls) 
		for i in range(no_of_results):
			results_method_id1 = results_method_id_info.results_method_id_ls[i]
			for j in range(i, no_of_results):
				results_method_id2 = results_method_id_info.results_method_id_ls[j]
				row_index = results_method_id_info.results_method_id2index.get(results_method_id1)
				col_index = results_method_id_info.results_method_id2index.get(results_method_id2)
				if results_method_id1<0 or results_method_id2<0:
					data_value = -3
				elif results_method_id1 in results_method_id2gene_set and results_method_id2 in results_method_id2gene_set:
					data_value = len(results_method_id2gene_set[results_method_id1]&results_method_id2gene_set[results_method_id2])
					if min_value==None:
						min_value = data_value
					elif data_value<min_value:
						min_value = data_value
					
					if max_value==None:
						max_value=data_value
					elif data_value>max_value:
						max_value =data_value
				else:
					continue
					#data_value = -1
				data_matrix[row_index, col_index] = data_value
				data_matrix[col_index, row_index] = data_value
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		sys.stderr.write("Done.\n")
		return return_data
	
	def run(self):	
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		results_method_id_info = self.getResultsMethodIDInfo(db)
		results_method_id2gene_set = self.getResultsMethodID2GeneSet(db, self.max_rank)
		rdata = self.getDataMatrix(results_method_id2gene_set, results_method_id_info)
		
		header = ['', ''] + results_method_id_info.results_method_id_label_ls
		strain_acc_list = results_method_id_info.results_method_id_label_ls
		category_list = results_method_id_info.results_method_id_ls
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
			im = drawMatrix(rdata.data_matrix, value2color_func, strain_acc_list,\
						strain_acc_list, with_grid=1, font=font)
			im = combineTwoImages(im, im_legend, font=font)
			im.save(self.fig_fname)
			
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = OutputOverlappingResultsByGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
