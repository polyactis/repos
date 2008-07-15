#!/usr/bin/env python
"""
Examples:
	GeneListRankTest.py -e 389 -l 1 -u yh
	
Description:
2008-07-14 program to do pvalue rank test based on a given candidate gene list.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList
from Results2DB_250k import Results2DB_250k
from GenomeBrowser import GenomeBrowser

class GeneListRankTest(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("results_method_id", 1, ): [None, '', 1, ''],\
							("min_distance", 1, int): [50000, '', 1, ''],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getChrPos2Pvalue(self, results_method_id):
		"""
		"""
		sys.stderr.write("Getting chrpos2pvalue ...")
		rm = ResultsMethod.get(results_method_id)
		
		genome_wide_result = GenomeBrowser.getGenomeWideResultFromFile(rm.filename, do_log10_transformation=True)
		chrpos2pvalue = {}
		for data_obj in genome_wide_result.data_obj_ls:
			chrpos2pvalue[(data_obj.chromosome, data_obj.position)] = data_obj.value
		sys.stderr.write("Done.\n")
		return chrpos2pvalue
	
	def getGeneID2hit(self, chrpos2pvalue, min_distance=50000):
		sys.stderr.write("Getting gene_id2hit ... \n")
		offset_index = 0
		block_size = 5000
		rows = SnpsContext.query.offset(offset_index).limit(block_size)
		gene_id2hit = {}
		counter = 0
		while rows.count()!=0:
			for row in rows:
				if row.disp_pos>=-min_distance and row.snp.end_position==None:	#it's a true SNP, not segment
					chrpos_key = (row.snp.chromosome, row.snp.position)
					pvalue = chrpos2pvalue.get(chrpos_key)
					if pvalue is None:	#no pvalue
						continue
					passingdata = PassingData(chr=row.snp.chromosome, pos=row.snp.position, pvalue=pvalue, snps_id=row.snp.id)
					if row.gene_id not in gene_id2hit:
						gene_id2hit[row.gene_id] = passingdata
					elif pvalue>gene_id2hit[row.gene_id].pvalue:	#pvalue is -log()
						gene_id2hit[row.gene_id] = passingdata
					counter += 1
				offset_index += 1
			sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))
			rows = SnpsContext.query.offset(offset_index).limit(block_size)
		
		sys.stderr.write("Done.\n")
		return gene_id2hit
	
	def getGeneList(self, list_type_id):
		sys.stderr.write("Getting gene_list ... ")
		rows = GeneList.query.filter_by(list_type_id=list_type_id)
		candidate_gene_list = []
		for row in rows:
			candidate_gene_list.append(row.gene_id)
		sys.stderr.write("Done.\n")
		return candidate_gene_list
	
	def prepareDataForRankTest(self, candidate_gene_list, gene_id2hit):
		sys.stderr.write("Preparing data for rank test ... ")
		candidate_gene_pvalue_list = []
		from sets import Set
		candidate_gene_set = Set(candidate_gene_list)
		candidate_gene_ls = []
		
		non_candidate_gene_ls = []
		non_candidate_gene_pvalue_list = []
		
		gene_id_ls = gene_id2hit.keys()
		gene_id_ls.sort()
		for gene_id in gene_id_ls:
			hit = gene_id2hit[gene_id]
			if gene_id in candidate_gene_set:
				candidate_gene_ls.append(gene_id)
				candidate_gene_pvalue_list.append(hit.pvalue)
			else:
				non_candidate_gene_ls.append(gene_id)
				non_candidate_gene_pvalue_list.append(hit.pvalue)
		sys.stderr.write("Done.\n")
		passingdata = PassingData(candidate_gene_ls=candidate_gene_ls, candidate_gene_pvalue_list=candidate_gene_pvalue_list,\
								non_candidate_gene_ls=non_candidate_gene_ls, non_candidate_gene_pvalue_list=non_candidate_gene_pvalue_list)
		return passingdata
	
	def run(self):
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		chrpos2pvalue = self.getChrPos2Pvalue(self.results_method_id)
		gene_id2hit = self.getGeneID2hit(chrpos2pvalue, self.min_distance)
		candidate_gene_list = self.getGeneList(self.list_type_id)
		passingdata = self.prepareDataForRankTest(candidate_gene_list, gene_id2hit)
		import rpy
		w_result = rpy.r.wilcox_test(passingdata.candidate_gene_pvalue_list, passingdata.non_candidate_gene_pvalue_list, conf_int=rpy.r.TRUE)
		print w_result['p.value']
		print w_result['statistic']
		print passingdata.candidate_gene_pvalue_list
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GeneListRankTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()