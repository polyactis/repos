#!/usr/bin/env python
"""
Examples:
	GeneListRankTest.py -e 389 -l 1 -u yh -c
	
	#debug, quick testing
	GeneListRankTest.py -e 389,190 -l 1 -u yh -b
	
Description:
	2008-07-14 program to do pvalue rank test based on a given candidate gene list.
	
	It verifies against the db several things:
	1. whether the results_method_id is available
	2. whether results_method_type_id is 1 (association)
	3. whether the same (results_method_id, list_type_id) pair has been in the candidate_gene_rank_sum_test_result table.
	
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
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, CandidateGeneRankSumTestResult
from Results2DB_250k import Results2DB_250k
from pymodule import getGenomeWideResultFromFile

from pymodule.db import TableClass
class SnpsContextWrapper(TableClass):
	chrpos2row_index_ls = None
	data_matrix = None
	no_of_rows = 0
	
	def addOneSnpsContext(self, snps_context):
		if self.data_matrix is None:
			self.data_matrix = []
		self.data_matrix.append([snps_context.snps_id, snps_context.disp_pos, snps_context.gene_id])
		self.no_of_rows += 1
		if self.chrpos2row_index_ls is None:
			self.chrpos2row_index_ls = {}
		chrpos_key = (snps_context.snp.chromosome, snps_context.snp.position)
		if chrpos_key not in self.chrpos2row_index_ls:
			self.chrpos2row_index_ls[chrpos_key] = []
		self.chrpos2row_index_ls[chrpos_key].append(self.no_of_rows-1)
	
	def returnGeneLs(self, chromosome, position):
		return_matrix = []
		chrpos_key = (chromosome, position)
		row_index_ls = self.chrpos2row_index_ls.get(chrpos_key)
		if row_index_ls is not None:
			for row_index in row_index_ls:
				return_matrix.append(self.data_matrix[row_index])
		return return_matrix
		
	

class GeneListRankTest(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("results_method_id_ls", 1, ): [None, 'e', 1, 'comma-separated results_method_id list'],\
							("min_distance", 1, int): [50000, 'm', 1, ''],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							("output_fname", 0, ): [None, 'o', 1, ''],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		results_method_id_ls = self.results_method_id_ls.split(',')
		self.results_method_id_ls = map(int, results_method_id_ls)
		self.results_method_id_ls.sort()
			
	def constructDataStruc(self, min_distance=50000):
		"""
		2008-07-16
		"""
		sys.stderr.write("Loading data structure from db ... \n")
		offset_index = 0
		block_size = 5000
		rows = SnpsContext.query.offset(offset_index).limit(block_size)
		data_struc = SnpsContextWrapper()
		counter = 0
		while rows.count()!=0:
			for row in rows:
				if row.disp_pos>=-min_distance and row.snp.end_position==None:	#it's a true SNP, not segment					
					data_struc.addOneSnpsContext(row)
					counter += 1
				offset_index += 1
			sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))
			if self.debug and offset_index > 1000:
				break
			rows = SnpsContext.query.offset(offset_index).limit(block_size)
		sys.stderr.write("Done.\n")
		return data_struc
	
	def getChrPos2Pvalue(self, results_method_id):
		"""
		2008-07-16
			deprecated
		"""
		sys.stderr.write("Getting chrpos2pvalue ...")
		rm = ResultsMethod.get(results_method_id)
		
		genome_wide_result = getGenomeWideResultFromFile(rm.filename, do_log10_transformation=True)
		chrpos2pvalue = {}
		for data_obj in genome_wide_result.data_obj_ls:
			chrpos2pvalue[(data_obj.chromosome, data_obj.position)] = data_obj.value
		sys.stderr.write("Done.\n")
		return chrpos2pvalue
	
	def getGeneID2hit(self, rm, snps_context_wrapper):
		"""
		2008-07-17
			no do_log10_transformation
		2008-07-16
			reverse the order of 1st-read results file, 2nd-read db.snps_context
		"""
		sys.stderr.write("Getting gene_id2hit ... \n")
		
		genome_wide_result = getGenomeWideResultFromFile(rm.filename)
		
		gene_id2hit = {}
		counter = 0
		for data_obj in genome_wide_result.data_obj_ls:
			pvalue = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				passingdata = PassingData(chr=data_obj.chromosome, pos=data_obj.position, pvalue=pvalue, disp_pos=disp_pos, snps_id=snps_id)
				if gene_id not in gene_id2hit:
					gene_id2hit[gene_id] = passingdata
				elif pvalue>gene_id2hit[gene_id].pvalue:	#pvalue is -log()
					gene_id2hit[gene_id] = passingdata
			counter += 1
			#sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))		
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
	
	def output_gene_id2hit(self, gene_id2hit, output_fname):
		sys.stderr.write("Outputting gene_id2hit ... ")
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['gene-id', 'chromosome', 'position', 'snps-id', 'disp_pos', 'pvalue', 'rank'])
		gene_id_ls = gene_id2hit.keys()
		gene_id_ls.sort()
		pvalue_ls = [gene_id2hit[gene_id].pvalue for gene_id in gene_id_ls]
		import rpy
		rank_ls = rpy.r.rank(pvalue_ls)
		for i in range(len(gene_id_ls)):
			gene_id = gene_id_ls[i]
			rank = rank_ls[i]
			hit = gene_id2hit[gene_id]
			one_row = [gene_id, hit.chr, hit.pos, hit.snps_id, hit.disp_pos, hit.pvalue, rank]
			writer.writerow(one_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def run_wilcox_test(self, results_method_id, snps_context_wrapper, list_type_id):
		"""
		2008-07-17
			split out as a standalone function so that MpiGeneListRankTest.py could call it more easily.
		"""
		if self.debug:
			sys.stderr.write("Running wilcox test ... ")
		rm = ResultsMethod.get(results_method_id)
		if not rm:
			sys.stderr.write("No results method available for results_method_id=%s.\n"%results_method_id)
			return None
		if rm.results_method_type_id!=1:
			sys.stderr.write("skip non-association results. results_method_type_id=%s, results_method_id=%s.\n"%(rm.results_method_type_id, results_method_id))
			return None
		db_results = CandidateGeneRankSumTestResult.query.filter_by(results_method_id=results_method_id).filter_by(list_type_id=self.list_type_id)
		if db_results.count()>0:	#done before
			db_result = db_results.first()
			sys.stderr.write("It's done already. id=%s, results_method_id=%s, list_type_id=%s, pvalue=%s, statistic=%s.\n"%\
							(db_result.id, db_result.results_method_id, db_result.list_type_id, db_result.pvalue, db_result.statistic))
			return None
		
		try:
			gene_id2hit = self.getGeneID2hit(rm, snps_context_wrapper)
			#if getattr(self, 'output_fname', None):
			#	self.output_gene_id2hit(gene_id2hit, self.output_fname)
			candidate_gene_list = self.getGeneList(list_type_id)
			passingdata = self.prepareDataForRankTest(candidate_gene_list, gene_id2hit)
			import rpy
			w_result = rpy.r.wilcox_test(passingdata.candidate_gene_pvalue_list, passingdata.non_candidate_gene_pvalue_list, conf_int=rpy.r.TRUE)
		except:
			traceback.print_exc()
			print sys.exc_info()
			return None
		candidate_gene_rank_sum_test_result = CandidateGeneRankSumTestResult(list_type_id=list_type_id, statistic=w_result['statistic']['W'],\
																			pvalue=w_result['p.value'])
		candidate_gene_rank_sum_test_result.results_method = rm
		if self.debug:
			sys.stderr.write("Done.\n")
		return candidate_gene_rank_sum_test_result
	
	def run(self):
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		session = db.session
		if self.commit:
			session.begin()
		#chrpos2pvalue = self.getChrPos2Pvalue(self.results_method_id)
		#gene_id2hit = self.getGeneID2hit(chrpos2pvalue, self.min_distance)
		snps_context_wrapper = self.constructDataStruc(self.min_distance)
		
		if getattr(self, 'output_fname', None):
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			writer.writerow(['results_method_id', 'list_type_id', 'wilcox.test.pvalue', 'statistic'])
		else:
			writer = None
		for results_method_id in self.results_method_id_ls:
			candidate_gene_rank_sum_test_result = self.run_wilcox_test(results_method_id, snps_context_wrapper, self.list_type_id)
			if candidate_gene_rank_sum_test_result is not None:
				row = [results_method_id, self.list_type_id, candidate_gene_rank_sum_test_result.pvalue, candidate_gene_rank_sum_test_result.statistic]
				print row
				if writer:
					writer.writerow(row)
				session.save(candidate_gene_rank_sum_test_result)
				session.flush()
		if self.commit:
			session.flush()
			session.commit()
		#print passingdata.candidate_gene_pvalue_list
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GeneListRankTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()