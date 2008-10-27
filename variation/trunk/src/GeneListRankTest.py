#!/usr/bin/env python
"""
Examples:
	GeneListRankTest.py -e 389 -l 1 -u yh -c
	
	#debug, quick testing
	GeneListRankTest.py -e 389,190 -l 1 -u yh -b
	
	#apart from doing rank test, pickle the snps_context_wrapper into a file.
	GeneListRankTest.py -e 389 -l 1 -u yh -c -s /tmp/snps_context
	
	#
	GeneListRankTest.py -e 389 -l1 -u yh -j1 -y2 -s ./mnt2/panfs/250k/snps_context_g0_m20000  -p yh324 -b
Description:
	2008-10-11
		test_type=4,5,6 is same as test_type=1,2,3 but allow one SNP to be assigned to both candidate and non-candidate gene group
	
	2008-10-09 the permutation rank sum test (test_type=2,3) preserves the order of SNPs within each chromosome
		1. shuffle chromosome order(not for test_type=3) and put all SNPs into a list in order
		2. shift all SNPs in that list by a random number
		3. get rank sums of scores from fixed candidate gene SNP indices that
	
	2008-07-14 program to do pvalue rank test based on a given candidate gene list.
	
	It verifies against the db several things:
	1. whether the results_id is available
	2. whether results_method_type_id is 1 (association)
	3. whether the same (results_id, list_type_id) pair has been in the candidate_gene_rank_sum_test_result table.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, importNumericArray
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, CandidateGeneRankSumTestResult, \
	ResultsByGene, CandidateGeneRankSumTestResultMethod
from Results2DB_250k import Results2DB_250k
from pymodule import getGenomeWideResultFromFile
from sets import Set
from pymodule.db import TableClass
import random

num = importNumericArray()
numpy = num

class SnpsContextWrapper(object):
	def __init__(self, get_closest=False):
		self.get_closest = get_closest
		self.chrpos2snps_id = {}
		self.snps_id2left_or_right2gene_id_ls = {}
		self.data_matrix = None
		self.no_of_rows = 0
	
	def addOneSnpsContext(self, snps_context):
		"""
		if self.data_matrix is None:
			self.data_matrix = []
		self.data_matrix.append([snps_context.snps_id, snps_context.disp_pos, snps_context.gene_id])
		"""
		snps_id = snps_context.snps_id
		if snps_id not in self.snps_id2left_or_right2gene_id_ls:
			self.snps_id2left_or_right2gene_id_ls[snps_id] = {}
		left_or_right = snps_context.left_or_right
		if left_or_right not in self.snps_id2left_or_right2gene_id_ls[snps_id]:
			self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right] = [(snps_context.disp_pos, snps_context.gene_id)]
			self.no_of_rows += 1
		elif left_or_right=='touch':
			self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right].append((snps_context.disp_pos, snps_context.gene_id))
			self.no_of_rows += 1
		elif self.get_closest:	#it's left or right. and get closest
			if abs(snps_context.disp_pos)==abs(self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right][0][0]):
				self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right].append((snps_context.disp_pos, snps_context.gene_id))
				self.no_of_rows += 1
			elif abs(snps_context.disp_pos)<abs(self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right][0][0]):
				self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right] = [(snps_context.disp_pos, snps_context.gene_id)]
		else:
			self.snps_id2left_or_right2gene_id_ls[snps_id][left_or_right].append((snps_context.disp_pos, snps_context.gene_id))
			self.no_of_rows += 1
		chrpos_key = (snps_context.snp.chromosome, snps_context.snp.position)
		if chrpos_key not in self.chrpos2snps_id:
			self.chrpos2snps_id[chrpos_key] = snps_context.snps_id
		
		#self.chrpos2row_index_ls[chrpos_key].append(self.no_of_rows-1)
	
	def returnGeneLs(self, chromosome, position):
		return_matrix = []
		chrpos_key = (chromosome, position)
		snps_id = self.chrpos2snps_id.get(chrpos_key)
		if snps_id is not None:
			if self.get_closest and 'touch' in self.snps_id2left_or_right2gene_id_ls[snps_id]:	#'touch' is closest, forget about all others if there's 'touch'
				for disp_pos, gene_id in self.snps_id2left_or_right2gene_id_ls[snps_id]['touch']:
					return_matrix.append([snps_id, disp_pos, gene_id])
			else:
				for left_or_right, gene_id_ls in self.snps_id2left_or_right2gene_id_ls[snps_id].iteritems():
					for disp_pos, gene_id in gene_id_ls:
						return_matrix.append([snps_id, disp_pos, gene_id])
		return return_matrix


class GeneListRankTest(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("results_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated results_by_gene id list, like 1,3-7'],\
							("min_distance", 1, int): [50000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 0, float): [0, 'n', 1, 'minimum Minor Allele Frequency. deprecated.'],\
							('min_sample_size', 0, int): [5, 'i', 1, 'minimum size for both candidate and non-candidate sets to do wilcox.test'],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("output_fname", 0, ): [None, 'o', 1, 'To store rank test results into this file as a backup version of db'],\
							("snps_context_picklef", 0, ): [None, 's', 1, 'given the option, if the file does not exist yet, to store a pickled snps_context_wrapper into it, min_distance and flag get_closest will be attached to the filename. If the file exists, load snps_context_wrapper out of it.'],\
							("results_type", 1, int): [1, 'w', 1, 'which type of results. 1; ResultsMethod, 2: ResultsByGene'],\
							("test_type_id", 1, int): [1, 'y', 1, 'which type of rank sum test. 1: r.wilcox.test() 2: loop-permutation. 3: loop-permutation with chromosome order kept. 4,5,6 are their counterparts which allow_two_sample_overlapping.'],\
							("no_of_permutations", 1, int): [40000, '', 1, 'no of permutations to carry out'],\
							("no_of_min_breaks", 1, int): [30, '', 1, 'no of minimum times that rank_sum_stat_perm>=rank_sum_stat to break away. if 0, no breaking'],\
							('null_distribution_type_id', 0, int):[1, 'C', 1, 'Type of null distribution. 1=original, 2=permutation, 3=random gene list. in db table null_distribution_type'],\
							("allow_two_sample_overlapping", 1, int): [0, '', 0, 'whether to allow one SNP to be assigned to both candidate and non-candidate gene group'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	#("allow_two_sample_overlapping", 1, int): [0, 'x', 0, 'whether to allow one SNP to be assigned to both candidate and non-candidate gene group'],\
	def __init__(self,  **keywords):
		"""
		2008-10-09
			add option results_type, test_type, no_of_permutations, permutation_type
		2008-07-24
			split results_id_ls if it exists, to accomodate MpiGeneListRankTest which removed this option
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if getattr(self, 'results_id_ls', None):
			self.results_id_ls = getListOutOfStr(self.results_id_ls, data_type=int)
			self.results_id_ls.sort()
		else:
			self.results_id_ls = []
			
	def constructDataStruc(self, min_distance=50000, get_closest=0):
		"""
		2008-08-28
			in debug mode, increase the break-off offset_index from 10000 to 50000
		2008-08-14
			add get_closest
		2008-07-16
		"""
		sys.stderr.write("Loading data structure from db ... \n")
		offset_index = 0
		block_size = 10000
		rows = SnpsContext.query.offset(offset_index).limit(block_size)
		data_struc = SnpsContextWrapper(get_closest)
		counter = 0
		while rows.count()!=0:
			for row in rows:
				if row.left_or_right=='touch':
					data_struc.addOneSnpsContext(row)
					counter += 1
				elif abs(row.disp_pos)<=min_distance and row.snp.end_position==None:	#it's a true SNP, not segment					
					data_struc.addOneSnpsContext(row)
					counter += 1
				offset_index += 1
			sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))
			if self.debug and offset_index > 50000:
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
	
	def getResultMethodContent(cls, rm, results_directory=None, min_MAF=0.1, construct_chr_pos2index=False, pdata=None):
		"""
		2008-10-23
			if pdata doesn't have construct_chr_pos2index defined. otherwise, pdata overrides the option.
		2008-10-22
			before deciding do_log10_transformation based on analysis_method, try to get it from pdata
		2008-10-21
			add pdata to conceal the passing of chr_pos2index to getGenomeWideResultFromFile()
		2008-10-15
			cache the genome_wide_result under cls.genome_wide_result, if cls.genome_wide_result.results_id==rm.id, directly return that.
		2008-09-24
			add option construct_chr_pos2index
		2008-09-16
			if result_fname is not a file, return None
		2008-09-15
			use field smaller_score_more_significant from ResultsMethod to set do_log10_transformation
		2008-08-13
			split from getGeneID2MostSignificantHit()
		"""
		genome_wide_result = getattr(cls, 'genome_wide_result', None)
		do_log10_transformation = getattr(pdata, 'do_log10_transformation', None)
		if genome_wide_result is not None and genome_wide_result.results_id==rm.id:
			return genome_wide_result
		
		if rm.analysis_method_id==13:
			sys.stderr.write("Analysis method id=%s is not supported.\n"%rm.analysis_method_id)
			return None
		
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		
		if do_log10_transformation is None:
			#based on the analysis method id, whether do -log() or not. it'll affect the later step of taking maximum pvalue out of SNPs associated with one gene
			if hasattr(rm, 'analysis_method'):
					if rm.analysis_method.smaller_score_more_significant==1:
						do_log10_transformation = True
					else:
						do_log10_transformation = False
			else:
				return None
		if pdata is None:
			pdata = PassingData()
		pdata.min_MAF = min_MAF
		if not hasattr(pdata, 'construct_chr_pos2index'):	#if pdata doesn't have construct_chr_pos2index defined. otherwise, pdata overrides the option.
			pdata.construct_chr_pos2index = construct_chr_pos2index
		if os.path.isfile(result_fname):
			genome_wide_result = getGenomeWideResultFromFile(result_fname, do_log10_transformation=do_log10_transformation, pdata=pdata)
			genome_wide_result.results_id = rm.id
		else:
			sys.stderr.write("Skip. %s doesn't exist.\n"%result_fname)
			genome_wide_result = None
		cls.genome_wide_result = genome_wide_result
		return genome_wide_result
	
	getResultMethodContent = classmethod(getResultMethodContent)
	
	def getGeneID2MostSignificantHit(self, rm, snps_context_wrapper, results_directory=None, min_MAF=0.1):
		"""
		2008-09-16
			if genome_wide_result is None, return None
		2008-08-13
			add min_MAF
		2008-07-21
			based on the analysis method id, whether do -log() or not. it'll affect the later step of taking maximum pvalue out of SNPs associated with one gene
		2008-07-17
			add results_directory in case of results in a different directory
		2008-07-17
			no do_log10_transformation
		2008-07-16
			reverse the order of 1st-read results file, 2nd-read db.snps_context
		"""
		sys.stderr.write("Getting gene_id2hit ... \n")
		genome_wide_result = self.getResultMethodContent(rm, results_directory, min_MAF)
		if genome_wide_result is None:
			return None
		gene_id2hit = {}
		counter = 0
		for data_obj in genome_wide_result.data_obj_ls:
			score = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				passingdata = PassingData(chr=data_obj.chromosome, pos=data_obj.position, score=score, disp_pos=disp_pos, snps_id=snps_id, comment=None)
				if gene_id not in gene_id2hit:
					gene_id2hit[gene_id] = passingdata
				elif score>gene_id2hit[gene_id].score:	#score is -log()
					gene_id2hit[gene_id] = passingdata
			counter += 1
			#sys.stderr.write("%s%s\t%s"%('\x08'*40, offset_index, counter))		
		sys.stderr.write("Done.\n")
		return gene_id2hit
	
	def getGeneID2MostSignificantHitFromSNPPairFile(self, rm, snps_context_wrapper, results_directory=None, min_MAF=0.1):
		"""
		2008-09-16
			similar to getGeneID2MostSignificantHit()
			but the results is in a different format. output of MpiIntraGeneSNPPairAsso.py
		"""
		sys.stderr.write("Getting gene_id2hit ... \n")
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		if not os.path.isfile(result_fname):
			return None
		
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		gene_id2hit = {}
		counter = 0
		for row in reader:
			gene_id = int(row[col_name2index['gene1_id']])
			score = -math.log(float(row[col_name2index['pvalue']]))
			snp1_id = row[col_name2index['snp1_id']]
			snp2_id = row[col_name2index['snp2_id']]
			bool_type = row[col_name2index['bool_type']]
			count1 = int(row[col_name2index['count1']])
			count2 = int(row[col_name2index['count2']])
			maf = float(min(count1, count2))/(count1+count2)
			if maf>=min_MAF:
				snps_id = snp1_id
				if snp2_id:
					snps_id += '|%s'%snp2_id
				passingdata = PassingData(score=score, snps_id=snps_id, disp_pos=None, comment='bool_type:%s'%bool_type)
				if gene_id not in gene_id2hit:
					gene_id2hit[gene_id] = passingdata
				elif score>gene_id2hit[gene_id].score:	#score is -log()
					gene_id2hit[gene_id] = passingdata
			counter += 1
			if self.report and counter%10000==0:
				sys.stderr.write("%s%s"%('\x08'*100, counter))
		sys.stderr.write("Done.\n")
		return gene_id2hit
		
	def getGeneList(cls, list_type_id):
		sys.stderr.write("Getting gene_list ... ")
		rows = GeneList.query.filter_by(list_type_id=list_type_id)
		candidate_gene_list = []
		for row in rows:
			candidate_gene_list.append(row.gene_id)
		sys.stderr.write("%s genes. Done.\n"%(len(candidate_gene_list)))
		return candidate_gene_list
	getGeneList = classmethod(getGeneList)
	
	def prepareDataForRankTestGivenGeneID2Hit(self, candidate_gene_list, gene_id2hit):
		sys.stderr.write("Preparing data for rank test ... ")
		candidate_gene_score_list = []
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
				#candidate_gene_ls.append(gene_id)
				candidate_gene_pvalue_list.append(hit.score)
			else:
				#non_candidate_gene_ls.append(gene_id)
				non_candidate_gene_pvalue_list.append(hit.score)
		sys.stderr.write("Done.\n")
		passingdata = PassingData(candidate_gene_ls=candidate_gene_ls, candidate_gene_pvalue_list=candidate_gene_pvalue_list,\
								non_candidate_gene_ls=non_candidate_gene_ls, non_candidate_gene_pvalue_list=non_candidate_gene_pvalue_list)
		return passingdata
	
	def prepareDataForRankTest(self, rm, snps_context_wrapper, candidate_gene_list, results_directory=None, min_MAF=0.1):
		"""
		2008-08-13
			for each gene, don't take the most significant hit
		"""
		sys.stderr.write("Preparing data for rank test ... ")
		candidate_gene_set = Set(candidate_gene_list)
		candidate_gene_ls = []
		candidate_gene_pvalue_list = []
		non_candidate_gene_ls = []
		non_candidate_gene_pvalue_list = []
		
		genome_wide_result = self.getResultMethodContent(rm, results_directory, min_MAF)
		
		counter = 0
		for data_obj in genome_wide_result.data_obj_ls:
			score = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					candidate_gene_pvalue_list.append(score)
				else:
					non_candidate_gene_pvalue_list.append(score)
		passingdata = PassingData(candidate_gene_ls=candidate_gene_ls, candidate_gene_pvalue_list=candidate_gene_pvalue_list,\
								non_candidate_gene_ls=non_candidate_gene_ls, non_candidate_gene_pvalue_list=non_candidate_gene_pvalue_list)
		sys.stderr.write("Done.\n")
		return passingdata
	
	def getTopResultsByGene(self, rbg, param_data):
		"""
		2008-10-24
			handle param_data.candidate_gene_set. if it's present, no_of_top_lines is relative to that gene set
		2008-10-02
			add analysis_method id to each ResultsByGene entry
		2008-09-30
			function to read a certain number of genes out of a ResultsByGene file.
			param_data has property results_directory, no_of_top_lines
		"""
		if self.report:
			sys.stderr.write("Getting results_by_gene ... ")
		if param_data.results_directory:	#given a directory where all results are.
			result_fname = os.path.join(param_data.results_directory, os.path.basename(rbg.filename))
		else:
			result_fname = rbg.filename
		if not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		no_of_lines_to_read = getattr(param_data, 'no_of_top_lines', None)
		candidate_gene_set = getattr(param_data, 'candidate_gene_set', None)	#2008-10-24
		return_data_ls = []
		for row in reader:
			gene_id = int(row[col_name2index['gene_id']])
			if candidate_gene_set and gene_id not in candidate_gene_set:	#2008-10-24
				continue
			score = float(row[col_name2index['score']])
			snps_id = int(row[col_name2index['snps_id']])
			disp_pos = int(row[col_name2index['disp_pos']])
			comment = row[col_name2index['comment']]
			entry = PassingData(gene_id=gene_id, score=score, snps_id=snps_id, disp_pos=disp_pos, comment=comment, \
							analysis_method=rbg.results_method.analysis_method_id)
			return_data_ls.append(entry)
			counter += 1
			if no_of_lines_to_read is not None and counter>=no_of_lines_to_read:
				break
		del reader
		if self.report:
			sys.stderr.write("Done.\n")
		return return_data_ls
	
	def prepareDataForRankTestFromResultsByGene(self, rm, param_data):
		"""
		2008-09-16
			if result_fname doesn't exist, return None
		2008-09-16
			add function to only read a number of top lines
		2008-09-16
			file associated with a results_by_gene entry is gene-based. in it, one gene only appears once and is ordered by its score.
			just get the gene_id and test whether it's candidate or not.
		"""
		sys.stderr.write("Preparing data for rank test ... ")
		candidate_gene_set = Set(param_data.candidate_gene_list)
		candidate_gene_ls = []
		candidate_gene_pvalue_list = []
		non_candidate_gene_ls = []
		non_candidate_gene_pvalue_list = []
		if param_data.results_directory:	#given a directory where all results are.
			result_fname = os.path.join(param_data.results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		if not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		no_of_lines_to_read = getattr(param_data, 'no_of_top_lines', None)
		for row in reader:
			gene_id = int(row[col_name2index['gene_id']])
			score = float(row[col_name2index['score']])
			if gene_id in candidate_gene_set:
				candidate_gene_ls.append(gene_id)
				candidate_gene_pvalue_list.append(score)
			else:
				non_candidate_gene_ls.append(gene_id)
				non_candidate_gene_pvalue_list.append(score)
			counter += 1
			if no_of_lines_to_read is not None and counter>=no_of_lines_to_read:
				break
		del reader
		passingdata = PassingData(candidate_gene_ls=candidate_gene_ls, candidate_gene_pvalue_list=candidate_gene_pvalue_list,\
								non_candidate_gene_ls=non_candidate_gene_ls, non_candidate_gene_pvalue_list=non_candidate_gene_pvalue_list)
		sys.stderr.write("Done.\n")
		return passingdata
	
	def convertChrIndex2GenomeIndex(self, chr_index_ls, chr2cumu_no_of_snps, chr2no_of_snps):
		"""
		2008-10-11
			convert (chr, index) into genome_index
		"""
		index_ls = []
		for chr, chr_index in chr_index_ls:
			index_ls.append(chr2cumu_no_of_snps[chr]-chr2no_of_snps[chr]+chr_index)	#take SNPs on all previous chromosomes into account
		return index_ls
	
	def get_candidate_gene_snp_index_ls(self, candidate_gene_set, chr_pos_ls, snps_context_wrapper):
		"""
		2008-10-22
			moved from CheckCandidateGeneRank.py
		2008-10-21
			get index of all associated snps, which are supposedly stored in a giant list
		"""
		if self.debug:
			sys.stderr.write("Getting candidate_gene_snp_index_ls ...")
		candidate_gene_snp_index_ls = []
		non_candidate_gene_snp_index_ls = []
		
		for index in range(len(chr_pos_ls)):
			chr, pos = chr_pos_ls[index]
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					assign_snp_candidate_gene = 1
					break
			if assign_snp_candidate_gene:
				candidate_gene_snp_index_ls.append(index)
		candidate_gene_snp_index_ls = numpy.array(candidate_gene_snp_index_ls, numpy.int)
		if self.debug:
			sys.stderr.write("Done.\n")
		return candidate_gene_snp_index_ls
	
	def prepareDataForPermutationRankTest(self, rm, snps_context_wrapper, param_data):
		"""
		2008-10-26
			add another way to take top snps, based on the min_score cutoff.
		2008-10-23
			if no_of_top_lines is not defined in param_data, fetch the data_obj from gwr.data_obj_ls (preserving chromosome,position order)
			add more structures to return
		2008-10-22 pass param_data directly to getResultMethodContent()
		2008-10-17
			handle option need_the_value from param_data
		2008-10-15
			try to take no_of_snps from param_data, the snps are taken in order from starting_rank (get from param_data or just 1)
		2008-10-09
			prepare data for a permutation-based rank sum test
			the permutation preserves the order of SNPs within each chromosome
			1. shuffle chromosome order and put all SNPs into a list in order
			2. shift all SNPs in that list by a random number
			3. get rank sums of scores from fixed candidate gene SNP indices that 
			
		"""
		if self.debug:
			sys.stderr.write("Preparing data for permutation rank test ... ")
		candidate_gene_set = getattr(param_data, 'candidate_gene_set', Set(param_data.candidate_gene_list))
		
		chr2rank_ls = {}
		candidate_gene_snp_chr_pos_ls = []
		non_candidate_gene_snp_chr_pos_ls = []
		candidate_gene_snp_index_ls = []
		non_candidate_gene_snp_index_ls = []
		candidate_gene_snp_rank_ls = []
		non_candidate_gene_snp_rank_ls = []
		candidate_gene_snp_value_ls = []
		non_candidate_gene_snp_value_ls = []
		#chr2no_of_snps = {}
		#chr2cumu_no_of_snps = {}
		chr_pos_ls = []	#to store a list of (chr,position)
		
		#input file is in chromosome,position order
		param_data.construct_chr_pos2index = True	#always need chr_pos2index from now on
		genome_wide_result = self.getResultMethodContent(rm, param_data.results_directory, param_data.min_MAF, pdata=param_data)
		if genome_wide_result is None:
			return None
		score_ls = [data_obj.value for data_obj in genome_wide_result.data_obj_ls]
		import rpy
		rank_ls = rpy.r.rank(score_ls)	#rpy.rank also exists!!
		no_of_total_snps = len(score_ls)
		starting_rank = getattr(param_data, 'starting_rank', 1)	#2008-10-15 rank starts from 1
		no_of_snps = getattr(param_data, 'no_of_top_snps', no_of_total_snps)
		need_the_value = getattr(param_data, 'need_the_value', 0)	#get the pvalues/scores as well
		need_chr_pos_ls = getattr(param_data, 'need_chr_pos_ls', 0)	#get the pvalues/scores as well
		min_score = getattr(param_data, 'min_score', None)	#2008-10-25
		
		if min_score is not None:
			no_of_snps = no_of_total_snps
		else:
			ranks_to_be_checked = range(starting_rank, starting_rank+int(no_of_snps))	#2008-10-15 rank starts from 1
		rank_sum_stat = 0
		
		score_range = [None, None]
		for i in range(int(no_of_snps)):	#due to the duality(both score and rank gap) of MpiTopSNPTest.py's rank_gap, no_of_snps becomes float.
			if min_score is not None:	#2008-10-25
				data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)	#i starts from 0, rank starts from 1.
				data_obj_index = genome_wide_result.get_data_obj_index_given_rank(i+1)
				rank = rank_ls[data_obj_index]	#watch this, the index is not i
				if data_obj.value<min_score:	#score is below threshold, forget about it
					break
			elif no_of_snps!=no_of_total_snps:	#get top ranked snps
				rank_to_be_checked = ranks_to_be_checked[i]
				data_obj = genome_wide_result.get_data_obj_at_given_rank(rank_to_be_checked)
				data_obj_index = genome_wide_result.get_data_obj_index_given_rank(rank_to_be_checked)
				rank = rank_ls[data_obj_index]	#watch this, the index is not i
				
			else:	#get whole genome, not top rank thingy. loop in (chromosome, position) order
				data_obj = genome_wide_result.data_obj_ls[i]
				rank = rank_ls[i]
			
			chr_pos = (data_obj.chromosome, data_obj.position)
			if score_range[0] is None or data_obj.value>score_range[0]:	#maximum score
				score_range[0] = data_obj.value
				
			if score_range[1] is None or data_obj.value<score_range[1]:
				score_range[1] = data_obj.value
			
			chr = data_obj.chromosome
			if chr not in chr2rank_ls:
				chr2rank_ls[chr] = []
			chr2rank_ls[chr].append(rank)
			
			if need_chr_pos_ls:
				chr_pos_ls.append(chr_pos)
			
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					assign_snp_candidate_gene = 1
					if not param_data.allow_two_sample_overlapping:	#early break only if two samples are NOT allowed to be overlapping
						break
				else:
					assign_snp_non_candidate_gene = 1
					"""
					#2008-10-26, this condition is buggy, enabling this would cause SNPs with candidate genes nearby prematurely assigned to non-candidate category
					if not param_data.allow_two_sample_overlapping:	#early break only if two samples are NOT allowed to be overlapping
						break
					"""
				if assign_snp_candidate_gene==1 and assign_snp_non_candidate_gene==1:	#both are set to 1, no need to check first genes around
					break
			if assign_snp_candidate_gene:
				candidate_gene_snp_index_ls.append(genome_wide_result.chr_pos2index[chr_pos])
				candidate_gene_snp_chr_pos_ls.append(chr_pos)
				rank_sum_stat += rank
				candidate_gene_snp_rank_ls.append(rank)
				if need_the_value:
					candidate_gene_snp_value_ls.append(data_obj.value)
			
			#condition to assign this snp to non candidate gene
			assign_snp_non_candidate_gene = (not param_data.allow_two_sample_overlapping and not assign_snp_candidate_gene) or \
				(param_data.allow_two_sample_overlapping and assign_snp_non_candidate_gene)
			
			"""
			if not param_data.allow_two_sample_overlapping:	#if two samples are NOT allowed to be overlapping, only check if assign_snp_candidate_gene is set or not
				if not assign_snp_candidate_gene:
					non_candidate_gene_snp_chr_index_ls.append((chr, chr2no_of_snps[chr]-1))
					non_candidate_gene_snp_rank_ls.append(rank)
			"""
			if assign_snp_non_candidate_gene:	#two samples allowed overlapping, then have to check the other tag
				non_candidate_gene_snp_index_ls.append(genome_wide_result.chr_pos2index[chr_pos])
				non_candidate_gene_snp_chr_pos_ls.append(chr_pos)
				non_candidate_gene_snp_rank_ls.append(rank)
				if need_the_value:
					non_candidate_gene_snp_value_ls.append(data_obj.value)
		from common import get_chr_id2cumu_size
		#chr2cumu_no_of_snps, chr_gap, chr_id_ls = get_chr_id2cumu_size(genome_wide_result.chr2no_of_snps, chr_gap=0)
		
		#candidate_gene_snp_index_ls = self.convertChrIndex2GenomeIndex(candidate_gene_snp_chr_index_ls, chr2cumu_no_of_snps, chr2no_of_snps)
		#non_candidate_gene_snp_index_ls = self.convertChrIndex2GenomeIndex(non_candidate_gene_snp_chr_index_ls, chr2cumu_no_of_snps, chr2no_of_snps)
		
		#turn them into numpy arrays
		candidate_gene_snp_index_ls = num.array(candidate_gene_snp_index_ls, num.int)
		non_candidate_gene_snp_index_ls = num.array(non_candidate_gene_snp_index_ls, num.int)
		
		total_chr_pos_ar = None
		if need_chr_pos_ls:
			chr_pos_ls = num.array(chr_pos_ls)
			total_chr_pos_ls = genome_wide_result.chr_pos2index.keys()
			total_chr_pos_ls.sort()
			total_chr_pos_ar = num.array(total_chr_pos_ls)
		passingdata = PassingData(candidate_gene_snp_index_ls=candidate_gene_snp_index_ls,\
								non_candidate_gene_snp_index_ls=non_candidate_gene_snp_index_ls,\
								candidate_gene_snp_chr_pos_ls=candidate_gene_snp_chr_pos_ls,
								non_candidate_gene_snp_chr_pos_ls=non_candidate_gene_snp_chr_pos_ls,
								rank_sum_stat=rank_sum_stat,
								chr2rank_ls=chr2rank_ls, 
								chr2no_of_snps=genome_wide_result.chr2no_of_snps,\
								no_of_total_snps=no_of_total_snps,\
								no_of_snps=no_of_snps,\
								candidate_gene_snp_rank_ls=candidate_gene_snp_rank_ls,\
								non_candidate_gene_snp_rank_ls=non_candidate_gene_snp_rank_ls,\
								candidate_gene_snp_value_ls=candidate_gene_snp_value_ls,\
								non_candidate_gene_snp_value_ls=non_candidate_gene_snp_value_ls,\
								score_range=score_range,\
								chr_pos_ls=chr_pos_ls,\
								total_chr_pos_ar=total_chr_pos_ar)
		del genome_wide_result
		if self.debug:
			sys.stderr.write("Done.\n")
		return passingdata
	
	def getPermutationRankSumStat(self, chr2rank_ls, candidate_gene_snp_index_ls, non_candidate_gene_snp_index_ls, no_of_snps, \
								chr2no_of_snps, permutation_type=1):
		"""
		2008-10-11
			add option non_candidate_gene_snp_index_ls
			if overlapping is found between candidate_gene_snp_index_ls & non_candidate_gene_snp_index_ls, need to re-rank
		2008-10-09
			permutation_type
				0: no permutation
				1: shuffle chromosome order and shift all SNPs
				2: keep chromosome order and shift all SNPs
			
			the permutation preserves the order of SNPs within each chromosome
			1. shuffle chromosome order and put all SNPs into a list in order
			2. shift all SNPs in that list by a random number
			3. get rank sums of scores from fixed candidate gene SNP indices that 
		"""
		if self.debug:
			sys.stderr.write("Getting permutation rank sum stat ...")
		chr_id_ls = chr2rank_ls.keys()
		if permutation_type==0:
			chr_id_ls.sort()
		elif permutation_type==1:
			random.shuffle(chr_id_ls)
		elif permutation_type==2:
			chr_id_ls.sort()
		whole_rank_ar = num.zeros([no_of_snps], num.int)
		snp_start_index = 0
		for chr_id in chr_id_ls:
			snp_stop_index = snp_start_index + chr2no_of_snps[chr_id]
			whole_rank_ar[snp_start_index:snp_stop_index] = chr2rank_ls[chr_id]
			snp_start_index = snp_stop_index
		
		if permutation_type!=0:
			shift = random.randint(1, no_of_snps)
			candidate_gene_snp_index_ls_perm = (candidate_gene_snp_index_ls+shift)%no_of_snps	#modulo to recycle
			non_candidate_gene_snp_index_ls_perm = (non_candidate_gene_snp_index_ls+shift)%no_of_snps
		else:
			candidate_gene_snp_index_ls_perm = candidate_gene_snp_index_ls
			non_candidate_gene_snp_index_ls_perm = non_candidate_gene_snp_index_ls
		no_of_snps_in_candidate_genes = len(candidate_gene_snp_index_ls)
		no_of_snps_in_non_candidate_genes = len(non_candidate_gene_snp_index_ls)
		if no_of_snps_in_candidate_genes+no_of_snps_in_non_candidate_genes>no_of_snps:
			#two samples are allowed to have overlapping SNPs
			#have to rerank the whole thing to get new ranks
			import rpy
			rank_ls = num.hstack((whole_rank_ar[candidate_gene_snp_index_ls_perm], whole_rank_ar[non_candidate_gene_snp_index_ls_perm]))
			rank_ls = rpy.r.rank(rank_ls)	#rpy.rank also exists!!
			rank_sum_stat = sum(rank_ls[:no_of_snps_in_candidate_genes])
		else:
			rank_sum_stat = sum(whole_rank_ar[candidate_gene_snp_index_ls_perm])
		if self.debug:
			sys.stderr.write("Done.\n")
		return rank_sum_stat
	
	def getPermutationRankSumPvalue(self, chr2rank_ls, candidate_gene_snp_index_ls, non_candidate_gene_snp_index_ls, rank_sum_stat, no_of_snps, \
								chr2no_of_snps, permutation_type=1, no_of_permutations=10000, no_of_min_breaks=25):
		"""
		2008-10-11
			re-calculate rank_sum_stat (of the sample before permutation) due to possibility of overlapping in candidate_gene_snp_index_ls and non_candidate_gene_snp_index_ls
				,which causes re-ranking is needed in getPermutationRankSumStat()
			add option non_candidate_gene_snp_index_ls
		2008-10-09
			do permutatioins in a smarter way. break away early if no_of_hits>= no_of_min_breaks
		2008-10-09
		
		"""
		if self.debug:
			sys.stderr.write("Getting permutation rank sum pvalue ...")
		rank_sum_stat = self.getPermutationRankSumStat(chr2rank_ls, candidate_gene_snp_index_ls, non_candidate_gene_snp_index_ls,
								no_of_snps, chr2no_of_snps, permutation_type=0)
		i = 0
		no_of_hits = 0
		while i<no_of_permutations:
			rank_sum_stat_perm = self.getPermutationRankSumStat(chr2rank_ls, candidate_gene_snp_index_ls, non_candidate_gene_snp_index_ls,
								no_of_snps, chr2no_of_snps, permutation_type=permutation_type)
			if rank_sum_stat_perm>=rank_sum_stat:
				no_of_hits += 1
			i+=1
			if no_of_min_breaks>0 and no_of_hits>=no_of_min_breaks:	#if no_of_min_breaks<=0, no smart breaking
				break
		pvalue = no_of_hits/float(i)
		if self.debug:
			sys.stderr.write("%s/%s tests in total. Done.\n"%(no_of_hits, i))
		return rank_sum_stat, pvalue
	
	def output_gene_id2hit(self, gene_id2hit, output_fname):
		sys.stderr.write("Outputting gene_id2hit ... ")
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['gene-id', 'chromosome', 'position', 'snps-id', 'disp_pos', 'score', 'rank'])
		gene_id_ls = gene_id2hit.keys()
		gene_id_ls.sort()
		pvalue_ls = [gene_id2hit[gene_id].score for gene_id in gene_id_ls]
		import rpy
		rank_ls = rpy.r.rank(pvalue_ls)
		for i in range(len(gene_id_ls)):
			gene_id = gene_id_ls[i]
			rank = rank_ls[i]
			hit = gene_id2hit[gene_id]
			one_row = [gene_id, hit.chr, hit.pos, hit.snps_id, hit.disp_pos, hit.score, rank]
			writer.writerow(one_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def dealWithCandidateGeneList(self, list_type_id, return_set=False):
		"""
		2008-10-23
			add option return_set
		2008-10-15
			to make caching candidate gene list possible
		"""
		if list_type_id not in self.list_type_id2candidate_gene_list_info:	#internal cache
			candidate_gene_list = self.getGeneList(list_type_id)
			self.list_type_id2candidate_gene_list_info[list_type_id] = PassingData(candidate_gene_list=candidate_gene_list, candidate_gene_set=Set(candidate_gene_list))
		if return_set:
			return self.list_type_id2candidate_gene_list_info[list_type_id].candidate_gene_set
		else:
			return self.list_type_id2candidate_gene_list_info[list_type_id].candidate_gene_list
		
	
	list_type_id2candidate_gene_list_info = {}
	
	def run_wilcox_test(self, pd):
		"""
		2008-10-15
			more test_types
		2008-08-20
			use pd to summarize all parameters
		2008-08-15
			fix option max_pvalue_per_gene, to have the choice to take the most significant SNP associated with one gene or not
		2008-08-14
			add min_MAF
			just prepareDataForRankTest().
		2008-07-24
			fix a bug. still using self.list_type_id, rather than just list_type_id in one line
		2008-07-17
			split out as a standalone function so that MpiGeneListRankTest.py could call it more easily.
		"""
		if self.debug:
			sys.stderr.write("Running wilcox test ... ")
		if pd.results_type==1:
			ResultsClass = ResultsMethod
			TestResultClass = CandidateGeneRankSumTestResultMethod
			rm = ResultsClass.get(pd.results_id)
			min_distance = pd.min_distance
			min_MAF = pd.min_MAF
			get_closest = pd.get_closest
		elif pd.results_type==2:
			ResultsClass = ResultsByGene
			TestResultClass = CandidateGeneRankSumTestResult
			rm = ResultsClass.get(pd.results_id)
			min_distance = rm.min_distance
			min_MAF = rm.min_MAF
			get_closest = rm.get_closest
		else:
			sys.stderr.write("Invalid results type : %s.\n"%pd.results_type)
			return None
		
		if not rm:
			sys.stderr.write("No results method available for results_id=%s.\n"%pd.results_id)
			return None
		db_results = TestResultClass.query.filter_by(results_id=pd.results_id).\
			filter_by(list_type_id=pd.list_type_id).filter_by(min_distance=pd.min_distance).\
			filter_by(min_MAF=pd.min_MAF).filter_by(get_closest=pd.get_closest).filter_by(test_type=pd.test_type)
		if db_results.count()>0:	#done before
			db_result = db_results.first()
			sys.stderr.write("It's done already. id=%s, results_id=%s, list_type_id=%s, pvalue=%s, statistic=%s.\n"%\
							(db_result.id, db_result.results_id, db_result.list_type_id, db_result.pvalue, db_result.statistic))
			return None
		try:
			import rpy
			candidate_gene_list = self.dealWithCandidateGeneList(pd.list_type_id)	#internal cache
			if pd.test_type>3:
				allow_two_sample_overlapping = 1
			else:
				allow_two_sample_overlapping = 0
			param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, \
									min_MAF=pd.min_MAF, allow_two_sample_overlapping=allow_two_sample_overlapping)
			if pd.results_type==2:
				passingdata = self.prepareDataForRankTestFromResultsByGene(rm, param_data)
				candidate_sample_size = len(passingdata.candidate_gene_pvalue_list)
				non_candidate_sample_size = len(passingdata.non_candidate_gene_pvalue_list)
				if candidate_sample_size>=pd.min_sample_size and non_candidate_sample_size>=pd.min_sample_size:	#2008-08-14
					w_result = rpy.r.wilcox_test(passingdata.candidate_gene_pvalue_list, passingdata.non_candidate_gene_pvalue_list, alternative='greater', conf_int=rpy.r.TRUE)
				else:
					sys.stderr.write("Ignore. sample size less than %s. %s vs %s.\n"%(pd.min_sample_size, candidate_sample_size, non_candidate_sample_size))
					return None
				statistic=w_result['statistic']['W']
				pvalue=w_result['p.value']
			elif pd.results_type==1:	#for ResultsMethod
				permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
				candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
				non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
				not_enough_sample = 0
				if candidate_sample_size<pd.min_sample_size or non_candidate_sample_size<pd.min_sample_size:
					sys.stderr.write("Ignore. sample size less than %s. %s vs %s.\n"%(pd.min_sample_size, candidate_sample_size, non_candidate_sample_size))
					return None
				if pd.test_type%3==1:
					w_result = rpy.r.wilcox_test(permData.candidate_gene_snp_rank_ls, permData.non_candidate_gene_snp_rank_ls, alternative='greater')
					statistic=w_result['statistic']['W']
					pvalue=w_result['p.value']
				elif pd.test_type%3==2 or pd.test_type%3==0:
					if pd.test_type%3==2:
						pd.permutation_type = 1
					elif pd.test_type%3==0:
						pd.permutation_type = 2
					rank_sum_stat, pvalue = self.getPermutationRankSumPvalue(permData.chr2rank_ls, permData.candidate_gene_snp_index_ls, permData.non_candidate_gene_snp_index_ls,\
												permData.rank_sum_stat, \
												permData.no_of_snps, permData.chr2no_of_snps, permutation_type=pd.permutation_type, \
												no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks)
					statistic = rank_sum_stat-(candidate_sample_size*(candidate_sample_size-1))/2.
				else:
					sys.stderr.write("Test_type %s not supported.\n"%(pd.test_type))
					return None
			else:
				sys.stderr.write("Results_type %s not supported.\n"%(pd.results_type))
				return None
		except:
			sys.stderr.write("Exception happened for results_id=%s, list_type_id=%s.\n"%(pd.results_id, pd.list_type_id))
			traceback.print_exc()
			sys.stderr.write('%s.\n'%repr(sys.exc_info()))
			return None
		candidate_gene_rank_sum_test_result = TestResultClass(list_type_id=pd.list_type_id, statistic=statistic,\
																			pvalue=pvalue)
		candidate_gene_rank_sum_test_result.results_id = pd.results_id
		candidate_gene_rank_sum_test_result.min_distance = min_distance
		candidate_gene_rank_sum_test_result.min_MAF = min_MAF
		candidate_gene_rank_sum_test_result.get_closest = get_closest
		candidate_gene_rank_sum_test_result.candidate_sample_size = candidate_sample_size
		candidate_gene_rank_sum_test_result.non_candidate_sample_size = non_candidate_sample_size
		candidate_gene_rank_sum_test_result.test_type = pd.test_type
		candidate_gene_rank_sum_test_result.max_pvalue_per_gene = 0	#always set to 0 in order to be compatible with previous approaches
		if self.debug:
			sys.stderr.write("Done.\n")
		return candidate_gene_rank_sum_test_result
	
	def dealWithSnpsContextWrapper(self, snps_context_picklef, min_distance, get_closest):
		"""
		2008-09-08
			split out of run()
			
			--constructDataStruc()
		"""
		sys.stderr.write("Dealing with snps_context_wrapper ...")
		if snps_context_picklef:
			if os.path.isfile(snps_context_picklef):	#if this file is already there, suggest to un-pickle it.
				picklef = open(snps_context_picklef)
				snps_context_wrapper = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle snps_context_wrapper into it
				snps_context_wrapper = self.constructDataStruc(min_distance, get_closest)
				#2008-09-07 pickle the snps_context_wrapper object
				picklef = open('%s_g%s_m%s'%(snps_context_picklef, get_closest, min_distance), 'w')
				cPickle.dump(snps_context_wrapper, picklef, -1)
				picklef.close()
		else:
			snps_context_wrapper = self.constructDataStruc(min_distance, get_closest)
		sys.stderr.write("Done.\n")
		return snps_context_wrapper
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		
		if getattr(self, 'output_fname', None):
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			header_row = []
			for column in CandidateGeneRankSumTestResult.c.keys():
				header_row.append(column)
			writer.writerow(header_row)
		else:
			writer = None
		pd = PassingData(list_type_id=self.list_type_id, snps_context_wrapper=snps_context_wrapper,\
							results_directory=self.results_directory, \
							min_MAF=self.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance, \
							min_sample_size=self.min_sample_size, test_type=self.test_type, \
							results_type=self.results_type, no_of_permutations=self.no_of_permutations,\
							no_of_min_breaks=self.no_of_min_breaks)
		for results_id in self.results_id_ls:
			pd.results_id = results_id
			candidate_gene_rank_sum_test_result = self.run_wilcox_test(pd)
			if candidate_gene_rank_sum_test_result is not None:
				row = []
				for column in candidate_gene_rank_sum_test_result.c.keys():
					row.append(getattr(candidate_gene_rank_sum_test_result, column))
					print '%s: %s'%(column, row[-1])
				if writer:
					writer.writerow(row)
				#session.save(candidate_gene_rank_sum_test_result)
				if self.commit:
					session.flush()
		#if self.commit:
		#	session.flush()
		#	session.commit()
		#print passingdata.candidate_gene_pvalue_list
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GeneListRankTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
