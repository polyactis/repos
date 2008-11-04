#!/usr/bin/env python
"""
Examples:
	#debug, quick testing
	TopSNPTest.py -e 389 -l 1 -u yh  -o /tmp/top_snp_test.out -b
	
	TopSNPTest.py -e 389 -l 1 -u yh  -o /tmp/top_snp_test.out
	
	#
	TopSNPTest.py -e 389 -l 28 -u yh -f 300 -s ./mnt2/panfs/250k/snps_context_g0_m20000 -m 20000 -y 15
	
Description:
	2008-08-20 program to do hypergeometric test on a number of genes from top SNPs based on a given candidate gene list.
	
	It verifies against the db several things:
	1. whether the results_method_id is available
	2. whether results_method_type_id is 1 (association)
	3. whether the same (results_method_id, list_type_id) pair has been in the top_snp_test table.
	
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
import Stock_250kDB
from Stock_250kDB import Snps, SnpsContext, ResultsMethod, GeneList, CandidateGeneRankSumTestResult, CandidateGeneTopSNPTest, ResultsByGene
from Results2DB_250k import Results2DB_250k
from pymodule import getGenomeWideResultFromFile
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
from sets import Set
from heapq import heappush, heappop, heapreplace
from common import get_total_gene_ls
import rpy, random, numpy

class TopSNPTest(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict.update({('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.']})
	option_default_dict.update({('gene_table', 1, ): ['genome.gene', 'a', 1, 'to get the number of total genes from database, which table.']})
	option_default_dict.update({('no_of_top_snps', 1, int): [200, 'f', 1, 'how many number of top snps based on score or -log(pvalue).']})
	option_default_dict.update({("test_type_id", 1, int): [15, 'y', 1, 'which type of tests. check db table analysis_method. likely be 14,15 etc.']})
	option_default_dict.update({('min_score', 0, float): [None, 'M', 1, 'alternative way to get top snps. if specified, no_of_top_snps is ignored.']})
	
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		GeneListRankTest.__init__(self, **keywords)
	
	def getNoOfTotalGenes(self, db, gene_table='genome.gene', tax_id=3702):
		"""
		2008-10-21
			add following condition to sql
				no mitochondrial chromosome, chromosome has to be known,
		"""
		if self.debug:
			sys.stderr.write("Getting no of total genes ... ")
		rows = db.metadata.bind.execute("select count(gene_id) as count from %s where tax_id=%s and chromosome is not null and chromosome!='MT'"%(gene_table, tax_id))
		row = rows.fetchone()
		if self.debug:
			sys.stderr.write("Done.\n")
		return row.count
	
	def prepareDataForHGTest(self, rm, snps_context_wrapper, candidate_gene_list, results_directory, min_MAF, no_of_top_snps):
		"""
		2008-08-20
		"""
		sys.stderr.write("Preparing data for HG test ... ")
		genome_wide_result = self.getResultMethodContent(rm, results_directory, min_MAF)
		genome_wide_result.data_obj_ls.sort()	#in value descending order. each SNP object has a defined method for comparison based on its value
		genome_wide_result.data_obj_ls.reverse()
		candidate_gene_set = Set(candidate_gene_list)
		candidate_gene_in_top_set = Set([])
		non_candidate_gene_in_top_set = Set([])
		for i in range(no_of_top_snps):
			data_obj = genome_wide_result.data_obj_ls[i]
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					candidate_gene_in_top_set.add(gene_id)
				else:
					non_candidate_gene_in_top_set.add(gene_id)
		passingdata = PassingData(candidate_gene_in_top_set=candidate_gene_in_top_set, non_candidate_gene_in_top_set=non_candidate_gene_in_top_set)
		sys.stderr.write("Done.\n")
		return passingdata
	
	def prepareDataForHGTest_SNPPair(self, rm, snps_context_wrapper, candidate_gene_list, results_directory, min_MAF, no_of_top_snps):
		"""
		2008-09-08
			a different format in the results file, which is output of MpiIntraGeneSNPPairAsso.py
		"""
		sys.stderr.write("Preparing data for HG test from SNP pair results ... ")
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		reader = csv.reader(open(result_fname), delimiter='\t')
		reader.next()
		heap_ls = []
		for row in reader:
			gene_id = int(row[1])
			pvalue = float(row[5])
			count1 = int(row[-2])
			count2 = int(row[-1])
			maf = float(min(count1, count2))/(count1+count2)
			if maf>=min_MAF:
				if len(heap_ls)<no_of_top_snps:
					heappush(heap_ls, (-pvalue, gene_id))
				else:
					heapreplace(heap_ls, (-pvalue, gene_id))
		candidate_gene_set = Set(candidate_gene_list)
		candidate_gene_in_top_set = Set([])
		non_candidate_gene_in_top_set = Set([])
		while heap_ls:
			pvalue_gene_id = heappop(heap_ls)
			gene_id = pvalue_gene_id[1]
			if gene_id in candidate_gene_set:
				candidate_gene_in_top_set.add(gene_id)
			else:
				non_candidate_gene_in_top_set.add(gene_id)
		passingdata = PassingData(candidate_gene_in_top_set=candidate_gene_in_top_set, non_candidate_gene_in_top_set=non_candidate_gene_in_top_set)
		sys.stderr.write("Done.\n")
		return passingdata
	
	def dealWithNoOfSNPsAssociatedWithCandidateGeneList(self, list_type_id, rm, pd):
		"""
		2008-10-15
			it's designed for hypergeometric test to get the size of the black urn.
			get the number of snps associated with the whole candidate_gene_list.
		"""
		candidate_gene_list = self.dealWithCandidateGeneList(list_type_id)	#internal cache
		if getattr(self.list_type_id2candidate_gene_list_info[list_type_id], 'no_of_snps_associated', None) is None:
			"""
			if pd.test_type>3:
				allow_two_sample_overlapping = 1
			else:
				allow_two_sample_overlapping = 0
			"""
			param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, \
										min_MAF=pd.min_MAF, allow_two_sample_overlapping=pd.allow_two_sample_overlapping)	#starting_rank=1, go through all SNPs
			permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
			candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_associated=candidate_sample_size
			self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_unassociated=non_candidate_sample_size			
		return self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_associated
	
	def get_looped_chr_pos_ls(self, candidate_gene_snp_index_ls, no_of_snps, total_chr_pos_ar, shift=None):
		"""
		2008-10-30
			add option shift
		2008-10-22
			total_chr_pos_ar is of numpy array type
		"""
		if shift is None:
			shift = random.randint(1, no_of_snps)
		candidate_gene_snp_index_ls_perm = (candidate_gene_snp_index_ls+shift)%no_of_snps	#modulo to recycle
		looped_chr_pos_ar = total_chr_pos_ar[candidate_gene_snp_index_ls_perm]
		return looped_chr_pos_ar
		
	def get_enrichment_pvalue_by_gw_looping(self, candidate_sample_size, top_snp_index_ls, candidate_gene_set, \
										snps_context_wrapper, \
										no_of_total_snps, total_chr_pos_ar=None, no_of_permutations=20000, no_of_min_breaks=30):
		"""
		2008-10-30
		2008-10-22
			get enrichment pvalue by genome-wide looping of SNP positions. a permutation to preserve LD.
		"""
		if self.debug:
			sys.stderr.write("Getting enrichment pvalue by gw-looping ... ")
		i = 0
		no_of_hits = 0
		while i<no_of_permutations:
			looped_chr_pos_ls = self.get_looped_chr_pos_ls(top_snp_index_ls, no_of_total_snps, total_chr_pos_ar)
			looped_candidate_gene_snp_index_ls = self.get_candidate_gene_snp_index_ls(candidate_gene_set, \
																					looped_chr_pos_ls, snps_context_wrapper)
			new_candidate_sample_size = len(looped_candidate_gene_snp_index_ls)
			if new_candidate_sample_size>=candidate_sample_size:	#pvalue = Prob(X>=candidate_sample_size)
				no_of_hits += 1
			i+=1
			if no_of_min_breaks>0 and no_of_hits>=no_of_min_breaks:	#if no_of_min_breaks<=0, no smart breaking
				break
		pvalue = no_of_hits/float(i)
		return_data = PassingData(pvalue=pvalue, no_of_tests=i, no_of_tests_passed=no_of_hits)
		if self.debug:
			sys.stderr.write("%s/%s tests in total. Done.\n"%(no_of_hits, i))
		return return_data
	
	def get_enrichment_pvalue_by_random_gene_list(self, sample_pvalue, total_gene_id_ls, candidate_gene_set, \
												total_chr_pos_ls, snps_context_wrapper, top_snp_chr_pos_ls, n,k,\
												no_of_permutations=20000, no_of_min_breaks=30):
		"""
		2008-10-22
		"""
		if self.debug:
			sys.stderr.write("Getting enrichment pvalue by random gene list ... ")
		i = 0
		no_of_hits = 0
		no_of_candidate_genes = len(candidate_gene_set)
		no_of_total_snps = len(total_chr_pos_ls)
		no_of_top_snps = len(top_snp_chr_pos_ls)
		while i<no_of_permutations:
			random_candidate_gene_set = Set(random.sample(total_gene_id_ls, no_of_candidate_genes))
			random_candidate_gene_snp_gw_index_ls = self.get_candidate_gene_snp_index_ls(random_candidate_gene_set, total_chr_pos_ls, snps_context_wrapper)
			random_candidate_gene_snp_sample_index_ls = self.get_candidate_gene_snp_index_ls(random_candidate_gene_set, top_snp_chr_pos_ls, snps_context_wrapper)
			x = len(random_candidate_gene_snp_sample_index_ls)
			m = len(random_candidate_gene_snp_gw_index_ls)
			n = no_of_total_snps - m
			k = no_of_top_snps
			new_sample_pvalue = rpy.r.phyper(x-1,m,n,k, lower_tail = rpy.r.FALSE)
			if new_sample_pvalue<=sample_pvalue:	#watch: pvalue = Prob(X<=sample_pvalue). chance of getting more significant (smaller) pvalues
				no_of_hits += 1
			i+=1
			if no_of_min_breaks>0 and no_of_hits>=no_of_min_breaks:	#if no_of_min_breaks<=0, no smart breaking
				break
		pvalue = no_of_hits/float(i)
		return_data = PassingData(pvalue=pvalue, no_of_tests=i, no_of_tests_passed=no_of_hits)
		if self.debug:
			sys.stderr.write("%s/%s tests in total. Done.\n"%(no_of_hits, i))
		return return_data
	
	def TestResultClassInDB(self, TestResultClass, results_id, list_type_id,starting_rank, type_id, min_distance, \
						no_of_top_snps, min_score=None):
		"""
		2008-10-30
			split out from runHGTest()
			check if a result is in db or not
		"""
		query = TestResultClass.query.filter_by(results_id=results_id).\
				filter_by(list_type_id=list_type_id).\
				filter_by(starting_rank=starting_rank).\
				filter_by(type_id=type_id).\
				filter_by(min_distance=min_distance)
		#2008-10-27 better to use integer no_of_top_snps to match db entry. float min_score is inaccurate.
		if min_score is not None:	#if min_score is present, check db differently. pd.no_of_top_snps doesn't correspond to SNPs whose score is above min_score.
			db_results = query.filter(TestResultClass.min_score>=min_score-0.0001).filter(TestResultClass.min_score<=min_score+0.0001)
		else:
			db_results = query.filter_by(no_of_top_snps=int(no_of_top_snps))
		return db_results
	
	def returnResultFromDB(self, TestResultClass, results_id, list_type_id,starting_rank, type_id, min_distance, \
						no_of_top_snps, min_score=None):
		"""
		2008-10-30
		"""
		result = None
		db_results = self.TestResultClassInDB(TestResultClass, results_id, list_type_id,starting_rank, type_id, min_distance, \
											no_of_top_snps, min_score)
		
		if db_results.count()==1:	#done before
			result = db_results.first()
		elif db_results.count()>1:
			result = db_results.first()
			sys.stderr.write("More than 1 result matched with results_id=%s, list_type_id=%s, starting_rank=%s,\
							type_id=%s, min_distance=%s, no_of_top_snps=%s, min_score=%s. First (id=%s) taken.\n"%\
							(results_id, list_type_id, starting_rank,
							type_id, min_distance, no_of_top_snps, min_score, result.id))
		return result
	
	def getTestResult(self, session, rm, TestResultClass, pd):
		"""
		2008-10-30
			get one instance of CandidateGeneTopSNPTestRM for CandidateGeneTopSNPTestRMNullData
			
			if it's in db, fetch it
			it not, generate one
		"""
		starting_rank = getattr(pd, 'starting_rank', 1)	#from which rank to look down for top snps
		commit = getattr(pd, 'commit', 0)	#2008-10-30 save objects right away
		need_permData = getattr(pd, 'need_permData', 0)	#2008-10-30 a flag indicating whether permData is needed
		
		return_data = PassingData(result=None, permData=None)
		
		candidate_gene_set = self.dealWithCandidateGeneList(pd.list_type_id, return_set=True)	#internal cache
		param_data = PassingData(results_directory=pd.results_directory, candidate_gene_set=candidate_gene_set, \
									min_MAF=pd.min_MAF, allow_two_sample_overlapping=pd.allow_two_sample_overlapping,\
									no_of_top_lines=pd.no_of_top_snps, no_of_top_snps=pd.no_of_top_snps, \
									starting_rank=starting_rank, \
									min_score=pd.min_score,\
									need_chr_pos_ls=1, )	#need_chr_pos_ls is for null_distribution_type_id=2/3
		
		result = None
		if pd.type_id:	#the type is already in database. now check if the same top snp test has been done or not
			result = self.returnResultFromDB(TestResultClass, pd.results_id, pd.list_type_id, starting_rank,
											pd.type_id, pd.min_distance, pd.no_of_top_snps, pd.min_score)
		if result is None:	#2008-10-30 it's still None, need to generate the observed data	
			permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
			if permData is None:
				if self.debug:
					sys.stderr.write("No permData from prepareDataForPermutationRankTest().\n")
				return return_data
			
			score_range = permData.score_range
			candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			if candidate_sample_size<pd.min_sample_size or non_candidate_sample_size<pd.min_sample_size:	#don't look at how many non-candidates are (different from whole genome rank test)
				if self.debug:
					sys.stderr.write("Ignore. sample size less than %s. %s vs %s.\n"%(pd.min_sample_size, candidate_sample_size, non_candidate_sample_size))
				return_data.permData = permData
				return return_data
			if pd.type_id and pd.min_score is not None:
				#2008-10-28 check db again because no_of_top_snps=candidate_sample_size+non_candidate_sample_size.
				#sometimes, different min_score cutoff gives same no_of_top_snps.
				result = self.returnResultFromDB(TestResultClass, pd.results_id, pd.list_type_id, starting_rank,
												pd.type_id, pd.min_distance, candidate_sample_size+non_candidate_sample_size)
			if result is None:
				import rpy
				
				candidate_gw_size = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
				x = candidate_sample_size
				m = candidate_gw_size
				n = permData.no_of_total_snps - m
				k = candidate_sample_size + non_candidate_sample_size
				
				
				result = TestResultClass(list_type_id=pd.list_type_id)	#2008-10-30 no pvalue
				result.pvalue = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
				result.results_id = pd.results_id
				result.candidate_sample_size = candidate_sample_size
				result.non_candidate_sample_size = non_candidate_sample_size
				result.min_distance = pd.min_distance
				result.candidate_gw_size = candidate_gw_size
				result.non_candidate_gw_size = permData.no_of_total_snps - candidate_gw_size
				result.no_of_top_snps = candidate_sample_size + non_candidate_sample_size
				result.starting_rank = starting_rank
				result.type_id = pd.type_id
				if score_range:
					result.max_score = score_range[0]
					if pd.min_score is not None:
						result.min_score = pd.min_score
					else:
						result.min_score = score_range[1]
				session.save(result)
				if commit:
					try:
						session.flush()
					except:
						session.expunge(result)
						sys.stderr.write("Exception happened for results_method_id=%s, list_type_id=%s.\n"%(result.results_id, result.list_type_id))
						for column in result.c.keys():
							sys.stderr.write("\t%s=%s.\n"%(column, getattr(result, column)))
						traceback.print_exc()
						sys.stderr.write('%s.\n'%repr(sys.exc_info()))
						result = None
		elif need_permData:
			permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
		else:
			permData = None
		return_data = PassingData(result=result, permData=permData)
		return return_data
	
	def runEnrichmentTestToGetNullData(self, session, pd):
		"""
		2008-10-30
			run enrichment test, also to get NULL data based on either null distribution
		"""
		if self.debug:
			sys.stderr.write("Running GW-looping test on results_id=%s, list_type_id=%s, no_of_top_snps=%s, type_id=%s, min_score=%s, ... "%\
							(getattr(pd, 'results_id',-1), getattr(pd, 'list_type_id', -1), getattr(pd, 'no_of_top_snps', -1),\
							getattr(pd, 'type_id', -1), getattr(pd, 'min_score', -1)))
		
		ResultsClass = ResultsMethod
		TestResultClass = Stock_250kDB.CandidateGeneTopSNPTestRM
		rm = ResultsClass.get(pd.results_id)
		min_distance = pd.min_distance
		min_MAF = pd.min_MAF
		get_closest = pd.get_closest
		
		no_of_top_snps_ls = getattr(pd, 'no_of_top_snps_ls', [])
		min_score_ls = getattr(pd, 'min_score_ls', [])
		if no_of_top_snps_ls:
			cutoff_ls = no_of_top_snps_ls
			cutoff_type = 1
		else:
			cutoff_ls = min_score_ls
			cutoff_type = 2
		commit = getattr(pd, 'commit', 0)	#2008-10-30 save objects right away
		
		if not rm:
			sys.stderr.write("No results available for results_id=%s.\n"%pd.results_id)
			return None
		
		
		
		candidate_gene_set = self.dealWithCandidateGeneList(pd.list_type_id, return_set=True)	#internal cache
		no_of_candidate_genes = len(candidate_gene_set)
		no_of_total_snps = None	#same for all tests from the same rm
		null_data_ls = []
		if pd.null_distribution_type_id==2 or pd.null_distribution_type_id==3:
			pd.need_permData = 1	#need to permData in getTestResult() even when the result is directly found in the database
			for i in range(pd.no_of_permutations):
				if pd.null_distribution_type_id==2:
					if no_of_total_snps is None:	#need to get it from the file
						if cutoff_type==1:
							pd.no_of_top_snps = cutoff_ls[0]
						elif cutoff_type==2:
							pd.min_score = cutoff_ls[0]
						return_data = self.getTestResult(session, rm, TestResultClass, pd)
						if return_data.permData:	#permData could be None
							no_of_total_snps = return_data.permData.no_of_total_snps
						else:
							if self.debug:
								sys.stderr.write("Warning: No permData from getTestResult(). aborted.\n")
							break
					shift = random.randint(1, no_of_total_snps)
					run_no = shift	#use this to link all NULL data under different no_of_top_snps together
				else:
					if no_of_candidate_genes>len(pd.total_gene_id_ls):
						if self.debug:
							sys.stderr.write("no_of_candidate_genes %s is bigger than no_of_total_genes, %s.\n"%\
											(no_of_candidate_genes, len(pd.total_gene_id_ls)))
						break
					random_candidate_gene_ls = random.sample(pd.total_gene_id_ls, no_of_candidate_genes)
					run_no = sum(random_candidate_gene_ls)%1000000	#take the sum of all gene ids and modulo to make it under 1 million. very little chance any two random gene lists have identical this number.
					random_candidate_gene_set = Set(random_candidate_gene_ls)
					random_candidate_gene_snp_gw_index_ls = None	#set it to None before every permutation
				for cutoff in cutoff_ls:
					if cutoff_type==1:
						pd.no_of_top_snps = cutoff
					elif cutoff_type==2:
						pd.min_score = cutoff
					return_data = self.getTestResult(session, rm, TestResultClass, pd)
					result = return_data.result
					permData = return_data.permData
					if return_data.result:
						"""
						#2008-11-04	doesn't care repeating run_no. the chance is small, however, this incurs a huge load on db server.
						rows = Stock_250kDB.TopSNPTestRMNullData.query.\
									filter_by(observed_id=result.id).\
									filter_by(run_no=run_no).\
									filter_by(null_distribution_type_id=pd.null_distribution_type_id)
						if rows.count()>0:
							if self.debug:
								sys.stderr.write("null data for observed_id=%s, run_no=%s, null_distribution_type_id=%s already in db.\n"%\
												(result.id, run_no, pd.null_distribution_type_id))
							continue
						"""
						if len(result.null_data_ls)>pd.no_of_permutations:	#skip if it's too many already
							continue
						
						if pd.null_distribution_type_id==2:
							top_snp_index_ls = numpy.hstack((permData.candidate_gene_snp_index_ls, permData.non_candidate_gene_snp_index_ls))
											
							looped_chr_pos_ls = self.get_looped_chr_pos_ls(top_snp_index_ls, permData.no_of_total_snps, permData.total_chr_pos_ar, \
																	shift=shift)
							looped_candidate_gene_snp_index_ls = self.get_candidate_gene_snp_index_ls(candidate_gene_set, \
																						looped_chr_pos_ls, \
																						pd.snps_context_wrapper)
							new_candidate_sample_size = len(looped_candidate_gene_snp_index_ls)
							new_candidate_gw_size = result.candidate_gw_size	#same as observed
						else:
							top_snp_chr_pos_ls = permData.candidate_gene_snp_chr_pos_ls + permData.non_candidate_gene_snp_chr_pos_ls
							
							if random_candidate_gene_snp_gw_index_ls is None:	#2008-10-31 if it's None, generate it. same for every simulation
								random_candidate_gene_snp_gw_index_ls = self.get_candidate_gene_snp_index_ls(random_candidate_gene_set,\
																				permData.total_chr_pos_ar, pd.snps_context_wrapper)
							random_candidate_gene_snp_sample_index_ls = self.get_candidate_gene_snp_index_ls(random_candidate_gene_set, \
																				top_snp_chr_pos_ls, pd.snps_context_wrapper)
							new_candidate_sample_size = len(random_candidate_gene_snp_sample_index_ls)
							new_candidate_gw_size = len(random_candidate_gene_snp_gw_index_ls)
						null_data = Stock_250kDB.TopSNPTestRMNullData(observed=result,\
																	candidate_sample_size=new_candidate_sample_size,\
																	candidate_gw_size=new_candidate_gw_size,\
																	run_no=run_no,\
																	null_distribution_type_id=pd.null_distribution_type_id)
						session.save(null_data)
						if commit:
							session.flush()
		elif pd.null_distribution_type_id==1:
			for cutoff in cutoff_ls:
				if cutoff_type==1:
					pd.no_of_top_snps = cutoff
				elif cutoff_type==2:
					pd.min_score = cutoff
				return_data = self.getTestResult(session, rm, TestResultClass, pd)
		
		else:
			sys.stderr.write("null_distribution_type %s not supported.\n"%(pd.null_distribution_type_id))
			return None
		if self.debug:
			sys.stderr.write("Done.\n")
		return None
		
	def runHGTest(self, pd):
		"""
		2008-10-30
			use TestResultClassInDB() to check if result has been stored in the db or not.
			if pd.min_score is not None, run TestResultClassInDB() again.
		2008-10-26
			handle the situation that uses min_score to get top snps to do enrichment test
		2008-10-23
			'test_type' = 'test_type_id'
			test_type_id is now formally put into db table analysis_method. 15 is functional now.
			add code to handle null_distribution_type_id=1,2,3
		2008-10-15
			add test_type 1,2,3 and more
		2008-09-18
			when checking the db whether a top snp test result is there or not, add no_of_top_snps as one more filter
		2008-09-16
			robust programming, check if prepareDataForRankTestFromResultsByGene() returns anything before doing stat test
		2008-09-16
			now get results from table results_by_gene, instead of results_method.
			results are directly linked to gene. no_of_top_snps is same as the number of top genes/lines.
		2008-09-09
			call prepareDataForHGTest_SNPPair if analysis_method_id==13
		2008-08-20
		"""
		if self.debug:
			sys.stderr.write("Running hypergeometric test on results_id=%s, list_type_id=%s, no_of_top_snps=%s, type_id=%s, min_score=%s, ... "%\
							(getattr(pd, 'results_id',-1), getattr(pd, 'list_type_id', -1), getattr(pd, 'no_of_top_snps', -1),\
							getattr(pd, 'type_id', -1), getattr(pd, 'min_score', -1)))
		if pd.results_type==1:
			ResultsClass = ResultsMethod
			TestResultClass = Stock_250kDB.CandidateGeneTopSNPTestRM
			rm = ResultsClass.get(pd.results_id)
			min_distance = pd.min_distance
			min_MAF = pd.min_MAF
			get_closest = pd.get_closest
		elif pd.results_type==2:
			ResultsClass = ResultsByGene
			TestResultClass = CandidateGeneTopSNPTest
			rm = ResultsClass.get(pd.results_id)
			min_distance = rm.min_distance
			min_MAF = rm.min_MAF
			get_closest = rm.get_closest
		else:
			sys.stderr.write("Invalid results type : %s.\n"%pd.results_type)
			return None
		
		if not rm:
			sys.stderr.write("No results available for results_id=%s.\n"%pd.results_id)
			return None
		starting_rank = getattr(pd, 'starting_rank', 1)	#from which rank to look down for top snps
		if pd.type_id:	#the type is already in database. now check if the same top snp test has been done or not
			db_results = self.TestResultClassInDB(TestResultClass, pd.results_id, pd.list_type_id, starting_rank,
												pd.type_id, pd.min_distance, pd.no_of_top_snps, pd.min_score)
			
			if db_results.count()>0:	#done before
				db_result = db_results.first()
				if self.debug:
					sys.stderr.write("It's done already. id=%s, results_id=%s, list_type_id=%s, no_of_top_snps=%s, pvalue=%s.\n"%\
								(db_result.id, db_result.results_id, db_result.list_type_id, db_result.no_of_top_snps, db_result.pvalue))
				return None
		
		score_range = None	#to be added as comment in dbs		
		m = None
		n = None
		no_of_tests_passed = None
		no_of_tests = None
		try:
			import rpy
			candidate_gene_list = self.dealWithCandidateGeneList(pd.list_type_id)	#internal cache
			"""
			if rm.analysis_method_id==13:
				passingdata = self.prepareDataForHGTest_SNPPair(rm, pd.snps_context_wrapper, candidate_gene_list, pd.results_directory, pd.min_MAF, pd.no_of_top_snps)
			else:
				passingdata = self.prepareDataForHGTest(rm, pd.snps_context_wrapper, candidate_gene_list, pd.results_directory, pd.min_MAF, pd.no_of_top_snps)
			"""
			param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, \
									min_MAF=pd.min_MAF, allow_two_sample_overlapping=pd.allow_two_sample_overlapping,\
									no_of_top_lines=pd.no_of_top_snps, no_of_top_snps=pd.no_of_top_snps, \
									starting_rank=starting_rank, \
									min_score=pd.min_score,\
									need_chr_pos_ls=1, )	#need_chr_pos_ls is for null_distribution_type_id=2/3
			#param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, no_of_top_lines=pd.no_of_top_snps)
			if pd.results_type==2:
				passingdata = self.prepareDataForRankTestFromResultsByGene(rm, param_data)
				if passingdata is None:
					sys.stderr.write("No data got from this result (id=%s, filename=%s).\n"%(rm.id, rm.filename))
					return None
				
				x = len(passingdata.candidate_gene_ls)
				m = len(candidate_gene_list)
				n = pd.no_of_total_genes - m
				k = x + len(passingdata.non_candidate_gene_ls)
				pvalue = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
				candidate_sample_size = x
				non_candidate_sample_size = len(passingdata.non_candidate_gene_ls)
			elif pd.results_type==1:	#for ResultsMethod
				permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
				if permData is None:
					if self.debug:
						sys.stderr.write("No permData from prepareDataForPermutationRankTest().\n")
					return None
				
				score_range = permData.score_range
				candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
				non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
				
				if pd.min_score is not None:
					#2008-10-28 check db again because no_of_top_snps=candidate_sample_size+non_candidate_sample_size.
					#sometimes, different min_score cutoff gives same no_of_top_snps.
					db_results = self.TestResultClassInDB(TestResultClass, pd.results_id, pd.list_type_id, starting_rank,
													pd.type_id, pd.min_distance, candidate_sample_size+non_candidate_sample_size)
					if db_results.count()>0:	#done before
						db_result = db_results.first()
						if self.debug:
							sys.stderr.write("It's done already. id=%s, results_id=%s, list_type_id=%s, no_of_top_snps=%s, pvalue=%s.\n"%\
										(db_result.id, db_result.results_id, db_result.list_type_id, db_result.no_of_top_snps, db_result.pvalue))
						return None
				
				not_enough_sample = 0
				if candidate_sample_size<pd.min_sample_size:	#don't look at how many non-candidates are (different from whole genome rank test)
					if self.debug:
						sys.stderr.write("Ignore. sample size less than %s. %s vs %s.\n"%(pd.min_sample_size, candidate_sample_size, non_candidate_sample_size))
					return None
				if pd.test_type_id==14 and pd.null_distribution_type_id==1:	#defunct
					w_result = rpy.r.wilcox_test(permData.candidate_gene_snp_rank_ls, permData.non_candidate_gene_snp_rank_ls, alternative='greater')
					statistic=w_result['statistic']['W']
					pvalue=w_result['p.value']
				elif pd.test_type_id==14 and pd.null_distribution_type_id==2:	#defunct
					pd.permutation_type = 1
					rank_sum_stat, pvalue = self.getPermutationRankSumPvalue(permData.chr2rank_ls, permData.candidate_gene_snp_index_ls, permData.non_candidate_gene_snp_index_ls,\
														permData.rank_sum_stat, \
													permData.no_of_snps, permData.chr2no_of_snps, permutation_type=pd.permutation_type, \
													no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks)
					statistic = rank_sum_stat-(candidate_sample_size*(candidate_sample_size-1))/2.
				elif pd.test_type_id==15:
					if pd.null_distribution_type_id==1:
						x = candidate_sample_size
						m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
						n = permData.no_of_total_snps - m
						k = candidate_sample_size + non_candidate_sample_size
						pvalue = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
					elif pd.null_distribution_type_id==2:
						m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
						n = permData.no_of_total_snps - m
						candidate_gene_set = self.dealWithCandidateGeneList(pd.list_type_id, return_set=True)	#internal cache
						top_snp_index_ls = numpy.hstack((permData.candidate_gene_snp_index_ls, permData.non_candidate_gene_snp_index_ls))
						return_data = self.get_enrichment_pvalue_by_gw_looping(candidate_sample_size, top_snp_index_ls, candidate_gene_set, \
										pd.snps_context_wrapper, permData.no_of_total_snps, total_chr_pos_ar=permData.total_chr_pos_ar, \
										no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks)
						pvalue = return_data.pvalue
						no_of_tests = return_data.no_of_tests
						no_of_tests_passed = return_data.no_of_tests_passed
					elif pd.null_distribution_type_id==3:
						x = candidate_sample_size
						m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
						n = permData.no_of_total_snps - m
						k = candidate_sample_size + non_candidate_sample_size
						sample_pvalue = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
						candidate_gene_set = self.dealWithCandidateGeneList(pd.list_type_id, return_set=True)	#internal cache
						top_snp_chr_pos_ls = permData.candidate_gene_snp_chr_pos_ls + permData.non_candidate_gene_snp_chr_pos_ls
						return_data = self.get_enrichment_pvalue_by_random_gene_list(sample_pvalue, pd.total_gene_id_ls, candidate_gene_set, \
												permData.total_chr_pos_ar, pd.snps_context_wrapper, top_snp_chr_pos_ls, n,k,\
												no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks)
						pvalue = return_data.pvalue
						no_of_tests = return_data.no_of_tests
						no_of_tests_passed = return_data.no_of_tests_passed
					else:
						sys.stderr.write("null_distribution_type %s not supported.\n"%(pd.null_distribution_type_id))
						return None
				else:
					sys.stderr.write("Test_type %s not supported.\n"%(pd.test_type))
					return None
		except:
			sys.stderr.write("Exception happened for results_method_id=%s, list_type_id=%s.\n"%(pd.results_id, pd.list_type_id))
			traceback.print_exc()
			sys.stderr.write('%s.\n'%repr(sys.exc_info()))
			return None
		result = TestResultClass(list_type_id=pd.list_type_id, pvalue=pvalue)
		result.results_id = pd.results_id
		result.candidate_sample_size = candidate_sample_size
		result.non_candidate_sample_size = non_candidate_sample_size
		#result.no_of_top_candidate_genes = x
		#result.no_of_top_genes = k
		result.min_distance = pd.min_distance
		result.no_of_tests_passed = no_of_tests_passed
		result.no_of_tests = no_of_tests
		result.candidate_gw_size = m
		result.non_candidate_gw_size = n
		result.no_of_top_snps = candidate_sample_size + non_candidate_sample_size
		result.starting_rank = starting_rank
		result.type_id = pd.type_id
		if score_range:
			result.max_score = score_range[0]
			if pd.min_score is not None:
				result.min_score = pd.min_score
			else:
				result.min_score = score_range[1]
		if self.debug:
			sys.stderr.write("Done.\n")
		return result
	
	def getTopSNPTestType(self, get_closest, min_MAF, allow_two_sample_overlapping, results_type,\
				test_type_id, null_distribution_type_id):
		"""
		2008-10-30
			null_distribution_type_id in CandidateGeneTopSNPTestRMType doesn't matter anymore.
			set it to 1
		2008-10-26
			min_distance is removed from CandidateGeneTopSNPTestRMType.
		2008-10-16
			check whcih TopSNPTest type this is, create one if it doesn't exist in db
		"""
		if self.debug:
			sys.stderr.write("Getting  CandidateGeneTopSNPTestRMType ...")
		rows = Stock_250kDB.CandidateGeneTopSNPTestRMType.query.\
				filter_by(get_closest =get_closest).\
				filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF>=min_MAF-0.0001).filter(Stock_250kDB.CandidateGeneTopSNPTestRMType.min_MAF<=min_MAF+0.0001).\
				filter_by(allow_two_sample_overlapping = allow_two_sample_overlapping).filter_by(results_type=results_type).\
				filter_by(test_type_id=test_type_id).\
				filter_by(null_distribution_type_id=1)
		if rows.count()>0:
			_type = rows.first()
		else:
			_type = Stock_250kDB.CandidateGeneTopSNPTestRMType(get_closest =get_closest,\
										min_MAF = min_MAF,\
										allow_two_sample_overlapping = allow_two_sample_overlapping, \
										results_type=results_type,\
										test_type_id=test_type_id,\
										null_distribution_type_id=1)
		if self.debug:
			sys.stderr.write("Done.\n")
		return _type
	
	
	def run(self):
		"""
		2008-08-19
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		total_gene_id_ls = get_total_gene_ls(db.metadata.bind)
		no_of_total_genes = len(total_gene_id_ls)
		#no_of_total_genes = self.getNoOfTotalGenes(db, self.gene_table, self.tax_id)
		
		#if self.commit:
		#	session.begin()
		_type = self.getTopSNPTestType(self.get_closest, self.min_MAF, \
									self.allow_two_sample_overlapping, self.results_type,\
									self.test_type_id, self.null_distribution_type_id)
		
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		pd = PassingData(list_type_id=self.list_type_id,
						snps_context_wrapper=snps_context_wrapper, \
						no_of_total_genes=no_of_total_genes,
						results_directory=self.results_directory, \
						min_MAF=self.min_MAF,
						get_closest=self.get_closest,
						min_distance=self.min_distance,\
						no_of_top_snps=self.no_of_top_snps,
						min_sample_size=self.min_sample_size,
						test_type_id=self.test_type_id, \
						results_type=self.results_type,
						no_of_permutations=self.no_of_permutations,\
						no_of_min_breaks=self.no_of_min_breaks,
						type_id=_type.id,\
						null_distribution_type_id=self.null_distribution_type_id,\
						allow_two_sample_overlapping=self.allow_two_sample_overlapping,
						total_gene_id_ls=total_gene_id_ls,\
						min_score=self.min_score,
						commit=self.commit)
		if getattr(self, 'output_fname', None):
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			header_row = []
			for column in CandidateGeneTopSNPTest.c.keys():
				header_row.append(column)
			writer.writerow(header_row)
		else:
			writer = None
		
		#2008-10-31 setting up list accordingly
		if self.min_score:
			pd.min_score_ls = [self.min_score]
		else:
			pd.no_of_top_snps_ls = [self.no_of_top_snps]
		for results_id in self.results_id_ls:
			pd.results_id = results_id
			self.runEnrichmentTestToGetNullData(session, pd)
			"""
			result = self.runHGTest(pd)
			if result is not None:
				result.type = _type	#assign the type here
				row = []
				for column in result.c.keys():
					row.append(getattr(result, column))
					print '%s: %s'%(column, row[-1])
				if writer:
					writer.writerow(row)
				session.save(result)
				if self.commit:
					session.flush()
			"""
			
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = TopSNPTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
