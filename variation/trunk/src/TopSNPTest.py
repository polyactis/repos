#!/usr/bin/env python
"""
Examples:
	#debug, quick testing
	TopSNPTest.py -e 389 -l 1 -u yh  -o /tmp/top_snp_test.out -b
	
	TopSNPTest.py -e 389 -l 1 -u yh  -o /tmp/top_snp_test.out
	
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


class TopSNPTest(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict.update({('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.']})
	option_default_dict.update({('gene_table', 1, ): ['genome.gene', 'a', 1, 'to get the number of total genes from database, which table.']})
	option_default_dict.update({('no_of_top_snps', 1, int): [50, 'f', 1, 'how many number of top snps based on score or -log(pvalue).']})
	option_default_dict.update({("test_type", 1, int): [1, 'y', 1, 'which type of tests. 1: r.wilcox.test() 2: loop-permutation. 3: hypergeometric. 4,5,6 are their counterparts which allow_two_sample_overlapping.']})
	
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		GeneListRankTest.__init__(self, **keywords)
	
	def getNoOfTotalGenes(self, db, gene_table='genome.gene', tax_id=3702):
		if self.debug:
			sys.stderr.write("Getting no of total genes ... ")
		rows = db.metadata.bind.execute("select count(gene_id) as count from %s where tax_id=%s"%(gene_table, tax_id))
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
			if pd.test_type>3:
				allow_two_sample_overlapping = 1
			else:
				allow_two_sample_overlapping = 0
			param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, \
										min_MAF=pd.min_MAF, allow_two_sample_overlapping=allow_two_sample_overlapping)	#starting_rank=1, go through all SNPs
			permData = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, param_data)
			candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
			non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
			self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_associated=candidate_sample_size
			self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_unassociated=non_candidate_sample_size			
		return self.list_type_id2candidate_gene_list_info[list_type_id].no_of_snps_associated
	
	def runHGTest(self, pd):
		"""
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
			sys.stderr.write("Running hypergeometric test ... ")
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
			sys.stderr.write("No results available for results_id=%s.\n"%pd.results_method_id)
			return None
		starting_rank = getattr(pd, 'starting_rank', 1)	#from which rank to look down for top snps
		db_results = TestResultClass.query.filter_by(results_id=pd.results_id).filter_by(list_type_id=pd.list_type_id).\
			filter_by(no_of_top_snps=pd.no_of_top_snps).filter_by(min_distance=pd.min_distance).\
			filter_by(min_MAF=pd.min_MAF).filter_by(get_closest=pd.get_closest).filter_by(test_type=pd.test_type).\
			filter_by(starting_rank=starting_rank)
		if db_results.count()>0:	#done before
			db_result = db_results.first()
			sys.stderr.write("It's done already. id=%s, results_id=%s, list_type_id=%s, no_of_top_snps=%s, pvalue=%s.\n"%\
							(db_result.id, db_result.results_id, db_result.list_type_id, db_result.no_of_top_snps, db_result.pvalue))
			return None
		
		try:
			import rpy
			candidate_gene_list = self.dealWithCandidateGeneList(pd.list_type_id)	#internal cache
			"""
			if rm.analysis_method_id==13:
				passingdata = self.prepareDataForHGTest_SNPPair(rm, pd.snps_context_wrapper, candidate_gene_list, pd.results_directory, pd.min_MAF, pd.no_of_top_snps)
			else:
				passingdata = self.prepareDataForHGTest(rm, pd.snps_context_wrapper, candidate_gene_list, pd.results_directory, pd.min_MAF, pd.no_of_top_snps)
			"""
			if pd.test_type>3:
				allow_two_sample_overlapping = 1
			else:
				allow_two_sample_overlapping = 0
			param_data = PassingData(results_directory=pd.results_directory, candidate_gene_list=candidate_gene_list, \
									min_MAF=pd.min_MAF, allow_two_sample_overlapping=allow_two_sample_overlapping,\
									no_of_top_lines=pd.no_of_top_snps, no_of_top_snps=pd.no_of_top_snps, \
									starting_rank=starting_rank)
			
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
				candidate_sample_size = len(permData.candidate_gene_snp_rank_ls)
				non_candidate_sample_size = len(permData.non_candidate_gene_snp_rank_ls)
				not_enough_sample = 0
				if candidate_sample_size<pd.min_sample_size:	#don't look at how many non-candidates are (different from whole genome rank test)
					sys.stderr.write("Ignore. sample size less than %s. %s vs %s.\n"%(pd.min_sample_size, candidate_sample_size, non_candidate_sample_size))
					return None
				if pd.test_type%3==1:
					w_result = rpy.r.wilcox_test(permData.candidate_gene_snp_rank_ls, permData.non_candidate_gene_snp_rank_ls, alternative='greater')
					statistic=w_result['statistic']['W']
					pvalue=w_result['p.value']
				elif pd.test_type%3==2:
					pd.permutation_type = 1
					rank_sum_stat, pvalue = self.getPermutationRankSumPvalue(permData.chr2rank_ls, permData.candidate_gene_snp_index_ls, permData.non_candidate_gene_snp_index_ls,\
														permData.rank_sum_stat, \
													permData.no_of_snps, permData.chr2no_of_snps, permutation_type=pd.permutation_type, \
													no_of_permutations=pd.no_of_permutations, no_of_min_breaks=pd.no_of_min_breaks)
					statistic = rank_sum_stat-(candidate_sample_size*(candidate_sample_size-1))/2.
				elif pd.test_type%3==0:
					x = candidate_sample_size
					m = self.dealWithNoOfSNPsAssociatedWithCandidateGeneList(pd.list_type_id, rm, pd)	#cache is internally going on
					n = permData.no_of_total_snps - m
					k = pd.no_of_top_snps
					pvalue = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
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
		result.min_distance = min_distance
		result.min_MAF = min_MAF
		result.get_closest = get_closest
		result.candidate_sample_size = candidate_sample_size
		result.non_candidate_sample_size = non_candidate_sample_size
		#result.no_of_top_candidate_genes = x
		#result.no_of_top_genes = k
		result.no_of_top_snps = pd.no_of_top_snps
		result.test_type = pd.test_type
		result.starting_rank = starting_rank
		if self.debug:
			sys.stderr.write("Done.\n")
		return result
	
	def run(self):
		"""
		2008-08-19
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		no_of_total_genes = self.getNoOfTotalGenes(db, self.gene_table, self.tax_id)
		#if self.commit:
		#	session.begin()
		
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		pd = PassingData(list_type_id=self.list_type_id, snps_context_wrapper=snps_context_wrapper, \
							no_of_total_genes=no_of_total_genes, results_directory=self.results_directory, \
							min_MAF=self.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance,\
							no_of_top_snps=self.no_of_top_snps, min_sample_size=self.min_sample_size, test_type=self.test_type, \
							results_type=self.results_type, no_of_permutations=self.no_of_permutations,\
							no_of_min_breaks=self.no_of_min_breaks)
		if getattr(self, 'output_fname', None):
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			header_row = []
			for column in CandidateGeneTopSNPTest.c.keys():
				header_row.append(column)
			writer.writerow(header_row)
		else:
			writer = None
			
		for results_id in self.results_id_ls:
			pd.results_id = results_id
			result = self.runHGTest(pd)
			if result is not None:
				row = []
				for column in result.c.keys():
					row.append(getattr(result, column))
					print '%s: %s'%(column, row[-1])
				if writer:
					writer.writerow(row)
				#session.save(result)
				if self.commit:
					session.flush()
			
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = TopSNPTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
