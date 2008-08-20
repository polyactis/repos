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
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, CandidateGeneRankSumTestResult, CandidateGeneTopSNPTest
from Results2DB_250k import Results2DB_250k
from pymodule import getGenomeWideResultFromFile
from GeneListRankTest import GeneListRankTest
from sets import Set



class TopSNPTest(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict.update({('tax_id', 1, int): [3702, 'x', 1, 'to get the number of total genes from database, which species.']})
	option_default_dict.update({('gene_table', 1, ): ['genome.gene', 'a', 1, 'to get the number of total genes from database, which table.']})
	option_default_dict.update({('no_of_top_snps', 1, int): [50, 'f', 1, 'how many number of top snps based on score or -log(pvalue).']})
	option_default_dict.pop(('max_pvalue_per_gene', 0, int))
	
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
		genome_wide_result.data_obj_ls.sort(reverse=True)	#in value descending order. each SNP object has a defined method for comparison based on its value
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
		
	def runHGTest(self, pd):
		"""
		2008-08-20
		"""
		if self.debug:
			sys.stderr.write("Running hypergeometric test ... ")
		rm = ResultsMethod.get(pd.results_method_id)
		if not rm:
			sys.stderr.write("No results method available for results_method_id=%s.\n"%pd.results_method_id)
			return None
		if rm.results_method_type_id!=1:
			sys.stderr.write("skip non-association results. results_method_type_id=%s, results_method_id=%s.\n"%(rm.results_method_type_id, pd.results_method_id))
			return None
		db_results = CandidateGeneTopSNPTest.query.filter_by(results_method_id=pd.results_method_id).filter_by(list_type_id=pd.list_type_id).filter_by(min_distance=pd.min_distance).filter_by(min_MAF=pd.min_MAF).filter_by(get_closest=pd.get_closest)
		if db_results.count()>0:	#done before
			db_result = db_results.first()
			sys.stderr.write("It's done already. id=%s, results_method_id=%s, list_type_id=%s, pvalue=%s.\n"%\
							(db_result.id, db_result.results_method_id, db_result.list_type_id, db_result.pvalue))
			return None
		try:
			candidate_gene_list = self.getGeneList(pd.list_type_id)
			passingdata = self.prepareDataForHGTest(rm, pd.snps_context_wrapper, candidate_gene_list, pd.results_directory, pd.min_MAF, pd.no_of_top_snps)
			import rpy
			x = len(passingdata.candidate_gene_in_top_set)
			m = len(candidate_gene_list)
			n = pd.no_of_total_genes - m
			k = x + len(passingdata.non_candidate_gene_in_top_set)
			p_value = rpy.r.phyper(x-1,m,n,k,lower_tail = rpy.r.FALSE)
		except:
			sys.stderr.write("Exception happened for results_method_id=%s, list_type_id=%s.\n"%(pd.results_method_id, pd.list_type_id))
			traceback.print_exc()
			sys.stderr.write('%s.\n'%repr(sys.exc_info()))
			return None
		result = CandidateGeneTopSNPTest(list_type_id=pd.list_type_id, pvalue=p_value)
		result.results_method_id = pd.results_method_id
		result.min_distance = pd.min_distance
		result.min_MAF = pd.min_MAF
		result.get_closest = pd.get_closest
		result.no_of_top_candidate_genes = x
		result.no_of_top_genes = k
		result.no_of_top_snps = pd.no_of_top_snps
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
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		session = db.session
		no_of_total_genes = self.getNoOfTotalGenes(db, self.gene_table, self.tax_id)
		#if self.commit:
		#	session.begin()
		
		snps_context_wrapper = self.constructDataStruc(self.min_distance, self.get_closest)
		pd = PassingData(snps_context_wrapper=snps_context_wrapper, list_type_id=self.list_type_id, \
							no_of_total_genes=no_of_total_genes, results_directory=self.results_directory, \
							min_MAF=self.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance,\
							no_of_top_snps=self.no_of_top_snps)
		if getattr(self, 'output_fname', None):
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			header_row = []
			for column in CandidateGeneTopSNPTest.c.keys():
				header_row.append(column)
			writer.writerow(header_row)
		else:
			writer = None
			
		for results_method_id in self.results_method_id_ls:
			pd.results_method_id = results_method_id
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