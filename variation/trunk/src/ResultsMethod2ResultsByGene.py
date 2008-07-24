#!/usr/bin/env python
"""
Examples:
	ResultsMethod2ResultsByGene.py -e 758 -u yh -c
	
Description:
	program to pull one results_method from db and convert its SNP-based score from a file into gene-based score.
	The results go into table results_by_gene.
	
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
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, ResultsByGene
from Results2DB_250k import Results2DB_250k
from pymodule import getGenomeWideResultFromFile
from GeneListRankTest import GeneListRankTest

class ResultsMethod2ResultsByGene(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict
	option_default_dict.pop(("list_type_id", 1, int))
	
	def __init__(self,  **keywords):
		GeneListRankTest.__init__(self, **keywords)
	
	def saveResultsByGene(self, session, rm, snps_context_wrapper, results_directory=None, commit=False):
		"""
		2008-07-19
		"""
		sys.stderr.write("Saving ResultsByGene ... \n")
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rm.filename))
		else:
			result_fname = rm.filename
		genome_wide_result = getGenomeWideResultFromFile(result_fname)
		
		gene_id2hit = {}
		counter = 0
		score_ls = []
		id_ls = []
		for data_obj in genome_wide_result.data_obj_ls:
			pvalue = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				rbg = ResultsByGene(gene_id=gene_id, snps_id=snps_id, disp_pos=disp_pos, score = pvalue)
				rbg.results_method = rm
				session.save(rbg)
				if commit:
					session.flush()
				score_ls.append(rbg.score)
				id_ls.append(rbg.id)
				
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		
		sys.stderr.write("%s%s\n"%('\x08'*40, counter))
		"""
		#put rank into the db
		import rpy
		rank_ls = rpy.r.rank(score_ls)
		for i in range(len(id_ls)):
			rbg_id = id_ls[i]
			rank = rank_ls[i]
			rbg = ResultsByGene.get(rbg_id)
			rbg.rank = rank
			session.save_or_update(rbg)
			session.flush()
		"""
		session.clear()
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2008-07-17
		"""
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		session = db.session
		snps_context_wrapper = self.constructDataStruc(self.min_distance)
		for results_method_id in self.results_method_id_ls:
			rm = ResultsMethod.get(results_method_id)
			if not rm:
				sys.stderr.write("No results method available for results_method_id=%s.\n"%results_method_id)
				continue
			self.saveResultsByGene(session, rm, snps_context_wrapper, commit=self.commit)
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultsMethod2ResultsByGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()