#!/usr/bin/env python
"""
Examples:
	ResultsMethod2ResultsByGene.py -e 758 -u yh -c
	
	#input is all results with call_method_id=17
	ResultsMethod2ResultsByGene.py -l 17 -u yh -c
	
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
from pymodule.db import formReadmeObj
import Stock_250kDB
from TopSNPTest import TopSNPTest

class ResultsMethod2ResultsByGene(TopSNPTest):
	__doc__ = __doc__
	option_default_dict = TopSNPTest.option_default_dict.copy()
	option_default_dict.pop(("results_method_id_ls", 1, ))
	option_default_dict.update({("results_method_id_ls", 0, ): [None, 'e', 1, 'comma-separated results_method_id list']})
	option_default_dict.update({('call_method_id', 0, int):[0, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.']})
	option_default_dict.pop(("list_type_id", 1, int))
	option_default_dict.pop(('tax_id', 1, int))
	option_default_dict.pop(('min_sample_size', 0, int))
	
	def __init__(self,  **keywords):
		"""
		2008-08-27
			inherit from TopSNPTest
		"""
		TopSNPTest.__init__(self, **keywords)
	
	def saveResultsByGene(self, session, rm, snps_context_wrapper, param_data):
		"""
		2008-08-27
			derive the rank from the SNP in genome_wide_result
			ResultsByGene has one more column, readme_id.
		2008-07-19
		"""
		sys.stderr.write("Saving ResultsByGene ... \n")
		genome_wide_result = self.getResultMethodContent(rm, param_data.results_directory, param_data.min_MAF)
		genome_wide_result.data_obj_ls.sort()	#in value descending order. each SNP object has a defined method for comparison based on its value
		genome_wide_result.data_obj_ls.reverse()
		counter = 0
		id_ls = []
		score_ls = []
		for i in range(param_data.no_of_top_snps):
			data_obj = genome_wide_result.data_obj_ls[i]
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				rbg = Stock_250kDB.ResultsByGene(gene_id=gene_id, snps_id=snps_id, disp_pos=disp_pos, score=data_obj.value, rank=i+1)
				rbg.results_method = rm
				rbg.readme = param_data.readme
				session.save(rbg)
				if param_data.commit:
					session.flush()
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
	
	def getResultsMethodIDLs(self, call_method_id):
		"""
		2008-08-27
			get all results method id given call_method_id
		"""
		sys.stderr.write("Getting all results method ids ...")
		i = 0
		block_size = 5000
		if call_method_id!=0:
			query = Stock_250kDB.ResultsMethod.query.filter_by(results_method_type_id=1).filter_by(call_method_id=call_method_id)
		else:
			query = Stock_250kDB.ResultsMethod.query.filter_by(results_method_type_id=1)
		rows = query.offset(i).limit(block_size)
		results_method_id_ls = []
		while rows.count()!=0:
			for row in rows:
				results_method_id_ls.append(row.id)
				i += 1
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("%s results.\n"%(len(results_method_id_ls)))
		return results_method_id_ls
	
	def run(self):
		"""
		2008-07-17
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		snps_context_wrapper = self.constructDataStruc(self.min_distance, self.get_closest)
		param_data = PassingData()
		param_data.results_directory = self.results_directory
		param_data.commit = self.commit
		param_data.min_MAF = self.min_MAF
		param_data.no_of_top_snps = self.no_of_top_snps
		
		readme = formReadmeObj(sys.argv, self.ad, Stock_250kDB.README)
		session.save(readme)
		param_data.readme = readme
		
		if not self.results_method_id_ls:
			self.results_method_id_ls = self.getResultsMethodIDLs(self.call_method_id)
		
		for results_method_id in self.results_method_id_ls:
			rm = Stock_250kDB.ResultsMethod.get(results_method_id)
			if not rm:
				sys.stderr.write("No results method available for results_method_id=%s.\n"%results_method_id)
				continue
			self.saveResultsByGene(session, rm, snps_context_wrapper, param_data)
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultsMethod2ResultsByGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()