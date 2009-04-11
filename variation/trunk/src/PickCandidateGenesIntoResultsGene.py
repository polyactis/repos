#!/usr/bin/env python
"""
Examples:
	#phenotype_method_id_ls=1-7,39, analysis_method_id_ls=1,5,6,7, list_type_id_ls=28,64
	PickCandidateGenesIntoResultsGene.py -A 1-7,39 -e 1,5,6,7 -f 5000 -j 17 -l 28,64 -s ./mnt2/panfs/250k/snps_context_g0_m0 -m 0 -c 

	#change the results_directory when running on cluster
	PickCandidateGenesIntoResultsGene.py -t ~/panfs/db/results/type_1/ -A 9-13,32-38,65-74 -e 1,5,6,7 -f 5000 -j 17 -l 24,29,30 -s ~/panfs/250k/snps_context_g1_m50000 -g -m 50000 -c 

Description:
	2008-10-28	Pick top snps from results that are associated with candidate genes (under certain association rule)
		and put them into table ResultsGene.
	
	Ignore options: min_sample_size, no_of_min_breaks, no_of_permutations, test_type_id, results_type, tax_id.
	
	2009-1-8 list_type_id 0 = ensemble of genes.
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
from pymodule import getGenomeWideResultFromFile
from GeneListRankTest import SnpsContextWrapper
from CheckCandidateGeneRank import CheckCandidateGeneRank
from MpiTopSNPTest import MpiTopSNPTest
from sets import Set
from heapq import heappush, heappop, heapreplace
from common import get_total_gene_ls
import rpy, random, numpy
import Stock_250kDB

class PickCandidateGenesIntoResultsGene(MpiTopSNPTest):
	__doc__ = __doc__
	option_default_dict = MpiTopSNPTest.option_default_dict.copy()
	#option_default_dict.update({('call_method_id', 1, int): [17, '', 1, 'which call method']})
	def __init__(self,  **keywords):
		"""
		2008-11-11
			inherit from MpiTopSNPTest
			use (analysis_method_id_ls, phenotype_method_id_ls, call_method_id) to pick results_method
		2008-10-28
		"""
		MpiTopSNPTest.__init__(self, **keywords)
	
	def pick_candidate_genes(self, pd):
		"""
		2009-4-10
			rank is now part of the whole unique constraint in table Stock_250kDB.ResultsGene.
		2008-11-12
			if entry already exists in ResultsGene, make sure it doesn't have the current type associated with
			upon failure of 'session.flush()', expunge the current entry and report error. avoid failure of the whole program.
		2008-10-28
		"""
		sys.stderr.write("Picking candidates genes from results (id=%s) for list_type (id=%s) ..."%(pd.results_id, pd.list_type_id))
		rm = Stock_250kDB.ResultsMethod.get(pd.results_id)
		session = getattr(pd, 'session', None)
		if session is None:
			sys.stderr.write("session is None. no db connection.\n")
			return None
		if pd.list_type_id==0:
			total_gene_id_ls = get_total_gene_ls(Stock_250kDB.ResultsMethod.table.bind)
			candidate_gene_set = Set(total_gene_id_ls)
		else:
			candidate_gene_set = self.dealWithCandidateGeneList(pd.list_type_id, return_set=True)	#internal cache
		pd.construct_data_obj_id2index = False	#default in getResultMethodContent is True
		pd.construct_chr_pos2index = False	#no need for this as well
		pd.need_candidate_association = True
		pd.need_snp_index = False
		pd.candidate_gene_set = candidate_gene_set
		return_data = self.prepareDataForPermutationRankTest(rm, pd.snps_context_wrapper, pd)
		counter = 0
		for snps_id, disp_pos, gene_id, score, rank in return_data.candidate_association_ls:
			rows = Stock_250kDB.ResultsGene.query.filter_by(snps_id=snps_id).\
				filter_by(gene_id=gene_id).\
				filter_by(results_id=pd.results_id).filter_by(rank=rank)
			if rows.count()==1:
				row = rows.first()
				already_in_db = 0
				if pd.type.id:	#2008-11-12
					for hist_type in row.types:
						if hist_type.id==pd.type.id:
							already_in_db = 1
							break
				if not already_in_db:	#2008-11-12
					row.types.append(pd.type)
					session.save_or_update(row)
					counter += 1
			elif rows.count()>1:
				sys.stderr.write("Error: more than 1 db entries with snps_id=%s, gene_id=%s, results_id=%s.\n"%\
								(snps_id, gene_id, results_id))
				continue
			else:
				row = Stock_250kDB.ResultsGene(snps_id=snps_id, gene_id=gene_id, disp_pos=disp_pos,\
											results_id=pd.results_id, score=score, rank=rank)
				row.types.append(pd.type)
				session.save_or_update(row)
				counter += 1
			if pd.commit:
				try:	#2008-11-12 don't wanna db failure to bog down the whole program
					session.flush()
				except:
					session.expunge(row)
					for column in row.c.keys():
						sys.stderr.write("\t%s=%s.\n"%(column, getattr(row, column)))
					traceback.print_exc()
					sys.stderr.write('%s.\n'%repr(sys.exc_info()))
		sys.stderr.write("%s entries saved. Done.\n"%counter)
	
	def run(self):
		"""
		2008-10-28
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		hist_type = CheckCandidateGeneRank.getHistType(self.call_method_id, self.min_distance, self.get_closest, self.min_MAF, \
									self.allow_two_sample_overlapping, self.results_type, self.null_distribution_type_id)
		
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		
		param_obj = PassingData(call_method_id=self.call_method_id, \
								analysis_method_id=getattr(self, 'analysis_method_id', None),\
								analysis_method_id_ls=getattr(self, 'analysis_method_id_ls', None),\
								phenotype_method_id_ls=getattr(self, 'phenotype_method_id_ls', None),\
								list_type_id_ls=self.list_type_id_ls, \
								results_type=self.results_type,\
								no_check_gene_list=True)
		params_ls = self.generate_params(param_obj)
		
		pd = PassingData(snps_context_wrapper=snps_context_wrapper, \
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
						type=hist_type,\
						null_distribution_type_id=self.null_distribution_type_id,\
						allow_two_sample_overlapping=self.allow_two_sample_overlapping,
						min_score=self.min_score,
						session=session,\
						commit=self.commit)
		
		for results_id, list_type_id in params_ls:
			pd.list_type_id = list_type_id
			pd.results_id = results_id
			self.pick_candidate_genes(pd)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PickCandidateGenesIntoResultsGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()