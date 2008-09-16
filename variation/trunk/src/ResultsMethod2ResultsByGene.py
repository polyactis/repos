#!/usr/bin/env python
"""
Examples:
	ResultsMethod2ResultsByGene.py -e 758 -u yh -c
	
	#input is all results with call_method_id=17
	ResultsMethod2ResultsByGene.py -l 17 -u yh -c
	
	#run the program on the cluster, one results_method=2077, change the db input/output directory, given snps_context pickle file
	ResultsMethod2ResultsByGene.py -e 2077 -o ~/panfs/db/results_by_gene/ -c -u yh -s ~/panfs/250k/snps_context_g0_m20000 -m 20000 -t ~/panfs/db/results/type_1/
	
	#ditto, but all results_method entries from call_method_id=17 and analysis_method_id=1
	ResultsMethod2ResultsByGene.py -o ~/panfs/db/results_by_gene/ -c -u yh -s ~/panfs/250k/snps_context_g0_m20000 -m 20000 -t ~/panfs/db/results/type_1/ -p yh324 -l 17 -a 1
	
Description:
	program to pull one results_method from db and convert its SNP-based score from a file into gene-based score.
	The record goes into table results_by_gene while the actual results go to a file.
	
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
from heapq import heappush, heappop, heapreplace
from GeneListRankTest import SnpsContextWrapper

class ResultsMethod2ResultsByGene(TopSNPTest):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_distance", 1, int): [50000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency'],\
							('call_method_id', 0, int):[0, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[0, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							("results_method_id_ls", 0, ): [None, 'e', 1, 'comma-separated results_method_id list'],\
							('input_db_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							('output_db_directory', 0, ):[None, 'o', 1, 'The file system directory corresponding to table results_by_gene. Supply this to overwrite the default. Database records still use the default_output_db_directory.'],\
							('default_output_db_directory', 0, ):['/Network/Data/250k/db/results_by_gene/', 'f', 1, 'The file system directory corresponding to table results_by_gene. It goes into database record. Usually no need to change this one..'],\
							("snps_context_picklef", 0, ): [None, 's', 1, 'given the option, if the file does not exist yet, to store a pickled snps_context_wrapper into it, min_distance and flag get_closest will be attached to the filename. If the file exists, load snps_context_wrapper out of it.'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-27
			inherit from TopSNPTest
		"""
		TopSNPTest.__init__(self, **keywords)
	
	def sortGeneIDBasedOnScore(self, gene_id2hit):
		"""
		2008-09-16
			use a max heap to sort gene_id so that most significant gene would come out first when heappop
		"""
		if self.report or self.debug:
			sys.stderr.write("Sorting gene id based on score ... \n")
		heap_ls = []
		for gene_id, hit in gene_id2hit.iteritems():
			heappush(heap_ls, (-hit.score, gene_id))	#score is already properly transformed. bigger it is , the more significant. heap is max heap. heappop would pop the smallest item.
		if self.report or self.debug:
			sys.stderr.write("Done.\n")
		return heap_ls
		
		
	def saveResultsByGene(self, session, rm, param_data):
		"""
		2008-09-16
			add some checkup (if database already has this entry, if the result file exist...) before starting
		2008-09-16
			overhaul, table results_by_gene only stores the link to a file, which stores max score for all genes, ordered by their score.
		2008-08-27
			derive the rank from the SNP in genome_wide_result
			ResultsByGene has one more column, readme_id.
		2008-07-19
		"""
		sys.stderr.write("Saving ResultsByGene ... \n")
		#check if it's in db already or not
		rbg = Stock_250kDB.ResultsByGene.query.filter_by(results_method_id=rm.id).filter_by(min_distance=param_data.min_distance).\
			filter_by(get_closest=param_data.get_closest)
		if rbg.count()>0:
			rbg = rbg.first()
			sys.stderr.write("Skip. An entry exists in results_by_gene (id=%s) with same (results_method_id, min_distance, get_closest)=(%s, %s, %s).\n"\
							%(rbg.id, rbg.results_method_id, rbg.min_distance, rbg.get_closest))
			return
		
		gene_id2hit = self.getGeneID2MostSignificantHit(rm, param_data.snps_context_wrapper, param_data.results_directory, param_data.min_MAF)
		if gene_id2hit is None:
			sys.stderr.write("Skip. gene_id2hit is None.\n")
			return
		gene_id_heap_ls = self.sortGeneIDBasedOnScore(gene_id2hit)
		rbg = Stock_250kDB.ResultsByGene(short_name='%s_m%s_g%s_by_gene'%(rm.short_name, param_data.min_distance, param_data.get_closest),\
										min_distance = param_data.min_distance, get_closest=param_data.get_closest, min_MAF=param_data.min_MAF)
		rbg.results_method = rm
		session.save(rbg)
		session.flush()
		rbg.filename = os.path.join(param_data.default_output_db_directory, '%s.tsv'%rbg.id)
		session.save_or_update(rbg)
		session.flush()
		
		if param_data.output_db_directory:
			output_fname = os.path.join(param_data.output_db_directory, '%s.tsv'%rbg.id)
		else:
			output_fname = rbg.filename
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['gene_id', 'score', 'snps_id', 'disp_pos', 'comment']
		writer.writerow(header)
		while len(gene_id_heap_ls)>0:
			gene_id = heappop(gene_id_heap_ls)[1]
			hit = gene_id2hit[gene_id]
			row = [gene_id, hit.score, hit.snps_id, hit.disp_pos, None]
			writer.writerow(row)
		del writer
		sys.stderr.write("Done.\n")
	
	def getResultsMethodIDLs(self, pdata):
		"""
		2008-08-27
			get all results method id given call_method_id
		"""
		sys.stderr.write("Getting all results method ids ...")
		i = 0
		block_size = 5000
		query = Stock_250kDB.ResultsMethod.query.filter_by(results_method_type_id=1)
		if pdata.call_method_id!=0:
			query = query.filter_by(call_method_id=pdata.call_method_id)
		if pdata.analysis_method_id!=0:
			query = query.filter_by(analysis_method_id=pdata.analysis_method_id)
		
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
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		
		param_data = PassingData()
		param_data.results_directory = self.input_db_directory
		param_data.default_output_db_directory = self.default_output_db_directory
		param_data.output_db_directory = self.output_db_directory
		param_data.commit = self.commit
		param_data.min_MAF = self.min_MAF
		param_data.min_distance = self.min_distance
		param_data.get_closest = self.get_closest
		param_data.snps_context_wrapper = snps_context_wrapper
		
		if not self.results_method_id_ls:
			pdata = PassingData(call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id)
			self.results_method_id_ls = self.getResultsMethodIDLs(pdata)
		
		for results_method_id in self.results_method_id_ls:
			rm = Stock_250kDB.ResultsMethod.get(results_method_id)
			if not rm:
				sys.stderr.write("No results method available for results_method_id=%s.\n"%results_method_id)
				continue
			self.saveResultsByGene(session, rm, param_data)
		
		if self.commit:
			session.commit()
			session.clear()
		else:
			session.rollback()
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ResultsMethod2ResultsByGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()