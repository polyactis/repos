#!/usr/bin/env python
"""

Examples:
	CheckCandidateGeneRank.py -e 1-5 -l 1  -o /tmp/hist_of_results_by_gene_candidate_score_rank
	
Description:
	2008-09-28
	program to check histogram of ranks of candidate genes vs those of non-candidate genes in results_by_gene.
	it also draws histogram of scores of results_by_gene.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
#	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
#	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr
import Stock_250kDB
from Stock_250kDB import ResultsByGene, ResultsMethod
from sets import Set
from GeneListRankTest import GeneListRankTest
#from sqlalchemy.orm import join
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6
import pylab


class CheckCandidateGeneRank(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("phenotype_id_ls", 1, ): [None, 'e', 1, 'comma/dash-separated phenotype_method id list, like 1,3-7'],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency. deprecated.'],\
							('min_sample_size', 0, int): [5, 'i', 1, 'minimum size for both candidate and non-candidate sets to do wilcox.test'],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("output_dir", 1, ): [None, 'o', 1, 'directory to store output'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	
	def __init__(self,  **keywords):
		"""
		2008-07-24
			split results_id_ls if it exists, to accomodate MpiGeneListRankTest which removed this option
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.phenotype_id_ls = getListOutOfStr(self.phenotype_id_ls, data_type=int)
	
	def getResultsIDLs(self, db, phenotype_id_ls, min_distance, get_closest, min_MAF):
		"""
		2008-09-28
		"""
		sys.stderr.write("Getting results_id_ls ...")
		phenotype_id2results_id_ls = {}
		rows = ResultsByGene.query.filter_by(min_distance=min_distance).filter_by(get_closest=get_closest).\
			filter(ResultsByGene.results_method.has(ResultsMethod.phenotype_method_id.in_(phenotype_id_ls))).\
			filter(ResultsByGene.min_MAF>=min_MAF-0.0001).filter(ResultsByGene.min_MAF<=min_MAF+0.0001)
			#.order_by(ResultsMethod.analysis_method_id)
		counter = 0
		for row in rows:
			phenotype_id = row.results_method.phenotype_method_id
			if phenotype_id not in phenotype_id2results_id_ls:
				phenotype_id2results_id_ls[phenotype_id] = []
			phenotype_id2results_id_ls[phenotype_id].append((row.results_method.analysis_method_id, row.id))
			counter += 1
		
		#sort results_id_ls for each phenotype according to analysis_method_id
		for phenotype_id, results_id_ls in phenotype_id2results_id_ls.iteritems():
			results_id_ls.sort()	#sorted by analysis_method_id
			phenotype_id2results_id_ls[phenotype_id] = [row[1] for row in results_id_ls]
		sys.stderr.write("%s results. Done.\n"%(counter))
		return phenotype_id2results_id_ls
	
	def getScoreRank(self, rbg, candidate_gene_set, results_directory):
		"""
		2008-09-28
		"""
		sys.stderr.write("Getting score & rank list ...")
		if results_directory:	#given a directory where all results are.
			result_fname = os.path.join(results_directory, os.path.basename(rbg.filename))
		else:
			result_fname = rbg.filename
		if not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		candidate_score_ls = []
		non_candidate_score_ls = []
		candidate_rank_ls = []
		non_candidate_rank_ls = []
		for row in reader:
			gene_id = int(row[col_name2index['gene_id']])
			score = float(row[col_name2index['score']])
			if gene_id in candidate_gene_set:
				candidate_score_ls.append(score)
				candidate_rank_ls.append(counter)
			else:
				non_candidate_score_ls.append(score)
				non_candidate_rank_ls.append(counter)
			counter += 1
		del reader
		analysis_method = Stock_250kDB.AnalysisMethod.get(rbg.results_method.analysis_method_id)
		
		score_rank_data = PassingData(candidate_score_ls=candidate_score_ls, candidate_rank_ls=candidate_rank_ls,\
								non_candidate_score_ls=non_candidate_score_ls, non_candidate_rank_ls=non_candidate_rank_ls,\
								analysis_method=analysis_method)
		
		sys.stderr.write("Done.\n")
		return score_rank_data
	
	def plotSubHistogram(self, candidate_data_ls, non_candidate_data_ls, which_figure, sub_title, xlabel):
		"""
		2008-09-28
		"""
		pylab.subplot(2,2,which_figure)
		pylab.title(sub_title)
		pylab.xlabel(xlabel)
		hist_patch_ls = []
		legend_ls = []
		no_of_bins = min(100, int(len(non_candidate_data_ls)/5))
		h1 = pylab.hist(non_candidate_data_ls, no_of_bins, alpha=0.2, normed=1)
		hist_patch_ls.append(h1[2][0])
		legend_ls.append('non-candidate gene')
		no_of_bins = min(50, int(len(candidate_data_ls)/5))
		h2 = pylab.hist(candidate_data_ls, no_of_bins, alpha=0.2, normed=1, facecolor='r')
		hist_patch_ls.append(h2[2][0])
		legend_ls.append('candidate gene')
		pylab.legend(hist_patch_ls, legend_ls, )
	
	def plotHistForOnePhenotype(self, phenotype_method, list_type, score_rank_data_ls, output_dir, data_type='score'):
		"""
		2008-09-28
		"""
		sys.stderr.write("Drawing histogram ...")
		pylab.clf()
		title = '%s_%s_%s_%s'%(phenotype_method.id, phenotype_method.short_name, list_type.id, list_type.short_name)
		pylab.title(title)
		
		output_fname_prefix = os.path.join(output_dir, title.replace('/', '_'))
		output_fname_prefix = '%s_%s'%(output_fname_prefix, data_type)
		for i in range(len(score_rank_data_ls)):
			score_rank_data = score_rank_data_ls[i]
			if data_type=='score':
				candidate_data_ls = score_rank_data.candidate_score_ls
				non_candidate_data_ls = score_rank_data.non_candidate_score_ls
			else:
				candidate_data_ls = score_rank_data.candidate_rank_ls
				non_candidate_data_ls = score_rank_data.non_candidate_rank_ls
			sub_title = score_rank_data.analysis_method.short_name
			xlabel = data_type
			self.plotSubHistogram(candidate_data_ls, non_candidate_data_ls, i+1, sub_title, xlabel)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		phenotype_id2results_id_ls = self.getResultsIDLs(db, self.phenotype_id_ls, self.min_distance, self.get_closest, self.min_MAF)
		candidate_gene_list = self.getGeneList(self.list_type_id)
		candidate_gene_set = Set(candidate_gene_list)
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		list_type = Stock_250kDB.GeneListType.get(self.list_type_id)
		for phenotype_id, results_id_ls in phenotype_id2results_id_ls.iteritems():
			phenotype_method = Stock_250kDB.PhenotypeMethod.get(phenotype_id)

			score_rank_data_ls = []
			for results_id in results_id_ls:
				rbg = Stock_250kDB.ResultsByGene.get(results_id)
				score_rank_data = self.getScoreRank(rbg, candidate_gene_set, self.results_directory)
				score_rank_data_ls.append(score_rank_data)
			self.plotHistForOnePhenotype(phenotype_method, list_type, score_rank_data_ls, self.output_dir, data_type='score')
			self.plotHistForOnePhenotype(phenotype_method, list_type, score_rank_data_ls, self.output_dir, data_type='rank')
			
			
			
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CheckCandidateGeneRank
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()