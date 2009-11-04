#!/usr/bin/env python
"""

Examples:
	CalculateOverlappingStat.py -m 32 -e 1-4 -u yh -c
	
Description:
	2009-11-3 calculate the number of overlapping SNPs among different subsets of analysis methods for one phenotype.
	It fetches SNP-based associations from table results, which if is empty or not sufficient data, it will invoke ResultsMethod2Results.rm2result() to fill in.
	Program would only go ahead if the relevant stat is NOT in table AssociationOverlappingStat.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback, subprocess
from pymodule import figureOutDelimiter, getListOutOfStr
from pymodule.algorithm import listSubsets
import Stock_250kDB
from ResultsMethod2Results import ResultsMethod2Results
from DrawSNPRegion import DrawSNPRegion

class CalculateOverlappingStat(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('phenotype_method_id_ls', 0, ): ['', 'e', 1, 'comma or dash connected phenotype method ids', ],\
							('call_method_id', 1, int): [32, 'm', 1, 'which call method', ],\
							('no_of_top_snps', 1, int): [1000, 'n', 1, 'Number of top SNPs from each association result.', ],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None, then use the one given by db.'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.'],\
							('run_type', 0, int):[1, 'y', 1, 'Run type (useless)']}
	
	def __init__(self, **keywords):
		"""
		2009-10-28
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		if self.phenotype_method_id_ls is not None:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
			
	def getPhenotype2AnalysisMethods(self, call_method_id, phenotype_method_id_ls=None):
		"""
		2009-11-2
		"""
		sys.stderr.write("Getting Phenotype2AnalysisMethods ...")
		phenotype_method_id2analysis_method_id_ls = {}
		query = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id)
		if phenotype_method_id_ls:
			query = query.filter(Stock_250kDB.ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls))
		for row in query:
			phenotype_id = row.phenotype_method_id
			if phenotype_id not in phenotype_method_id2analysis_method_id_ls:
				phenotype_method_id2analysis_method_id_ls[phenotype_id] = []
			phenotype_method_id2analysis_method_id_ls[phenotype_id].append(row.analysis_method_id)
		no_of_associations_ls = map(len, phenotype_method_id2analysis_method_id_ls.values())
		no_of_associations_per_phenotype = sum(no_of_associations_ls)/float(len(no_of_associations_ls))
		sys.stderr.write("%s phenotypes, %s associations per phenotype. Done.\n"%\
						(len(phenotype_method_id2analysis_method_id_ls), no_of_associations_per_phenotype))
		return phenotype_method_id2analysis_method_id_ls
	
	def getAssociationOverlappingType(self, session, analysis_method_id_ls):
		"""
		2009-11-2
		"""
		sys.stderr.write("Getting AssociationOverlappingType for %s ..."%(repr(analysis_method_id_ls)))
		analysis_method_id_ls.sort()
		analysis_method_id_ls_in_str = map(str, analysis_method_id_ls)
		type_short_name = ','.join(analysis_method_id_ls_in_str)
		
		overlapping_type = Stock_250kDB.AssociationOverlappingType.get_by(short_name=type_short_name)
		if overlapping_type is None:
			overlapping_type = Stock_250kDB.AssociationOverlappingType(short_name=type_short_name)
			for analysis_method_id in analysis_method_id_ls:
				am = Stock_250kDB.AnalysisMethod.get(analysis_method_id)
				overlapping_type.analysis_method_ls.append(am)
			session.save(overlapping_type)
			session.flush()
		return overlapping_type
	
	results_id2snp_id_set = {}
	snp_info = None
	def calculateOverlappingStatForOneCombo(self, db, phenotype_method_id, call_method_id, analysis_method_id_ls, \
										no_of_top_snps=1000, association_overlapping_type=None, commit=False, \
										results_directory=None):
		"""
		2009-11-2
		"""
		sys.stderr.write("Calculating overlapping stat for phenotype %s and combo %s ...\n"%(phenotype_method_id, \
																						repr(analysis_method_id_ls),))
		session = db.session
		snp_id_set_ls = []
		for analysis_method_id in analysis_method_id_ls:
			rm = Stock_250kDB.ResultsMethod.query.filter_by(phenotype_method_id=phenotype_method_id).\
					filter_by(call_method_id=call_method_id).filter_by(analysis_method_id=analysis_method_id).first()
			if rm.id in self.results_id2snp_id_set:
				snp_id_set = self.results_id2snp_id_set.get(rm.id)
			else:
				association_entries = Stock_250kDB.Results.query.filter_by(results_id=rm.id).\
						filter(Stock_250kDB.Results.rank<=no_of_top_snps)
				no_of_association_entries = association_entries.count()
				if no_of_association_entries<no_of_top_snps:
					min_rank = no_of_association_entries+1
					max_rank = no_of_top_snps
					if self.snp_info is None:
						self.snp_info = DrawSNPRegion.getSNPInfo(db)
					ResultsMethod2Results.rm2result(session, rm, self.snp_info, min_rank=min_rank, max_rank=max_rank, \
												commit=commit, results_directory=results_directory)
					association_entries = Stock_250kDB.Results.query.filter_by(results_id=rm.id).\
							filter(Stock_250kDB.Results.rank<=no_of_top_snps)
				no_of_association_entries = association_entries.count()
				if no_of_association_entries!=no_of_top_snps:
					sys.stderr.write("Error: The number of SNPs %s from Result %s (analysis_method_id %s) doesn't match the no_of_top_snps %s.\n"%(no_of_association_entries, rm.id, rm.analysis_method_id, no_of_top_snps))
					return
				snp_id_set = set()
				for entry in association_entries:
					snp_id_set.add(entry.snps_id)
				self.results_id2snp_id_set[rm.id] = snp_id_set
			snp_id_set_ls.append(snp_id_set)
		overlapping_snp_id_set = snp_id_set_ls[0]
		
		for i in range(1, len(snp_id_set_ls)):
			snp_id_set = snp_id_set_ls[i]
			overlapping_snp_id_set = overlapping_snp_id_set&snp_id_set
		no_of_overlapping_snps = len(overlapping_snp_id_set)
		
		entry = Stock_250kDB.AssociationOverlappingStat(phenotype_method_id=phenotype_method_id, call_method_id=call_method_id, \
										no_of_top_snps=no_of_top_snps, no_of_overlapping_snps=no_of_overlapping_snps)
		entry.overlapping_type = association_overlapping_type
		session.save(entry)
		session.flush()
		sys.stderr.write("%s overlapping SNPs out of %s results. Done.\n"%(no_of_overlapping_snps, len(snp_id_set_ls)))
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  			password=self.db_passwd, hostname=self.hostname, database=self.dbname, 
				   			schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		phenotype_method_id2analysis_method_id_ls = self.getPhenotype2AnalysisMethods(self.call_method_id, self.phenotype_method_id_ls)
		
		for phenotype_method_id, analysis_method_id_ls in phenotype_method_id2analysis_method_id_ls.iteritems():
			subset_ls = listSubsets(analysis_method_id_ls)
			for subset in subset_ls:
				if len(subset)>1:	# more than 1 analysis methods included
					#subset.sort()
					association_overlapping_type = self.getAssociationOverlappingType(session, subset)
					db_entry = Stock_250kDB.AssociationOverlappingStat.query.filter_by(phenotype_method_id=phenotype_method_id).\
							filter_by(call_method_id=self.call_method_id).\
							filter_by(overlapping_type_id=association_overlapping_type.id).\
							filter_by(no_of_top_snps=self.no_of_top_snps).first()
					if db_entry:
						sys.stderr.write("Overlapping stat for call_method_id=%s, phenotype_method_id=%s, type %s, no_of_top_snps=%s already in db.\n"%\
										(self.call_method_id, phenotype_method_id, association_overlapping_type.id, self.no_of_top_snps))
						
					else:
						self.calculateOverlappingStatForOneCombo(db, phenotype_method_id, self.call_method_id, \
															subset, \
															association_overlapping_type=association_overlapping_type, \
															no_of_top_snps=self.no_of_top_snps, \
															commit=self.commit,\
															results_directory=self.results_directory)
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CalculateOverlappingStat
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
