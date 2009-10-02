#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results_by_gene/ -o ~/top_snp_test.out -s 100 -c -l 17

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -s 100 -b -c
	
	#use null_distribution_type 2 to do enrichment test among top 500 snps
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results/type_1/ -j 17 -f 500 -q10 -s ~/panfs/250k/snps_context_g0_m1000 -m 1000 -y 15 -C 2 -c
	
	#use min_score cutoff rather than no_of_top_snps, set a negative rank_gap, a new stop_rank. only results whose analysis_method_id=1 or 7.
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results/type_1/  -j 17 -f 200 -q80 -s ~/panfs/250k/snps_context_g0_m20000 -m 20000 -y 15 -C 3 -c -M 8 -E -0.2 -F 2 -e 1,7
	
	#ditto but with score cutoff from low to high
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results/type_1/  -j 17 -f 200 -q80 -s ~/panfs/250k/snps_context_g0_m20000 -m 20000 -y 15 -C 3 -c -M 2 -E 0.2 -F 8 -e 1,7
	
Description:
	MPI version TopSNPTest.py. No need to specify list_type_id and results_method_id_ls. Automatically calculates
	for all combinations of results_method_id and list_type_id, skipping the ones that have been done.
	
	2008-10-26
	Two major ways to specify how to do top SNP test.
		1. (no_of_top_snps, rank_gap, stop_rank) 
		2. (min_score, rank_gap, stop_rank). in this scenario, rank_gap, stop_rank become score_gap and stop_score.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, traceback
import Numeric, cPickle
from Scientific import MPI
from pymodule.MPIwrapper import mpi_synchronize, MPIwrapper
from pymodule import PassingData, getListOutOfStr
from TopSNPTest import TopSNPTest
import Stock_250kDB
from Stock_250kDB import Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, CandidateGeneTopSNPTest
from sets import Set
from MpiGeneListRankTest import MpiGeneListRankTest
from GeneListRankTest import SnpsContextWrapper
from common import get_total_gene_ls

class MpiTopSNPTest(TopSNPTest, MpiGeneListRankTest, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = TopSNPTest.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[100, 'q', 1, 'How many results one computing node should handle.']})
	option_default_dict.update({('call_method_id', 0, int):[0, 'j', 1, 'Restrict results based on this call_method. Default is no such restriction.']})
	option_default_dict.update({('analysis_method_id_ls', 0, ):[None, 'e', 1, 'Restrict results based on this set of analysis_method. Default is no such restriction. i.e., 1,7']})
	option_default_dict.update({("list_type_id_ls", 0, ): ['1-3,6,8,28,51,64,65,68,71,76,129', 'l', 1, 'comma/dash-separated list of gene list type ids. Each id has to encompass>=10 genes.']})
	option_default_dict.update({("phenotype_method_id_ls", 0, ): ['1-7,39-61,80-82', 'A', 1, 'Restrict results based on this set of phenotype_method ids.']})
	option_default_dict.pop(("list_type_id", 1, int))	#already poped in MpiGeneListRankTest
	option_default_dict.pop(("results_id_ls", 1, ))
	option_default_dict.update({('rank_gap', 1, float): [200, 'E', 1, 'the number of SNPs added onto previous no_of_top_snps to look for enrichment. negative gap could also be used. in that case, stop_rank means minimum score allowed.']})
	option_default_dict.update({('stop_rank', 1, float): [10000, 'F', 1, 'program iterates over i different types of i*no_of_top_snps until i*no_of_top_snps > this number. maximum rank allowed.']})
	option_default_dict.update({('alter_hostname', 1, ):['banyan.usc.edu', '', 1, 'host for non-output nodes to connect, since they only query and not save objects. this host can be a slave.']})
	option_default_dict.update({('store_null_data', 0, int):[0, 'S', 0, 'whether to store the NULL/permutation results for the enrichment test in the database']})
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		TopSNPTest.__init__(self, **keywords)
		self.list_type_id_ls = getListOutOfStr(self.list_type_id_ls, data_type=int)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	def addCutoffToParamsLs(self, cutoff_ls, params_ls):
		"""
		2008-11-04
			followon for generate_params() of MpiGeneListRankTest.py to add cutoff's into params_ls
		"""
		sys.stderr.write("Adding cutoff to params_ls ...")
		new_params_ls = []
		for params_tup in params_ls:
			for cutoff in cutoff_ls:
				new_params_ls.append(list(params_tup)+[cutoff])
		
		sys.stderr.write("%s parameters. Done.\n"%(len(new_params_ls)))
		return new_params_ls
	
	def generate_cutoff_ls(self, no_of_top_snps, min_score, rank_gap, stop_rank):
		"""
		2008-11-04
			split out of computing_node_handler. common for every computing run
		
			#construct min_score_ls & no_of_top_snps_ls
		"""
		sys.stderr.write("Generating cutoff_ls ...")
		if rank_gap<0:
			stop_marker = -stop_rank
		else:
			stop_marker = stop_rank
		current_marker = stop_marker - 1
		if min_score is None:	#this is to generate cutoffs for no_of_top_snps_ls
			if rank_gap<0:
				current_marker = min(current_marker, -no_of_top_snps)	#both in min() are negative
			else:
				current_marker = min(current_marker, no_of_top_snps)	#both in min() are positive
		cutoff_ls = []
		i = 0
		while current_marker<stop_marker:	#add one more layer to look at certain top genes
			if min_score is not None:
				current_marker = min_score +i*rank_gap
			else:
				t = min(max(0, int(abs(current_marker)/4000)), 6)	#enlarge the rank_gap as the rank goes further, every 4000 is a window.
				#t must be bigger than 0 and less than 6
				
				coeffcient = math.pow(2,t)/(rank_gap/50.)	#it kind of makes rank_gap useless. always assume rank_gap=50
				#sys.stderr.write("i=%s, t=%s, coeffcient=%s, rank_gap=%s, current_marker=%s, no_of_top_snps=%s.\n"%(i, t, coeffcient, rank_gap, current_marker, no_of_top_snps))
				current_marker = int(current_marker + coeffcient*rank_gap)
			cutoff_ls.append(current_marker)
			
			if rank_gap<0:
				current_marker = -current_marker
			else:
				current_marker = current_marker
			i += 1
		if self.debug:
			sys.stderr.write('%s\n'%repr(cutoff_ls))
		sys.stderr.write("%s cutoffs. Done.\n"%(len(cutoff_ls)))
		return cutoff_ls
		
	def computing_node_handler(self, communicator, data, comp_param_obj):
		"""
		2009-1-22
			deal with option self.store_null_data
		2008-11-12
			turn runHGTest() back into life
			turn off runEnrichmentTestToGetNullData()
		2008-10-31
			runEnrichmentTestToGetNullData() is gonna get data at all different no_of_top_snps's or min_score's
		2008-10-26
			handle (min_score, rank_gap, stop_rank)
			handle scenario that rank_gap is negative and so the parameters tried are descending.
		2008-08-20
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		null_data_ls = []
		pd = PassingData(snps_context_wrapper=comp_param_obj.snps_context_wrapper,\
							no_of_total_genes=comp_param_obj.no_of_total_genes, \
							results_directory=comp_param_obj.results_directory, \
							min_MAF=comp_param_obj.min_MAF, \
							get_closest=self.get_closest, 
							min_distance=self.min_distance, \
							no_of_top_snps=self.no_of_top_snps, #2008-10-25 no_of_top_snps is useless. overwritten later
							min_sample_size=self.min_sample_size,
							test_type_id=self.test_type_id, \
							results_type=self.results_type, 
							no_of_permutations=self.no_of_permutations,\
							no_of_min_breaks=self.no_of_min_breaks,
							type_id=comp_param_obj.type_id,\
							null_distribution_type_id=self.null_distribution_type_id,\
							allow_two_sample_overlapping=self.allow_two_sample_overlapping,
							total_gene_id_ls=comp_param_obj.total_gene_id_ls,\
							min_score=self.min_score,
							commit=self.commit)	#2008-10-25 min_score is useless. overwritten later
		#2008-10-25
		#if rank_gap is negative, stop_marker means the minimum cutoff
		#if rank_gap is positive, stop_marker means the maximum cutoff
		#both signs have to be swapped in the case of negative rank_gap
		"""	
		if self.rank_gap<0:
			stop_marker = -self.stop_rank
		else:
			stop_marker = self.stop_rank
		
		for results_id, list_type_id in data:
			if self.debug:
				sys.stderr.write("working on results_id=%s, list_type_id=%s, type_id=%s .\n"%(results_id, list_type_id, pd.type_id))
			i = 0
			#reset it to zero!!
			if self.rank_gap<0:	#has to be less than -self.stop_rank in order to pass first round. because stop_marker=-stop_rank when rank_gap<0.
				current_marker = stop_marker - 1
			else:
				current_marker = stop_marker -1
			
			while current_marker<stop_marker:	#add one more layer to look at certain top genes
				if self.min_score is not None:
					current_marker = self.min_score +i*self.rank_gap
					pd.min_score = current_marker
				else:
					current_marker = self.no_of_top_snps + i*self.rank_gap
					pd.no_of_top_snps = current_marker
				
				if self.rank_gap<0:
					current_marker = -current_marker
				else:
					current_marker = current_marker
				
				pd.results_id = results_id
				pd.list_type_id = list_type_id
				if self.debug:
					sys.stderr.write("working on results_id=%s, list_type_id=%s, current_marker=%s.\n"%\
									(pd.results_id, pd.list_type_id, current_marker))
				i += 1
				result = self.runHGTest(pd)
				if result is not None:
					result_ls.append(result)
		"""
		
		pd.commit = 0	#commit once afterwards. commit runtime would render 'Lock wait timeout exceeded; try restarting transaction'
		for results_id, list_type_id, cutoff in data:
			if self.debug:
				sys.stderr.write("working on results_id=%s, list_type_id=%s, type_id=%s, cutoff %s.\n"%(results_id, list_type_id, pd.type_id, cutoff))
			pd.results_id = results_id
			pd.list_type_id = list_type_id
			if self.min_score:
				pd.min_score_ls = [cutoff]
				pd.min_score = cutoff
			else:
				pd.no_of_top_snps_ls = [cutoff]
				pd.no_of_top_snps = cutoff
			if self.store_null_data:
				return_data = self.runEnrichmentTestToGetNullData(comp_param_obj.session, pd)
			else:
				return_data = self.runHGTest(pd)
			if return_data:
				result_ls += return_data.result_ls
				null_data_ls += return_data.null_data_ls
		
		#if self.commit:
		#	comp_param_obj.session.flush()
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return_data = PassingData(result_ls=result_ls, null_data_ls=null_data_ls)
		return return_data
	
	"""
	2008-10-15
		mothball output_node_handler() cuz it's exactly same as MpiGeneListRankTest.py
		
		2008-08-20
			derived from MpiGeneListRankTest.output_node_handler
	def output_node_handler(self, communicator, parameter_list, data):
		
		writer, session, commit, TestResultClass= parameter_list
		table_obj_ls = cPickle.loads(data)
		for table_obj in table_obj_ls:
			row = []
			
			result = TestResultClass()
			#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
			#can't save table_obj because it's associated with a different db thread
			for column in table_obj.c.keys():
				row.append(getattr(table_obj, column))
				setattr(result, column, getattr(table_obj, column))
			if writer:
				writer.writerow(row)
			session.save(result)
			if commit:
				session.flush()
	"""
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2008-11-04
			save both results and null_data(only for runEnrichmentTestToGetNullData()
		"""
		return_data = cPickle.loads(data)
		#must save result_ls first, to allow new results have ids assigned
		self.sub_output(output_param_obj, return_data.result_ls, output_param_obj.TestResultClass)
		
		self.saveNullData(output_param_obj, return_data.null_data_ls, output_param_obj.TestResultClass)
	
	#2009-9-14 report the progress of how many were saved
	saved_null_data_count = 0
	unsaved_null_data_count = 0
	
	def saveNullData(self, output_param_obj, table_obj_ls, TestResultClass):
		"""
		2008-11-04
		
			similar to sub_output() but has to handle unsaved table_obj.observed
		"""
		commit = output_param_obj.commit
		for table_obj in table_obj_ls:
			row = []
			
			result = Stock_250kDB.TopSNPTestRMNullData()
			#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
			#can't save table_obj because it's associated with a different db thread
			for column in table_obj.c.keys():
				if output_param_obj.writer:	#2008-10-30 append only when writer is not None
					row.append(getattr(table_obj, column))
				setattr(result, column, getattr(table_obj, column))
			
			if output_param_obj.writer:
				output_param_obj.writer.writerow(row)
			
			if table_obj.observed.id is not None:	#handle the previously unsaved TestResultClass'es
				result.observed_id = table_obj.observed.id
			else:
				new_observed = self.returnResultFromDB(TestResultClass, table_obj.observed.results_id, table_obj.observed.list_type_id,
													table_obj.observed.starting_rank,\
											table_obj.observed.type_id, table_obj.observed.min_distance, table_obj.observed.no_of_top_snps,\
											table_obj.observed.min_score)
				if new_observed:
					result.observed_id = new_observed.id
				else:
					sys.stderr.write("TestResultClass with results_id=%s, list_type_id=%s, starting_rank=%s,\
							type_id=%s, min_distance=%s, no_of_top_snps=%s, min_score=%s) is not found for this null data.\n"%\
							(table_obj.observed.results_id, table_obj.observed.list_type_id, table_obj.observed.starting_rank,
							table_obj.observed.type_id, table_obj.observed.min_distance, table_obj.observed.no_of_top_snps, \
							table_obj.observed.min_score))
					continue	#skip if TestResultClass is not found
			
			output_param_obj.session.save(result)
			if commit:
				try:
					output_param_obj.session.flush()
					self.saved_null_data_count += 1
				except:
					#2008-10-30 remove it from memory. otherwise, next flush() will try on this old object again.
					output_param_obj.session.expunge(result)
					#output_param_obj.session.delete(result)
					sys.stderr.write("Exception happened when saving null data with observed_id=%s, run_no=%s, null_distribution_type_id=%s.\n"%\
									(getattr(result, 'observed_id', None),\
									getattr(result, 'run_no', None),\
									getattr(result, 'null_distribution_type_id', None)))
					for column in table_obj.c.keys():
						sys.stderr.write("\t%s=%s.\n"%(column, getattr(table_obj, column)))
					traceback.print_exc()
					sys.stderr.write('%s.\n'%repr(sys.exc_info()))
					self.unsaved_null_data_count += 1
		sys.stderr.write("%s(%s null data saved, %s unsaved).\n"%(self.saved_null_data_count+self.unsaved_null_data_count, self.saved_null_data_count,\
																 self.unsaved_null_data_count))  
	
	def run(self):
		"""
		2008-08-20
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		#2008-10-30 comment out because computing node is gonna save the stuff itself.
		if node_rank!=output_node_rank:		#to reduce the number of connections/queries to the master
			self.hostname = self.alter_hostname
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		if self.results_type==1:
			ResultsClass = Stock_250kDB.ResultsMethod
			TestResultClass = Stock_250kDB.CandidateGeneTopSNPTestRM
		elif self.results_type==2:
			ResultsClass = Stock_250kDB.ResultsByGene
			TestResultClass = Stock_250kDB.CandidateGeneTopSNPTest
		elif self.results_type==3:
			ResultsClass = Stock_250kDB.ResultsMethod
			TestResultClass = Stock_250kDB.CandidateGeneTopSNPTestRG
		else:
			sys.stderr.write("Invalid results type : %s.\n"%pd.results_type)
			sys.exit(3)
		
		if node_rank == 0:
			pdata_for_computing = PassingData()
			pdata_for_computing.total_gene_id_ls = get_total_gene_ls(db.metadata.bind)
			pdata_for_computing.no_of_total_genes = len(pdata_for_computing.total_gene_id_ls)
			param_obj = PassingData(call_method_id=self.call_method_id, \
								analysis_method_id=getattr(self, 'analysis_method_id', None),\
								analysis_method_id_ls=getattr(self, 'analysis_method_id_ls', None),\
								phenotype_method_id_ls=getattr(self, 'phenotype_method_id_ls', None),\
								list_type_id_ls=self.list_type_id_ls, \
								results_type=self.results_type)
			params_ls = self.generate_params(param_obj, self.min_no_of_genes)
			cutoff_ls = self.generate_cutoff_ls(self.no_of_top_snps, self.min_score, self.rank_gap, self.stop_rank)
			params_ls = self.addCutoffToParamsLs(cutoff_ls, params_ls)
			
			pdata_for_computing.snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
			if self.debug:
				params_ls = params_ls[:100]
			pdata_for_computing_pickle = cPickle.dumps(pdata_for_computing, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(pdata_for_computing_pickle, node, 0)
				sys.stderr.write(".\n")
			del pdata_for_computing_pickle
			del pdata_for_computing
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			data =  cPickle.loads(data)
			sys.stderr.write(".\n")
		else:
			pass
		
		_type = self.getTopSNPTestType(self.get_closest, self.min_MAF, \
										self.allow_two_sample_overlapping, self.results_type,\
										self.test_type_id, self.null_distribution_type_id)
		self.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls]
			self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			comp_param_obj = PassingData(snps_context_wrapper=data.snps_context_wrapper, \
												results_directory=self.results_directory, min_MAF=self.min_MAF,\
												no_of_total_genes=data.no_of_total_genes, \
												total_gene_id_ls=data.total_gene_id_ls,\
												type_id=_type.id,	#_type is placeholder. output_node decides on this.
												session=session)
			self.computing_node(comp_param_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				header_row = []
				for column in TestResultClass.c.keys():
					header_row.append(column)
				writer.writerow(header_row)
			else:
				writer = None
			
			output_param_obj = PassingData(writer=writer, session=session, commit=self.commit, TestResultClass=TestResultClass,
										_type=_type)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del writer		
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiTopSNPTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
