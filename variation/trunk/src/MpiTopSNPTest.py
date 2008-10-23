#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results_by_gene/ -o ~/top_snp_test.out -s 100 -c -l 17

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -s 100 -b -c
	
	#use null_distribution_type 2 to do enrichment test among top 500 snps
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results/type_1/ -j 17 -f 500 -q10 -s ~/panfs/250k/snps_context_g0_m1000 -m 1000 -y 15 -C 2 -c
	
Description:
	MPI version TopSNPTest.py. No need to specify list_type_id and results_method_id_ls. Automatically calculates
	for all combinations of results_method_id and list_type_id, skipping the ones that have been done.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
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
	option_default_dict.update({('analysis_method_id', 0, int):[0, '', 1, 'Restrict results based on this analysis_method. Default is no such restriction.']})
	option_default_dict.update({("list_type_id_ls", 0, ): [None, 'l', 1, 'comma/dash-separated list of gene list type ids. ids not present in db will be filtered out. Each id has to encompass>=10 genes.']})
	option_default_dict.pop(("list_type_id", 1, int))	#already poped in MpiGeneListRankTest
	option_default_dict.pop(("results_id_ls", 1, ))
	option_default_dict.update({('starting_rank_gap', 1, int): [50, '', 1, 'the gap between the rank of the leading snp in this window and that of the next window. deprecated for now.']})
	option_default_dict.update({('stop_rank', 1, int): [5000, '', 1, 'program iterates over i different types of i*no_of_top_snps until i*no_of_top_snps > this number.']})
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		TopSNPTest.__init__(self, **keywords)
		self.list_type_id_ls = getListOutOfStr(self.list_type_id_ls, data_type=int)
	
	def computing_node_handler(self, communicator, data, comp_param_obj):
		"""
		2008-08-20
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		pd = PassingData(snps_context_wrapper=comp_param_obj.snps_context_wrapper,\
							no_of_total_genes=comp_param_obj.no_of_total_genes, \
							results_directory=comp_param_obj.results_directory, \
							min_MAF=comp_param_obj.min_MAF, \
							get_closest=self.get_closest, 
							min_distance=self.min_distance, \
							no_of_top_snps=self.no_of_top_snps, 
							min_sample_size=self.min_sample_size, 
							test_type_id=self.test_type_id, \
							results_type=self.results_type, 
							no_of_permutations=self.no_of_permutations,\
							no_of_min_breaks=self.no_of_min_breaks,
							type_id=comp_param_obj.type_id,\
							null_distribution_type_id=self.null_distribution_type_id,\
							allow_two_sample_overlapping=self.allow_two_sample_overlapping,
							total_gene_id_ls=comp_param_obj.total_gene_id_ls)
		for results_id, list_type_id in data:
			i = 0
			while pd.no_of_top_snps<self.stop_rank:	#add one more layer to look at certain top genes
				i += 1
				pd.no_of_top_snps = self.no_of_top_snps*i
				pd.results_id = results_id
				pd.list_type_id = list_type_id
				result = self.runHGTest(pd)
				if result is not None:
					result_ls.append(result)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
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
	
	def run(self):
		"""
		2008-08-20
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
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
		
		if node_rank == 0:
			pdata_for_computing = PassingData()
			pdata_for_computing.snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
			pdata_for_computing.total_gene_id_ls = get_total_gene_ls(db.metadata.bind)
			pdata_for_computing.no_of_total_genes = len(pdata_for_computing.total_gene_id_ls)
			param_obj = PassingData(call_method_id=self.call_method_id, analysis_method_id=getattr(self, 'analysis_method_id', None),\
								list_type_id_ls=self.list_type_id_ls, results_type=self.results_type)
			params_ls = self.generate_params(param_obj)
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
		
		self.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls]
			self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			_type = self.getTopSNPTestType(self.min_distance, self.get_closest, self.min_MAF, \
										self.allow_two_sample_overlapping, self.results_type,\
										self.test_type_id, self.null_distribution_type_id)
			comp_param_obj = PassingData(snps_context_wrapper=data.snps_context_wrapper, \
												results_directory=self.results_directory, min_MAF=self.min_MAF,\
												no_of_total_genes=data.no_of_total_genes, \
												total_gene_id_ls=data.total_gene_id_ls,\
												type_id=_type.id)	#_type is placeholder. output_node decides on this.
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
			
			_type = self.getTopSNPTestType(self.min_distance, self.get_closest, self.min_MAF, \
										self.allow_two_sample_overlapping, self.results_type,\
										self.test_type_id, self.null_distribution_type_id)
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
