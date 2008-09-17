#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -t ~/panfs/db/results_by_gene/ -o ~/top_snp_test.out -s 100 -c -l 17

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiTopSNPTest.py -u yh -p passw**d -s 100 -b -c
	
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
from pymodule import PassingData
from TopSNPTest import TopSNPTest
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, CandidateGeneTopSNPTest
from sets import Set
from MpiGeneListRankTest import MpiGeneListRankTest
from GeneListRankTest import SnpsContextWrapper

class MpiTopSNPTest(TopSNPTest, MpiGeneListRankTest, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = TopSNPTest.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[200, 's', 1, 'How many results one computing node should handle.']})
	option_default_dict.update({('call_method_id', 0, int):[0, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.']})
	option_default_dict.update({('analysis_method_id', 0, int):[0, '', 1, 'Restrict results based on this analysis_method. Default is no such restriction.']})
	option_default_dict.pop(("list_type_id", 1, int))	#already poped in MpiGeneListRankTest
	option_default_dict.pop(("results_id_ls", 1, ))
	
	def __init__(self,  **keywords):
		"""
		2008-08-20
		"""
		TopSNPTest.__init__(self, **keywords)
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-08-20
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		pd = PassingData(snps_context_wrapper=computing_parameter_obj.snps_context_wrapper,\
							no_of_total_genes=computing_parameter_obj.no_of_total_genes, results_directory=computing_parameter_obj.results_directory, \
							min_MAF=computing_parameter_obj.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance, \
							no_of_top_snps=self.no_of_top_snps)
		for results_method_id, list_type_id in data:
			pd.results_method_id = results_method_id
			pd.list_type_id = list_type_id
			result = self.runHGTest(pd)
			if result is not None:
				result_ls.append(result)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		2008-08-20
			derived from MpiGeneListRankTest.output_node_handler
		"""
		writer, session, commit = parameter_list
		table_obj_ls = cPickle.loads(data)
		for table_obj in table_obj_ls:
			row = []
			
			result = CandidateGeneTopSNPTest()
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
	
	def run(self):
		"""
		2008-08-20
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		if node_rank == 0:
			pdata_for_computing = PassingData()
			#pdata_for_computing.snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
			pdata_for_computing.no_of_total_genes = self.getNoOfTotalGenes(db, self.gene_table, self.tax_id)
			param_obj = PassingData(call_method_id=self.call_method_id, analysis_method_id=getattr(self, 'analysis_method_id', None))
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
			#snps_context_wrapper = data.snps_context_wrapper
			no_of_total_genes = data.no_of_total_genes
			del data
			sys.stderr.write(".\n")
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls]
			self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snps_context_wrapper=None, \
												results_directory=self.results_directory, min_MAF=self.min_MAF,\
												no_of_total_genes=no_of_total_genes)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				header_row = []
				for column in CandidateGeneTopSNPTest.c.keys():
					header_row.append(column)
				writer.writerow(header_row)
			else:
				writer = None
			
			parameter_list = [writer, session, self.commit]
			self.output_node(free_computing_nodes, parameter_list, self.output_node_handler)
			del writer		
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiTopSNPTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
