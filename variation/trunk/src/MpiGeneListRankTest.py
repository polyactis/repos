#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiGeneListRankTest.py -u yh -p passw**d -t ~/panfs/db/results/type_1/ -o ~/mpigene_list_rank_test.out -s 100 -c

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiGeneListRankTest.py -u yh -m 20000 -g -p passw**d -s 100 -b -c
	
Description:
	MPI version GeneListRankTest.py. No need to specify list_type_id and results_method_id_ls. Automatically calculates
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
from GeneListRankTest import GeneListRankTest
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, CandidateGeneRankSumTestResult
from sets import Set


class MpiGeneListRankTest(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[200, 's', 1, 'How many results one computing node should handle.']})
	option_default_dict.update({('call_method_id', 0, int):[0, 'l', 1, 'Restrict results based on this call_method. Default is no such restriction.']})
	option_default_dict.pop(("list_type_id", 1, int))
	option_default_dict.pop(("results_method_id_ls", 1, ))
	
	def __init__(self,  **keywords):
		GeneListRankTest.__init__(self, **keywords)
	
	def generate_params(self, call_method_id, min_no_of_genes=10):
		"""
		2008-08-19
			add call_method_id
		2008-08-15
			stop filtering if CandidateGeneRankSumTestResult has (results_method_id, list_type_id) combo
		2008-07-24
			only association results (results_method_type_id=1)
			only candidate gene lists with >min_no_of_genes genes
			skip ones that been done
		"""
		sys.stderr.write("Generating parameters ...")
		i = 0
		block_size = 5000
		if call_method_id!=0:
			query = ResultsMethod.query.filter_by(results_method_type_id=1).filter_by(call_method_id=call_method_id)
		else:
			query = ResultsMethod.query.filter_by(results_method_type_id=1)
		rows = query.offset(i).limit(block_size)
		results_method_id_ls = []
		while rows.count()!=0:
			for row in rows:
				results_method_id_ls.append(row.id)
				i += 1
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("%s results. "%(len(results_method_id_ls)))
		
		i = 0
		rows = GeneListType.query.offset(i).limit(block_size)
		list_type_id_ls = []
		while rows.count()!=0:
			for row in rows:
				if len(row.gene_list)>=min_no_of_genes:
					list_type_id_ls.append(row.id)
				i += 1
			rows = GeneListType.query.offset(i).limit(block_size)
		sys.stderr.write("%s candidate gene lists. "%(len(list_type_id_ls)))
		
		rm_id_lt_id_set = Set()
		"""
		i = 0
		rows = CandidateGeneRankSumTestResult.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				rm_id_lt_id_set.add((row.results_method_id, row.list_type_id))
				i += 1
			rows = CandidateGeneRankSumTestResult.query.offset(i).limit(block_size)
		sys.stderr.write("%s candidate gene rank sum test results. "%(len(rm_id_lt_id_set)))
		"""
		
		params_ls = []
		for results_method_id in results_method_id_ls:
			for list_type_id in list_type_id_ls:
				rm_id_lt_id = (results_method_id, list_type_id)
				if rm_id_lt_id not in rm_id_lt_id_set:
					params_ls.append(rm_id_lt_id)
		sys.stderr.write("generating params done.\n")
		return params_ls
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		params_ls = parameter_list[0]
		data_to_return = []
		for i in range(message_size):
			if len(params_ls)>0:
				one_parameter = params_ls.pop(0)
				data_to_return.append(one_parameter)
			else:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-08-20
			wrap all parameters into pd and pass it to run_wilcox_test
		2008-07-17
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		pd = PassingData(snps_context_wrapper=computing_parameter_obj.snps_context_wrapper,\
							results_directory=computing_parameter_obj.results_directory, max_pvalue_per_gene=self.max_pvalue_per_gene,\
							min_MAF=computing_parameter_obj.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance, \
							min_sample_size=self.min_sample_size)
		for results_method_id, list_type_id in data:
			pd.results_method_id = results_method_id
			pd.list_type_id = list_type_id
			result = self.run_wilcox_test(pd)
			if result is not None:
				result_ls.append(result)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		2008-05-19
			row is strain. col is snp. reversed due to utilization of SNPData
		05/14/2008
			flush outf and outf_avg
		05/12/2008
			common_var_name_ls
		"""
		writer, session, commit = parameter_list
		table_obj_ls = cPickle.loads(data)
		for table_obj in table_obj_ls:
			row = []
			
			candidate_gene_rank_sum_test_result = CandidateGeneRankSumTestResult()
			#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
			#can't save table_obj because it's associated with a different db thread
			for column in table_obj.c.keys():
				row.append(getattr(table_obj, column))
				setattr(candidate_gene_rank_sum_test_result, column, getattr(table_obj, column))
			if writer:
				writer.writerow(row)
			session.save(candidate_gene_rank_sum_test_result)
			if commit:
				session.flush()
	
	def run(self):
		"""
		2008-07-17
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		if node_rank in free_computing_node_set:	#to reduce the number of connections on papaya
			self.hostname = 'banyan.usc.edu'
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		
		if node_rank == 0:
			snps_context_wrapper = self.constructDataStruc(self.min_distance, self.get_closest)
			params_ls = self.generate_params(self.call_method_id)
			if self.debug:
				params_ls = params_ls[:100]
			snps_context_wrapper_pickle = cPickle.dumps(snps_context_wrapper, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(snps_context_wrapper_pickle, node, 0)
				sys.stderr.write(".\n")
			del snps_context_wrapper_pickle
			del snps_context_wrapper
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			snps_context_wrapper =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		else:
			pass
		
		mw = MPIwrapper(self.communicator, debug=self.debug, report=self.report)
		mw.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls]
			mw.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snps_context_wrapper=snps_context_wrapper, \
												results_directory=self.results_directory, min_MAF=self.min_MAF)
			mw.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				header_row = []
				for column in CandidateGeneRankSumTestResult.c.keys():
					header_row.append(column)
				writer.writerow(header_row)
			else:
				writer = None
			
			parameter_list = [writer, session, self.commit]
			mw.output_node(free_computing_nodes, parameter_list, self.output_node_handler)
			del writer		
		mw.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiGeneListRankTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
