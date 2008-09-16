#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiRM2ResultsByGene.py -o ~/panfs/db/results_by_gene/ -c -u yh -s ~/panfs/250k/snps_context_g0_m20000 -m 20000 -t ~/panfs/db/results/type_1/ -p passw**d -l 17 -a 4

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiRM2ResultsByGene.py.py ...
	
Description:
	MPI version ResultsMethod2ResultsByGene.py.
	
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
from ResultsMethod2ResultsByGene import ResultsMethod2ResultsByGene
import Stock_250kDB
from sets import Set
from GeneListRankTest import SnpsContextWrapper

class MpiRM2ResultsByGene(ResultsMethod2ResultsByGene, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = ResultsMethod2ResultsByGene.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[20, 'i', 1, 'How many results one computing node should handle.']})
	
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def computing_node_handler(self, communicator, data, param_data):
		"""
		2008-09-16
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		results_method_id_ls = cPickle.loads(data)
		result_ls = []
		counter = 0
		for results_method_id in results_method_id_ls:
			rm = Stock_250kDB.ResultsMethod.get(results_method_id)
			if not rm:
				sys.stderr.write("No results method available for results_method_id=%s.\n"%results_method_id)
				continue
			self.saveResultsByGene(param_data.session, rm, param_data)
			counter += 1
		
		if param_data.commit:
			param_data.session.commit()
			param_data.session.clear()
			param_data.session.begin()	#restart the transaction
		else:
			param_data.session.rollback()
			param_data.session.begin()
		sys.stderr.write("Node no.%s done with %s jobs.\n"%(node_rank, counter))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-09-16
			do nothing
		"""
		pass
	
	def run(self):
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		if node_rank == 0:
			snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
			if not self.results_method_id_ls:
				pdata = PassingData(call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id)
				self.results_method_id_ls = self.getResultsMethodIDLs(pdata)
			
			snps_context_wrapper_pickle = cPickle.dumps(snps_context_wrapper, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(snps_context_wrapper_pickle, node, 0)
				sys.stderr.write(".\n")
			del snps_context_wrapper_pickle, snps_context_wrapper
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			snps_context_wrapper =  cPickle.loads(data)
			del data
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(params_ls=self.results_method_id_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.input_node(param_obj, free_computing_nodes, input_handler=self.input_fetch_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			param_data = PassingData(session=session)
			param_data.results_directory = self.input_db_directory
			param_data.default_output_db_directory = self.default_output_db_directory
			param_data.output_db_directory = self.output_db_directory
			param_data.commit = self.commit
			param_data.min_MAF = self.min_MAF
			param_data.min_distance = self.min_distance
			param_data.get_closest = self.get_closest
			param_data.snps_context_wrapper = snps_context_wrapper
			self.computing_node(param_data, self.computing_node_handler)
		else:
			param_obj = PassingData()
			self.output_node(free_computing_nodes, param_obj, self.output_node_handler)
		self.synchronize()	#to avoid some node early exits


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiRM2ResultsByGene
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()