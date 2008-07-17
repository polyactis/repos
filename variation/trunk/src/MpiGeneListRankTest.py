#!/usr/bin/env mpipython
"""

Examples:
	#test parallel run on desktop, using Strain X SNP format
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiGeneListRankTest.py -c -e 389,190 -l 1 -u yh -b -p passw**d
	
Description:
	MPI version GeneListRankTest.py
	
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
from GeneListRankTest import GeneListRankTest
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, CandidateGeneRankSumTestResult

class MpiGeneListRankTest(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict
	option_default_dict.update({('message_size', 1, int):[1, 's', 1, 'How many results one computing node should handle.']})
	
	def __init__(self,  **keywords):
		GeneListRankTest.__init__(self, **keywords)
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		results_method_id_ls = parameter_list[0]
		data_to_return = []
		for i in range(message_size):
			if len(results_method_id_ls)>0:
				results_method_id = results_method_id_ls.pop(0)
				data_to_return.append(results_method_id)
			else:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2007-03-07
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		for results_method_id in data:
			result = self.run_wilcox_test(results_method_id, computing_parameter_obj.snps_context_wrapper, computing_parameter_obj.list_type_id)
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
			row = [results_method_id, table_obj.list_type_id, table_obj.pvalue, table_obj.statistic]
			if writer:
				writer.writerow(row)
			if commit:
				session.begin()
			session.save(table_obj)
			session.flush()
			if commit:
				session.commit()
	
	def run(self):
		"""
		2008-07-17
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		session = db.session
		
		if node_rank == 0:
			snps_context_wrapper = self.constructDataStruc(self.min_distance)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				data_pickle = cPickle.dumps(snps_context_wrapper, -1)
				self.communicator.send(data_pickle, node, 0)
				sys.stderr.write("good")
				del data_pickle
				sys.stderr.write(".\n")
			del snps_context_wrapper
		elif node_rank in free_computing_nodes:
			sys.stderr.write("receiving ...")
			data, source, tag = self.communicator.receiveString(0, 0)
			sys.stderr.write("got.\n")
			snps_context_wrapper =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		else:
			pass
		
		mw = MPIwrapper(self.communicator, debug=self.debug, report=self.report)
		mw.synchronize()
		if node_rank == 0:
			parameter_list = [self.results_method_id_ls]
			mw.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_nodes:
			computing_parameter_obj = PassingData(snps_context_wrapper=snps_context_wrapper, list_type_id=self.list_type_id)
			mw.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				writer.writerow(['results_method_id', 'list_type_id', 'wilcox.test.pvalue', 'statistic'])
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