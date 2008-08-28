#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiQC149CrossMatch.py -m 4 -u yh -p passw**d -s 1000 -c

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
Description:
	MPI version QC_149_cross_match.py. Message_size determines how many pairs one computing node should handle.
	
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
import matplotlib; matplotlib.use("Agg")	#to disable interactive. so that "import networkx as nx" in variation.src.__init__ won't cause "RuntimeError: could not create GdkCursor object"
from variation.src import StockDB
from QC_149_cross_match import QC_149_cross_match
from sets import Set
from pymodule.db import formReadmeObj

class MpiQC149CrossMatch(QC_149_cross_match):
	__doc__ = __doc__
	option_default_dict = QC_149_cross_match.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[200, 's', 1, 'How many results one computing node should handle.']})
	
	def __init__(self,  **keywords):
		QC_149_cross_match.__init__(self, **keywords)
	
	def input_node(self, communicator, parameter_list, free_computing_nodes, message_size, report=0):
		"""
		2008-08-28
		"""
		node_rank = communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		twoSNPData = parameter_list[0]
		output_node_rank = parameter_list[1]
		counter = 0
		row_id_pair_list = []
		for row_id1 in twoSNPData.SNPData1.row_id_ls:
			for row_id2 in twoSNPData.SNPData2.row_id_ls:
				row_id_pair_list.append((row_id1, row_id2))
				if len(row_id_pair_list)==message_size:
					communicator.send("1", output_node_rank, 1)	#WATCH: tag is 1, to the output_node.
					free_computing_node, source, tag = communicator.receiveString(output_node_rank, 2)
					#WATCH: tag is 2, from the output_node
					data_pickle = cPickle.dumps(row_id_pair_list, -1)
					communicator.send(data_pickle, int(free_computing_node),0)	#WATCH: int()
					row_id_pair_list = []	#clear the list
					if report:
						sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
					counter += 1
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-08-28
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		twoSNPData = computing_parameter_obj.twoSNPData
		QC_method_id = computing_parameter_obj.QC_method_id
		for row_id1, row_id2 in data:
			NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs = twoSNPData.cmpOneRow(row_id1, row_id2)
			#the 2nd position in the row-id1 tuple is strain id
			if QC_method_id==4:	#the 2nd position in the row-id2 tuple is strain id
				target_id = row_id2[1]
			else:
				target_id = row_id2
			qc_cross_match = PassingData(strainid=row_id1[1], target_id=target_id, mismatch_rate=mismatch_rate, \
												no_of_mismatches=no_of_mismatches, no_of_non_NA_pairs=no_of_non_NA_pairs)
			result_ls.append(qc_cross_match)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		2008-08-28
		"""
		session, commit, QC_method_id, readme = parameter_list
		table_obj_ls = cPickle.loads(data)
		for table_obj in table_obj_ls:			
			qc_cross_match = StockDB.QCCrossMatch()
			#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
			#can't save table_obj because it's associated with a different db thread
			for column in qc_cross_match.c.keys():	#it's qc_cross_match, not table_obj because table_obj is not linked to db.
				setattr(qc_cross_match, column, getattr(table_obj, column, None))
			qc_cross_match.qc_method_id = QC_method_id
			qc_cross_match.readme = readme
			session.save(qc_cross_match)
			if commit:
				session.flush()
	
	def run(self):
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank not in free_computing_node_set:	#computing nodes don't need db connection
			db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
								password=self.db_passwd, hostname=self.hostname, database=self.dbname)
			db.setup(create_tables=False)
			session = db.session
		
		if node_rank == 0:
			twoSNPData = self.prepareTwoSNPData(db)
			twoSNPData_pickle = cPickle.dumps(twoSNPData, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(twoSNPData_pickle, node, 0)
				sys.stderr.write(".\n")
			del twoSNPData_pickle
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			twoSNPData =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		else:
			readme = formReadmeObj(sys.argv, self.ad, StockDB.README)
			session.save(readme)
		
		mw = MPIwrapper(self.communicator, debug=self.debug, report=self.report)
		mw.synchronize()
		if node_rank == 0:
			parameter_list = [twoSNPData, output_node_rank]
			self.input_node(self.communicator, parameter_list, free_computing_nodes, message_size=self.message_size, report=self.report)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(twoSNPData=twoSNPData, QC_method_id=self.QC_method_id)
			mw.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			parameter_list = [session, self.commit, self.QC_method_id, readme]
			mw.output_node(free_computing_nodes, parameter_list, self.output_node_handler)
		mw.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiQC149CrossMatch
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()