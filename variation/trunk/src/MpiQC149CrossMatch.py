#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiQC149CrossMatch.py -m 4 -u yh -p passw**d -s 1000 -c

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
Description:
	MPI version QC_149_cross_match.py. Message_size determines how many pairs one computing node should handle.

	If it's 149 self-cross-match (QC_method_id=4) and output_fname is specified, it'll do half pairwise calculation.
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
	option_default_dict.update({('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputted into it and NOT to the db.']})
	option_default_dict.update({('message_size', 1, int):[200, 's', 1, 'How many results one computing node should handle.']})
	
	def __init__(self,  **keywords):
		QC_149_cross_match.__init__(self, **keywords)
	
	def inputNodeHandler(self, param_obj, row_id_pair_list):
		"""
		2008-08-28
			split out of input_node()
		"""
		param_obj.communicator.send("1", param_obj.output_node_rank, 1)	#WATCH: tag is 1, to the output_node.
		free_computing_node, source, tag = param_obj.communicator.receiveString(param_obj.output_node_rank, 2)
		#WATCH: tag is 2, from the output_node
		data_pickle = cPickle.dumps(row_id_pair_list, -1)
		param_obj.communicator.send(data_pickle, int(free_computing_node),0)	#WATCH: int()
		if param_obj.report:
			sys.stderr.write("block %s sent to %s.\n"%(param_obj.counter, free_computing_node))
		param_obj.counter += 1
	
	def input_node(self, param_obj, free_computing_nodes, message_size):
		"""
		2008-09-06
			found a bug in this program
			after the double-for loop, row_id_pair_list can still have some entries in it if the last one's size !=message_size. send it out.
		2008-08-28
		"""
		node_rank = param_obj.communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		twoSNPData = param_obj.twoSNPData
		output_node_rank = param_obj.output_node_rank
		
		param_obj.counter = 0
		row_id_pair_list = []
		if param_obj.QC_method_id==4 and param_obj.output_fname:	#if it's 149 self-cross-match and output goes to file
			for i in range(len(twoSNPData.SNPData1.row_id_ls)):
				for j in range(i, len(twoSNPData.SNPData1.row_id_ls)):	#SNPData2 is same data as SNPData1
					row_id1 = twoSNPData.SNPData1.row_id_ls[i]
					row_id2 = twoSNPData.SNPData1.row_id_ls[j]
					row_id_pair_list.append((row_id1, row_id2))
					if len(row_id_pair_list)==message_size:
						self.inputNodeHandler(param_obj, row_id_pair_list)
						row_id_pair_list = []	#clear the list
		else:
			for row_id1 in twoSNPData.SNPData1.row_id_ls:
				for row_id2 in twoSNPData.SNPData2.row_id_ls:
					row_id_pair_list.append((row_id1, row_id2))
					if len(row_id_pair_list)==message_size:
						self.inputNodeHandler(param_obj, row_id_pair_list)
						row_id_pair_list = []	#clear the list
		if len(row_id_pair_list)>0:	#don't forget the last batch
			self.inputNodeHandler(param_obj, row_id_pair_list)
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			param_obj.communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-09-10
			add source_id to PassingData
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
			qc_cross_match = PassingData(source_id=row_id1[0], strainid=row_id1[1], target_id=target_id, mismatch_rate=mismatch_rate, \
												no_of_mismatches=no_of_mismatches, no_of_non_NA_pairs=no_of_non_NA_pairs)
			result_ls.append(qc_cross_match)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-09-10
			write header in this function. the output columns are directly guessed from the data passed from computing_node.
		2008-08-28
			add functionality to output into file
		2008-08-28
		"""
		#writer, session, commit, QC_method_id, readme = parameter_list
		table_obj_ls = cPickle.loads(data)
		for table_obj in table_obj_ls:
			if param_obj.writer:
				if not param_obj.is_header_written:
					header_row = []
					for column_name in table_obj.__dict__:
						header_row.append(column_name)
					header_row.append('qc_method_id')
					param_obj.writer.writerow(header_row)
					param_obj.is_header_written = True
				row = []
				for column in table_obj.__dict__:	#table_obj has one id that QCCrossMatch table doesn't, source_id
					row.append(getattr(table_obj, column, None))
				row.append(param_obj.QC_method_id)
				param_obj.writer.writerow(row)
			else:
				qc_cross_match = StockDB.QCCrossMatch()
				#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
				#can't save table_obj because it's associated with a different db thread
				for column in qc_cross_match.c.keys():	#it's qc_cross_match, not table_obj because table_obj is not linked to db.
					setattr(qc_cross_match, column, getattr(table_obj, column, None))
				qc_cross_match.qc_method_id = param_obj.QC_method_id
				qc_cross_match.readme = param_obj.readme
				param_obj.session.save(qc_cross_match)
				if param_obj.commit:
					param_obj.session.flush()
	
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
			param_obj = PassingData(communicator=self.communicator, twoSNPData=twoSNPData, output_node_rank=output_node_rank, \
								QC_method_id=self.QC_method_id, output_fname=getattr(self, 'output_fname', None), report=self.report)
			self.input_node(param_obj, free_computing_nodes, self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(twoSNPData=twoSNPData, QC_method_id=self.QC_method_id)
			mw.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			else:
				writer = None
			param_obj = PassingData(writer=writer, session=session, commit=self.commit, QC_method_id=self.QC_method_id, readme=readme, is_header_written=False)
			mw.output_node(free_computing_nodes, param_obj, self.output_node_handler)
			del writer
		mw.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiQC149CrossMatch
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
