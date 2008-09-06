#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiLD.py -i /tmp/SNPmatrix.tsv -o /tmp/LD.tsv -s 1000

	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
Description:
	2008-09-05 a MPI program to calculate LD of an input SNP matrix. The input format is StrainXSNP matrix, either tab or comma-delimited.
	Apart from 1st row as header, 1st column as strain id, 2nd column is affiliating information for each strain.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, SNPData, read_data
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper

num = importNumericArray()

class MpiLD(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('input_fname',1, ): [None, 'i', 1, 'a file containing StrainXSNP matrix.'],\
							("output_fname", 1, ): [None, 'o', 1, 'Filename to store data matrix'],\
							('message_size', 1, int):[1000, 's', 1, 'How many results one computing node should handle.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def input_node(self, param_obj, free_computing_nodes, message_size):
		"""
		2008-09-05
			similar to MpiQC149CrossMatch.py's input_node()
		"""
		node_rank = self.communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		snpData = param_obj.snpData
		output_node_rank = param_obj.output_node_rank
		
		param_obj.counter = 0
		id_pair_list = []
		no_of_cols = len(snpData.col_id_ls)
		for i in range(no_of_cols):
			for j in range(i+1, no_of_cols):
				col_id1 = snpData.col_id_ls[i]
				col_id2 = snpData.col_id_ls[j]
				id_pair_list.append((col_id1, col_id2))
				if len(id_pair_list)==message_size:
					self.inputNodeHandler(param_obj, id_pair_list)
					id_pair_list = []	#clear the list
		
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			self.communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-08-28
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		snpData = computing_parameter_obj.snpData
		for col_id1, col_id2 in data:
			LD_data = snpData.calLD(col_id1, col_id2)
			if LD_data is not None:
				result_ls.append(LD_data)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-08-28
			add functionality to output into file
		2008-08-28
		"""
		writer = param_obj.writer
		result_ls = cPickle.loads(data)
		for LD_data in result_ls:
			if writer:
				row = []
				if not param_obj.is_header_written:
					header_row = []
					for column_name in LD_data.__dict__:
						if column_name=='allele_freq':
							header_row.append('allele1_freq')
							header_row.append('allele2_freq')
						elif column_name=='snp_pair_ls':
							header_row.append('snp1')
							header_row.append('snp2')
						else:
							header_row.append(column_name)
					writer.writerow(header_row)
					param_obj.is_header_written = True
				
				for key,value in LD_data.__dict__.iteritems():
					if key=='allele_freq':
						row.append(value[0])
						row.append(value[1])
					elif key=='snp_pair_ls':
						row.append(value[0])
						row.append(value[1])
					else:
						row.append(value)
				writer.writerow(row)
	
	def run(self):
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		"""
		if node_rank!=output_node_rank:
			header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix)	#category_list is not used to facilitate row-id matching
		"""
		if node_rank == 0:
			header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix)	#category_list is not used to facilitate row-id matching
			snpData_pickle = cPickle.dumps(snpData, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(snpData_pickle, node, 0)
				sys.stderr.write(".\n")
			del snpData_pickle
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			snpData =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(snpData=snpData, output_node_rank=output_node_rank, report=self.report)
			self.input_node(param_obj, free_computing_nodes, self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snpData=snpData)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				#header_row = ['snp1_id', 'snp2_id', 'r2', 'D', "D'", "no_of_pairs"]
				#writer.writerow(header_row)
			else:
				writer = None
			param_obj = PassingData(writer=writer, is_header_written=False)
			self.output_node(free_computing_nodes, param_obj, self.output_node_handler)
			del writer
		self.synchronize()	#to avoid some node early exits
	

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiLD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()