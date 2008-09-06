#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiLD.py -i /tmp/SNPmatrix.tsv -o /tmp/LD.tsv -s 100

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
							('block_size', 1, int):[1000, 's', 1, 'square of it is the number (or half of that) of LDs each computing node handles. Imagine a SNPXSNP LD matrix. block_size is the dimension of the block each node handles.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def generate_params(self, no_of_snps, block_size=1000):
		"""
		2008-09-06
			cut the SNPXSNP LD matrix into blocks and put the two dimensions into params_ls
		"""
		sys.stderr.write("Generating parameters ...")
		params_ls = []
		no_of_blocks = int(no_of_snps/block_size)+1
		for i in range(no_of_blocks+1):
			min_index1 = i*block_size
			if min_index1>=no_of_snps:	#out of bound ignore
				continue
			stop1 = min((i+1)*block_size, no_of_snps)
			
			for j in range(i, no_of_blocks+1):	#the 2nd index is equal or bigger than the 1st index
				min_index2 = j*block_size
				if min_index2>=no_of_snps:	#out of bound
					continue
				stop2 = min((j+1)*block_size, no_of_snps)
				params_ls.append([(min_index1,stop1), (min_index2,stop2)])
		sys.stderr.write("Done.\n")
		return params_ls
				
	def input_handler(self, param_obj, message_size=1, report=0):
		"""
		2008-09-06
		"""
		if param_obj.report:
			sys.stderr.write("Fetching stuff...\n")
		params_ls = param_obj.params_ls
		data_to_return = []
		for i in range(message_size):
			if len(params_ls)>0:
				one_parameter = params_ls.pop(0)
				data_to_return.append(one_parameter)
				param_obj.counter += 1
			else:
				break
		if param_obj.report:
			sys.stderr.write("Fetching done at counter=%s.\n"%(param_obj.counter))
		return data_to_return
	
	def _input_node(self, param_obj, free_computing_nodes, message_size):
		"""
		2008-09-06
			deprecated
			use generate_params() and call MPIwrapper's input_node()
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
		
		if id_pair_list:	#don't forget the last batch
			self.inputNodeHandler(param_obj, id_pair_list)
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			self.communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_node_handler(self, communicator, data, computing_parameter_obj):
		"""
		2008-09-06
			data from input_node is changed
		2008-09-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		snpData = computing_parameter_obj.snpData
		for col1_range, col2_range in data:
			min_index1, stop1 = col1_range
			min_index2, stop2 = col2_range
			for i in range(min_index1, stop1):
				for j in range(max(i+1, min_index2), stop2):	#the lower bound of j is the bigger one of i+1 and min_index2
					col_id1 = snpData.col_id_ls[i]
					col_id2 = snpData.col_id_ls[j]
					LD_data = snpData.calLD(col_id1, col_id2)
					if LD_data is not None:
						result_ls.append(LD_data)
						
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-09-05
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
			params_ls = self.generate_params(len(snpData.col_id_ls), self.block_size)
			del snpData
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			snpData =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(params_ls=params_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.input_node(param_obj, free_computing_nodes, input_handler=self.input_handler, message_size=1)
			#self.input_node(param_obj, free_computing_nodes, self.message_size)
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