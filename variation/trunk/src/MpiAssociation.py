#!/usr/bin/env mpipython
"""

Examples:
	mpiexec ~/script/variation/src/MpiAssociation.py -i ~/panfs/250k/dataset/call_method_17_test.tsv -p ~/panfs/250k/phenotype.tsv -o ~/panfs/250k/association_results/call_method_17_test_y3.tsv -y3 -w 187-190,210-221

Description:
	MPI version of Association.py, doesn't support test_type==5 or test_type==6.
	
	class to do association test on SNP data. option 'test_type' decides which test to run.
	
	Input genotype file format is Strain X SNP format (Yu's format, Output by DB_250k2data.py Or Output250KSNPs.py + ConvertBjarniSNPFormat2Yu.py).
	Input phenotype file format is Strain X phenotype format (Output by OutputPhenotype.py). 
	
	It requires a minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis or linear model test.
	
	For all methods, it will automatically match strains in two files.
		NO worry for missing/extra data in either input file.
	
	All methods iterate through phenotypes given by '-w' except that Method "5: LM two phenotypes with PCs" takes two phenotypes from '-w'.
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import csv, numpy, traceback
from pymodule import read_data, ProcessOptions, PassingData, SNPData, getListOutOfStr

import rpy, cPickle
from sets import Set
from Association import Association
from pymodule.MPIwrapper import MPIwrapper
from Scientific import MPI

class MpiAssociation(Association, MPIwrapper):
	__doc__ = __doc__
	# 2009-5-15 add debug here because self.debug set in __init__() is not visible in classmethod thru cls.debug
	debug = 0
	option_default_dict = Association.option_default_dict.copy()
	def __init__(self, **keywords):
		"""
		2009-3-20
		"""
		Association.__init__(self, **keywords)
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2009-3-20
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		
		initData = param_obj.initData
		min_data_point = param_obj.min_data_point
		which_PC_index_ls = param_obj.which_PC_index_ls
		which_phenotype = data
		environment_matrix = None
		gene_environ_interaction = False
		results = self.run_whole_matrix[param_obj.test_type](initData.snpData.data_matrix, initData.phenData.data_matrix[:, which_phenotype], \
													min_data_point, PC_matrix=initData.PC_matrix, \
													which_PC_index_ls=which_PC_index_ls, environment_matrix=environment_matrix, \
													gene_environ_interaction=gene_environ_interaction)
		result_ls.append([which_phenotype, results])
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2009-3-20
			doesn't support self.test_type==5 or self.test_type==6
		"""
		result_ls = cPickle.loads(data)
		which_phenotype, results = result_ls[0]
		phenotype_name = param_obj.phenotypeColIDLs[which_phenotype]
		phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
		output_fname='%s_pheno_%s.tsv'%(os.path.splitext(param_obj.output_fname)[0], phenotype_name)	#make up a new name corresponding to this phenotype
		
		output_results_func = self.output_results.get(param_obj.test_type)
		if output_results_func is None:
			output_results_func = self.output_lm_results
		output_results_func(results, param_obj.snpDataColIDLs, output_fname, param_obj.minus_log_pvalue)
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			initData = self.readInData(self.phenotype_fname, self.input_fname, self.eigen_vector_fname, self.phenotype_method_id_ls, self.test_type, self.report)
			initDataPickle = cPickle.dumps(initData, -1)
			which_phenotype_ls = initData.which_phenotype_ls
			
			#send the output node some data
			outputNodeData = PassingData(phenotypeColIDLs=initData.phenData.col_id_ls, snpDataColIDLs=initData.snpData.col_id_ls)
			outputNodeDataPickle = cPickle.dumps(outputNodeData, -1)
			self.communicator.send(outputNodeDataPickle, output_node_rank, 0)
			del outputNodeDataPickle
			
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(initDataPickle, node, 0)
				sys.stderr.write(".\n")
			del initDataPickle
			
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			initData =  cPickle.loads(data)
			del data
		else:
			data, source, tag = self.communicator.receiveString(0, 0)
			outputNodeData = cPickle.loads(data)
			phenotypeColIDLs = outputNodeData.phenotypeColIDLs
			snpDataColIDLs = outputNodeData.snpDataColIDLs
			
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = which_phenotype_ls)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(initData=initData, min_data_point=self.min_data_point, \
												which_PC_index_ls=self.which_PC_index_ls, \
												test_type=self.test_type)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			param_obj = PassingData(output_fname=self.output_fname, test_type=self.test_type, minus_log_pvalue=self.minus_log_pvalue,\
								phenotypeColIDLs=phenotypeColIDLs, snpDataColIDLs=snpDataColIDLs)
			self.output_node(free_computing_nodes, param_obj, self.output_node_handler)
		self.synchronize()	#to avoid some node early exits
	
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiAssociation
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()