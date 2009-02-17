#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster on all Flowering phenotypes with Flowering times genes
	mpiexec ~/script/variation/src/MpiInterGeneSNPPairAsso.py -g ~/panfs/250k/ft_gene.tsv -i ~/panfs/250k/call_method_17_binary.tsv -n ~/panfs/250k/snps_context_g0_m5000_g0_m5000  -p ~/panfs/250k/phenotype.tsv -o ~/panfs/250k/boolean_snp_pair_ft_gene/ -w 1-7,39-61,80-82
	
	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
	#to test and debug a portion of code (mainly input_node's stuff)
	python ~/script/variation/src/MpiInterGeneSNPPairAsso.py -i ... -b
	
	
	
Description:
	2008-11-25
	Program to do inter-gene SNP pair association (vs MpiIntraGeneSNPPairAsso.py).
	
	2009-2-8 'gene_id_fname' contains genes to pair with each other. at least 1 column, which is gene id.
	
	Input is StrainXSNP matrix with SNP 0 or 1 or <0(=NA). Be careful in choosing the block_size which determines the lower bound of #tests each computing node handles.
		The input node tells the computing node which genes it should work on. A computing node might have to deal with lots of tests regardless of the block_size, if one gene has lots of SNPs in it with #tests far exceeding the lower bound.
		The amount each computing node handles is multiplied by the number of phenotypes. Different genes harbor different no of SNPs.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import cPickle
from pymodule import PassingData, importNumericArray, SNPData, read_data, getListOutOfStr, figureOutDelimiter
from pymodule.SNP import NA_set
from sets import Set
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper
from Kruskal_Wallis import Kruskal_Wallis
from GeneListRankTest import SnpsContextWrapper

from MpiIntraGeneSNPPairAsso import MpiIntraGeneSNPPairAsso, num, bool_type2merge_oper
from PlotGroupOfSNPs import PlotGroupOfSNPs

class MpiInterGeneSNPPairAsso(MpiIntraGeneSNPPairAsso):
	__doc__ = __doc__
	option_default_dict = MpiIntraGeneSNPPairAsso.option_default_dict.copy()
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		if self.phenotype_method_id_ls is not None:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	def generate_params(self, gene_id_fname, pdata, block_size=1000, **keywords):
		"""
		2009-2-12
			use yield to become a generator
		2008-11-25
			read gene ids from gene_id_fname and generate pairs among them
			each node handle a number of pairs, depending on many SNP pairs it incurs
		"""
		#sys.stderr.write("Generating parameters ...")
		#params_ls = []
		no_of_phenotypes = len(pdata.phenotype_index_ls)
		start_index = 0	#for each computing node: the index of gene >= start_index
		#no_of_genes = len(pdata.gene_id2snps_id_ls)
		no_of_tests_per_node = 0
		
		reader = csv.reader(open(gene_id_fname), delimiter=figureOutDelimiter(gene_id_fname))
		gene_id_ls = []
		for row in reader:
			gene_id = int(row[0])
			gene_id_ls.append(gene_id)
		del reader
		
		no_of_genes = len(gene_id_ls)
		gene_id_pairs_for_each_node = []
		for i in range(no_of_genes):
			gene1_id = gene_id_ls[i]
			n1 = len(pdata.gene_id2snps_id_ls[gene1_id])	#no_of_snps_of_this_gene
			#for gene2_id in pdata.gene_id2snps_id_ls:	#2009-2-8 another setting: gene_id_fname vs all genes
			for j in range(i+1, no_of_genes):
				gene2_id = gene_id_ls[j]
				n2 = len(pdata.gene_id2snps_id_ls[gene2_id])
				est_no_of_tests = (n1*n2)*no_of_phenotypes	#this is the upper bound for the number of tests for each gene on a computing node. data missing would make the number smaller.
				no_of_tests_per_node += est_no_of_tests
				gene_id_pairs_for_each_node.append((gene1_id, gene2_id))
				if no_of_tests_per_node>=block_size:
					yield gene_id_pairs_for_each_node
					#reset gene_id_pairs_for_each_node
					gene_id_pairs_for_each_node = []
					no_of_tests_per_node = 0	#reset this to 0
				
		#pick up the last gene_id_pairs_for_each_node if it's not empty
		if gene_id_pairs_for_each_node:
			yield gene_id_pairs_for_each_node
			#no need to cleanup because this is the end of loop
		#sys.stderr.write("%s params. Done.\n"%(len(params_ls)))
		#return params_ls
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2008-11-25
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		snpData = param_obj.snpData
		no_of_phenotypes = param_obj.phenData.data_matrix.shape[1]
		for gene_id_pair_ls in data:
			for gene_id_pair in gene_id_pair_ls:
				gene1_id, gene2_id = gene_id_pair
				snps_id_ls1 = param_obj.gene_id2snps_id_ls[gene1_id]
				snps_id_ls2 = param_obj.gene_id2snps_id_ls[gene2_id]
				for snp1_id in snps_id_ls1:
					snp1_index = snpData.col_id2col_index.get(snp1_id)
					if snp1_index is None:	#snp1_id not in input matrix
						continue
					for snp2_id in snps_id_ls2:
						snp2_index = snpData.col_id2col_index.get(snp2_id)
						if snp2_index is None or snp2_index==snp1_index:	#snp2_id not in input matrix or two SNPs are same
							continue
						for phenotype_index in param_obj.phenotype_index_ls:
							if phenotype_index>=no_of_phenotypes:	#skip if out of bound
								continue
							phenotype_ls = param_obj.phenData.data_matrix[:, phenotype_index]
							
							pdata_ls = self.test_type2test[param_obj.test_type](snp1_id, gene1_id, snpData.data_matrix[:,snp1_index], \
																		snp2_id, gene2_id, snpData.data_matrix[:,snp2_index], \
																		phenotype_index, phenotype_ls, \
																		min_data_point=param_obj.min_data_point)
							result_ls += pdata_ls
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def run(self):
		"""
		2008-11-25
		"""
		if self.debug:
			#for one-node testing purpose
			import pdb
			pdb.set_trace()
			sys.exit(2)
		
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			dstruc = self.inputNodePrepare()
			params_ls = dstruc.params_ls
			
			#send the output node the phenotype_label_ls
			self.communicator.send(dstruc.output_node_data_pickle, output_node_rank, 0)
			del dstruc.output_node_data_pickle
			
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(dstruc.snpData_pickle, node, 0)
				self.communicator.send(dstruc.other_data_pickle, node, 0)
				sys.stderr.write(".\n")
			del dstruc
		elif node_rank in free_computing_node_set:
			data, source, tag = self.communicator.receiveString(0, 0)
			snpData =  cPickle.loads(data)
			del data
			data, source, tag = self.communicator.receiveString(0, 0)
			other_data = cPickle.loads(data)
			del data
			self.phenotype_index_ls = other_data.phenotype_index_ls
		else:
			data, source, tag = self.communicator.receiveString(0, 0)
			output_node_data_pickle = cPickle.loads(data)
			phenotype_label_ls = output_node_data_pickle.phenotype_label_ls
			self.phenotype_index_ls = output_node_data_pickle.phenotype_index_ls
			
		self.synchronize()
		if node_rank == 0:
			param_obj = PassingData(params_ls=params_ls, output_node_rank=output_node_rank, report=self.report, counter=0)
			self.inputNode(param_obj, free_computing_nodes, param_generator = params_ls)
			#self.input_node(param_obj, free_computing_nodes, input_handler=self.input_fetch_handler, message_size=1)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snpData=snpData, gene_id_ls=other_data.gene_id_ls, \
												gene_id2snps_id_ls=other_data.gene_id2snps_id_ls, phenData=other_data.phenData,
												phenotype_index_ls=self.phenotype_index_ls, min_data_point=self.min_data_point,
												test_type=self.test_type)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			self.general_output_node(self.output_dir, self.phenotype_index_ls, phenotype_label_ls, free_computing_nodes)
		self.synchronize()	#to avoid some node early exits
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiInterGeneSNPPairAsso
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()