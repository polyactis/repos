#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster on all Flowering phenotypes with Flowering times genes
	mpiexec ~/script/variation/src/MpiInterSNPPairAsso.py -g ~/panfs/250k/ft_gene.tsv -i ~/panfs/250k/call_method_17_binary.tsv -n ~/panfs/250k/snps_context_g0_m5000_g0_m5000  -p ~/panfs/250k/phenotype.tsv -o ~/panfs/250k/boolean_snp_pair_ft_gene/ -w 1-7,39-61,80-82
	
	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
	#to test and debug a portion of code (mainly input_node's stuff)
	python ~/script/variation/src/MpiInterSNPPairAsso.py -i ... -b
	
	
	
Description:
	2009-1-22
	Program to do SNP pair (SNPs around genes from gene_id_fname versus all SNPs) association. 
	
	1. Boolean Mode: Each SNP pair is tested by 5 different boolean operations
		(AND, Inhibition, Inhibition, XOR, OR).
	2. Conventional Linear Model interaction detection. y = b + SNP1 + SNP2 + SNP1xSNP2 + e
	
	'gene_id_fname' contains genes from which SNPs are selected to pair with all SNPs. if not given, pairwise between top 10K snps for each phenotype.
	
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
from GeneListRankTest import SnpsContextWrapper, GeneListRankTest

from MpiIntraGeneSNPPairAsso import MpiIntraGeneSNPPairAsso, num, bool_type2merge_oper
from PlotGroupOfSNPs import PlotGroupOfSNPs
from DrawSNPRegion import DrawSNPRegion
import Stock_250kDB

class MpiInterSNPPairAsso(MpiIntraGeneSNPPairAsso):
	__doc__ = __doc__
	option_default_dict = MpiIntraGeneSNPPairAsso.option_default_dict.copy()
	option_default_dict.update({('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ]})
	option_default_dict.update({('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ]})
	option_default_dict.update({('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ]})
	option_default_dict.update({('schema', 0, ): [None, 'k', 1, 'database schema name', ]})
	option_default_dict.update({('db_user', 1, ): [None, 'u', 1, 'database username', ]})
	option_default_dict.update({('db_passwd', 1, ): [None, 'a', 1, 'database password', ]})
	#option_default_dict.update({('gene_id_fname', 1, ): [None, 'g', 1, 'A file with gene id on each line. Try to detect interactions of SNPs around these genes vs all SNPs.']})
	#option_default_dict.update({})
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		if self.phenotype_method_id_ls is not None:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	def getTopNumberSNPs(self, call_method_id, phenotype_method_id, analysis_method_id=1, results_directory=None, min_MAC=7,\
						no_of_top_snps=10000):
		"""
		2009-3-18
			if data_obj from genome_wide_result is None, ignore it
		2009-2-16
			get a certain number of top SNPs from a result (according to call_method_id, phenotype_method_id, analysis_method_id)
		"""
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).first()
		result_id = rm.id
		param_data = PassingData(min_MAC=min_MAC)
		genome_wide_result = GeneListRankTest.getResultMethodContent(rm, results_directory=results_directory, min_MAF=0., \
																	pdata=param_data)
		
		chr_pos_ls = []
		for i in range(no_of_top_snps):
			data_obj = genome_wide_result.get_data_obj_at_given_rank(i+1)
			if data_obj is not None:
				chr_pos_ls.append((data_obj.chromosome, data_obj.position))
		chr_pos_ls.sort()
		return chr_pos_ls
	
	def generate_params(self, gene_id_fname, pdata, block_size=1000, **keywords):
		"""
		2009-2-18
			if gene_id_fname is given and is a file:
				yield (gene1_id, snp_start_index, snp_stop_index)
			else:
				yield (phenotype_index, snp_start_index1, snp_stop_index1, snp_start_index2, snp_stop_index2)
		2009-2-12
			use yield to become a generator
			called by inputNodePrepare()
		2009-1-22
		"""
		no_of_phenotypes = len(pdata.phenotype_index_ls)
		start_index = 0	#for each computing node: the index of gene >= start_index
		#no_of_genes = len(pdata.gene_id2snps_id_ls)
		no_of_tests_per_node = 0
		
		if gene_id_fname and os.path.isfile(gene_id_fname):
			reader = csv.reader(open(gene_id_fname), delimiter=figureOutDelimiter(gene_id_fname))
			gene_id_ls = []
			for row in reader:
				gene_id = int(row[0])
				gene_id_ls.append(gene_id)
			del reader
			no_of_genes = len(gene_id_ls)
			no_of_total_snps = len(pdata.snp_info.chr_pos_ls)
			for i in range(no_of_genes):
				gene1_id = gene_id_ls[i]
				n1 = len(pdata.gene_id2snps_id_ls[gene1_id])	#no_of_snps_of_this_gene
				snp_start_index = 0
				while snp_start_index < no_of_total_snps:
					no_of_snps_to_consider = block_size/(n1*no_of_phenotypes)
					snp_stop_index = snp_start_index+no_of_snps_to_consider
					if snp_stop_index > no_of_total_snps:
						snp_stop_index = no_of_total_snps
					yield (gene1_id, snp_start_index, snp_stop_index)
					snp_start_index += no_of_snps_to_consider
		else:
			#no gene_id_fname. pairwise among all SNPs
			no_of_snps_to_consider = int(math.sqrt(block_size))
			no_of_total_snps = len(pdata.snp_info.chr_pos_ls)
			for phenotype_index in pdata.phenotype_index_ls:
				for snp_start_index1 in range(0, no_of_total_snps, no_of_snps_to_consider):
					snp_stop_index1 = min(no_of_total_snps, snp_start_index1+no_of_snps_to_consider)
					for snp_start_index2 in range(snp_start_index1, no_of_total_snps, no_of_snps_to_consider):
						snp_stop_index2 = min(no_of_total_snps, snp_start_index2+no_of_snps_to_consider)
						yield (phenotype_index, snp_start_index1, snp_stop_index1, snp_start_index2, snp_stop_index2)
	
	def computeGeneVsSNP(self, gene_id, snp_start_index, snp_stop_index, param_obj):
		"""
		2009-2-16
			refactored out of computing_node_handler()
			
			handle the pairwise computation between SNPs from one gene and SNPs defined by (snp_start_index, snp_stop_index)
		"""
		result_ls = []
		snpData = param_obj.snpData
		no_of_phenotypes = param_obj.phenData.data_matrix.shape[1]		
		snps_id_ls1 = param_obj.gene_id2snps_id_ls[gene_id]
		snp_chr_pos_ls = param_obj.snp_info.chr_pos_ls[snp_start_index:snp_stop_index]
		for snp1_id in snps_id_ls1:
			snp1_index = snpData.col_id2col_index.get(snp1_id)
			if snp1_index is None:	#snp1_id not in input matrix
				continue
			for chr_pos in snp_chr_pos_ls:
				snp2_id = '%s_%s'%(chr_pos[0], chr_pos[1])
				snp2_index = snpData.col_id2col_index.get(snp2_id)
				if snp2_index is None or snp2_index==snp1_index:	#snp2_id not in input matrix, or two SNPs are same
					continue
				for phenotype_index in param_obj.phenotype_index_ls:
					if phenotype_index>=no_of_phenotypes:	#skip if out of bound
						continue
					phenotype_ls = param_obj.phenData.data_matrix[:, phenotype_index]
					
					pdata_ls = self.test_type2test[param_obj.test_type](snp1_id, gene_id, snpData.data_matrix[:,snp1_index], \
																snp2_id, None, snpData.data_matrix[:,snp2_index], \
																phenotype_index, phenotype_ls, \
																min_data_point=param_obj.min_data_point)
					result_ls += pdata_ls
		return result_ls
	
	def computeSNPVsSNP(self, phenotype_index, snp_start_index1, snp_stop_index1, snp_start_index2, snp_stop_index2, param_obj):
		"""
		2009-2-16
			in contrast to computeGeneVsSNP()
			this is pairwise between [snp_start_index1:snp_stop_index1] and [snp_start_index2:snp_stop_index2].
		"""
		result_ls = []
		no_of_phenotypes = param_obj.phenData.data_matrix.shape[1]
		if phenotype_index>=no_of_phenotypes:	#skip if out of bound
			return result_ls
		phenotype_ls = param_obj.phenData.data_matrix[:, phenotype_index]
		snpData = param_obj.snpData
		phenotype_method_id = param_obj.phenData.phenotype_method_id_ls[phenotype_index]
		which_chr_pos_ls = param_obj.snp_info.phenotype_method_id2chr_pos_ls[phenotype_method_id]
		snp_chr_pos_ls1 = which_chr_pos_ls[snp_start_index1:snp_stop_index1]
		snp_chr_pos_ls2 = which_chr_pos_ls[snp_start_index2:snp_stop_index2]
		is_two_chr_pos_ls_same = 0
		if snp_start_index1==snp_start_index2 and snp_stop_index1==snp_stop_index2:
			is_two_chr_pos_ls_same = 1
		
		for i in range(len(snp_chr_pos_ls1)):
			chr_pos1 = snp_chr_pos_ls1[i]
			snp1_id = '%s_%s'%(chr_pos1[0], chr_pos1[1])
			snp1_index = snpData.col_id2col_index.get(snp1_id)
			if snp1_index is None:	#snp1_id not in input matrix
				continue
			if is_two_chr_pos_ls_same:	#to avoid redundancy
				chr_pos_ls = snp_chr_pos_ls2[i+1:]
			else:
				chr_pos_ls = snp_chr_pos_ls2
			for chr_pos in chr_pos_ls:
				snp2_id = '%s_%s'%(chr_pos[0], chr_pos[1])
				snp2_index = snpData.col_id2col_index.get(snp2_id)
				if snp2_index is None or snp2_index==snp1_index:	#snp2_id not in input matrix, or two SNPs are same
					continue
				pdata_ls = self.test_type2test[param_obj.test_type](snp1_id, None, snpData.data_matrix[:,snp1_index], \
																snp2_id, None, snpData.data_matrix[:,snp2_index], \
																phenotype_index, phenotype_ls, \
																min_data_point=param_obj.min_data_point)
				result_ls += pdata_ls
		return result_ls
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2009-2-16
			call computeGeneVsSNP() or computeSNPVsSNP() according to data passed from input node
		2009-1-22 modified from MpiInterGeneSNPPairAsso.computing_node_handler()
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		
		if len(data)==3:
			gene_id, snp_start_index, snp_stop_index = data
			result_ls = self.computeGeneVsSNP(gene_id, snp_start_index, snp_stop_index, param_obj)
		elif len(data)==5:
			phenotype_index, snp_start_index1, snp_stop_index1, snp_start_index2, snp_stop_index2 = data
			result_ls = self.computeSNPVsSNP(phenotype_index, snp_start_index1, snp_stop_index1, snp_start_index2, snp_stop_index2, param_obj)
		else:
			result_ls = []
		
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def run(self):
		"""
		2009-1-22
		"""
		if self.debug:
			#for one-node testing purpose
			import pdb
			pdb.set_trace()
			sys.exit()
		
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			#2009-1-22
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				  				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			
			#snp_info = DrawSNPRegion.getSNPInfo(db)	#2009-2-12 turn this off
			
			#2009-2-16 get top snps for each phenotype_method_id.
			snp_info = PassingData(phenotype_method_id2chr_pos_ls={})
			for phenotype_method_id in self.phenotype_method_id_ls:
				chr_pos_ls = self.getTopNumberSNPs(self.call_method_id, phenotype_method_id, analysis_method_id=1, \
							results_directory=self.results_directory, min_MAC=7, no_of_top_snps=self.no_of_top_snps)
				snp_info.phenotype_method_id2chr_pos_ls[phenotype_method_id] = chr_pos_ls
				if not hasattr(snp_info, 'chr_pos_ls'):	#this is from a random phenotype (1st one). chr_pos_ls is for generate_params(). The actual SNPs in it dont' matter.
					snp_info.chr_pos_ls = chr_pos_ls
			
			dstruc = self.inputNodePrepare(snp_info)
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
			#self.input_node(param_obj, free_computing_nodes, self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snpData=snpData, gene_id_ls=other_data.gene_id_ls, \
												gene_id2snps_id_ls=other_data.gene_id2snps_id_ls, phenData=other_data.phenData,
												phenotype_index_ls=self.phenotype_index_ls, min_data_point=self.min_data_point,\
												snp_info=other_data.snp_info,\
												test_type=self.test_type)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			self.general_output_node(self.output_dir, self.phenotype_index_ls, phenotype_label_ls, free_computing_nodes)
		self.synchronize()	#to avoid some node early exits
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiInterSNPPairAsso
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()