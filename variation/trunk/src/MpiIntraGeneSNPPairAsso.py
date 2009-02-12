#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiIntraGeneSNPPairAsso.py -i /mnt/nfs/250k/call_method_17_binary.tsv -n /mnt/nfs/250k/snps_context_wrapper_g1_m50000  -p /mnt/nfs/250k/phenotypes_transformed_publishable_v2.tsv -o /mnt/nfs/250k/boolean_snp_pair/

	#test parallel run on desktop
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/...
	
	#to test and debug a portion of code (mainly input_node's stuff)
	python ~/script/variation/src/MpiIntraGeneSNPPairAsso.py -i ... -b
	
Description:
	2008-09-07
	Program to do intra-gene SNP pair association. two models:
	
		1. Each SNP pair is tested by 5 different boolean operations (AND, Inhibition, Inhibition, XOR, OR).
	
		2. Conventional Linear Model interaction detection. y = b + SNP1 + SNP2 + SNP1xSNP2 + e
	
	Input is StrainXSNP matrix with SNP 0 or 1 or <0(=NA). Be careful in choosing the block_size which determines the lower bound of #tests each computing node handles.
		The input node tells the computing node which genes it should work on. A computing node might have to deal with lots of tests regardless of the block_size, if one gene has lots of SNPs in it with #tests far exceeding the lower bound.
		The amount each computing node handles is multiplied by the number of phenotypes. Different genes harbor different no of SNPs.
	
	"gene_id_fname" is given to restrict the analysis only those genes. at least 1 column with NCBI gene id.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, numpy
import cPickle
from pymodule import PassingData, importNumericArray, SNPData, read_data, getListOutOfStr, figureOutDelimiter
from pymodule.SNP import NA_set
from sets import Set
from Scientific import MPI
from pymodule.MPIwrapper import MPIwrapper
from Kruskal_Wallis import Kruskal_Wallis
from Association import Association
from GeneListRankTest import SnpsContextWrapper
from PlotGroupOfSNPs import PlotGroupOfSNPs

num = importNumericArray()

bool_type2merge_oper = {1: num.array([[0,0],[0,1]], num.int8),	#AND
					2: num.array([[0,0],[1,0]], num.int8),	#Inhibition
					4: num.array([[0,1],[0,0]], num.int8),	#Inhibition
					6: num.array([[0,1],[1,0]], num.int8),	#XOR
					7: num.array([[0,1],[1,1]], num.int8)}	#OR

class TwoSNPAssociation(object):
	"""
	2009-2-10
		class containing methods to do two-SNP association
		
		refactored out of MpiIntraGeneSNPPairAsso
	"""
	def mergeTwoGenotypeLs(cls, genotype_ls1, genotype_ls2, merge_oper_matrix):
		"""
		2008-09-07
			
		"""
		genotype_ls = []
		no_of_genotypes = len(genotype_ls1)
		for i in range(no_of_genotypes):
			genotype1 = genotype_ls1[i]
			genotype2 = genotype_ls2[i]
			if genotype1>=0 and genotype2>=0:	#assume it's binary. negative value is regarded as NA.
				new_genotype = merge_oper_matrix[genotype1][genotype2] + 1	#+1 to avoid 0 being treated as NA in Kruskal_Wallis._kruskal_wallis
			else:
				new_genotype = 0	#0=NA for Kruskal_Wallis._kruskal_wallis
			genotype_ls.append(new_genotype)
		return genotype_ls
	mergeTwoGenotypeLs = classmethod(mergeTwoGenotypeLs)
	
	def KWOnBooleanSNP(cls, snp1_id, gene1_id, genotype_ls1, snp2_id, gene2_id, genotype_ls2, phenotype_index, phenotype_ls, \
					min_data_point=3):
		"""
		2009-2-8
			refactor out of computing_node_handler()
			kruskal wallis on a boolean-merged SNP
		"""
		return_ls = []
		for bool_type in bool_type2merge_oper:
			merge_oper_matrix = bool_type2merge_oper[bool_type]
			genotype_ls = cls.mergeTwoGenotypeLs(genotype_ls1, genotype_ls2, merge_oper_matrix)
			pdata = Kruskal_Wallis._kruskal_wallis(genotype_ls, phenotype_ls, min_data_point)
			if pdata:
				pdata = PassingData(snp1_id=snp1_id, gene1_id=gene1_id, snp2_id=snp2_id, gene2_id=gene2_id, pvalue=pdata.pvalue,\
								count1=pdata.count_ls[0], count2=pdata.count_ls[1], bool_type=bool_type, phenotype_index=phenotype_index)
				return_ls.append(pdata)
		return return_ls
	KWOnBooleanSNP  = classmethod(KWOnBooleanSNP)
	
	def IntLMOnTwoSNPs(cls, snp1_id, gene1_id, genotype_ls1, snp2_id, gene2_id, genotype_ls2, phenotype_index, phenotype_ls, \
					min_data_point=3):
		"""
		2009-2-8
			interaction detection linear model
			y = b + SNP1xSNP2 + SNP1 + SNP2 + e
			interaction is the 1st term. therefore the pvalue directly returned is also for this term.
		"""
		return_ls = []
		genotype_ls1 = numpy.resize(genotype_ls1, [len(genotype_ls1),1])
		genotype_ls2 = numpy.resize(genotype_ls2, [len(genotype_ls2),1])
		snp_int_matrix = genotype_ls1*genotype_ls2
		genotype_ls = numpy.hstack((snp_int_matrix, genotype_ls1, genotype_ls2))	#interaction variable is the 1st position
		
		pdata = Association.linear_model(genotype_ls, phenotype_ls, min_data_point, snp_index=snp1_id+snp2_id)
		
		if pdata:
			pdata = PassingData(snp1_id=snp1_id, gene1_id=gene1_id, snp2_id=snp2_id, gene2_id=gene2_id, pvalue=pdata.pvalue,\
							count1=pdata.count_ls[0], count2=pdata.count_ls[1], bool_type=None, phenotype_index=phenotype_index,\
							var_perc=pdata.var_perc, coeff_list=pdata.coeff_list, coeff_p_value_list=pdata.coeff_p_value_list)
			return_ls.append(pdata)
		return return_ls
	
	IntLMOnTwoSNPs = classmethod(IntLMOnTwoSNPs)

class MpiIntraGeneSNPPairAsso(MPIwrapper):
	__doc__ = __doc__
	option_default_dict = {('snps_context_fname',1, ): [None, 'n', 1, 'a file containing a pickled snps_context_wrapper. outputted by GeneListRankTest.constructDataStruc()'],\
							('input_fname',1, ): [None, 'i', 1, 'a file containing StrainXSNP matrix. must be binary matrix.'],\
							("output_dir", 1, ): [None, 'o', 1, 'directory to store output association results. one phenotype, one file.'],\
							('phenotype_fname', 1, ): [None, 'p', 1, 'phenotype file. Same format as input_fname. but replace the data matrix with phenotype data.', ],\
							('min_data_point', 1, int): [3, 'm', 1, 'minimum number of ecotypes for either alleles of a single SNP to be eligible for kruskal wallis test'],\
							('phenotype_method_id_ls', 0, ): [None, 'w', 1, 'which phenotypes to work on. a comma-separated list phenotype_method ids in the phenotype file. Check db Table phenotype_method. Default is to take all.',],\
							('block_size', 1, int):[10000, 's', 1, '~Maximum number of tests each computing node is gonna handle. The computing node loops over all phenotypes, test all pairwise SNPs within a gene.'],\
							('test_type', 1, int): [1, 'y', 1, 'Which type of test to do. 1:Kruskal_Wallis on a SNP boolean-merged from two SNPs, 2:linear model(y=b + SNP1 + SNP2 + SNP1xSNP2 + e)'],\
							('gene_id_fname', 0, ): [None, 'g', 1, 'A file with gene id on each line.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-11-25
			change option phenotype_index_ls to phenotype_method_id_ls to pick phenotypes. phenotype_method_id is more intuitive.
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if self.phenotype_method_id_ls is not None:
			self.phenotype_method_id_ls = getListOutOfStr(self.phenotype_method_id_ls, data_type=int)
	
	test_type2test = {1: TwoSNPAssociation.KWOnBooleanSNP,
					2: TwoSNPAssociation.IntLMOnTwoSNPs}
	def get_gene_id2snps_id_ls(self, snps_context_wrapper):
		"""
		2009-2-8
			SNPs that are within genes are not exclusive to those genes anymore.
			gene-SNP association is totally based on the snps_context_wrapper
		2008-09-07
		"""
		sys.stderr.write("Getting gene_id2snps_id_ls  ...")
		gene_id2snps_id_ls = {}
		for chrpos_key, snps_id in snps_context_wrapper.chrpos2snps_id.iteritems():
			chr_pos_label = '%s_%s'%(chrpos_key[0], chrpos_key[1])
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chrpos_key[0], chrpos_key[1])
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id not in gene_id2snps_id_ls:
					gene_id2snps_id_ls[gene_id] = []
				gene_id2snps_id_ls[gene_id].append(chr_pos_label)
			"""
			if 'touch' in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id]:	#'touch' is closest, forget about all others if there's 'touch'
				for disp_pos, gene_id in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id]['touch']:
					if gene_id not in gene_id2snps_id_ls:
						gene_id2snps_id_ls[gene_id] = []
					gene_id2snps_id_ls[gene_id].append(chr_pos_label)
			else:
				for left_or_right, gene_id_ls in snps_context_wrapper.snps_id2left_or_right2gene_id_ls[snps_id].iteritems():
					for disp_pos, gene_id in gene_id_ls:
						if gene_id not in gene_id2snps_id_ls:
							gene_id2snps_id_ls[gene_id] = []
						gene_id2snps_id_ls[gene_id].append(chr_pos_label)
			"""
		sys.stderr.write("Done.\n")
		return gene_id2snps_id_ls
	
	def generate_params(self, gene_id_fname, pdata, block_size=1000):
		"""
		2009-2-12
			use yield make this function a generator
		2009-2-9
			add argument gene_id_fname to restrict analysis on genes from it
		2008-09-09
			estimate the number of tests each gene would encompass, and decide how many genes should be included in a set to send out
		2008-09-06
			each node handles a certain number of genes. identified by the index of the 1st gene and the index of the last gene.
		"""
		#sys.stderr.write("Generating parameters ...")
		#params_ls = []
		no_of_phenotypes = len(pdata.phenotype_index_ls)
		start_index = 0	#for each computing node: the index of gene >= start_index
		no_of_tests_per_node = 0
		
		#2009-2-9
		if gene_id_fname and os.path.isfile(gene_id_fname):
			reader = csv.reader(open(gene_id_fname), delimiter=figureOutDelimiter(gene_id_fname))
			gene_id_ls = []
			for row in reader:
				gene_id = int(row[0])
				gene_id_ls.append(gene_id)
			del reader
			pdata.gene_id_ls = gene_id_ls	#replace pdata's gene_id_ls
		
		no_of_genes = len(pdata.gene_id_ls)
		
		for i in range(no_of_genes):
			gene_id = pdata.gene_id_ls[i]
			n = len(pdata.gene_id2snps_id_ls[gene_id])	#no_of_snps_of_this_gene
			est_no_of_tests = (n*(n-1)*5/2.0 + n)*no_of_phenotypes	#this is the upper bound for the number of tests for each gene on a computing node. data missing would make the number smaller.
			no_of_tests_per_node += est_no_of_tests
			if no_of_tests_per_node>=block_size:
				yield (start_index, i+1)	#the computing node is gonna handle genes from pdata.gene_id_ls[start_index] to pdata.gene_id_ls[i]
				#reset the starting pointer to the index of the next gene
				start_index = i+1
				no_of_tests_per_node = 0	#reset this to 0
			elif i==no_of_genes-1:	#this is the last gene, have to include them
				yield (start_index, i+1)
				#no need to cleanup because this is the end of loop
		#sys.stderr.write("%s params. Done.\n"%(len(params_ls)))
		#return params_ls
	
	def inputNodePrepare(self, snp_info=None):
		"""
		2009-2-11
			refactored out of run()
		"""
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, \
						data_matrix=data_matrix, turn_into_array=1)	#category_list is not used to facilitate row-id matching
		
		picklef = open(self.snps_context_fname)
		snps_context_wrapper = cPickle.load(picklef)
		del picklef
		gene_id2snps_id_ls = self.get_gene_id2snps_id_ls(snps_context_wrapper)
		del snps_context_wrapper
		gene_id_ls = gene_id2snps_id_ls.keys()
		gene_id_ls.sort()
		
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, data_matrix=data_matrix_phen)
		phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, phenData.row_id_ls, phenData.data_matrix)
		
		self.phenotype_index_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(phenData, Set(self.phenotype_method_id_ls))
		
		if not self.phenotype_index_ls:
			self.phenotype_index_ls = range(len(phenData.col_id_ls))
		pdata = PassingData(gene_id_ls=gene_id_ls, gene_id2snps_id_ls=gene_id2snps_id_ls, \
						phenotype_index_ls=self.phenotype_index_ls, snp_info=snp_info)
		params_ls = self.generate_params(self.gene_id_fname, pdata, self.block_size)
		
		other_data = PassingData(gene_id2snps_id_ls=gene_id2snps_id_ls, gene_id_ls=pdata.gene_id_ls, phenData=phenData, \
								phenotype_index_ls=self.phenotype_index_ls, snp_info=snp_info)
		other_data_pickle = cPickle.dumps(other_data, -1)
		del other_data
		
		output_node_data = PassingData(phenotype_label_ls=phenData.col_id_ls, \
								phenotype_index_ls=self.phenotype_index_ls)
		output_node_data_pickle = cPickle.dumps(output_node_data, -1)
		
		snpData_pickle = cPickle.dumps(snpData, -1)
		del snpData, data_matrix
		return_data = PassingData(snpData_pickle=snpData_pickle, other_data_pickle=other_data_pickle,\
								output_node_data_pickle=output_node_data_pickle, params_ls=params_ls)
		return return_data
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2008-09-07
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		snpData = param_obj.snpData
		no_of_phenotypes = param_obj.phenData.data_matrix.shape[1]
		min_index1, stop1 = data
		for gene_id in param_obj.gene_id_ls[min_index1: stop1]:
			snps_id_ls = param_obj.gene_id2snps_id_ls[gene_id]
			no_of_snps = len(snps_id_ls)
			for phenotype_index in param_obj.phenotype_index_ls:
				if phenotype_index>=no_of_phenotypes:	#skip if out of bound
					continue
				phenotype_ls = param_obj.phenData.data_matrix[:, phenotype_index]
				for i in range(no_of_snps):
					snp1_id = snps_id_ls[i]
					#single SNP test
					snp1_index = snpData.col_id2col_index.get(snp1_id)
					if snp1_index is None:	#snp1_id not in input matrix
						continue
					"""
					pdata = Kruskal_Wallis._kruskal_wallis(snpData.data_matrix[:,snp1_index]+1, phenotype_ls, param_obj.min_data_point)	#watch, +1 to avoid 0 being treated as NA in Kruskal_Wallis._kruskal_wallis
					if pdata:
						pdata = PassingData(snp1_id=snp1_id, gene1_id=gene_id, snp2_id=None, gene2_id=None, pvalue=pdata.pvalue, count1=pdata.count_ls[0],\
										count2=pdata.count_ls[1], bool_type=0, phenotype_index=phenotype_index)	#single SNP
						result_ls.append(pdata)
					"""
					for j in range(i+1, no_of_snps):
						snp2_id = snps_id_ls[j]
						snp2_index = snpData.col_id2col_index.get(snp2_id)
						if snp2_index is None or snp2_index==snp1_index:	#snp2_id not in input matrix
							continue
						pdata_ls = self.test_type2test[param_obj.test_type](snp1_id, gene_id, snpData.data_matrix[:,snp1_index], \
																		snp2_id, None, snpData.data_matrix[:,snp2_index], \
																		phenotype_index, phenotype_ls, \
																		min_data_point=param_obj.min_data_point)
						result_ls += pdata_ls
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def general_output_node(self, output_dir, phenotype_index_ls, phenotype_label_ls, free_computing_nodes):
		"""
		2009-2-8
			general strategy for output node to do while it's computing
			
			refactored out of run()
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		writer_dict = {}
		header_row = ['snp1_id', 'gene1_id', 'snp2_id', 'gene2_id', 'bool_type', 'pvalue', 'count1', 'count2', 'var_perc', 'coeff_list', 'coeff_p_value_list']
		for phenotype_index in phenotype_index_ls:
			phenotype_label = phenotype_label_ls[phenotype_index]
			phenotype_label = phenotype_label.replace('/', '_')	#'/' is taken as folder separator
			output_fname = os.path.join(output_dir, 'SNPpair_%s.tsv'%phenotype_label)
			writer = csv.writer(open(output_fname, 'w'), lineterminator='\n', delimiter='\t')
			writer.writerow(header_row)
			writer_dict[phenotype_index] = writer
		param_obj = PassingData(writer_dict=writer_dict, header_row=header_row)
		self.output_node(free_computing_nodes, param_obj, self.output_node_handler)
		del writer_dict
		
	def output_node_handler(self, communicator, param_obj, data):
		"""
		2008-09-05
		"""
		writer_dict = param_obj.writer_dict
		result_ls = cPickle.loads(data)
		for pdata in result_ls:
			writer = writer_dict[pdata.phenotype_index]
			row = []
			for col_name in param_obj.header_row:
				row.append(getattr(pdata, col_name, None))
			writer.writerow(row)
	
	def run(self):
		"""
		2008-09-06
		"""
		if self.debug:
			#for one-node testing purpose
			import pdb
			pdb.set_trace()
			header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname)
			snpData = SNPData(header=header, strain_acc_list=strain_acc_list, \
							data_matrix=data_matrix, turn_into_array=1)	#category_list is not used to facilitate row-id matching
			
			picklef = open(self.snps_context_fname)
			snps_context_wrapper = cPickle.load(picklef)
			del picklef
			gene_id2snps_id_ls = self.get_gene_id2snps_id_ls(snps_context_wrapper)
			gene_id_ls = gene_id2snps_id_ls.keys()
			gene_id_ls.sort()
			
			header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(self.phenotype_fname, turn_into_integer=0)
			phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list, data_matrix=data_matrix_phen)	#row label is that of the SNP matrix, because the phenotype matrix is gonna be re-ordered in that way
			phenData.data_matrix = Kruskal_Wallis.get_phenotype_matrix_in_data_matrix_order(snpData.row_id_ls, phenData.row_id_ls, phenData.data_matrix)
			
			other_data = PassingData(gene_id2snps_id_ls=gene_id2snps_id_ls, gene_id_ls=gene_id_ls, phenData=phenData)
			other_data_pickle = cPickle.dumps(other_data, -1)
			phenotype_label_ls_pickle = cPickle.dumps(phenData.col_id_ls, -1)
			snpData_pickle = cPickle.dumps(snpData, -1)
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
	main_class = MpiIntraGeneSNPPairAsso
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()