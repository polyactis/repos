#!/usr/bin/env mpipython
"""

Examples:
	mpiexec MpiQCCall.py -i /mnt/nfs/NPUTE_data/input/250k_method_3.csv -p /mnt/nfs/NPUTE_data/input/perlgen.csv -f /mnt/nfs/NPUTE_data/input/2010_149_384.csv -o /tmp/param_qc.out
	
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython ~/script/variation/src/MpiQCCall.py
	
	#test
	mpirun -np 5 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiQCCall.py -i /mnt/nfs/NPUTE_data/input/250K_m3_70_n1000.csv -p /mnt/nfs/NPUTE_data/input/perlgen.csv -f /mnt/nfs/NPUTE_data/input/2010_149_384.csv -o /tmp/param_qc.out -n 20 -v 0.3 -a 0.4 -w 0.2

Description:
	a parallel program to do QC on genotype calls before/after imputation under various parameter settings.
	
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
#from FilterStrainSNPMatrix import FilterStrainSNPMatrix	#for read_data()
from common import nt2number,number2nt
import dataParsers, FilterAccessions, FilterSnps, MergeSnpsData
from variation.genotyping.NPUTE.SNPData import SNPData as NPUTESNPData
from variation.genotyping.NPUTE.NPUTE import imputeData
from variation.src.QC_250k import SNPData, TwoSNPData
from pymodule import PassingData
import copy

class MpiQCCall(object):
	__doc__ = __doc__
	option_default_dict = {('input_fname', 1, ): ['', 'i', 1, '250k call file'],\
							('callProbFile', 0, ): ['', 't', 1, '250k probability file'],\
							('fname_2010_149_384', 1, ): ['', 'f', 1, ''],\
							('fname_perlegen', 1, ): ['', 'p', 1, ''],\
							('output_fname', 1, ): [None, 'o', 1, '', ],\
							('min_call_probability_ls', 1, ): ["0.6,0.7,0.75,0.8,0.85,0.9,0.95,0.975,0.982", 'y', 1, 'minimum probability for a call to be non-NA if there is a 3rd column for probability.', ],\
							('max_call_mismatch_rate_ls', 1, ): ["0.05,0.1,0.15,0.20,0.25,0.3", 'x', 1, 'maximum mismatch rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_call_NA_rate_ls', 1, ): ["0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8", 'a', 1, 'maximum NA rate of an array call_info entry. used to exclude bad arrays.'],\
							('max_snp_mismatch_rate_ls', 1, ): ["0.05,0.1,0.15,0.2,0.25,0.3", 'w', 1, 'maximum snp error rate, used to exclude bad SNPs', ],\
							('max_snp_NA_rate_ls', 1, ): ["0.1,0.2,0.3,0.4,0.5,0.6,0.7", 'v', 1, 'maximum snp NA rate, used to exclude SNPs with too many NAs', ],\
							('npute_window_size_ls', 1): ["10,20,30,40,50",],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-05-08
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def generate_parameters(self, parameter_names, parameter_depth=2):
		"""
		parameter_depth is not used now
		"""
		sys.stderr.write( "Generating parameter settings ...")
		param_d = PassingData()
		for parameter_name in parameter_names:
			parameter_value = getattr(self, parameter_name)
			parameter_value = parameter_value.split(',')
			parameter_value = map(float, parameter_value)
			setattr(self, parameter_name, parameter_value)
		
		#figure out call probability from input_fname
		import re
		call_prob_pattern = re.compile(r'_(\d+)\.csv')
		call_prob_p_result = call_prob_pattern.search(self.input_fname)
		if call_prob_p_result:
			min_call_probability = float(call_prob_p_result.groups()[0])
		else:
			min_call_probability = -1
		
		
		#only 1st 4, last 2 passed to computing node
		parameters = []
		for max_call_mismatch_rate in getattr(self, parameter_names[0]):
			for max_call_NA_rate in getattr(self, parameter_names[1]):
				for max_snp_mismatch_rate in getattr(self, parameter_names[2]):
					parameters.append([min_call_probability, max_call_mismatch_rate, max_call_NA_rate, max_snp_mismatch_rate])
		
		param_d.parameters = parameters
		param_d.max_snp_NA_rate_ls = self.max_snp_NA_rate_ls
		param_d.npute_window_size_ls = self.npute_window_size_ls
		sys.stderr.write("Done.\n")
		return param_d
	
	def input_handler(self, data, message_size, report=0):
		"""
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		param_d = data
		if param_d.index<len(param_d.parameters):
			data_to_return = param_d.parameters[param_d.index]
			param_d.index += 1
		else:
			data_to_return = None
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def computing_node_handler(self, communicator, data, init_data):
		"""
		2007-03-07
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result = []
		min_call_probability, max_call_mismatch_rate, max_call_NA_rate, max_snp_mismatch_rate = data
		
		#for i in range(0,len(init_data.snpsd_250k)):
		#	sys.stderr.write("Fraction converted ="%(init_data.snpsd_250k[i].convertBadCallsToNA(min_call_probability)))
		
		FilterAccessions.filterByError(init_data.snpsd_250k, init_data.snpsd_2010_149_384, max_call_mismatch_rate, withArrayIds=1)
		FilterAccessions.filterByNA(init_data.snpsd_250k, max_call_NA_rate, withArrayIds=1)
		
		FilterSnps.filterByError(init_data.snpsd_250k, init_data.snpsd_perlegen, max_snp_mismatch_rate)
		
		result = []
		for max_snp_NA_rate in init_data.param_d.max_snp_NA_rate_ls:
			snpsd_250k_tmp = copy.copy(init_data.snpsd_250k)
			FilterSnps.filterByNA(snpsd_250k_tmp, max_snp_NA_rate)
			for npute_window_size in init_data.param_d.npute_window_size_ls:
				MergeSnpsData.merge(snpsd_250k_tmp, init_data.snpsd_perlegen, unionType=0, priority=2)
				FilterSnps.filterMonomorphic(snpsd_250k_tmp)
				
				chr_pos_ls = []
				data_matrix_250k_0 = []
				data_matrix_250k = []
				chr_pos_ls_ref = []
				data_matrix_ref = []
				for i in range(len(snpsd_250k_tmp)):
					data_matrix_250k_0 += snpsd_250k_tmp[i].snps
					
					npute_data_struc = NPUTESNPData(inFile=snpsd_250k_tmp[i], input_NA_char='NA', input_file_format=4, lower_case_for_imputation=0)
					imputeData(npute_data_struc, int(npute_window_size))
					
					this_chr_pos_ls = [(i+1, pos) for pos in snpsd_250k_tmp[i].positions]	#chromosome is i+1
					chr_pos_ls += this_chr_pos_ls
					this_chr_pos_ls = [(i+1, pos) for pos in init_data.snpsd_2010_149_384[i].positions]	#chromosome is i+1
					chr_pos_ls_ref += this_chr_pos_ls
					data_matrix_250k += npute_data_struc.snps
					data_matrix_ref += init_data.snpsd_2010_149_384[i].snps
				
				snpData0 = SNPData(header=['', '']+snpsd_250k_tmp[0].accessions,\
					strain_acc_list=chr_pos_ls, category_list=chr_pos_ls, data_matrix=data_matrix_250k_0)
				del data_matrix_250k_0
				
				snpData1 = SNPData(header=['', '']+snpsd_250k_tmp[0].accessions,\
					strain_acc_list=chr_pos_ls, category_list=chr_pos_ls, data_matrix=data_matrix_250k)
				del data_matrix_250k
				
				snpData2 = SNPData(header=['', '']+init_data.snpsd_2010_149_384[0].accessions, \
					strain_acc_list=chr_pos_ls_ref, category_list= chr_pos_ls_ref, data_matrix=data_matrix_ref)
				del data_matrix_ref
				
				qcdata = PassingData()
				
				twoSNPData0 = TwoSNPData(SNPData1=snpData0, SNPData2=snpData2, \
								row_id_matching_type=2)
				qcdata.row_id2NA_mismatch_rate0 = twoSNPData0.cmp_row_wise()
				qcdata.col_id2NA_mismatch_rate0 = twoSNPData0.cmp_col_wise()
				del twoSNPData0, snpData0
				
				twoSNPData1 = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2, \
								row_id_matching_type=2)
				qcdata.row_id2NA_mismatch_rate1 = twoSNPData1.cmp_row_wise()
				qcdata.col_id2NA_mismatch_rate1 = twoSNPData1.cmp_col_wise()
				del twoSNPData1, snpData1, snpData2
				
				qcdata.min_call_probability = min_call_probability
				qcdata.max_call_mismatch_rate = max_call_mismatch_rate
				qcdata.max_call_NA_rate = max_call_NA_rate
				qcdata.max_snp_mismatch_rate = max_snp_mismatch_rate
				qcdata.max_snp_NA_rate = max_snp_NA_rate
				qcdata.npute_window_size = npute_window_size
				result.append(qcdata)
		sys.stderr.write("Node no.%s done with %s QC.\n"%(node_rank, len(result)))
		return result
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		"""
		writer = parameter_list[0]
		result = cPickle.loads(data)
		id2NA_mismatch_rate_name_ls = ['row_id2NA_mismatch_rate0', 'col_id2NA_mismatch_rate0', 'row_id2NA_mismatch_rate1', 'col_id2NA_mismatch_rate1']
		for qcdata in result:
			common_ls = [qcdata.min_call_probability, qcdata.max_call_mismatch_rate, qcdata.max_call_NA_rate,\
			qcdata.max_snp_mismatch_rate, qcdata.max_snp_NA_rate, qcdata.npute_window_size]
			for id2NA_mismatch_rate_name in id2NA_mismatch_rate_name_ls:
				if not hasattr(qcdata, id2NA_mismatch_rate_name):
					continue
				id2NA_mismatch_rate = getattr(qcdata, id2NA_mismatch_rate_name)
				row_id_ls = id2NA_mismatch_rate.keys()
				row_id_ls.sort()	#try to keep them in call_info_id order
				if id2NA_mismatch_rate_name[:3]=='row':
					type_of_id = 'snp'
				else:
					type_of_id = 'strain'
				if id2NA_mismatch_rate_name[-1]=='0':
					type_of_id += '(before imputation)'
				else:
					type_of_id += '(after imputation)'
				for row_id in row_id_ls:
					NA_mismatch_ls = id2NA_mismatch_rate[row_id]
					writer.writerow(['%s %s'%(type_of_id, repr(row_id))] + common_ls + NA_mismatch_ls)
	
	def run(self):
		"""
		2008-05-09
			(rank==0)
				--get_chr_start_ls()
			elif free_computing_nodes:
				-- (receive data)
			
			--mpi_synchronize()
			
			(rank==0)
				--input_node()
					--input_handler()
			elif free_computing_nodes:
				--computing_node()
					--computing_node_handler()
						--identify_ancestry_with_min_jumps()
							--initialize_score_trace_matrix()
								--is_child_heterozygous_SNP_compatible_with_parents()
							(for loop)
								--identify_ancestry_of_one_chr_with_DP()
									--is_child_heterozygous_SNP_compatible_with_parents()
							--trace()
								--recursive_trace()
			else:
				--output_node()
					--output_node_handler()
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		self.parameter_names = ['max_call_mismatch_rate_ls', 'max_call_NA_rate_ls', \
			'max_snp_mismatch_rate_ls', 'max_snp_NA_rate_ls', 'npute_window_size_ls']
		data_to_pickle_name_ls = ['snpsd_250k', 'snpsd_2010_149_384', 'snpsd_perlegen', 'param_d']
		if node_rank == 0:
			init_data = PassingData()
			init_data.snpsd_250k = dataParsers.parseCSVData(self.input_fname, withArrayIds=True)
			init_data.snpsd_2010_149_384 = dataParsers.parseCSVData(self.fname_2010_149_384, deliminator=',')
			init_data.snpsd_perlegen = dataParsers.parseCSVData(self.fname_perlegen)
			param_d = self.generate_parameters(self.parameter_names)
			init_data.param_d = param_d
			
			sys.stderr.write("passing data to nodes from %s ..."%node_rank)
			for node in free_computing_nodes:	#send it to the computing_node
				for data_to_pickle_name in data_to_pickle_name_ls:
					data_pickle = cPickle.dumps(getattr(init_data, data_to_pickle_name), -1)
					self.communicator.send(data_pickle, node, 0)
					del data_pickle
			del init_data
			sys.stderr.write("Done.\n")
		elif node_rank in free_computing_nodes:
			init_data = PassingData()
			for data_to_pickle_name in data_to_pickle_name_ls:
				data, source, tag = self.communicator.receiveString(0, 0)
				setattr(init_data, data_to_pickle_name, cPickle.loads(data))
				del data
		else:
			pass
		
		mpi_synchronize(self.communicator)
		
		mw = MPIwrapper(self.communicator)
		if node_rank == 0:
			param_d.index = 0
			mw.input_node(param_d, free_computing_nodes, report=self.report, input_handler=self.input_handler)
		elif node_rank in free_computing_nodes:
			parameter_list = init_data
			mw.computing_node(parameter_list, self.computing_node_handler, report=self.report)
		else:
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			header = ['strain/snp', 'min_call_probability'] + self.parameter_names + ['NA_rate', 'mismatch_rate', 'no_of_NAs', 'no_of_totals', 'no_of_mismatches', 'no_of_non_NA_pairs']
			writer.writerow(header)
			parameter_list = [writer]
			mw.output_node(free_computing_nodes, parameter_list, self.output_node_handler, self.report)
			del writer

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiQCCall
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	
