#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	chr=1;mpiexec ~/script/variation/src/MpiGADA.py -u yh -p passw**d -i ~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr$chr.tsv -o ...

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiGADA.py -u yh -p passw**d -M 5-7 -T 2,4,8,12 -a 0.2,0.5,0.8,1.0 ...
	
Description:
	MPI version RunGADA.py. Automatically calculates all combinations of aAlpha, TBackElim, MinSegLen.
	All output files have the same prefix as output_fname_prefix, followed by paramter combination.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, traceback
import cPickle, traceback
from Scientific import MPI
from pymodule.MPIwrapper import mpi_synchronize, MPIwrapper
from pymodule import PassingData, getListOutOfStr
from RunGADA import RunGADA
import Stock_250kDB
from sets import Set


class MpiGADA(RunGADA, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = RunGADA.option_default_dict.copy()
	option_default_dict.pop(('aAlpha', 1, float))
	option_default_dict.pop(('TBackElim', 1, float))
	option_default_dict.pop(('MinSegLen', 1, int))
	option_default_dict.pop(('output_fname', 1, ))
	option_default_dict.update({('output_fname_prefix', 1, ): ['', 'o', 1, 'output files all have this prefix + parameter combo. Non-existent folder would be created.', ]})
	option_default_dict.update({('aAlpha_ls', 1, ):["0.5", 'A', 1, 'a coma-separated list of alphas to be tried']})
	option_default_dict.update({('TBackElim_ls', 1, ):["4", 'T', 1, 'a list of (amp1-amp2)/stddev in GADA']})
	option_default_dict.update({('MinSegLen_ls', 1, ):["5,10", 'M', 1, 'a coma-separated list of minimum no of probes to call a segment in GADA']})
	option_default_dict.update({('alter_hostname', 1, ):['banyan.usc.edu', '', 1, 'host for non-output nodes to connect, since they only query and not save objects. this host can be a slave.']})
	option_default_dict.update({('message_size', 1, int):[1, 'q', 1, 'How many parameter combos one computing node should handle.']})
	
	def __init__(self,  **keywords):
		RunGADA.__init__(self, **keywords)
		self.aAlpha_ls = getListOutOfStr(self.aAlpha_ls, data_type=float)
		self.TBackElim_ls = getListOutOfStr(self.TBackElim_ls, data_type=float)
		self.MinSegLen_ls = getListOutOfStr(self.MinSegLen_ls, data_type=int)
	
	def generate_params(cls, aAlpha_ls, TBackElim_ls, MinSegLen_ls, array_id_ls, param_obj=None):
		"""
		2009-10-27
		"""
		sys.stderr.write("Generating parameters ...")
		params_ls = []
		for aAlpha in aAlpha_ls:
			for TBackElim in TBackElim_ls:
				for MinSegLen in MinSegLen_ls:
					for array_id in array_id_ls:
						params_ls.append((aAlpha, TBackElim, MinSegLen, array_id))
		sys.stderr.write(" %s params generated.\n"%(len(params_ls)))
		return params_ls
	
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		2009-10-27
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		params_ls, array_id2col_index, data_matrix = parameter_list[:3]
		data_to_return = []
		for i in range(message_size):
			if len(params_ls)>0:
				one_parameter = params_ls.pop(0)
				array_id = one_parameter[-1]
				col_index = array_id2col_index[array_id]
				one_data = [one_parameter, data_matrix[:, col_index]]
				data_to_return.append(one_data)
			else:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2009-10-28
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		for one_data in data:
			one_parameter, intensity_ls = one_data[:2]
			aAlpha, TBackElim, MinSegLen, array_id = one_parameter[:4]
			tmp_input_fname = '%s_%s_A%sT%sM%s'%(param_obj.tmp_input_fname, array_id, aAlpha, TBackElim, MinSegLen)
			input_file = self.prepareGADAinput(intensity_ls, tmp_input_fname, IO_thru_file=self.IO_thru_file)
			tmp_output_fname = '%s_%s_A%sT%sM%s'%(param_obj.tmp_output_fname, array_id, aAlpha, TBackElim, MinSegLen)
			GADA_output = self._GADA(self.GADA_path, input_file, tmp_output_fname, aAlpha, \
									TBackElim, MinSegLen)
			
			if GADA_output:
				result_ls.append((one_parameter, GADA_output))
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	array_id2array = {}
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2009-10-28
		"""
		result_ls = cPickle.loads(data)
		counter = 0
		for result in result_ls:
			one_parameter, GADA_output = result
			aAlpha, TBackElim, MinSegLen, array_id = one_parameter[:4]
			if array_id in self.array_id2array:
				array = self.array_id2array.get(array_id)
			else:
				array = Stock_250kDB.ArrayInfo.get(array_id)
				self.array_id2array[array_id] = array
			
			if array:
				ecotype_id = array.maternal_ecotype_id
			else:
				ecotype_id = ''
			param_combo = (aAlpha, TBackElim, MinSegLen)
			if param_combo not in output_param_obj.param_combo2writer:
				output_fname_prefix = os.path.splitext(output_param_obj.output_fname_prefix)[0]
				output_fname = '%s_A%sT%sM%s.tsv'%(output_fname_prefix, aAlpha, TBackElim, MinSegLen)
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				self.output_header(writer)
				output_param_obj.param_combo2writer[param_combo] = writer
				output_param_obj.param_combo2array_id2no_of_segments[param_combo] = {}
			writer = output_param_obj.param_combo2writer[param_combo]
			no_of_segments = self.output_GADA_output(writer, GADA_output, array_id, ecotype_id, output_param_obj.chr_pos_ls, \
													output_param_obj.probe_id_ls)
			output_param_obj.param_combo2array_id2no_of_segments[param_combo][array_id] = no_of_segments
			counter += 1
		sys.stderr.write("%s GADA results were outputted.\n"%counter)
	
	def outputSegmentStat(self, param_combo2array_id2no_of_segments):
		"""
		2009-10-28
			output no_of_segments_per_array for each parameter combo
		"""
		param_combo_ls = param_combo2array_id2no_of_segments.keys()
		param_combo_ls.sort()
		for param_combo in param_combo_ls:
			array_id2no_of_segments = param_combo2array_id2no_of_segments.get(param_combo)
			no_of_arrays = len(array_id2no_of_segments)
			no_of_segments_ls = array_id2no_of_segments.values()
			no_of_segments_per_array = sum(no_of_segments_ls)/float(no_of_arrays)
			sys.stderr.write("Param-combo (a, T, M) %s: %s segments per array.\n"%(repr(param_combo), no_of_segments_per_array))
	
	def run(self):
		"""
		2009-10-28
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		
		if node_rank == 0:
			data_matrix, probe_id_ls, chr_pos_ls, header = self.get_input(self.input_fname)
			array_id_ls = header[1:-2]
			array_id2col_index = {}
			for i in range(len(array_id_ls)):
				array_id_ls[i] = int(array_id_ls[i])
				array_id2col_index[array_id_ls[i]] = i
			if self.which_array_id_ls:	# only deal with selected arrays
				array_id_ls = self.which_array_id_ls
			
			passing_to_output_node_obj = PassingData(probe_id_ls=probe_id_ls, chr_pos_ls=chr_pos_ls, array_id_ls=array_id_ls)
		 	passing_to_output_node_obj_pickle = cPickle.dumps(passing_to_output_node_obj, -1)
		 	sys.stderr.write("Passing data to output node %s from %s ... "%(output_node_rank, node_rank,))
			self.communicator.send(passing_to_output_node_obj_pickle, output_node_rank, 0)
			sys.stderr.write(".\n")
			del passing_to_output_node_obj_pickle, passing_to_output_node_obj 
			
			params_ls = self.generate_params(self.aAlpha_ls, self.TBackElim_ls, self.MinSegLen_ls, array_id_ls)
			if self.debug:
				params_ls = params_ls[:100]
			
		elif node_rank in free_computing_node_set:			
			pass
		else:
			db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
			db.setup(create_tables=False)
			session = db.session
			
			data, source, tag = self.communicator.receiveString(0, 0)
			passing_to_output_node_obj =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
		
		self.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls, array_id2col_index, data_matrix]
			self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			if self.IO_thru_file:
				tmp_input_dir = os.path.split(self.tmp_input_fname)[0]
				tmp_output_dir = os.path.split(self.tmp_output_fname)[0]
				try:
					if not os.path.isdir(tmp_input_dir):
						os.makedirs(tmp_input_dir)
					if not os.path.isdir(tmp_output_dir):
						os.makedirs(tmp_output_dir)
				except:
					 traceback.print_exc()
					 sys.stderr.write('Except type: %s\n'%repr(sys.exc_info()))
			computing_parameter_obj = PassingData(tmp_input_fname=self.tmp_input_fname, tmp_output_fname=self.tmp_output_fname)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			output_dir = os.path.split(self.output_fname_prefix)[0]
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			param_combo2writer = {}
			param_combo2array_id2no_of_segments = {}
			output_param_obj = PassingData(param_combo2writer=param_combo2writer, output_fname_prefix=self.output_fname_prefix,\
										chr_pos_ls = passing_to_output_node_obj.chr_pos_ls,\
										probe_id_ls = passing_to_output_node_obj.probe_id_ls,\
										param_combo2array_id2no_of_segments = param_combo2array_id2no_of_segments)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del param_combo2writer
			self.outputSegmentStat(param_combo2array_id2no_of_segments)
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiGADA
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
