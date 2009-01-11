#!/usr/bin/env mpipython
"""

Examples:
	#run it on hpc-cmb cluster
	mpiexec ~/script/variation/src/MpiGeneListRankTest.py -u yh -p passw**d -t ~/panfs/db/results_by_gene/ -o ~/mpigene_list_rank_test.out -s 100 -c

	#test parallel run on desktop
	mpirun -np 3 -machinefile  /tmp/hostfile /usr/bin/mpipython  ~/script/variation/src/MpiGeneListRankTest.py -u yh -m 20000 -g -p passw**d -s 100 -b -c
	
Description:
	MPI version GeneListRankTest.py. No need to specify list_type_id and results_method_id_ls. Automatically calculates
	for all combinations of results_method_id and list_type_id, skipping the ones that have been done.
	
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math, traceback
import Numeric, cPickle
from Scientific import MPI
from pymodule.MPIwrapper import mpi_synchronize, MPIwrapper
from pymodule import PassingData, getListOutOfStr
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
from Stock_250kDB import Stock_250kDB, Snps, SnpsContext, ResultsMethod, GeneList, GeneListType, \
	CandidateGeneRankSumTestResult, ResultsByGene, CandidateGeneRankSumTestResultMethod
from sets import Set


class MpiGeneListRankTest(GeneListRankTest, MPIwrapper):
	__doc__ = __doc__
	option_default_dict = GeneListRankTest.option_default_dict.copy()
	option_default_dict.update({('message_size', 1, int):[200, 'q', 1, 'How many results one computing node should handle.']})
	option_default_dict.update({('call_method_id', 0, int):[0, 'j', 1, 'Restrict results based on this call_method. Default is no such restriction.']})
	option_default_dict.update({('analysis_method_id_ls', 0, ):[None, '', 1, 'Restrict results based on this analysis_method. Default is no such restriction.']})
	option_default_dict.update({("list_type_id_ls", 0, ): [None, 'l', 1, 'comma/dash-separated list of gene list type ids. ids not present in db will be filtered out. Each id has to encompass>=10 genes.']})
	option_default_dict.update({('alter_hostname', 1, ):['banyan.usc.edu', '', 1, 'host for non-output nodes to connect, since they only query and not save objects. this host can be a slave.']})
	option_default_dict.pop(("list_type_id", 1, int))
	option_default_dict.pop(("results_id_ls", 1, ))
	
	def __init__(self,  **keywords):
		GeneListRankTest.__init__(self, **keywords)
		self.list_type_id_ls = getListOutOfStr(self.list_type_id_ls, data_type=int)
		self.analysis_method_id_ls = getListOutOfStr(self.analysis_method_id_ls, data_type=int)
		
	def generate_params(cls, param_obj, min_no_of_genes=10):
		"""
		2009-1-11
			handle param_obj.no_check_gene_list
				if it's there and true, list_type_id_ls = param_obj.list_type_id_ls. no db check to make sure each gene list has enough genes in it.
				added to let PickCandidateGenesIntoResultsGene.py deal with list_type_id=0 (all genes)
		2008-11-13
			add "results_type==3" same as "results_type==1"
			
			if results_type is not supported, report and return []
		2008-11-08
			become a classmethod
		2008-10-26
			restrict results via param_obj.analysis_method_id_ls  and param_obj.phenotype_method_id_ls
		2008-10-10
			depends on results_type, decide which ( ResultsMethod or ResultsByGene) to get data
			also which (CandidateGeneRankSumTestResult or CandidateGeneRankSumTestResultMethod) to store results
		2008-09-26
			deal with the situation that list_type_id_ls is given.
		2008-09-16
			modify it to get result ids from ResultsByGene
		2008-09-10
			add results_method filtering by analysis_method_id
		2008-08-19
			add call_method_id
		2008-08-15
			stop filtering if CandidateGeneRankSumTestResult has (results_method_id, list_type_id) combo
		2008-07-24
			only association results (results_method_type_id=1)
			only candidate gene lists with >min_no_of_genes genes
			skip ones that been done
		"""
		sys.stderr.write("Generating parameters ...")
		i = 0
		block_size = 5000
		if param_obj.results_type==1 or  param_obj.results_type==3:	#1 and 3 are same ResultsMethod class
			query = ResultsMethod.query
			if param_obj.call_method_id!=0:
				query = query.filter_by(call_method_id=param_obj.call_method_id)
			if hasattr(param_obj, 'analysis_method_id_ls') and param_obj.analysis_method_id_ls:
				query = query.filter(ResultsMethod.analysis_method_id.in_(param_obj.analysis_method_id_ls))
			if hasattr(param_obj, 'phenotype_method_id_ls') and param_obj.phenotype_method_id_ls:
				query = query.filter(ResultsMethod.phenotype_method_id.in_(param_obj.phenotype_method_id_ls))
		elif param_obj.results_type==2:
			query = ResultsByGene.query
			if param_obj.call_method_id!=0:
				query = query.filter(ResultsByGene.results_method.has(call_method_id=param_obj.call_method_id))
			if hasattr(param_obj, 'analysis_method_id') and param_obj.analysis_method_id!=0 and param_obj.analysis_method_id is not None:
				query = query.filter(ResultsByGene.results_method.has(analysis_method_id=param_obj.analysis_method_id))
			if param_obj.analysis_method_id_ls:
				query = query.filter(ResultsByGene.results_method.has(ResultsMethod.analysis_method_id.in_(param_obj.analysis_method_id_ls)))
			if hasattr(param_obj, 'phenotype_method_id_ls') and param_obj.phenotype_method_id_ls:
				query = query.filter(ResultsByGene.results_method.has(ResultsMethod.phenotype_method_id.in_(param_obj.phenotype_method_id_ls)))
		else:
			sys.stderr.write("results_type %s not supported.\n"%results_type)
			return []
		
		rows = query.offset(i).limit(block_size)
		results_method_id_ls = []
		while rows.count()!=0:
			for row in rows:
				results_method_id_ls.append(row.id)
				i += 1
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("%s results. "%(len(results_method_id_ls)))
		
		#if self.debug:	#2008-10-25 temporary testing
		#	results_method_id_ls = [2095, 2079]
		
		if getattr(param_obj, 'no_check_gene_list', None) and getattr(param_obj, 'list_type_id_ls', None):
			list_type_id_ls = param_obj.list_type_id_ls
		else:
			list_type_id_ls = []
			if getattr(param_obj, 'list_type_id_ls', None):	#if list_type_id_ls is given, check whether each one exists in db and has minimum number of genes.
				for list_type_id in param_obj.list_type_id_ls:
					glt = GeneListType.get(list_type_id)
					if glt and len(glt.gene_list)>=min_no_of_genes:
						list_type_id_ls.append(list_type_id)
			else:
				i = 0
				rows = GeneListType.query.offset(i).limit(block_size)
				while rows.count()!=0:
					for row in rows:
						if len(row.gene_list)>=min_no_of_genes:
							list_type_id_ls.append(row.id)
						i += 1
					rows = GeneListType.query.offset(i).limit(block_size)
		sys.stderr.write("%s candidate gene lists. "%(len(list_type_id_ls)))
		
		rm_id_lt_id_set = Set()
		"""
		i = 0
		rows = CandidateGeneRankSumTestResult.query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				rm_id_lt_id_set.add((row.results_method_id, row.list_type_id))
				i += 1
			rows = CandidateGeneRankSumTestResult.query.offset(i).limit(block_size)
		sys.stderr.write("%s candidate gene rank sum test results. "%(len(rm_id_lt_id_set)))
		"""
		
		params_ls = []
		for results_method_id in results_method_id_ls:
			for list_type_id in list_type_id_ls:
				rm_id_lt_id = (results_method_id, list_type_id)
				if rm_id_lt_id not in rm_id_lt_id_set:
					params_ls.append(rm_id_lt_id)
		sys.stderr.write(" %s params generated.\n"%(len(params_ls)))
		return params_ls
	
	generate_params=classmethod(generate_params)
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		params_ls = parameter_list[0]
		data_to_return = []
		for i in range(message_size):
			if len(params_ls)>0:
				one_parameter = params_ls.pop(0)
				data_to_return.append(one_parameter)
			else:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return data_to_return
	
	def computing_node_handler(self, communicator, data, param_obj):
		"""
		2008-08-20
			wrap all parameters into pd and pass it to run_wilcox_test
		2008-07-17
		
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result_ls = []
		pd = PassingData(snps_context_wrapper=param_obj.snps_context_wrapper,\
							results_directory=param_obj.results_directory,\
							min_MAF=param_obj.min_MAF, get_closest=self.get_closest, min_distance=self.min_distance, \
							min_sample_size=self.min_sample_size, test_type=self.test_type, \
							results_type=self.results_type, no_of_permutations=self.no_of_permutations,\
							no_of_min_breaks=self.no_of_min_breaks)
		for results_method_id, list_type_id in data:
			pd.results_id = results_method_id
			pd.list_type_id = list_type_id
			result = self.run_wilcox_test(pd)
			if result is not None:
				result_ls.append(result)
		sys.stderr.write("Node no.%s done with %s results.\n"%(node_rank, len(result_ls)))
		return result_ls
	
	def output_node_handler(self, communicator, output_param_obj, data):
		"""
		2008-11-04
			refactor a portion to become sub_output()
		2008-10-30
			upon session.flush() failure, get into 'except ...' and remove the result from memory and report traceback without interrupting program
		2008-10-23
			replace parameter_list with output_param_obj 
			handle result.type
		2008-05-19
			row is strain. col is snp. reversed due to utilization of SNPData
		05/14/2008
			flush outf and outf_avg
		05/12/2008
			common_var_name_ls
		"""
		table_obj_ls = cPickle.loads(data)
		self.sub_output(output_param_obj, table_obj_ls, output_param_obj.TestResultClass)
	
	def sub_output(self, output_param_obj, table_obj_ls, TableClass):
		"""
		2008-11-04
		
			split out of output_node_handler() so that MpiTopSNPTest's output_node_handler can call this.
		"""
		commit = output_param_obj.commit
		for table_obj in table_obj_ls:
			result = self.returnResultFromDB(TableClass, table_obj.results_id, table_obj.list_type_id, table_obj.starting_rank,
											table_obj.type_id, table_obj.min_distance, table_obj.no_of_top_snps, table_obj.min_score)
			if result:
				if self.debug:
					sys.stderr.write("TestResultClass with results_id=%s, list_type_id=%s, starting_rank=%s,\
							type_id=%s, min_distance=%s, no_of_top_snps=%s, min_score=%s) already in db with id=%s. skip.\n"%\
							(table_obj.results_id, table_obj.list_type_id, table_obj.starting_rank,
							table_obj.type_id, table_obj.min_distance, table_obj.no_of_top_snps, \
							table_obj.min_score, result.id))
				continue
			
			row = []
			
			result = TableClass()
			#pass values from table_obj to this new candidate_gene_rank_sum_test_result.
			#can't save table_obj because it's associated with a different db thread
			for column in table_obj.c.keys():
				if output_param_obj.writer:	#2008-10-30 append only when writer is not None
					row.append(getattr(table_obj, column))
				setattr(result, column, getattr(table_obj, column))
			
			if hasattr(result, 'type_id') and result.type_id is None and getattr(output_param_obj, '_type', None):	#2008-10-23 assign _type
				if self.debug:
					sys.stderr.write("type_id not avaiable on the result got from computing node. assign type here.\n")
				result.type = output_param_obj._type
			
			if output_param_obj.writer:
				output_param_obj.writer.writerow(row)
			
			output_param_obj.session.save(result)
			if commit:
				try:
					output_param_obj.session.flush()
				except:
					#2008-10-30 remove it from memory. otherwise, next flush() will try on this old object again.
					output_param_obj.session.expunge(result)
					#output_param_obj.session.delete(result)
					sys.stderr.write("Exception happened for results_method_id=%s, list_type_id=%s.\n"%(getattr(result, 'results_id', None),\
																									getattr(result, 'list_type_id', None)))
					for column in table_obj.c.keys():
						sys.stderr.write("\t%s=%s.\n"%(column, getattr(table_obj, column)))
					traceback.print_exc()
					sys.stderr.write('%s.\n'%repr(sys.exc_info()))
	
	def run(self):
		"""
		2008-07-17
		"""
		self.communicator = MPI.world.duplicate()
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		free_computing_node_set = Set(free_computing_nodes)
		output_node_rank = self.communicator.size-1
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		if self.results_type==1:
			ResultsClass = ResultsMethod
			TestResultClass = CandidateGeneRankSumTestResultMethod
		elif self.results_type==2:
			ResultsClass = ResultsByGene
			TestResultClass = CandidateGeneRankSumTestResult
		else:
			sys.stderr.write("Error: Invalid results type : %s.\n"%pd.results_type)
			
		if node_rank == 0:
			#snps_context_wrapper = self.constructDataStruc(self.min_distance, self.get_closest)
			snps_context_wrapper = self.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
			param_obj = PassingData(call_method_id=self.call_method_id,
								analysis_method_id=getattr(self, 'analysis_method_id', None),\
								analysis_method_id_ls=getattr(self, 'analysis_method_id_ls', None),\
								list_type_id_ls=self.list_type_id_ls,
								results_type=self.results_type,\
								phenotype_method_id_ls=getattr(self, 'phenotype_method_id_ls', None))
			params_ls = self.generate_params(param_obj)
			if self.debug:
				params_ls = params_ls[:100]
			
			snps_context_wrapper_pickle = cPickle.dumps(snps_context_wrapper, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				sys.stderr.write("passing initial data to nodes from %s to %s ... "%(node_rank, node))
				self.communicator.send(snps_context_wrapper_pickle, node, 0)
				sys.stderr.write(".\n")
			del snps_context_wrapper_pickle
			del snps_context_wrapper
			
		elif node_rank in free_computing_node_set:			
			data, source, tag = self.communicator.receiveString(0, 0)
			snps_context_wrapper =  cPickle.loads(data)
			del data
			sys.stderr.write(".\n")
			
		else:
			pass
		
		self.synchronize()
		if node_rank == 0:
			parameter_list = [params_ls]
			self.input_node(parameter_list, free_computing_nodes, input_handler=self.input_handler, message_size=self.message_size)
		elif node_rank in free_computing_node_set:
			computing_parameter_obj = PassingData(snps_context_wrapper=snps_context_wrapper, \
												results_directory=self.results_directory, min_MAF=self.min_MAF)
			self.computing_node(computing_parameter_obj, self.computing_node_handler)
		else:
			if getattr(self, 'output_fname', None):
				writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
				header_row = []
				for column in TestResultClass.c.keys():
					header_row.append(column)
				writer.writerow(header_row)
			else:
				writer = None
			
			output_param_obj = PassingData(writer=writer, session=session, commit=self.commit, TestResultClass=TestResultClass,
										_type=None)
			self.output_node(free_computing_nodes, output_param_obj, self.output_node_handler)
			del writer		
		self.synchronize()	#to avoid some node early exits

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = MpiGeneListRankTest
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
