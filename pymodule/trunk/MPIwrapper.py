#!/usr/bin/env python
"""
2007-10-23
	module which wraps MPI
	
	functions originally copied from annot.bin.codense.common.
"""
import sys, os, cPickle

def mpi_synchronize(communicator):
	"""
	05-19-05
		copied from MpiBiclustering.py
	"""
	sys.stdout.flush()
	sys.stderr.flush()
	communicator.barrier()
	
def fetch_cluster_block(curs, message_size, report=0):
	"""
	10-20-05
		commonly used by input_node()
	"""
	if report:
		sys.stderr.write("Fetching stuff...\n")
	curs.execute("fetch 50 from crs")
	rows = curs.fetchall()
	prediction_ls = []
	string_length = 0
	while rows:
		for row in rows:
			prediction_ls.append(list(row))
			string_length += len(repr(row))	#the length to control MPI message size
		if string_length>=message_size:
			break
		curs.execute("fetch 50 from crs")
		rows = curs.fetchall()
	if report:
		sys.stderr.write("Fetching done.\n")
	return prediction_ls

def mpi_schedule_jobs(communicator, job_list, node_function, node_parameter_list, debug=0):
	"""
	05-19-05
		a universal scheduling function, the elements in job_list
		maybe string, integer, or something else, It's 'repr'ed before send.
		WARNING: NO -1 in job_list.
		So node_function should handle it as 'repr'ed.
		
		node_function()
			input: (element of the job_list, node_parameter_list).
			output: returns a value('repr'ed)
		
		
	"""
	from sets import Set
	node_returned_value_list = []
	node_rank = communicator.rank
	if node_rank == 0:
		sys.stderr.write("\tTotally, %d jobs to be scheduled.\n"%len(job_list))
		seed_utilized = Set()
		for node in range(1, communicator.size):
			if len(job_list)==0:	#if #nodes > #jobs, tell those nodes to break their listening loop.
				stop_signal = "-1"
				communicator.send(stop_signal, node, 0)	#no more jobs, stop that node,
				if debug:
					sys.stderr.write("node %s stopped.\n"%node)
			else:
				job = job_list.pop(0)	#the first item poped first.
				communicator.send(repr(job), node, 0)	#string format
				if debug:
					sys.stderr.write("node %s schedule a job, %s to %s\n"%(node_rank, repr(job), node))
				seed_utilized.add(node)
		
		received_value, source, tag = communicator.receiveString(None, None)	#listen
		while received_value:
			node_returned_value_list.append(received_value)
			if len(job_list) == 0:	#first check if there're still files left, otherwise pop(0) raises error.
				stop_signal = "-1"
				communicator.send(stop_signal, source, 0)	#no more jobs, stop that node,
				if debug:
					sys.stderr.write("node %s stopped.\n"%source)
				seed_utilized.remove(source)
				if len(seed_utilized) == 0:	#all seed used have finished their jobs
					break
			else:
				job = job_list.pop(0)
				communicator.send(repr(job), source, 0)	#string format,
				if debug:
					sys.stderr.write("node %s get one more job, %s\n"%(source, repr(job)) )
			received_value, source, tag = communicator.receiveString(None, None)	#listen
	else:
		received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0,
		while received_data:
			if received_data=="-1":	#stop signal
				if debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				sys.stderr.write("node %s working on %s...\n"%(node_rank, received_data))
				node_return_value = node_function(received_data, node_parameter_list)
				sys.stderr.write("node %s work on %s finished.\n"%(node_rank, received_data))
				communicator.send(repr(node_return_value), 0, node_rank)
				
			received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0
	
	return node_returned_value_list

class MPIwrapper(object):
	"""
	2008-05-08
		a general wrapper for mpi programs.
		
		functions copied from annot.bin.codense.common. made into a class.
	"""
	def __init__(self, communicator, debug=0, report=0):
		self.communicator = communicator
		self.debug = debug
		self.report = report
	
	"""
	10-21-05
		below output_node(), fetch_cluster_block(),input_node(), computing_node() are shared
		functions for Mpi Programs
	"""
	def synchronize(self):
		"""
		05-19-05
			copied from MpiBiclustering.py
		"""
		sys.stdout.flush()
		sys.stderr.flush()
		self.communicator.barrier()
	
	def input_query_send_handler(self, param_obj, data_to_send):
		"""
		2008-09-05
			to be called by input_node
			coming from variation/src/MpiQC149CrossMatch.py
			1. query the output_node to get a free computing node
			2. send the data to the computing node
		"""
		self.communicator.send("1", param_obj.output_node_rank, 1)	#WATCH: tag is 1, to the output_node.
		free_computing_node, source, tag = self.communicator.receiveString(param_obj.output_node_rank, 2)
		#WATCH: tag is 2, from the output_node
		data_pickle = cPickle.dumps(data_to_send, -1)
		self.communicator.send(data_pickle, int(free_computing_node),0)	#WATCH: int()
		if param_obj.report:
			sys.stderr.write("block %s sent to %s.\n"%(param_obj.counter, free_computing_node))
		param_obj.counter += 1
	
	def input_fetch_handler(self, param_obj, message_size=1, report=0):
		"""
		2008-09-06
			coming from variation/src/MpiLD.py
			param_obj has 3 variables, params_ls, counter, report(optional)
			it pops a message_size number of entries out of param_obj.params_ls and then return it.
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
		if getattr(param_obj, 'report', None):
			sys.stderr.write("Fetching done at counter=%s.\n"%(param_obj.counter))
		return data_to_return
	
	def output_node(self, free_computing_nodes, parameter_list, handler, type=1):
		"""
		10-20-05
		10-21-05
			handle the situation when free_computing_node list is exhausted
			add type to specify the communicating data type
			based on type, different ways to get stop_signal
		10-21-05
			a common output_node() function
			two jobs:
			1. give node 0 the free_computing_node
			2. handle the data from the computing_nodes
			
			handler(communicator, parameter_list, data)
		"""
		communicator = self.communicator
		report = self.report
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s ready to accept output...\n"%node_rank)
		if type==1:
			data, source, tag = communicator.receiveString(None, 1)
		else:
			data, source, tag, count = communicator.receive(type, None, 1)
		no_of_computing_nodes = len(free_computing_nodes)
		no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
		free_computing_node_request = 0	#10-21-05 Flag used to check whether node 0 is requesting
		request_node = 0	#default request_node is 0
		while 1:
			if source==0:	#10-19-05 the input_node is asking me for free computing_node WATCH: it's array
				if len(free_computing_nodes)==0:	#10-21-05, it's empty. Flag the request.
					free_computing_node_request = 1
					request_node = source	#record the request_node
				else:
					free_computing_node = free_computing_nodes.pop(0)
					communicator.send(str(free_computing_node), source, 2)	#WATCH tag is 2.
			elif type==1 and data=="-1":
				no_of_resting_nodes += 1
				if report:
					sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
				if no_of_resting_nodes==no_of_computing_nodes:	#WATCH: its' size-3
					break
					if report:
						sys.stderr.write("node %s output finished.\n"%node_rank)
			elif type!=1 and data.toscalar()==-1:	#10-21-05 for non-string type
				no_of_resting_nodes += 1
				if report:
					sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
				if no_of_resting_nodes==no_of_computing_nodes:	#WATCH: its' size-3
					break
					if report:
						sys.stderr.write("node %s output finished.\n"%node_rank)
			else:
				if free_computing_node_request==1:	#there's a request for free_computing_node.
					communicator.send(str(source), request_node, 2)	#WATCH tag is 2.
					free_computing_node_request = 0
				else:
					free_computing_nodes.append(source)	#append the free computing_node
				handler(communicator, parameter_list, data)
			if type==1:
				data, source, tag = communicator.receiveString(None, 1)
			else:
				data, source, tag, count = communicator.receive(type, None, 1)
		sys.stderr.write("Node no.%s output done.\n"%node_rank)		
	
	def input_node(self, curs, free_computing_nodes, message_size=None, input_handler=fetch_cluster_block):
		"""
		10-20-05
		10-22-05
			add input_handler and regard curs as parameter_list
		"""
		communicator = self.communicator
		node_rank = communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		data = input_handler(curs, message_size, self.report)
		counter = 0
		while data:
			communicator.send("1", communicator.size-1, 1)	#WATCH: tag is 1, to the output_node.
			free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)
				#WATCH: tag is 2, from the output_node
			data_pickle = cPickle.dumps(data, -1)
			communicator.send(data_pickle, int(free_computing_node), 0)	#WATCH: int()
			if self.report:
				sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
			data = input_handler(curs, message_size, self.report)
			counter += 1
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			if self.report:
				sys.stderr.write("exit signal sent to %s.\n"%(node))
			communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def inputNode(self, param_obj, free_computing_nodes, param_generator=None):
		"""
		2009-2-10
			version similar to input_node() but use a generator (argument param_generator)
		"""
		communicator = self.communicator
		node_rank = communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		counter = 0
		for param in param_generator:
			communicator.send("1", param_obj.output_node_rank, 1)	#WATCH: tag is 1, to the output_node.
			free_computing_node, source, tag = communicator.receiveString(param_obj.output_node_rank, 2)
				#WATCH: tag is 2, from the output_node
			data_pickle = cPickle.dumps(param, -1)
			communicator.send(data_pickle, int(free_computing_node), 0)	#WATCH: int()
			if self.report:
				sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
			counter += 1
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			if self.report:
				sys.stderr.write("exit signal sent to %s.\n"%(node))
			communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_cleanup_handler(cls, communicator):
		"""
		10-22-05
			default cleanup_handler
		"""
		#tell the last node to stop
		communicator.send("-1", communicator.size-1, 1)	#tag is 1
	
	computing_cleanup_handler = classmethod(computing_cleanup_handler)
	
	def computing_node(self, parameter_list, node_fire_handler, cleanup_handler=None):
		"""
		2008-09-08
			delete result and result_pickle right after they are useless to free up memory
		10-21-05
			0 is the node where the data is from
			size-1 is the node where output goes
		10-22-05
			make computing_cleanup_handler to be the default cleanup_handler
		"""
		report = self.report
		communicator = self.communicator
		node_rank = communicator.rank
		data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
		while 1:
			if data=="-1":
				if report:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				result = node_fire_handler(communicator, data, parameter_list)
				result_pickle = cPickle.dumps(result, -1)
				del result
				communicator.send(result_pickle, communicator.size-1, 1)
				del result_pickle
			data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
		if not cleanup_handler:
			MPIwrapper.computing_cleanup_handler(communicator)
