#!/usr/bin/env mpipython
"""
Usage: MpiTrioAncestryInference.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	(db connection by username 'yuhuang')
	-i ...,	input file
	-o ...,	output file
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-q ...,	message_size when transmitting pattern_sig_lists 100,000(default)
	-s ...,	queue size, 8,000,000(default, 1G Mem)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	MpiTrioAncestryInference.py -i justin_data_filtered.csv -o justin_data_filtered.trio_ancestry
	

Description:
	The input StrainSNP matrix is in integer format. Program to seek whether one strain
	is a recombinant of the other two.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, math
import Numeric, cPickle
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, computing_node
from FilterStrainSNPMatrix import FilterStrainSNPMatrix	#for read_data()
from Queue import Queue
from threading import Thread

class trio_generation_thread(Thread):
	"""
	2007-03-15
		this is dead.
	"""
	def __init__(self, communicator, trio_queue, no_of_strains, debug=0, report=0):
		Thread.__init__(self)
		self.communicator = communicator
		self.trio_queue = trio_queue
		self.no_of_strains = int(no_of_strains)
		self.debug = int(debug)
		self.report = int(report)
		
	def run(self):
		if communicator.rank==0:	#only for 1st node
			if self.report:
				sys.stderr.write("node %s thread 0 starts to generate trios...\n"%communicator.rank)
			for i in range(no_of_strains):
				for j in range(i+1, no_of_strains):
					for k in range(j+1, no_of_strains):
						trio_queue.put([i,j,k])
			trio_queue.put(-1)
		else:
			if self.report:
				sys.stderr.write("node %s thread 0 skips trio generation.\n"%communicator.rank)
		if self.report:
			sys.stderr.write("node %s thread 0 done.\n"%communicator.rank)

class MpiTrioAncestryInference:
	def __init__(self, communicator=None, hostname='dl324b-1', dbname='yhdb', schema='dbsnp', \
		input_fname=None, output_fname=None, snp_locus_table='snp_locus', message_size=100000, debug=0, report=0):
		"""
		2007-03-07
		"""
		self.communicator = communicator
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.snp_locus_table = snp_locus_table
		self.message_size = int(message_size)
		self.debug = int(debug)
		self.report = int(report)
	
	def input_node(self, communicator, parameter_list, free_computing_nodes, message_size, report=0):
		"""
		2007-03-09
			the queue is way too slow
		"""
		node_rank = communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		no_of_strains = parameter_list[0]
		counter = 0
		trio_list = []
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				for k in range(j+1, no_of_strains):
					trio_list.append([i,j,k])
					if len(trio_list)==message_size:
						communicator.send("1", communicator.size-1, 1)	#WATCH: tag is 1, to the output_node.
						free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)
						#WATCH: tag is 2, from the output_node
						data_pickle = cPickle.dumps(trio_list, -1)
						communicator.send(data_pickle,int(free_computing_node),0)	#WATCH: int()
						trio_list = []	#clear the list
						if report:
							sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
						counter += 1
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		2007-03-07
			defunct
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		trio_queue = parameter_list[0]
		block = []
		trio = trio_queue.get()
		while trio!=-1:
			block.append(trio)
			if len(block)>=message_size:
				break
			trio = trio_queue.get()
		if report:
			sys.stderr.write("Fetching done.\n")
		return block
	
	def get_chr_start_ls(self, curs, snp_acc_list, snp_locus_table='snp_locus'):
		"""
		2007-03-07
			snp_acc_list is already in order
		"""
		sys.stderr.write("Getting chr_start_ls ...")
		chr_start_ls = [0]
		old_chromosome = -1
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome from %s where acc='%s'"%(snp_locus_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome = rows[0][0]
			if old_chromosome == -1:	#1st encounter
				old_chromosome = chromosome
			elif chromosome!=old_chromosome:
				chr_start_ls.append(i)
				old_chromosome = chromosome
		chr_start_ls.append(len(snp_acc_list))	#the last one as a stop
		sys.stderr.write("Done.\n")
		return chr_start_ls
	
	def computing_node_handler(self, communicator, data, parameter_list):
		"""
		2007-03-07
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data_matrix, chr_start_ls, trio_arrangement_ls = parameter_list
		data = cPickle.loads(data)
		no_of_successes = 0
		no_of_trials = 0
		result = []
		for trio in data:
			for trio_arrangement in trio_arrangement_ls:
				row1_index = trio[trio_arrangement[0]]
				row2_index = trio[trio_arrangement[1]]
				row3_index = trio[trio_arrangement[2]]
				ancestry_ls, no_of_jumps = self.identify_ancestry_with_min_jumps(data_matrix[row1_index], data_matrix[row2_index], data_matrix[row3_index], chr_start_ls)
				no_of_trials += 1
				if ancestry_ls:
					no_of_successes += 1
					result.append([row1_index, row2_index, row3_index] + ancestry_ls + [no_of_jumps])
		sys.stderr.write("Node no.%s done with %s successes out of %s trials.\n"%(node_rank, no_of_successes, no_of_trials))
		return result
	
	def initialize__score_trace_matrix(self, row1, row2, row3, chr_start_ls):
		score_matrix = Numeric.zeros([2, len(row1)])
		trace_matrix = Numeric.zeros([2, len(row1)])
		for l in chr_start_ls[:-1]:
			no_of_incompatibles = 0
			if row1[l]==row3[l] or row1[l]==0 or row3[l]==0:
				score_matrix[0,l] = 0
			else:
				score_matrix[0,l] = -1
				no_of_incompatibles += 1
			if row2[l]==row3[l] or row2[l]==0 or row3[l]==0:
				score_matrix[1,l] = 0
			else:
				score_matrix[1,l] = -1
				no_of_incompatibles += 1
			trace_matrix[0,l] = -1
			trace_matrix[1,l] = -1
			if no_of_incompatibles == 2:	#this is a bad spot, the whole alignment breaks down
				score_matrix = [-1]
				trace_matrix = [-1]
				break
		return score_matrix, trace_matrix
	
	def identify_ancestry_of_one_chr_with_DP(self, row1, row2, row3, score_matrix, trace_matrix, chr_start, next_chr_start):
		"""
		2007-03-07
		"""
		is_identified = 1
		no_of_rows, no_of_cols = score_matrix.shape
		if self.debug:
			import pdb
			pdb.set_trace()
		for col_index in range(chr_start+1, next_chr_start):	#the index goes from chr_start+1 to next_chr_start-1
			valid_cells = []
			if row1[col_index]==row3[col_index] or row1[col_index]==0 or row3[col_index]==0:
				valid_cells.append(0)
			else:
				score_matrix[0, col_index] = col_index + 2	#set it maximum no of jumps
				trace_matrix[0, col_index] = -1	#-1 means incompatible
			if row2[col_index]==row3[col_index] or row2[col_index]==0 or row3[col_index]==0:
				valid_cells.append(1)
			else:
				score_matrix[1, col_index] = col_index + 2
				trace_matrix[1, col_index] = -1
			if valid_cells:	#this position has compatible snp type
				for row_index in valid_cells:
					min_no_of_jumps = col_index + 2	#this is like the maximum no of jumps + 1
					for prev_row_index in range(no_of_rows):	#check previous spot
						if score_matrix[prev_row_index, col_index-1]!=-1:
							candidate_min_no_of_jumps = score_matrix[prev_row_index, col_index-1] + 1 - int(row_index==prev_row_index)
							if candidate_min_no_of_jumps < min_no_of_jumps:
								min_no_of_jumps = candidate_min_no_of_jumps
								jump_src = prev_row_index
					score_matrix[row_index, col_index] = min_no_of_jumps
					trace_matrix[row_index, col_index] = jump_src
			else:
				is_identified = 0
				break
		return is_identified
	
	def trace(self, score_matrix, trace_matrix, chr_start_ls):
		"""
		2007-03-07
			do it chromosome by chromosome
		"""
		no_of_jumps = 0
		ancestry_ls = []
		if self.debug:
			import pdb
			pdb.set_trace()
		for i in range(len(chr_start_ls)-1):
			chr_stop = chr_start_ls[i+1] - 1
			if score_matrix[0, chr_stop]<score_matrix[1, chr_stop]:
				end_chr_ancestry = 0
			else:
				end_chr_ancestry = 1
			no_of_jumps += score_matrix[end_chr_ancestry, chr_stop]
			ancestry_ls += self.recursive_trace(trace_matrix, end_chr_ancestry, chr_stop, chr_start_ls[i])
			ancestry_ls.append(end_chr_ancestry)	#the last chr ancestry
			ancestry_ls.append('||')	#the separator
		return ancestry_ls, no_of_jumps
	
	def recursive_trace(self, trace_matrix, i, j, chr_start):
		"""
		2007-03-07
		"""
		chr_ancestry_ls = []
		if j > chr_start:
			prev_chr_ancestry = trace_matrix[i,j]
			chr_ancestry_ls = self.recursive_trace(trace_matrix, prev_chr_ancestry, j-1, chr_start)
			chr_ancestry_ls.append(prev_chr_ancestry)
		return chr_ancestry_ls
	
	def identify_ancestry_with_min_jumps(self, row1, row2, row3, chr_start_ls):
		ancestry_ls = []
		no_of_jumps = 0
		score_matrix, trace_matrix = self.initialize__score_trace_matrix(row1, row2, row3, chr_start_ls)
		if len(score_matrix) == 1:
			return ancestry_ls, no_of_jumps
		
		for i in range(len(chr_start_ls)-1):
			is_identified = self.identify_ancestry_of_one_chr_with_DP(row1, row2, row3, score_matrix, trace_matrix, chr_start_ls[i], chr_start_ls[i+1])
			if not is_identified:
				break
		if is_identified:
			ancestry_ls, no_of_jumps = self.trace(score_matrix, trace_matrix, chr_start_ls)
		return ancestry_ls, no_of_jumps
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		2007-03-08
		"""
		writer = parameter_list[0]
		prediction_ls = cPickle.loads(data)
		for row in prediction_ls:
			writer.writerow(row)
	
	def run(self):
		"""
		2007-03-08
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
							--initialize__score_trace_matrix()
							(for loop)
								--identify_ancestry_of_one_chr_with_DP()
							--trace()
								--recursive_trace()
			else:
				--output_node()
					--output_node_handler()
		"""
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		if node_rank == 0:
			FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
			snp_acc_list = header[2:]
			data_matrix = Numeric.array(data_matrix)
			no_of_strains = data_matrix.shape[0]
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema, password='123456', user='yuhuang')
			chr_start_ls = self.get_chr_start_ls(curs, snp_acc_list, self.snp_locus_table)
			
			chr_start_ls_pickle = cPickle.dumps(chr_start_ls, -1)	#-1 means use the highest protocol
			data_matrix_pickle = cPickle.dumps(data_matrix, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				self.communicator.send(chr_start_ls_pickle, node, 0)
				self.communicator.send(data_matrix_pickle, node, 0)
		elif node_rank in free_computing_nodes:
			data, source, tag = self.communicator.receiveString(0, 0)
			chr_start_ls = cPickle.loads(data)	#take the data
			data, source, tag = self.communicator.receiveString(0, 0)
			data_matrix = cPickle.loads(data)
		
		mpi_synchronize(self.communicator)
		
		if node_rank == 0:
			parameter_list = [no_of_strains]
			self.input_node(self.communicator, parameter_list, free_computing_nodes, self.message_size, \
				self.report)
		elif node_rank in free_computing_nodes:
			trio_arrangement_ls = [[0,1,2], [1,2,0], [2,0,1]]	#three different ways to pick the parent-set and the child
			parameter_list = [data_matrix, chr_start_ls, trio_arrangement_ls]
			computing_node(self.communicator, parameter_list, self.computing_node_handler, report=self.report)
		else:
			writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
			parameter_list = [writer]
			output_node(self.communicator, free_computing_nodes, parameter_list, self.output_node_handler, self.report)
			del writer
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:n:q:s:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'dbsnp'
	input_fname = None
	output_fname = None
	snp_locus_table = 'snp_locus'
	message_size = 100000
	queue_size = 8000000
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-q",):
			message_size = int(arg)
		elif opt in ("-s",):
			queue_size = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname and hostname and dbname and schema:
		communicator = MPI.world.duplicate()
		instance = MpiTrioAncestryInference(communicator, hostname, dbname, schema, \
			input_fname, output_fname, snp_locus_table, message_size, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)