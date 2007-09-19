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
	The input StrainSNP matrix is in integer format (either homozygous or heterozygous).
	Program to seek whether one strain is a recombinant of the other two (use dynamic
	programming to determine the minimum number of jumps=recombinations).
	
	Note: if it happens to be two consecutive heterozygous SNPs in the child, no jump is counted
	between them as it's not easy to decide whether recombination occurs or not.
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
from common import nt2number,number2nt, nt_number_matching_matrix
nt_number_matching_matrix = Numeric.array(nt_number_matching_matrix)

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
	
	def is_child_heterozygous_SNP_compatible_with_parents(self, SNP_parent_1, SNP_parent_2, SNP_child):
		"""
		2007-04-16
			SNP_child>4
		"""
		if SNP_child<5:	#homo or NA
			return -1
		SNP_child_1 = nt2number[number2nt[SNP_child][0]]
		SNP_child_2 = nt2number[number2nt[SNP_child][1]]
		if (nt_number_matching_matrix[SNP_child_1, SNP_parent_1] and nt_number_matching_matrix[SNP_child_2, SNP_parent_2]) or\
			(nt_number_matching_matrix[SNP_child_1, SNP_parent_2] and nt_number_matching_matrix[SNP_child_2, SNP_parent_1]):
			return 1
		else:
			return 0
	
	def initialize_score_trace_matrix(self, row1, row2, row3, chr_start_ls):
		"""
		2007-04-16
			score_matrix(i,j) records the number of minimum jumps(recombination events) from the chromosome starting or the last heterozygous call of the child while the ancestry of j-th child SNP is parent i(0 or 1).
				-1 means incompatible
			trace_matrix(i,j) records the ancestry of the (j-1)-th child SNP assuming the ancestry of j-th child SNP is parent i.
				-2 means the starting position
				-1 means incompatible
				0, 1 means the parent index
				2 means the child SNP is heterozygous and comes from either parent
		2007-04-16
			deal with heterozygous calls
		"""
		score_matrix = Numeric.zeros([2, len(row1)])
		trace_matrix = Numeric.zeros([2, len(row1)])
		for l in chr_start_ls[:-1]:
			no_of_incompatibles = 0
			if row3[l]>4:	#heterozygous
				if self.is_child_heterozygous_SNP_compatible_with_parents(row1[l], row2[l], row3[l])==1:
					score_matrix[0,l] = 0
					score_matrix[1,l] = 0
					trace_matrix[0, l] = 2
					trace_matrix[1, l] = 2
				else:
					score_matrix[0,l] = -1
					score_matrix[1,l] = -1
					trace_matrix[0, l] = -1
					trace_matrix[1, l] = -1
					no_of_incompatibles += 2
			else:
				if nt_number_matching_matrix[row1[l], row3[l]]:
					score_matrix[0,l] = 0
					trace_matrix[0,l] = -2
				else:
					score_matrix[0,l] = l+2	#maximum possible # of jumps +1
					trace_matrix[0,l] = -1
					no_of_incompatibles += 1
				if nt_number_matching_matrix[row2[l], row3[l]]:
					score_matrix[1,l] = 0
					trace_matrix[1,l] = -2
				else:
					score_matrix[1,l] = l+2	#maximum possible # of jumps +1
					trace_matrix[1,l] = -1
					no_of_incompatibles += 1
			if no_of_incompatibles == 2:	#this is a bad spot, the whole alignment breaks down
				score_matrix = [-1]
				trace_matrix = [-1]
				break
		return score_matrix, trace_matrix
	
	def identify_ancestry_of_one_chr_with_DP(self, row1, row2, row3, score_matrix, trace_matrix, chr_start, next_chr_start):
		"""
		2007-03-07
		2007-04-16
			modify to deal with heterozygous SNPs
		"""
		is_identified = 1
		fake_chr_start_ls = []
		no_of_rows, no_of_cols = score_matrix.shape
		if self.debug:
			import pdb
			pdb.set_trace()
		for col_index in range(chr_start+1, next_chr_start):	#the index goes from chr_start+1 to next_chr_start-1
			valid_cells = []
			if row3[col_index]>4:	#2007-04-16 the child SNP is heterozygous
				if self.is_child_heterozygous_SNP_compatible_with_parents(row1[col_index], row2[col_index], row3[col_index])==1:
					score_matrix[0, col_index] = 0
					score_matrix[1, col_index] = 0
					trace_matrix[0, col_index] = 2
					trace_matrix[1, col_index] = 2
					fake_chr_start_ls.append(col_index)
				else:
					is_identified = 0
					break
			else:	#the child is homozygous
				if nt_number_matching_matrix[row1[col_index], row3[col_index]]:
					valid_cells.append(0)
				else:
					score_matrix[0, col_index] = col_index+2	#maximum possible # of jumps +1
					trace_matrix[0, col_index] = -1	#-1 means incompatible
				if nt_number_matching_matrix[row2[col_index], row3[col_index]]:
					valid_cells.append(1)
				else:
					score_matrix[1, col_index] = col_index+2	#maximum possible # of jumps +1
					trace_matrix[1, col_index] = -1
				if valid_cells:	#this position has compatible snp type
					for row_index in valid_cells:
						min_no_of_jumps = col_index + 2	#this is like the maximum no of jumps + 1
						for prev_row_index in range(no_of_rows):	#check previous spot
							if trace_matrix[prev_row_index, col_index-1]!=-1:	#the previous position is compatible
								if trace_matrix[0, col_index-1]==2 and trace_matrix[1, col_index-1]==2:
									#the previous position is heterozygous
									candidate_min_no_of_jumps = 1	#homozygous after heterozygous is always 1 jump
								else:
									candidate_min_no_of_jumps = score_matrix[prev_row_index, col_index-1] + 1 - int(row_index==prev_row_index)
								if candidate_min_no_of_jumps < min_no_of_jumps:
									min_no_of_jumps = candidate_min_no_of_jumps
									jump_src = prev_row_index
						score_matrix[row_index, col_index] = min_no_of_jumps
						trace_matrix[row_index, col_index] = jump_src
				else:
					is_identified = 0
					break
		return is_identified, fake_chr_start_ls
	
	def trace(self, score_matrix, trace_matrix, chr_start_ls, flag_of_real_chr_start_ls):
		"""
		2007-03-07
			do it chromosome by chromosome
		2007-04-16
			deal with heterozygous (fake chromosome starting)
			add flag_of_real_chr_start_ls
		"""
		no_of_jumps = 0
		ancestry_ls = []
		if self.debug:
			import pdb
			pdb.set_trace()
		for i in range(len(chr_start_ls)-1):
			chr_start = chr_start_ls[i]
			chr_stop = chr_start_ls[i+1] - 1
			if chr_start==chr_stop:
				if trace_matrix[0, chr_start]==trace_matrix[1, chr_start]==2:
					ancestry_ls.append(2)
					if trace_matrix[0, chr_start-1]!=2 and trace_matrix[1, chr_start-1]!=2 and flag_of_real_chr_start_ls[i]==0:	#the previous position is not heterozygous and this is fake chromosome starting
						no_of_jumps += 1
				elif score_matrix[0, chr_stop]<score_matrix[1, chr_stop]:	#this is real chromosome start
					ancestry_ls.append(0)
				else:
					ancestry_ls.append(1)
			else:
				if score_matrix[0, chr_stop]<score_matrix[1, chr_stop]:
					end_chr_ancestry = 0
				else:
					end_chr_ancestry = 1
				no_of_jumps += score_matrix[end_chr_ancestry, chr_stop]
				chr_ancestry_ls = self.recursive_trace(trace_matrix, end_chr_ancestry, chr_stop, chr_start)
				if trace_matrix[0, chr_start]==trace_matrix[1, chr_start]==2:	#it's heterozygous
					ancestry_ls.append(2)
					ancestry_ls += chr_ancestry_ls[1:]
					if trace_matrix[0, chr_start-1]!=2 and trace_matrix[1, chr_start-1]!=2 and flag_of_real_chr_start_ls[i]==0:	#the previous position is not heterozygous and this is fake chromosome starting
						no_of_jumps += 1
				else:
					ancestry_ls += chr_ancestry_ls
				ancestry_ls.append(end_chr_ancestry)	#the last chr ancestry
			if flag_of_real_chr_start_ls[i+1]==1:	#if the next chromosome starting position is real chromosome starting
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
		"""
		2007-04-16
			regard heterozygous SNPs as fake chromosome starting position
		"""
		ancestry_ls = []
		no_of_jumps = 0
		score_matrix, trace_matrix = self.initialize_score_trace_matrix(row1, row2, row3, chr_start_ls)
		if len(score_matrix) == 1:
			return ancestry_ls, no_of_jumps
		fake_chr_start_ls = []
		flag_of_real_chr_start_ls = []	#flag corresponding to each position in fake_chr_start_ls
		for i in range(len(chr_start_ls)-1):
			is_identified, fake_chr_start_ls_of_this_chr = self.identify_ancestry_of_one_chr_with_DP(row1, row2, row3, score_matrix, trace_matrix, chr_start_ls[i], chr_start_ls[i+1])
			if is_identified:	#2007-04-16
				fake_chr_start_ls.append(chr_start_ls[i])
				flag_of_real_chr_start_ls.append(1)	#this is real
				for fake_chr_start in fake_chr_start_ls_of_this_chr:
					fake_chr_start_ls.append(fake_chr_start)
					flag_of_real_chr_start_ls.append(0)
			else:
				break
		if is_identified:
			fake_chr_start_ls.append(chr_start_ls[-1])	#2007-04-16 append the last chromosome starting position
			flag_of_real_chr_start_ls.append(0)	#the last one is placeholder, not real chromosome start
			ancestry_ls, no_of_jumps = self.trace(score_matrix, trace_matrix, fake_chr_start_ls, flag_of_real_chr_start_ls)
		return ancestry_ls, no_of_jumps
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		2007-03-08
		2007-09-17
			add strain_acc_list into the parameter_list
		"""
		writer, strain_acc_list = parameter_list
		prediction_ls = cPickle.loads(data)
		for row in prediction_ls:
			for i in range(3):
				row[i] = strain_acc_list[row[i]]
			writer.writerow(row)
	
	def run(self):
		"""
		2007-04-16
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
		node_rank = self.communicator.rank
		free_computing_nodes = range(1, self.communicator.size-1)	#exclude the 1st and last node
		if node_rank == 0:
			FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
			header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
			snp_acc_list = header[2:]
			data_matrix = Numeric.array(data_matrix)
			no_of_strains = data_matrix.shape[0]
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema, password='123456', user='yuhuang')
			
			#2007-09-17 send strain_acc_list to the output_node
			strain_acc_list_pickle = cPickle.dumps(strain_acc_list, -1)
			self.communicator.send(strain_acc_list_pickle, self.communicator.size-1, 0)
			
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
		else:
			data, source, tag = self.communicator.receiveString(0, 0)
			strain_acc_list = cPickle.loads(data)
		
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
			parameter_list = [writer, strain_acc_list]
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