#!/usr/bin/env python
"""
Usage: SelectStrains.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-i ...,	input file
	-o ...,	output file
	-n ...,	no_of_strains_to_be_removed, 300(default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	SelectStrains.py -i justin_data.csv -o justin_data_selected.csv -r
	

Description:
	Program to select strains based on distance matrix.
	The distance measure is hamming distance. Also conditional entropy measure is also hidden.
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
import Numeric as num
from sets import Set

class SelectStrains:
	def __init__(self, input_fname=None, output_fname=None, no_of_strains_to_be_removed=300, \
		debug=0, report=0):
		"""
		2007-02-27
		"""
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.no_of_strains_to_be_removed = int(no_of_strains_to_be_removed)
		self.debug = int(debug)
		self.report = int(report)
	
	def calculate_hamming_distance(self, row1, row2):
		no_of_valid_entries = 0.0
		no_of_different_entries = 0.0
		for i in range(len(row1)):
			if row1[i]!=0 and row2[i]!=0:	#ignore the NA
				no_of_valid_entries += 1
				if row1[i]!=row2[i]:
					no_of_different_entries += 1
		if no_of_valid_entries:
			hamming_dist =  no_of_different_entries/no_of_valid_entries
		else:
			hamming_dist = 0
		return hamming_dist
	
	def get_conditional_entropy(self, row1, row2):
		two_by_two_table = num.zeros((2,2), num.Float)
		for i in range(len(row1)):
			if row1[i]!=0 and row2[i]!=0:	#ignore the NA
				two_by_two_table[row1[i]-1, row2[i]-1] += 1	#-1 due to the fact that 0 is reserved for NA
		row1_marginal_sum = num.sum(two_by_two_table, 1)
		total_sum = num.sum(row1_marginal_sum)
		H = 0
		for i in range(two_by_two_table.shape[0]):
			for j in range(two_by_two_table.shape[1]):
				if two_by_two_table[i,j]:	#if this is 0, row1_marginal_sum entry is also 0
					cnd_prob = two_by_two_table[i,j]/row1_marginal_sum[i]
				else:
					cnd_prob = 0
				if cnd_prob!=0 and cnd_prob!=1:
					H += -two_by_two_table[i,j]*math.log(cnd_prob)
		H /= total_sum
		return H
	
	def compute_conditional_entropy_matrix(self, data_matrix, upper_triangle_only=1):
		sys.stderr.write("Computing conditional entropy matrix ... ")
		no_of_strains = data_matrix.shape[0]
		cnd_entropy_matrix = num.zeros([no_of_strains, no_of_strains], num.Float)
		for i in range(no_of_strains):
			cnd_entropy_matrix[i,i] = 0
			if upper_triangle_only:	#symmetric distance matrix
				for j in range(i+1, no_of_strains):
					cnd_entropy_matrix[j,i] = cnd_entropy_matrix[i,j] = self.calculate_hamming_distance(data_matrix[i], data_matrix[j])
			else:	#non-symmetric
				for j in range(no_of_strains):	#j conditioned on i
					cnd_entropy_matrix[i,j] = self.get_conditional_entropy(data_matrix[i], data_matrix[j])
		sys.stderr.write("done.\n")
		return cnd_entropy_matrix
	
	def recursive_remove_redundant_strains(self, cnd_entropy_matrix, cnd_entropy_sum_vector, no_of_strains_to_be_removed):
		strain_to_be_deleted_index_list = []
		if no_of_strains_to_be_removed>0:
			strain_to_be_deleted = num.argmin(cnd_entropy_sum_vector)
			strain_to_be_deleted_index_list.append(strain_to_be_deleted)
			for i in range(cnd_entropy_sum_vector.shape[0]):	#remove 
				cnd_entropy_sum_vector[i] -= cnd_entropy_matrix[strain_to_be_deleted,i]
			cnd_entropy_sum_vector[strain_to_be_deleted] = cnd_entropy_matrix.shape[0]	#set this deleted strain's cond entropy to maximum
			no_of_strains_to_be_removed -= 1
			if self.debug:
				print
				print 'strain_to_be_deleted:', strain_to_be_deleted
				print 'strain_to_be_deleted_index_list'
				print strain_to_be_deleted_index_list
				print 'cnd_entropy_sum_vector'
				print cnd_entropy_sum_vector
			strain_to_be_deleted_index_list += self.recursive_remove_redundant_strains(cnd_entropy_matrix, cnd_entropy_sum_vector, no_of_strains_to_be_removed)
		return strain_to_be_deleted_index_list
	
	def remove_redundant_strains(self, cnd_entropy_matrix, no_of_strains_to_be_removed=100):
		sys.stderr.write("Removing redundant strains ... ")
		cnd_entropy_sum_vector = num.sum(cnd_entropy_matrix)	#in axis 0, sum all different entries along axis 0(column by column)
		if self.debug:
			print
			print 'cnd_entropy_matrix'
			print cnd_entropy_matrix
			print 'cnd_entropy_sum_vector'
			print cnd_entropy_sum_vector
		strain_to_be_deleted_index_list = self.recursive_remove_redundant_strains(cnd_entropy_matrix, cnd_entropy_sum_vector, no_of_strains_to_be_removed)
		sys.stderr.write("done.\n")
		return strain_to_be_deleted_index_list
	
	
	def convert_matrix_to_binary_data_matrix(self, data_matrix):
		"""
		2007-02-26
		"""
		sys.stderr.write("Converting data matrix into binary form ... ")
		no_of_columns = len(data_matrix[0])
		no_of_rows = len(data_matrix)
		new_data_matrix = num.zeros([no_of_rows, no_of_columns])
		for i in range(no_of_columns):	#column by column
			SNP2index = {}	#start from 1, 0 is reserved for NA
			for j in range(no_of_rows):
				entry = data_matrix[j][i]
				if entry!='0':
					if entry not in SNP2index:
						SNP2index[entry] = len(SNP2index)+1
					new_data_matrix[j,i] = SNP2index[entry]
		if self.debug:
			print
			print 'new_data_matrix'
			print new_data_matrix
		sys.stderr.write("Done.\n")
		return new_data_matrix
	
	def read_data(self, input_fname):
		sys.stderr.write("Reading data ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		data_matrix = []
		strain_acc_list = []
		category_list = []
		for row in reader:
			strain_acc_list.append(row[0])
			category_list.append(row[1])
			data_matrix.append(row[2:])
		del reader
		sys.stderr.write("Done.\n")
		return header, strain_acc_list, category_list, data_matrix
	
	def write_data_matrix(self, data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=Set()):
		"""
		2007-02-27
		"""
		sys.stderr.write("Writing data_matrix ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(header)
		for i in range(len(data_matrix)):
			if i not in rows_to_be_tossed_out:
				new_row = [strain_acc_list[i], category_list[i]]
				new_row += list(data_matrix[i])
				writer.writerow(new_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-02-27
		"""
		header, strain_acc_list, category_list, data_matrix = self.read_data(self.input_fname)
		new_data_matrix = self.convert_matrix_to_binary_data_matrix(data_matrix)
		cnd_entropy_matrix = self.compute_conditional_entropy_matrix(new_data_matrix)
		strain_to_be_deleted_index_list = self.remove_redundant_strains(cnd_entropy_matrix, self.no_of_strains_to_be_removed)
		
		self.write_data_matrix(data_matrix, self.output_fname, header, strain_acc_list, category_list, Set(strain_to_be_deleted_index_list))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:n:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_fname = None
	output_fname = None
	no_of_strains_to_be_removed = 300
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-n",):
			no_of_strains_to_be_removed = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname:
		instance = SelectStrains(input_fname, output_fname, no_of_strains_to_be_removed, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)	
