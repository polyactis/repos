#!/usr/bin/env python
"""
Usage: RemoveBadSNPs.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-i ...,	input file
	-o ...,	output file
	-c ...,	log probability cutoff, -0.5(default)
	-a ...,	bits to toggle input/output StrainSNP matrix format, 0 is integer.
		1 is alphabet. 1st digit is for input. 2nd is for output. 00 (default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	RemoveBadSNPs.py -i justin_data.csv -o justin_data_filtered.csv
	
Description:
	The input StrainSNP matrix is in integer format.
	Remove SNPs with too low probability.
	For each SNP locus, its homo-hetero pattern follows a binomial-like distribution.
	calculating probability for each snp locus, lower probability, worse the snp locus

"""
from __init__ import *

class RemoveBadSNPs:
	def __init__(self, input_fname=None, output_fname=None, min_log_prob=-0.5, nt_alphabet_bits='00', debug=0, report=0):
		"""
		2007-04-10
		2007-05-14
			add nt_alphabet_bits
		"""
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.min_log_prob = float(min_log_prob)
		self.nt_alphabet_bits = nt_alphabet_bits
		self.debug = int(debug)
		self.report = int(report)
	
	def cal_strain_homo_perc_vector(self, data_matrix):
		sys.stderr.write("Calculating strain_homo_perc_vector...")
		no_of_strains, no_of_snps = data_matrix.shape
		strain_homo_perc_vector = num.zeros(no_of_strains, num.float)
		for i in range(no_of_strains):
			no_of_valid_calls = 0.0
			no_of_homo_calls = 0.0
			for j in range(no_of_snps):
				if data_matrix[i][j] != 0:
					no_of_valid_calls += 1
					if data_matrix[i][j] <=4:
						no_of_homo_calls +=1
			strain_homo_perc_vector[i] = no_of_homo_calls/no_of_valid_calls
		sys.stderr.write("Done\n")
		return strain_homo_perc_vector
	
	def cal_snp_locus_log_prob(self, data_matrix, strain_homo_perc_vector):
		"""
		2007-07-16
			add an exception handler when no_of_valid_calls=0
			if no_of_valid_calls==0, then all snp calls of each strain within the population are homozygous.
			or all are heterozygous. (for arabidopsis, the former mostly)
		"""
		sys.stderr.write("Calculating snp_locus_log_prob ...")
		no_of_strains, no_of_snps = data_matrix.shape
		snp_locus_log_prob = num.zeros(no_of_snps, num.float)
		for j in range(no_of_snps):
			no_of_valid_calls = 0.0
			for i in range(no_of_strains):
				if data_matrix[i,j]!=0 and strain_homo_perc_vector[i]!=0 and strain_homo_perc_vector[i]!=1:
					no_of_valid_calls += 1
					is_call_homo = int(data_matrix[i,j]<=4)
					snp_locus_log_prob[j] += is_call_homo*math.log(strain_homo_perc_vector[i])+(1-is_call_homo)*math.log(1-strain_homo_perc_vector[i])
			if no_of_valid_calls!=0:
				snp_locus_log_prob[j] /= no_of_valid_calls
		sys.stderr.write("Done\n")
		return snp_locus_log_prob
	
	def run(self):
		"""
		2007-04-30
		2007-05-14
			add nt_alphabet_bits
		"""
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, int(self.nt_alphabet_bits[0]))
		data_matrix = num.array(data_matrix)
		strain_homo_perc_vector = self.cal_strain_homo_perc_vector(data_matrix)
		snp_locus_log_prob = self.cal_snp_locus_log_prob(data_matrix, strain_homo_perc_vector)
		from sets import Set
		cols_to_be_tossed_out_set = Set()
		for i in range(len(snp_locus_log_prob)):
			if snp_locus_log_prob[i]<=min_log_prob:
				cols_to_be_tossed_out_set.add(i)
		print "%sSNPs removed:"%(len(cols_to_be_tossed_out_set))
		for col_index in cols_to_be_tossed_out_set:
			print '\t%s\t%s'%(col_index, header[2+col_index])
		write_data_matrix(data_matrix, self.output_fname, header, strain_acc_list, category_list, cols_to_be_tossed_out=cols_to_be_tossed_out_set, nt_alphabet=int(self.nt_alphabet_bits[1]))
		import pylab
		pylab.title("histogram of snp locus log probability")
		pylab.hist(snp_locus_log_prob, 20)
		pylab.show()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:c:a:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_fname = None
	output_fname = None
	min_log_prob = -0.5
	nt_alphabet_bits = '00'
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
		elif opt in ("-c",):
			min_log_prob = float(arg)
		elif opt in ("-a",):
			nt_alphabet_bits = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname:
		instance = RemoveBadSNPs(input_fname, output_fname, min_log_prob, nt_alphabet_bits,\
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)