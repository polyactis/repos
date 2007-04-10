#!/usr/bin/env python
"""
Usage: EstimateSelfingGeneration.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-i ...,	input file
	-o ...,	output file
	-m ...,	max_selfing_generation, 20(default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/EstimateSelfingGeneration.py -i data/justin_data_filtered_y.csv -o data/justin_data_filtered_y.estimate.selfing.generation.csv

Description:
	Estimate the number of selfing generations since the last outcrossing event.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv
import Numeric, copy, pylab
from common import nt2number, number2nt

class EstimateSelfingGeneration:
	def __init__(self, input_fname=None, output_fname=None, max_selfing_generation=20, debug=0, report=0):
		"""
		2007-04-10
		"""
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.max_selfing_generation = int(max_selfing_generation)
		self.debug = int(debug)
		self.report = int(report)
	
	def cal_locus_allele_prob_vector(self, data_matrix):
		sys.stderr.write("Calculating locus allele probability ...")
		no_of_strains, no_of_snps = data_matrix.shape
		locus_allele_prob_vector = Numeric.zeros([no_of_snps, 2], Numeric.Float)
		for i in range(no_of_snps):
			nt2counter = {}
			no_of_valid_calls = 0.0
			for j in range(no_of_strains):
				if data_matrix[j,i]!=0:	#not NA
					no_of_valid_calls += 2	#watch this is 2, cuz it's diploid
					nt_string = number2nt[data_matrix[j,i]]
					if len(nt_string)==1:	#double the nt if it's homozygous
						nt_string += nt_string
					for k in range(len(nt_string)):
						nt = nt_string[k]
						if nt not in nt2counter:
							nt2counter[nt] = 0
						nt2counter[nt] += 1
			nt_key_ls = nt2counter.keys()
			for j in range(len(nt_key_ls)):
				locus_allele_prob_vector[i,j] = nt2counter[nt_key_ls[j]]/no_of_valid_calls
		sys.stderr.write("Done.\n")
		return locus_allele_prob_vector
	
	def cal_locus_heterozygous_prob_vector(self, locus_allele_prob_vector):
		sys.stderr.write("Calculating locus heterozygous probability ...")
		no_of_snps, no_of_alleles = locus_allele_prob_vector.shape
		locus_heterozygous_prob_vector = Numeric.zeros(no_of_snps, Numeric.Float)
		if self.debug:
			import pdb
			pdb.set_trace()
		for i in range(no_of_snps):
			homo_prob = 0.0
			for j in range(no_of_alleles):
				homo_prob += locus_allele_prob_vector[i,j]*locus_allele_prob_vector[i,j]
			locus_heterozygous_prob_vector[i] = 1-homo_prob
		sys.stderr.write("Done.\n")
		return locus_heterozygous_prob_vector
	
	def cal_locus_heterozygous_prob_matrix(self, locus_heterozygous_prob_vector, max_selfing_generation):
		sys.stderr.write("Calculating locus heterozygous probability matrix ...")
		locus_heterozygous_prob_matrix = Numeric.zeros([max_selfing_generation, len(locus_heterozygous_prob_vector)], Numeric.Float)
		locus_heterozygous_prob_matrix[0] = locus_heterozygous_prob_vector
		lhp_vector = copy.deepcopy(locus_heterozygous_prob_vector)
		if self.debug:
			import pdb
			pdb.set_trace()
		for i in range(1, max_selfing_generation):
			lhp_vector = lhp_vector/2.0	#every selfing generation halves the heterozygous probability
			locus_heterozygous_prob_matrix[i] = lhp_vector
		sys.stderr.write("Done.\n")
		return locus_heterozygous_prob_matrix
	
	def cal_individual_selfing_generation_prob_vector(self, data_vector, locus_heterozygous_prob_vector, locus_heterozygous_prob_matrix):
		max_selfing_generation, no_of_snps = locus_heterozygous_prob_matrix.shape
		likelihood_vector = Numeric.ones(max_selfing_generation, Numeric.Float)
		for i in range(max_selfing_generation):
			for j in range(no_of_snps):
				if data_vector[j] != 0:	#not NA
					heterozygous_bit = int(data_vector[j]>4)	#over 4 is heterozygous
					if heterozygous_bit==1:
						locus_likelihood = locus_heterozygous_prob_matrix[i,j]
					else:
						locus_likelihood = 1-locus_heterozygous_prob_matrix[i,j]
					likelihood_vector[i] *= locus_likelihood
		return likelihood_vector
	
	def cal_selfing_generation_prob(self, data_matrix, locus_heterozygous_prob_vector, strain_acc_list, category_list, locus_heterozygous_prob_matrix, output_fname):
		sys.stderr.write("Selfing generation probability...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		no_of_strains, no_of_snps = data_matrix.shape
		selfing_generation_ls = []
		for i in range(no_of_strains):
			if self.debug:
				import pdb
				pdb.set_trace()
			selfing_generation_prob_vector = self.cal_individual_selfing_generation_prob_vector(data_matrix[i], locus_heterozygous_prob_vector, locus_heterozygous_prob_matrix)
			mle_selfing_generation = Numeric.argmax(selfing_generation_prob_vector)
			mle_selfing_generation_prob = selfing_generation_prob_vector[mle_selfing_generation]
			writer.writerow([strain_acc_list[i], category_list[i], mle_selfing_generation, mle_selfing_generation_prob])
			selfing_generation_ls.append(mle_selfing_generation)
		del writer
		sys.stderr.write("Done.\n")
		return selfing_generation_ls
	
	def run(self):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
		data_matrix = Numeric.array(data_matrix)
		locus_allele_prob_vector = self.cal_locus_allele_prob_vector(data_matrix)
		locus_heterozygous_prob_vector = self.cal_locus_heterozygous_prob_vector(locus_allele_prob_vector)
		locus_heterozygous_prob_matrix = self.cal_locus_heterozygous_prob_matrix(locus_heterozygous_prob_vector, self.max_selfing_generation)
		selfing_generation_ls = self.cal_selfing_generation_prob(data_matrix, locus_heterozygous_prob_vector, strain_acc_list, category_list, locus_heterozygous_prob_matrix, self.output_fname)
		
		import pylab
		pylab.clf()
		pylab.hist(selfing_generation_ls, 20)
		pylab.title("hist of selfing generations")
		pylab.show()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:m:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_fname = None
	output_fname = None
	max_selfing_generation = 20
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
		elif opt in ("-m",):
			max_selfing_generation = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname:
		instance = EstimateSelfingGeneration(input_fname, output_fname, max_selfing_generation,\
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)