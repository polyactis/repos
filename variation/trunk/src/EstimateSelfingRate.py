#!/usr/bin/env python
"""
Usage: EstimateSelfingRate.py [OPTIONS] -i INPUT_FILE

Option:
	-i ...,	input file
	-o ...,	output file (place holder right now)
	-y ...,	which estimating method, 1(default)
	-a ...,	bits to toggle input/output StrainSNP matrix format, 0 is integer.
		1 is alphabet. 1st digit is for input. 2nd is for output. 00 (default)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	EstimateSelfingRate.py -i /tmp/pop_25.t1.14 -y2 -b

Description:
	Program to estimate selfing rate. estimating method includes
	Jarne2006, Robertson2004

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
import Numeric, pylab, rpy

class s_estimate_result:
	def __init__(self):
		self.FIS_vector = None
		self.selfing_rate_vector = None
		self.FIS_std = None
		self.s_std = None
		self.s_of_avg_FIS = None

class EstimateSelfingRate:
	def __init__(self, input_fname=None, output_fname=None, which_method=1,\
		nt_alphabet_bits='00', debug=0, report=0):
		"""
		2007-08-13
		"""
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.which_method = int(which_method)
		self.nt_alphabet_bits = nt_alphabet_bits
		self.debug = int(debug)
		self.report = int(report)
		self.estimate_method = {1: self.estimate_Jarne2006,
			2: self.estimate_Robertson2004}
	
	def cal_observed_heterozygosity_vector(self, data_matrix):
		sys.stderr.write("Calculating observed heterozygosity ...")
		no_of_strains, no_of_snps = data_matrix.shape
		observed_heterozygosity_vector = Numeric.zeros(no_of_snps, Numeric.Float)
		for i in range(no_of_snps):
			observed_heterozygosity_vector[i] = float(sum(data_matrix[:,i]>4))/sum(data_matrix[:,i]>0)
		sys.stderr.write("Done.\n")
		return observed_heterozygosity_vector
	
	def cal_FIS_vector(self, locus_allele_prob_vector, observed_heterozygosity_vector):
		sys.stderr.write("Calculating FIS vector ...")
		FIS_vector = []
		for i in range(len(locus_allele_prob_vector)):
			heterozygosity_e = 2*locus_allele_prob_vector[i,0]*(1-locus_allele_prob_vector[i,0])
			if heterozygosity_e >0:
				FIS_vector.append(1-observed_heterozygosity_vector[i]/heterozygosity_e)
		sys.stderr.write("Done.\n")
		return FIS_vector
	
	def cal_selfing_rate_vector(self, FIS_vector):
		sys.stderr.write("Calculating selfing rate vector ...")
		self_func = lambda x: 2*x/(1+x)
		selfing_rate_vector = map(self_func, FIS_vector)
		sys.stderr.write("Done.\n")
		return selfing_rate_vector
	
	def estimate_Jarne2006(self, data_matrix):
		"""
		2007-08-14
			method based on Jarne2006
		"""
		sys.stderr.write("Jarne2006 metod ...\n")
		s_estimate_result_instance = s_estimate_result()
		
		from EstimateSelfingGeneration import EstimateSelfingGeneration
		EstimateSelfingGeneration_instance = EstimateSelfingGeneration()
		locus_allele_prob_vector = EstimateSelfingGeneration_instance.cal_locus_allele_prob_vector(data_matrix)
		observed_heterozygosity_vector = self.cal_observed_heterozygosity_vector(data_matrix)
		FIS_vector = self.cal_FIS_vector(locus_allele_prob_vector, observed_heterozygosity_vector)
		selfing_rate_vector = self.cal_selfing_rate_vector(FIS_vector)
		s_estimate_result_instance.FIS_vector = FIS_vector
		s_estimate_result_instance.selfing_rate_vector = selfing_rate_vector
		return s_estimate_result_instance
	
	def cal_locus_allele_count_vector(self, data_matrix):
		sys.stderr.write("Calculating locus allele count ...")
		no_of_strains, no_of_snps = data_matrix.shape
		locus_allele_count_vector = Numeric.zeros([no_of_snps, 3], Numeric.Float)
		for i in range(no_of_snps):
			nt_number2index = {}
			for j in range(no_of_strains):
				nt_number = data_matrix[j,i]
				if nt_number!=0:	#not NA
					if nt_number>4:	#heterozygous, N12
						locus_allele_count_vector[i,2] += 1
					else:	#homozygous, N11 or N22
						if nt_number not in nt_number2index:
							nt_number2index[nt_number] = len(nt_number2index)
						index = nt_number2index[nt_number]
						locus_allele_count_vector[i,index] += 1
		sys.stderr.write("Done.\n")
		return locus_allele_count_vector
	
	def cal_Haldane_FIS_estimator(self, locus_allele_count_vector):
		"""
		2007-08-14
			simplified formula by maxima
		
		                                           2
							n12
					(4 (n11 n22 - ----) + n12) N
							4
		                   f = ----------------------------
					(N - 1) (2 N - N1) N1

		"""
		sys.stderr.write("Calculating Haldane FIS estimator...")
		FIS_vector = []
		selfing_rate_vector = []
		for i in range(len(locus_allele_count_vector)):
			N11, N22, N12 = locus_allele_count_vector[i]
			N = N11+N22+N12
			if N==N11 or N==N22:
				continue
			else:
				N1 = 2*N11+N12
				int_tmp = N11*N22-N12*N12/4
				F = (4*int_tmp+N12)*N/((N-1)*(2*N-N1)*N1)
				FIS_vector.append(F)
				s = 2*F/(1+F)
				selfing_rate_vector.append(s)
		sys.stderr.write("Done.\n")
		return FIS_vector, selfing_rate_vector
	
	def estimate_Robertson2004(self, data_matrix):
		"""
		2007-08-14
			method based on Robertson2004
		"""
		sys.stderr.write("Robertson2004 metod ...\n")
		s_estimate_result_instance = s_estimate_result()
		locus_allele_count_vector = self.cal_locus_allele_count_vector(data_matrix)
		FIS_vector, selfing_rate_vector = self.cal_Haldane_FIS_estimator(locus_allele_count_vector)
		s_estimate_result_instance.FIS_vector = FIS_vector
		s_estimate_result_instance.selfing_rate_vector = selfing_rate_vector
		return s_estimate_result_instance
	
	def run(self):
		"""
		2007-08-13
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname, int(self.nt_alphabet_bits[0]))
		data_matrix = Numeric.array(data_matrix)
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		s_estimate_result_instance = self.estimate_method[self.which_method](data_matrix)

		avg_s = sum(s_estimate_result_instance.selfing_rate_vector)/len( s_estimate_result_instance.selfing_rate_vector)
		s_std = rpy.r.sd(s_estimate_result_instance.selfing_rate_vector)
		avg_FIS = sum(s_estimate_result_instance.FIS_vector)/len(s_estimate_result_instance.FIS_vector)
		s_of_avg_FIS = 2*avg_FIS/(1+avg_FIS)
		print 'selfing_rate_vector:', s_estimate_result_instance.selfing_rate_vector
		print 'avg_s:', avg_s
		print 'std of s:', s_std
		print 'avg_FIS:', avg_FIS
		print 's_of_avg_FIS:', s_of_avg_FIS
		pylab.hist(s_estimate_result_instance.selfing_rate_vector, 20)
		pylab.show()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:y:a:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_fname = None
	output_fname = None
	which_method = 1
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
		elif opt in ("-y",):
			which_method = int(arg)
		elif opt in ("-a",):
			nt_alphabet_bits = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname:
		instance = EstimateSelfingRate(input_fname, output_fname, which_method,\
			nt_alphabet_bits, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)