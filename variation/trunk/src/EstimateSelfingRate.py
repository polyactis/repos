#!/usr/bin/env python
"""
Usage: EstimateSelfingRate.py [OPTIONS] -i INPUT_FILE -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input file
	-o ...,	output file
	-y ...,	which estimating method, 1(default)
	-a ...,	bits to toggle input/output StrainSNP matrix format, 0 is integer.
		1 is alphabet. 1st digit is for input. 2nd is for output. 00 (default)
	-s ...,	selfing rate table (need to specify if commit is toggled. to store selfing rate in db)
	-c,	commit
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	EstimateSelfingRate.py -i /tmp/pop_25.t1.14 -o /tmp/pop_25.t1.14.y2 -y2 -b

Description:
	Program to estimate selfing rate. estimating method includes
	Jarne2006(method 1), Robertson1984(2), Weir1984(3), Nordborg1997(4),
	David2007(by g2, 5)
	
	NA-filtering has to be carried out on the data beforehand.

"""
from __init__ import *
import Numeric, pylab, rpy

class s_estimate_result:
	def __init__(self):
		self.FIS_vector = []
		self.selfing_rate_vector = []
		self.FIS_std = None
		self.avg_s = None
		self.std_s = None
		self.s_of_avg_FIS = None
		self.weir1984_multi_loci_s = None	#for Weir1984
		self.s_M_ls = []	#for Nordborg1997
		self.theta_ls = []
		self.theta_M_ls = []
		self.g2_David2007 = None	#for David2007
		self.s_g2_David2007 = None
		self.genotyping_error_rate_vector = None

class EstimateSelfingRate:
	__doc__ = __doc__
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, 'database name', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname', 1, ): ['', 'i', 1, 'SNP data formatted in populations, outputted by OutputPopulation.py'],\
							('output_fname', 1, ): [None, 'o', 1, ''],\
							('which_method', 1, int): [1, 'y', 1, 'which estimating method'],\
							('nt_alphabet_bits', 1, ): ['00', 'a', 1, 'bits to toggle input/output StrainSNP matrix format, 0 is integer. 1 is alphabet. 1st digit is for input. 2nd is for output'],\
							('selfing_rate_table', 1, ): [None, 's', 1, 'table to store estimates into db'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-06-02
			use ProcessOptions
		2007-08-13
		2007-08-30 add database options
		"""
		
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		self.estimate_method = {1: self.estimate_Jarne2006,
			2: self.estimate_Robertson2004,
			3: self.estimate_Weir1984,
			4: self.estimate_Nordborg1997,
			5: self.estimate_David2007_g2}
		self.method2table_entry = {1: ['avg_s_Jarne2006', 'std_s_Jarne2006'],
			2: ['avg_s_Robertson1984', 'std_s_Robertson1984'],
			3: ['avg_s_Weir1984', 'std_s_Weir1984'],
			4: ['avg_s_Nordborg1997', 'std_s_Nordborg1997']}
	
	def cal_observed_heterozygosity_vector(self, data_matrix):
		sys.stderr.write("Calculating observed heterozygosity ...")
		no_of_strains, no_of_snps = data_matrix.shape
		observed_heterozygosity_vector = Numeric.zeros(no_of_snps, Numeric.Float)
		for i in range(no_of_snps):
			observed_heterozygosity_vector[i] = float(sum(data_matrix[:,i]>4))/sum(data_matrix[:,i]>0) #2007-08-14 data filtered, no n=0
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
		sys.stderr.write("Jarne2006 method ...\n")
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
			method based on Robertson1984
		"""
		sys.stderr.write("Robertson1984 method ...\n")
		s_estimate_result_instance = s_estimate_result()
		locus_allele_count_vector = self.cal_locus_allele_count_vector(data_matrix)
		FIS_vector, selfing_rate_vector = self.cal_Haldane_FIS_estimator(locus_allele_count_vector)
		s_estimate_result_instance.FIS_vector = FIS_vector
		s_estimate_result_instance.selfing_rate_vector = selfing_rate_vector
		return s_estimate_result_instance
	
	def cal_Weir1984_b_c_vector(self, data_matrix):
		"""
		2007-08-14
			regard r=1 in formula 3 and 4
			
		"""
		sys.stderr.write("Calculating Weir1984 b-c vector ...")
		no_of_strains, no_of_snps = data_matrix.shape
		b_c_vector = []
		FIS_vector = []
		for i in range(no_of_snps):
			nt2counter = {}
			no_of_valid_calls = 0.0
			no_of_heterozygous_calls = 0.0
			for j in range(no_of_strains):
				nt_number = data_matrix[j,i]
				if nt_number!=0:	#not NA
					no_of_valid_calls += 2
					nt_string = number2nt[data_matrix[j,i]]
					if len(nt_string)==1:	#double the nt if it's homozygous
						nt_string += nt_string
					else:
						no_of_heterozygous_calls += 1
					for k in range(len(nt_string)):
						nt = nt_string[k]
						if nt not in nt2counter:
							nt2counter[nt] = 0
						nt2counter[nt] += 1
			nt_key_1 = nt2counter.keys()[0]
			p = nt2counter[nt_key_1]/no_of_valid_calls
			n = no_of_valid_calls/2
			h = no_of_heterozygous_calls/n
			b = n/(n-1)*(p*(1-p)-(2*n-1)/(4*n)*h)	#2007-08-14 data filtered, no n=0
			c = h/2
			b_c_vector.append((b,c))
			if (b+c)!=0:
				FIS_vector.append(b/(b+c))
		sys.stderr.write("Done.\n")
		return b_c_vector, FIS_vector
	
	def cal_weir1984_multi_loci_s(self, b_c_vector):
		"""
		2007-08-14
		"""
		sys.stderr.write("Calculating weir1984_multi_loci_s ...")
		numerator = 0.0
		denominator = 0.0
		for b,c in b_c_vector:
			numerator += b
			denominator += (b+c)
		weir1984_multi_loci_s = numerator/denominator
		sys.stderr.write("Done.\n")
		return weir1984_multi_loci_s
	
	def estimate_Weir1984(self, data_matrix):
		"""
		2007-08-14
		"""
		sys.stderr.write("Weir1984 method ...\n")
		s_estimate_result_instance = s_estimate_result()
		b_c_vector, FIS_vector = self.cal_Weir1984_b_c_vector(data_matrix)
		s_estimate_result_instance.FIS_vector = FIS_vector
		self_func = lambda x: 2*x/(1+x)
		s_estimate_result_instance.selfing_rate_vector = map(self_func, FIS_vector)
		s_estimate_result_instance.weir1984_multi_loci_s = self.cal_weir1984_multi_loci_s(b_c_vector)
		return s_estimate_result_instance
	
	def estimate_Nordborg1997(self, data_matrix):
		"""
		2007-08-14
		"""
		sys.stderr.write("Nordborg1997 method ...\n")
		no_of_strains, no_of_snps = data_matrix.shape
		s_estimate_result_instance = s_estimate_result()
		from Nordborg1997.Simulate import Estimate
		Estimate_instance = Estimate()
		for i in range(no_of_snps):
			individual_ls = []
			for j in range(no_of_strains):
				nt_number = data_matrix[j,i]
				if nt_number!=0:	#not NA
					nt_string = number2nt[nt_number]
					if len(nt_string)==1:	#double the nt if it's homozygous
						nt_number1 = nt_number2 = nt_number
					else:
						nt_number1 = nt2number[nt_string[0]]
						nt_number2 = nt2number[nt_string[1]]
					individual_ls.append((nt_number1, nt_number2))
			if individual_ls:
				Hw = Estimate_instance.estimate_Hw(individual_ls)
				Hb = Estimate_instance.estimate_Hb(individual_ls)
				if (2*Hb-Hw)!=1.0:
					s = Estimate_instance.estimate_s(Hw, Hb)
					theta = Estimate_instance.estimate_theta(Hw, Hb)
					s_M = Estimate_instance.estimate_s_Milligan1996(Hw, Hb)
					theta_M = Estimate_instance.estimate_theta_Milligan1996(Hw, Hb)
					s_estimate_result_instance.selfing_rate_vector.append(s)
					s_estimate_result_instance.theta_ls.append(theta)
					s_estimate_result_instance.s_M_ls.append(s_M)
					s_estimate_result_instance.theta_M_ls.append(theta_M)
		sys.stderr.write("Done.\n")
		return s_estimate_result_instance
	
	def estimate_David2007_g2(self, data_matrix):
		"""
		2007-08-15
			David2007
		"""
		sys.stderr.write("David2007 g2 method ...\n")
		no_of_strains, no_of_snps = data_matrix.shape
		s_estimate_result_instance = s_estimate_result()
		Hijk = Numeric.zeros([no_of_snps, no_of_snps], Numeric.Float)	#each entry is \sum_k Hik*Hjk
		Hijkl = Numeric.zeros([no_of_snps, no_of_snps], Numeric.Float)	#each entry is \sum_k \sum_{l!=k} Hik*Hjl
		Mij = Numeric.zeros([no_of_snps, no_of_snps], Numeric.Float)	#number of individuals with NA at locus both i an j
		Mi = Numeric.zeros(no_of_snps, Numeric.Float)	#number of individuals with NA at locus i
		for i in range(no_of_snps):
			for j in range(no_of_snps):
				if j==i:
					continue
				for k in range(no_of_strains):
					if data_matrix[k,i]==0:
						Mi[i] += 1	#this one is inflated no_of_snps-1 times due to the inner j loop
						if data_matrix[k,j]==0:
							Mij[i,j] += 1
					Hik = int(data_matrix[k,i]>4)	#NA=0 is treated as homozygous
					Hjk = int(data_matrix[k,j]>4)
					Hijk[i,j] += (Hik*Hjk)
					for l in range(no_of_strains):
						if l==k:
							continue
						Hik = int(data_matrix[k,i]>4)	#NA=0 is treated as homozygous
						Hjl = int(data_matrix[l,j]>4)
						Hijkl[i,j] += (Hik*Hjl)
		Mi = Mi/(no_of_snps-1)
		
		numerator = 0.0
		denominator = 0.0
		for i in range(no_of_snps):
			for j in range(no_of_snps):
				if j==i:
					continue
				numerator += 1.0/(no_of_strains-Mij[i,j])*Hijk[i,j]
				denominator += 1.0/(no_of_strains*(no_of_strains-1)-Mi[i]*Mi[j] + Mij[i,j])*Hijkl[i,j]
		if denominator>0.0:
			g2 = numerator/denominator
			s_estimate_result_instance.g2_David2007 = g2
		else:
			g2=None
		if g2:	#g2=0 causes ZeroDivisionError: float division in calculating s
			s = (1+5*g2-math.sqrt(1+10*g2+9*g2*g2))/(2*g2)
			FIS = s/(2-s)
			
			from EstimateSelfingGeneration import EstimateSelfingGeneration
			EstimateSelfingGeneration_instance = EstimateSelfingGeneration()
			locus_allele_prob_vector = EstimateSelfingGeneration_instance.cal_locus_allele_prob_vector(data_matrix)
			observed_heterozygosity_vector = self.cal_observed_heterozygosity_vector(data_matrix)
			FIS_vector = self.cal_FIS_vector(locus_allele_prob_vector, observed_heterozygosity_vector)
			
			genotyping_error_func = lambda x: (x-FIS)/(1-FIS)
			genotyping_error_rate_vector = map(genotyping_error_func, FIS_vector)
			
			s_estimate_result_instance.s_g2_David2007 = s
			s_estimate_result_instance.FIS_vector = FIS_vector
			s_estimate_result_instance.genotyping_error_rate_vector = genotyping_error_rate_vector
		
		return s_estimate_result_instance
	
	def spruce_a_list(self, writer, ls, variable_name):
		"""
		2007-08-30
			return avg, std
		"""
		avg = sum(ls)/float(len(ls))
		std = rpy.r.sd(ls)
		writer.writerow(['list of %s:'%variable_name, ls])
		writer.writerow(['avg_%s:'%variable_name, avg])
		writer.writerow(['std_%s:'%variable_name, std])
		return avg, std
	
	def write_result(self, output_fname, s_estimate_result_instance):
		"""
		2007-08-13
		2007-08-30
			assign value to s_estimate_result_instance.avg_s and std_s
		"""
		sys.stderr.write("Writing result to output file ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		if s_estimate_result_instance.selfing_rate_vector:
			avg_s, std_s = self.spruce_a_list(writer, s_estimate_result_instance.selfing_rate_vector, 's')
			s_estimate_result_instance.avg_s = avg_s
			s_estimate_result_instance.std_s = std_s
		if s_estimate_result_instance.FIS_vector:
			self.spruce_a_list(writer, s_estimate_result_instance.FIS_vector, 'FIS')
			avg_FIS = sum(s_estimate_result_instance.FIS_vector)/len(s_estimate_result_instance.FIS_vector)
			s_of_avg_FIS = 2*avg_FIS/(1+avg_FIS)
			writer.writerow(['s_of_avg_FIS:', s_of_avg_FIS])
		if s_estimate_result_instance.weir1984_multi_loci_s:
			writer.writerow(['weir1984_multi_loci_s:', s_estimate_result_instance.weir1984_multi_loci_s])
		if s_estimate_result_instance.theta_ls:
			self.spruce_a_list(writer, s_estimate_result_instance.theta_ls, 'theta')
		if s_estimate_result_instance.theta_M_ls:
			self.spruce_a_list(writer, s_estimate_result_instance.theta_M_ls, 'theta_M')
		if s_estimate_result_instance.g2_David2007:
			writer.writerow(['g2:', s_estimate_result_instance.g2_David2007])
		if s_estimate_result_instance.s_g2_David2007:
			writer.writerow(['s_g2:', s_estimate_result_instance.s_g2_David2007])
		if s_estimate_result_instance.genotyping_error_rate_vector:
			self.spruce_a_list(writer, s_estimate_result_instance.genotyping_error_rate_vector, 'genotyping error rate vector')
		del writer
		sys.stderr.write("Done.\n")
	
	def create_selfing_rate_table(self, curs, selfing_rate_table):
		"""
		2007-08-30
		"""
		sys.stderr.write("Creating selfing rate table ...")
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			popid	integer not null,\
			avg_s_Jarne2006	float,\
			std_s_Jarne2006	float,\
			avg_s_Robertson1984	float,\
			std_s_Robertson1984	float,\
			avg_s_Weir1984	float,\
			std_s_Weir1984	float,\
			weir1984_multi_loci_s	float,\
			avg_s_Nordborg1997	float,\
			std_s_Nordborg1997	float,\
			s_g2_David2007	float)engine=INNODB;"%selfing_rate_table)
		sys.stderr.write("Done.\n")
	
	def check_table_existent(self, curs, selfing_rate_table):
		"""
		2007-08-30
			check if table exists or not. mysql-specific command is used.
		"""
		curs.execute("show tables")
		rows = curs.fetchall()
		for row in rows:
			if row[0]==selfing_rate_table:
				return 1
		return 0
	
	def submit_to_table(self, curs, selfing_rate_table, s_estimate_result_instance, popid, method2table_entry, which_method):
		"""
		2007-08-30
		"""
		sys.stderr.write("Submitting selfing rate ...")
		if not self.check_table_existent(curs, selfing_rate_table):
			self.create_selfing_rate_table(curs, selfing_rate_table)
		#check whether popid has data recorded in the table already
		curs.execute("select * from %s where popid=%s"%(selfing_rate_table, popid))
		rows = curs.fetchall()
		if not rows:
			curs.execute("insert %s(popid) values(%s)"%(selfing_rate_table, popid))
		if s_estimate_result_instance.avg_s is not None and s_estimate_result_instance.std_s is not None:
			curs.execute("update %s set %s=%s, %s=%s where popid=%s"%(selfing_rate_table,\
				method2table_entry[which_method][0], s_estimate_result_instance.avg_s,\
				method2table_entry[which_method][1], s_estimate_result_instance.std_s,\
				popid))
		if s_estimate_result_instance.weir1984_multi_loci_s is not None:
			curs.execute("update %s set weir1984_multi_loci_s=%s where popid=%s"%(selfing_rate_table,s_estimate_result_instance.weir1984_multi_loci_s, popid))
		if s_estimate_result_instance.s_g2_David2007 is not None:
			curs.execute("update %s set s_g2_David2007=%s where popid=%s"%(selfing_rate_table,s_estimate_result_instance.s_g2_David2007, popid))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-08-13
		"""
		header, strain_acc_list, category_list, data_matrix = read_data(self.input_fname, int(self.nt_alphabet_bits[0]))
		data_matrix = Numeric.array(data_matrix)
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		s_estimate_result_instance = self.estimate_method[self.which_method](data_matrix)
		self.write_result(self.output_fname, s_estimate_result_instance)
		
		import re
		pop_number_pattern = re.compile('\d+$')	#the trailing number in the input_fname is population number
		if pop_number_pattern.search(self.input_fname):
			pop_number = pop_number_pattern.search(self.input_fname).group()
		else:
			pop_number = '00'
		if s_estimate_result_instance.selfing_rate_vector:
			pylab.title("histogram of selfing rate. pop %s"%pop_number)
			pylab.hist(s_estimate_result_instance.selfing_rate_vector, 20)
			pylab.savefig('%s.png'%self.output_fname)
		if self.commit:
			if pop_number=='00':
				sys.stderr.write("Can't infer pop_number from input_fname. Exit!\n")
				sys.exit(1)
			if not self.selfing_rate_table:
				sys.stderr.write("Need to specify selfing_rate_table. Exit!\n")
				sys.exit(1)
			import MySQLdb
			conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user = self.user, passwd = self.passwd)
			curs = conn.cursor()
			self.submit_to_table(curs, self.selfing_rate_table, s_estimate_result_instance, pop_number, self.method2table_entry, self.which_method)
			conn.commit()
			
			
if __name__ == '__main__':
	main_class = EstimateSelfingRate
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:y:a:s:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	input_fname = None
	output_fname = None
	which_method = 1
	nt_alphabet_bits = '00'
	selfing_rate_table = None
	commit = 0
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
		elif opt in ("-y",):
			which_method = int(arg)
		elif opt in ("-a",):
			nt_alphabet_bits = arg
		elif opt in ("-s",):
			selfing_rate_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_fname and output_fname:
		instance = EstimateSelfingRate(hostname, dbname, schema, input_fname, output_fname,\
			which_method, nt_alphabet_bits, selfing_rate_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""