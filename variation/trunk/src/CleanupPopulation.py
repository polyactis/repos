#!/usr/bin/env python
"""

Examples:
	CleanupPopulation.py -l popid2ecotypeid_50 -o popid2snpid_50  -c
	
Description:
	2007-07-16
	
	Clean up population by removing sequentially
	1. rows(strains) with too many NAs
	2. cols(snps) with too many NAs
	3. bad snps (too many bogus heterozygous calls)
	
"""
from __init__ import *

class CleanupPopulation:
	__doc__ = __doc__
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, 'database name', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_table', 1, ): ['calls', 'i', 1, 'SNP data table'],\
							('output_table', 1, ): [None, 'o', 1, 'Table to store population id versus snpid'],\
							('strain_info_table', 1, ): ['ecotype', 's', 1, 'ecotype table'],\
							('snp_locus_table', 1, ): ['snps', 'n', 1, 'SNP locus table'],\
							('population_table', 1, ): [None, 'l', 1, 'Table storing population id versus ecotypeid'],\
							('min_no_of_strains_per_pop', 1, int, ): [5, 'm', 1, ''],\
							('max_row_NA_max_col_NA_min_log_prob', 1, ): ['1,1,-0.5', 't', 1, 'min row NA ratio to toss out, min column NA ratio, min log probability'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-06-02
			use ProcessOptions
		2007-07-16
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		cutoff_ls = self.max_row_NA_max_col_NA_min_log_prob.split(',')
		self.row_cutoff, self.col_cutoff, self.min_log_prob = map(float, cutoff_ls)
	
	def create_sub_data_matrix(self, popid, data_matrix, ecotypeid_ls, strain_id2index):
		sys.stderr.write("\tCreating sub data_matrix for population %s for %s ecotypes ..."%(popid, len(ecotypeid_ls)))
		no_of_strains, no_of_snps = data_matrix.shape
		sub_data_matrix = []
		ecotypeid_ls.sort()
		new_ecotypeid_ls = []
		for i in range(len(ecotypeid_ls)):
			#import pdb
			#pdb.set_trace()
			if ecotypeid_ls[i] in strain_id2index:
				sub_data_matrix.append(data_matrix[strain_id2index[ecotypeid_ls[i]],])
				new_ecotypeid_ls.append(ecotypeid_ls[i])
		sub_data_matrix = num.array(sub_data_matrix)
		sys.stderr.write("Done.\n")
		return sub_data_matrix, new_ecotypeid_ls
	
	def cleanup_one_population(self, FilterStrainSNPMatrix_instance, RemoveBadSNPs_instance, data_matrix, ecotypeid_ls, snp_id_list, min_no_of_strains_per_pop, row_cutoff, col_cutoff, min_log_prob):
		sys.stderr.write("\tCleaning up ... \n")
		
		if row_cutoff<=1 and row_cutoff>=0:
			remove_rows_data = FilterStrainSNPMatrix_instance.remove_rows_with_too_many_NAs(data_matrix, row_cutoff)
			rows_with_too_many_NAs_set = remove_rows_data.rows_with_too_many_NAs_set
			strain_index2no_of_NAs = remove_rows_data.row_index2no_of_NAs
		else:
			rows_with_too_many_NAs_set = Set()
		
		if col_cutoff<=1 and col_cutoff>=0:
			remove_cols_data = FilterStrainSNPMatrix_instance.remove_cols_with_too_many_NAs(data_matrix, col_cutoff, rows_with_too_many_NAs_set)
			cols_with_too_many_NAs_set = remove_cols_data.cols_with_too_many_NAs_set
		else:
			cols_with_too_many_NAs_set = Set()
		
		strain_id_selected = []
		strain_index_selected = []
		for i in range(len(ecotypeid_ls)):
			if i not in rows_with_too_many_NAs_set:
				strain_index_selected.append(i)
				strain_id_selected.append(ecotypeid_ls[i])
		snp_index_selected = []
		for i in range(len(snp_id_list)):
			if i not in cols_with_too_many_NAs_set:
				snp_index_selected.append(i)
		
		if len(strain_index_selected)<min_no_of_strains_per_pop:
			#too few strains
			return [], []
		
		data_matrix = num.take(data_matrix, strain_index_selected)
		data_matrix = num.take(data_matrix, snp_index_selected, 1)
		
		strain_homo_perc_vector = RemoveBadSNPs_instance.cal_strain_homo_perc_vector(data_matrix)
		snp_locus_log_prob = RemoveBadSNPs_instance.cal_snp_locus_log_prob(data_matrix, strain_homo_perc_vector)
		#2nd selection for snp ids
		snp_id_selected = []
		no_of_bad_snps = 0
		for i in range(len(snp_locus_log_prob)):
			if snp_locus_log_prob[i]>min_log_prob:
				snp_id_selected.append(snp_id_list[snp_index_selected[i]])
			else:
				no_of_bad_snps += 1
		sys.stderr.write("\t%s bad snps.\n"%(no_of_bad_snps))
		sys.stderr.write("Done.\n")
		return strain_id_selected, snp_id_selected
	
	def create_popid2snpid_table(self, curs, output_table):
		"""
		2007-07-16
		"""
		sys.stderr.write("Creating table %s ..."%(output_table))
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			popid	integer not null,\
			snpid	integer)engine=INNODB"%output_table)
		sys.stderr.write("Done.\n")
	
	def mark_strain_id_selected(self, curs, popid2strain_id_snp_id_ls, population_table):
		"""
		2007-07-16
		2007-07-16
			add popid=%s one more condition to database update command
		"""
		sys.stderr.write("Marking strain_id selected ...")
		for popid, strain_id_snp_id_list in popid2strain_id_snp_id_ls.iteritems():
			strain_id_list = strain_id_snp_id_list[0]
			for ecotypeid in strain_id_list:
				curs.execute("update %s set selected=1 where ecotypeid=%s and popid=%s"%(population_table, ecotypeid, popid))
		sys.stderr.write("Done.\n")
	
	def submit_popid2snpid_list(self, curs, popid2strain_id_snp_id_ls, population_table, output_table):
		sys.stderr.write("Submitting popid2snpid_list ...")
		for popid, strain_id_snp_id_list in popid2strain_id_snp_id_ls.iteritems():
			strain_id_list, snp_id_list = strain_id_snp_id_list
			for snpid in snp_id_list:
				curs.execute("insert into %s(popid, snpid) values(%s, %s)"%\
				(output_table, popid, snpid))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		from dbSNP2data import dbSNP2data
		dbSNP2data_instance = dbSNP2data(user=self.user, passwd=self.passwd, output_fname='whatever')
		
		snp_id2index, snp_id_list, snp_id2info = dbSNP2data_instance.get_snp_id2index_m(curs, self.input_table, self.snp_locus_table)
		#strain_id2index, strain_id_list
		strain_id2index, strain_id_list, nativename2strain_id, strain_id2acc, strain_id2category  = dbSNP2data_instance.get_strain_id2index_m(curs, self.input_table, self.strain_info_table)
		#2008-06-02 stuff returned by get_strain_id2index_m is totally changed.
		ecotype_id2row_index = {}
		for strain_id, acc in strain_id2acc.iteritems():
			row_index = strain_id2index[strain_id]
			ecotype_id2row_index[acc] = row_index
		
		#strain_id2acc, strain_id2category = dbSNP2data_instance.get_strain_id_info_m(curs, strain_id_list, self.strain_info_table)
		snp_id2acc = dbSNP2data_instance.get_snp_id_info_m(curs, snp_id_list, self.snp_locus_table)
		data_matrix = dbSNP2data_instance.get_data_matrix_m(curs, strain_id2index, snp_id2index, nt2number, self.input_table, need_heterozygous_call=1)
		
		
		from OutputPopulation import OutputPopulation
		
		popid2ecotypeid_ls = OutputPopulation.get_popid2ecotypeid_ls(curs, self.population_table)
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		
		from RemoveBadSNPs import RemoveBadSNPs
		RemoveBadSNPs_instance = RemoveBadSNPs()
		popid2strain_id_snp_id_ls = {}
		for popid, ecotypeid_ls in popid2ecotypeid_ls.iteritems():
			if len(ecotypeid_ls)>=self.min_no_of_strains_per_pop:
				sys.stderr.write("Population %s\n"%popid)
				sub_data_matrix, new_ecotypeid_ls = self.create_sub_data_matrix(popid, data_matrix, ecotypeid_ls, ecotype_id2row_index)
				if len(new_ecotypeid_ls)>=self.min_no_of_strains_per_pop:
					sys.stderr.write("\tPopulation %s has %s strains\n"%(popid, len(new_ecotypeid_ls)))
					strain_id_selected, snp_id_selected = self.cleanup_one_population(FilterStrainSNPMatrix_instance, RemoveBadSNPs_instance, sub_data_matrix, new_ecotypeid_ls, snp_id_list, self.min_no_of_strains_per_pop, self.row_cutoff, self.col_cutoff, self.min_log_prob)
					if strain_id_selected and snp_id_selected:
						popid2strain_id_snp_id_ls[popid] = [strain_id_selected, snp_id_selected]
		
		if self.commit:
			self.create_popid2snpid_table(curs, self.output_table)
			self.mark_strain_id_selected(curs, popid2strain_id_snp_id_ls, self.population_table)
			self.submit_popid2snpid_list(curs, popid2strain_id_snp_id_ls, self.population_table, self.output_table)
			conn.commit()
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = CleanupPopulation
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:n:p:m:t:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	input_table = 'calls'
	output_table = ''
	strain_info_table = 'ecotype'
	snp_locus_table = 'snps'
	population_table = ''
	min_no_of_strains_per_pop = 5
	row_cutoff = 0.4
	col_cutoff = 0.4
	min_log_prob = -0.5
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
			input_table = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-p",):
			population_table = arg
		elif opt in ("-m",):
			min_no_of_strains_per_pop = int(arg)
		elif opt in ("-t",):
			cutoff_ls = arg.split(',')
			print cutoff_ls
			row_cutoff, col_cutoff, min_log_prob = map(float, cutoff_ls)
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1

	if input_table and output_table and hostname and dbname and schema and population_table:
		instance = CleanupPopulation(hostname, dbname, schema, input_table, output_table, \
			strain_info_table, snp_locus_table, population_table, \
			min_no_of_strains_per_pop, row_cutoff, col_cutoff, min_log_prob, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
	"""