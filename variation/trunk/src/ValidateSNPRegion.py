#!/usr/bin/env python
"""  
Examples:
	ValidateSNPRegion.py -i Desktop/ftRegions.csv -o /Network/Data/250k/tmp-yh/ftRegions_sig.tsv -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf -D /Network/Data/250k/tmp-yh/call_method_17_LD_pickle_dummy_n0.1_m20000 -u yh -p secret

Description:
	2008-10-02
	program to validate snp regions (from magnus) under some criteria. several options are useless, inherited from DrawSNPRegion.py.
	no time to cleanup.
	
	Input Filename contains with 4 bare columns: chr1, pos1, chr2, pos2.
"""



import sys, os, math, csv
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment

from DrawSNPRegion import SNPPassingData
from pymodule import PassingData, figureOutDelimiter
from sets import Set

class ValidateSNPRegion(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('LD_fname', 0, ): [None, 'L', 1, 'the file containing LD info, output of MpiLD.py', ],\
							("min_distance", 1, int): [5000, 'm', 1, 'minimum distance for a given input region. expand the input region if it falls short.'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
							("list_type_id", 0, int): [0, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("input_fname", 1, ): [None, 'i', 1, 'Filename which contains with 4 bare columns: chr1, pos1, chr2, pos2.'],\
							('output_fname', 1, ): [None, 'o', 1, 'output filename.'],\
							('call_method_id', 0, int):[17, '', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							("which_LD_statistic", 1, int): [2, 'w', 1, 'which LD_statistic to plot, 1=r2, 2=|D_prime|, 3=|D|'],\
							("LD_info_picklef", 0, ): ['/Network/Data/250k/tmp-yh/call_method_17_LD_pickle_dummy_n0.1_m20000', 'D', 1, 'given the option, If the file does not exist yet, store a pickled LD_info into it (min_MAF and min_distance will be attached to the filename). If the file exists, load LD_info out of it.'],\
							("gene_annotation_picklef", 0, ): ['/Network/Data/250k/tmp-yh/at_gene_model_pickelf', 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
							('min_margarita_rank',0, int): [50, '', 1, 'Minimum rank of a margarita score to be considered significant'],\
							('min_rf_rank',0, int): [25, '', 1, 'Minimum rank of a RandomForest score to be considered significant'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 1, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-2
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def findSNPsInRegion(self, snp_info, chromosome, start, stop, center_snp_position=None):
		"""
		2008-10-1
			called by plotSNPRegion()
			find SNPs in this region, if center_snp_position is not given, find one.
			similar to getSNPsAroundThisSNP()
		"""
		if self.report:
			sys.stderr.write("Get SNPs in this region ...")
		from DrawSNPRegion import SNPPassingData
		chr_pos_ls = []
		chr_pos2adjacent_window = {}
		j = 0
		midpoint = (start+stop)/2.
		if center_snp_position is None:
			_center_snp_position = start
		else:
			_center_snp_position = center_snp_position
		center_snp = SNPPassingData(chromosome=chromosome, position=_center_snp_position, snps_id=None)
		for i in range(start-1, stop+2):
			new_pos = i
			new_chr_pos = (chromosome, new_pos)
			if new_chr_pos in snp_info.chr_pos2index:
				if center_snp_position is None and abs(new_pos-midpoint)<abs(center_snp.position-midpoint):	#this SNP is closer to the center
					center_snp.position = new_pos
				chr_pos_ls.append(new_chr_pos)
				if j!=0:
					#add_mid_point(chr_pos_ls, chr_pos2adjacent_window)
					pass
				j += 1
		center_snp.snps_id = '%s_%s'%(center_snp.chromosome, center_snp.position)
		snp_region = PassingData(chr_pos_ls=chr_pos_ls, chr_pos2adjacent_window=chr_pos2adjacent_window, center_snp=center_snp)
		if self.report:
			sys.stderr.write("Done.\n")
		return snp_region
	
	def get_phenotype_id_ls_given_biology_category_id(self, biology_category_id):
		import Stock_250kDB
		rows = Stock_250kDB.PhenotypeMethod.query.filter_by(biology_category_id=biology_category_id)
		phenotype_id_ls = []
		for row in rows:
			phenotype_id_ls.append(row.id)
		return phenotype_id_ls
	
	
	def checkIfRegionSignificant(self, DrawSNPRegion_instance, snp_region, snp_info, analysis_method_id2gwr, value_criteria):
		is_significant = 0
		#analysis_method_id_set = Set(value_criteria.keys())
		for chr_pos in snp_region.chr_pos_ls:
			for analysis_method_id, criteria_type in value_criteria:
				criteria_value = value_criteria[(analysis_method_id, criteria_type)]
				gwr = analysis_method_id2gwr[analysis_method_id]
				data_obj = gwr.get_data_obj_by_chr_pos(chr_pos[0], chr_pos[1])
				if not data_obj:
					continue
				if criteria_type=='rank':
					data_obj_at_this_rank = gwr.get_data_obj_at_given_rank(criteria_value)
					if data_obj>=data_obj_at_this_rank:
						is_significant = 1
						break
				elif data_obj.value>=criteria_value:
					is_significant = 1
					break
			if is_significant==1:	#break if already found significant
				break
		del analysis_method_id2gwr
		return is_significant
	
	def get_snp_region_ls(self, ft_region_fname, snp_info, min_distance=5000):
		sys.stderr.write("Get all snp regions ...")
		delimiter = figureOutDelimiter(ft_region_fname)
		ft_region_reader = csv.reader(open(ft_region_fname, 'r'), delimiter=delimiter)		
		snp_region_ls = []
		for row in ft_region_reader:
			row = map(int, row)
			chr1, pos1, chr2, pos2 = row
			if pos2<pos1:
				pos1, pos2 = pos2, pos1
			span = abs(pos2-pos1)
			if span < min_distance*2:
				extra_span = (min_distance*2-span)/2
				pos1 = max(pos1 - extra_span, 1)
				pos2 = pos2+extra_span
			
			snp_region = self.findSNPsInRegion(snp_info, chr1, pos1, pos2, center_snp_position=None)
			snp_region_ls.append(snp_region)
		del ft_region_reader
		sys.stderr.write("Done.\n")
		return snp_region_ls
	
	def checkRegions(self, DrawSNPRegion_instance, snp_info, snp_region_ls, output_fname, value_criteria, biology_category_id=1):
		sys.stderr.write("Checking regions ...")
		outf = open(output_fname, 'w')
		writer = csv.writer(outf, delimiter='\t')
		header = ['phenotype_id', 'chromosome', 'start', 'stop']
		writer.writerow(header)
		
		ft_phenotype_id_ls = self.get_phenotype_id_ls_given_biology_category_id(biology_category_id)
		ft_phenotype_id_ls.sort()
		analysis_method_id_ls = [row[0] for row in value_criteria.keys()]
		for phenotype_id in ft_phenotype_id_ls:
			analysis_method_id2gwr = DrawSNPRegion_instance.getSimilarGWResultsGivenResultsByGene(phenotype_id, call_method_id=self.call_method_id, \
																							results_directory=self.results_directory, \
																							analysis_method_id_ls=analysis_method_id_ls)
			if not analysis_method_id2gwr:	#no gwr results, skip
				continue
			for snp_region in snp_region_ls:
				if self.checkIfRegionSignificant(DrawSNPRegion_instance, snp_region, snp_info, analysis_method_id2gwr, value_criteria):
					data_row = [phenotype_id, snp_region.center_snp.chromosome, snp_region.chr_pos_ls[0][1]-1, snp_region.chr_pos_ls[-1][1]+1]
					writer.writerow(data_row)
			outf.flush()
		del writer
		sys.stderr.write("Done.\n")

	def run(self):
		if self.debug==1:
			import pdb
			pdb.set_trace()
		
		from DrawSNPRegion import DrawSNPRegion
		DrawSNPRegion_instance = DrawSNPRegion(db_user=self.db_user, db_passwd=self.db_passwd, hostname=self.hostname, \
											database=self.dbname, input_fname='/tmp/dumb', output_dir='/tmp', debug=0)
		
		grand_dataStructure = DrawSNPRegion_instance.loadDataStructure(self.gene_annotation_picklef, self.LD_info_picklef, self.LD_fname, min_MAF=self.min_MAF, min_distance=20000, list_type_id=None)		
		snp_region_ls = self.get_snp_region_ls(self.input_fname, grand_dataStructure.snp_info, self.min_distance)
		value_criteria={(1, 'value'):8., (7, 'value'):6., (5, 'rank'):self.min_margarita_rank,(6, 'rank'):self.min_rf_rank}	#minimum threshold for different analysis methods
		self.checkRegions(DrawSNPRegion_instance, grand_dataStructure.snp_info, snp_region_ls, self.output_fname, value_criteria)
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ValidateSNPRegion
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()