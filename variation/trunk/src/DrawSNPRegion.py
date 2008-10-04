#!/usr/bin/env python
"""
Examples:
	#output SNP plots by ranks according to results_by_gene's id=260 (one particular phenotype).
	DrawSNPRegion.py -e 260 -l 28 -L /Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv  -o /Network/Data/250k/tmp-yh/snp_region/ -i phenotype_1_c10_f200
	
	#output SNP plots by ranks according to results_by_gene (analysis_method_id=7-Emma, call_method_id=17) (covering all phenotypes)
	DrawSNPRegion.py -l 28 -L /Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv -o /Network/Data/250k/tmp-yh/snp_region_all/ -i phenotype_1_c10_f200
	
	#take LD_info and gene_annotation from a file contains the pickled data structure
	DrawSNPRegion.py -i ./banyan_fs/tmp/GWA_res_FT.csv -l 28 -D /Network/Data/250k/tmp-yh/call_method_17_LD_m0.1_n0.1_m40000 -o /Network/Data/250k/tmp-yh/snp_region/  -j /tmp/at_gene_model_pickelf
	
Description:
	2008-09-24 program to draw pvalues, gene-models, LD around one SNP.
		Top panel is pvalues from all different methods. Margarita and RF's values will be normalized in range of KW.
		Middle panel is gene model. Displaying CDS, intron, strand.
		Bottom panel is LD of all pairwise SNPs in that range.
	
	The input file contains at least 3 columns, chromosome/chr, position/pos, phenotype_id/phenot_id. tab or comma-delimited.
		The 1st row is a header telling which column is which. No particular order is required. The name could either be full (chromosome) or short(chr).
		Output of FindTopSharedGenes.py conforms to this standard.
		
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr, GeneModel
import Stock_250kDB
from sets import Set
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment
from GeneListRankTest import GeneListRankTest	#GeneListRankTest.getGeneList()
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6
from matplotlib.patches import Polygon, CirclePolygon
import pylab
from pymodule.yh_matplotlib_artists import ExonIntronCollection
import ImageColor
import numpy
from pymodule.latex import outputMatrixInLatexTable
from pymodule.latex import outputFigureInLatex

class LD_statistic(object):
	"""
	2008-09-30
		a class to get label/name conveniently given which LD_statistic is chosen
	"""
	which_LD_statistic2name = {1:'r2', 2:'D_prime', 3:'D'}
	which_LD_statistic2label = {1:'r^2', 2:"|D^'|", 3:'D'}
	
	def get_name(cls, which_LD_statistic):
		return cls.which_LD_statistic2name.get(which_LD_statistic)
	get_name = classmethod(get_name)
	
	def get_label(cls, which_LD_statistic):
		return cls.which_LD_statistic2label.get(which_LD_statistic)
	get_label = classmethod(get_label)


class SNPPassingData(PassingData):
	chromosome = None
	position = None
	def __cmp__(self, other):
		"""
		2008-09-30
			define how to compare to itself
		"""
		return cmp((self.chromosome, self.position), (other.chromosome, other.position))


class GeneSNPPhenotypeAssoData(object):
	"""
	2008-10-02
		the SNP could also denote a region
	2008-10-02
		a data structure to hold gene-snp-phenotype association data got from input file
		iterator not fully tested, ...
	"""
	def __init__(self, **keywords):
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
		self.snp2phenotype_id_ls = {}
		self.gene_id2snp_data = {}
		self.gene_id2phenotype_id_set = {}
		
		self._gene_id_ls = None
		self.gene_index = -1
		self.iter_block = []
		self.phenotype_id2phenotype = {}	#caches for return_matrix_of_snp_descriptions()
	
	def return_phenotype_short_name(self, phenotype_id):
		"""
		2008-10-1
			just for convenience, for DrawSNPRegion.run()
		"""
		phenotype = self.phenotype_id2phenotype.get(phenotype_id)
		if not phenotype:
			phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_id)
			self.phenotype_id2phenotype[phenotype_id] = phenotype
		return phenotype.short_name
	
	def return_phenotype_id2snp_ls(self):
		phenotype_id2snp_ls = {}
		for snp, phenotype_id_set in self.snp2phenotype_id_ls.iteritems():
			for phenotype_id in phenotype_id_set:
				if phenotype_id not in phenotype_id2snp_ls:
					phenotype_id2snp_ls[phenotype_id] = []
				phenotype_id2snp_ls[phenotype_id].append(snp)
		self.phenotype_id2snp_ls = phenotype_id2snp_ls
		return phenotype_id2snp_ls
	
	def no_of_genes(self):
		return len(self.gene_id2snp_data)
	
	no_of_genes = property(no_of_genes)
	
	def no_of_snps(self):
		return len(self.snp2phenotype_id_ls)
	
	no_of_snps = property(no_of_snps)
	
	def gene_id_ls(self):
		if self._gene_id_ls == None:
			self._gene_id_ls = self.gene_id2snp_data.keys()
			self._gene_id_ls.sort()
		return self._gene_id_ls
	gene_id_ls = property(gene_id_ls)
	
	def addGeneSNPPhenotype(self, gene_id, snp, phenotype_id):
		this_snp = (snp.chromosome, snp.position, snp.snps_id, snp.stop)
		if gene_id is not None:
			if gene_id not in self.gene_id2snp_data:
				self.gene_id2snp_data[gene_id] = {}
			snp_data = self.gene_id2snp_data[gene_id]
			if this_snp not in snp_data:
				snp_data[this_snp] = Set()
			snp_data[this_snp].add(phenotype_id)
			
			if gene_id not in self.gene_id2phenotype_id_set:
				self.gene_id2phenotype_id_set[gene_id] = Set()
			self.gene_id2phenotype_id_set[gene_id].add(phenotype_id)
		
		if this_snp not in self.snp2phenotype_id_ls:
			self.snp2phenotype_id_ls[this_snp] = Set()
		self.snp2phenotype_id_ls[this_snp].add(phenotype_id)
	
	def get_no_of_phenotypes_given_gene_id(self, gene_id):
		"""
		"""
		phenotype_id_set = self.gene_id2phenotype_id_set.get(gene_id)
		if phenotype_id_set:
			return len(phenotype_id_set)
		else:
			return 0
	
	snp_desc_names = ['chromosome', 'position', 'phenotypes']
	def return_matrix_of_snp_descriptions(self, gene_id):
		"""
		2008-10-02
			return a 2-D list
		"""
		snp_data = self.gene_id2snp_data.get(gene_id)
		matrix_of_snp_descriptions = []
		if not snp_data:
			return matrix_of_snp_descriptions
		snp_ls = snp_data.keys()
		snp_ls.sort()
		for this_snp in snp_ls:
			snp_desc_ls = []
			for snp_desc_name in self.snp_desc_names:
				if snp_desc_name=='phenotypes':
					element_ls = []
					for phenotype_id in snp_data[this_snp]:
						phenotype = self.phenotype_id2phenotype.get(phenotype_id)
						if not phenotype:
							phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_id)
							self.phenotype_id2phenotype[phenotype_id] = phenotype
						element_ls.append('%s(id %s)'%(phenotype.short_name, phenotype_id))
					element = ', '.join(element_ls)
				elif snp_desc_name=='chromosome':
					element = this_snp[0]
				elif snp_desc_name=='position':
					element = this_snp[1]
				else:
					element = ''
				snp_desc_ls.append(element)
			matrix_of_snp_descriptions.append(snp_desc_ls)
		return matrix_of_snp_descriptions
	
	def return_snp_data_given_gene_id(self, gene_id):
		return self.gene_id2snp_data.get(gene_id)
	
	def __iter__(self):
		return self
	
	def next(self):
		"""
		2008-10-02
			iterator not fully tested, ...
		"""
		self.read()
		return self.iter_block
	
	def read(self):
		"""
		2008-10-02
			iterator not fully tested, ...
		"""
		self.gene_index += 1
		if self.gene_index<self.no_of_genes:
			gene_id = self.gene_id_ls[self.gene_index]
			self.iter_block = [gene_id, self.gene_id2snp_data[gene_id]]
		if self.gene_index==len(self.gene_id_ls):
			self.gene_index = -1
			raise StopIteration
		
class DrawSNPRegion(GeneListRankTest):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('LD_fname', 0, ): [None, 'L', 1, 'the file containing LD info, output of MpiLD.py', ],\
							("min_distance", 1, int): [20000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [0, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency.'],\
							("list_type_id", 0, int): [0, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							('results_directory', 0, ):[None, 't', 1, 'The results directory. Default is None. use the one given by db.'],\
							("input_fname", 1, ): [None, 'i', 1, 'Filename which contains at least 3 columns: chromosome, position, phenotype_id. '],\
							("output_dir", 1, ): [None, 'o', 1, 'directory to store all images'],\
							('call_method_id', 0, int):[17, '', 1, 'Restrict results based on this call_method. Default is no such restriction.'],\
							('analysis_method_id', 0, int):[7, 'a', 1, 'Restrict results based on this analysis_method. Default is no such restriction.'],\
							('no_of_top_hits', 1, int): [1000, 'f', 1, 'how many number of top hits based on score or -log(pvalue).'],\
							("which_LD_statistic", 1, int): [2, 'w', 1, 'which LD_statistic to plot, 1=r2, 2=|D_prime|, 3=|D|'],\
							("LD_info_picklef", 0, ): [None, 'D', 1, 'given the option, If the file does not exist yet, store a pickled LD_info into it (min_MAF and min_distance will be attached to the filename). If the file exists, load LD_info out of it.'],\
							("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
							("latex_output_fname", 0, ): [None, 'x', 1, 'a file to store the latex of gene description tables and figures'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 1, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-09-24
		"""
		GeneListRankTest.__init__(self, **keywords)
		#from pymodule import ProcessOptions
		#self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	
	def get_LD(self, LD_fname, min_MAF, min_gap=40000):
		"""
		2008-09-29
			add min_MAF, min_gap, which_LD_statistic
			min_gap is the distance allowed between two SNPs
			which_LD_statistic returns appropriate LD statistic
		2008-09-24
		"""
		sys.stderr.write("Reading in LD info from %s ...\n"%(LD_fname))
		snp_pair2r2 = {}
		reader = csv.reader(open(LD_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = map(int, snp1)
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = map(int, snp2)
			if snp1[0]==snp2[0] and abs(snp1[1]-snp2[1])<=min_gap:	#on the same chromosome, and less than a certain distance
				allele1_freq = float(row[col_name2index['allele1_freq']])
				allele2_freq = float(row[col_name2index['allele2_freq']])
				if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
					r2 = float(row[col_name2index['r2']])
					D_prime = float(row[col_name2index['D_prime']])
					D = float(row[col_name2index['D']])
					if snp1<snp2:
						snp_pair = (snp1[0], snp1[1], snp2[0], snp2[1])
					else:
						snp_pair = (snp2[0], snp2[1], snp1[0], snp1[1])
					snp_pair2r2[snp_pair] = PassingData(r2=r2, D_prime=D_prime, D=D)
			counter += 1
			if counter%100000==0:
				sys.stderr.write('%s\t%s'%('\x08'*100, counter))
				if self.debug>0:
					break
					pass
		LD_info = PassingData(snp_pair2r2=snp_pair2r2)
		sys.stderr.write("%s entries. Done.\n"%len(snp_pair2r2))
		return LD_info
	
	def dealLD_info(self, LD_info_picklef, LD_fname=None, min_MAF=0.1, min_distance=20000):
		"""
		2008-09-30
			if the LD_info_picklef does not exist yet, store a pickled LD_info into it. If the file exists, load LD_info out of it.
		"""
		sys.stderr.write("Dealing with LD_info ...")
		if LD_info_picklef:
			if os.path.isfile(LD_info_picklef):	#if this file is already there, suggest to un-pickle it.
				picklef = open(LD_info_picklef)
				LD_info = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle into it
				LD_info = self.get_LD(LD_fname, min_MAF, min_distance*2+100)	#min_gap is 2 X min_distance + a few more bases
				picklef = open('%s_n%s_m%s'%(LD_info_picklef, min_MAF, min_distance), 'w')
				cPickle.dump(LD_info, picklef, -1)
				picklef.close()
		else:
			LD_info = self.get_LD(LD_fname, min_MAF, min_distance*2+100)
		sys.stderr.write("Done.\n")
		return LD_info
	
	def getGeneAnnotation(self):
		"""
		2008-10-01
			substitute GenomeBrowser.get_gene_id2model() with GenomeDB.GenomeDatabase.get_gene_id2model()
			
		2008-09-24
		"""
		from transfac.src import GenomeDB
		db = GenomeDB.GenomeDatabase(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database='genome', schema=self.schema)
		db.setup(create_tables=False)
		gene_id2model, chr_id2gene_id_ls = db.get_gene_id2model(tax_id=3702)
		gene_annotation = PassingData()
		gene_annotation.gene_id2model = gene_id2model
		gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
		return gene_annotation
	
	def dealWithGeneAnnotation(self, gene_annotation_picklef):
		"""
		2008-10-01
			similar to dealLD_info()
		"""
		sys.stderr.write("Dealing with gene_annotation ...")
		if gene_annotation_picklef:
			if os.path.isfile(gene_annotation_picklef):	#if this file is already there, suggest to un-pickle it.
				picklef = open(gene_annotation_picklef)
				gene_annotation = cPickle.load(picklef)
				del picklef
			else:	#if the file doesn't exist, but the filename is given, pickle into it
				gene_annotation = self.getGeneAnnotation()	#min_gap is 2 X min_distance + a few more bases
				picklef = open(gene_annotation_picklef, 'w')
				cPickle.dump(gene_annotation, picklef, -1)
				picklef.close()
		else:
			gene_annotation = self.getGeneAnnotation()
		sys.stderr.write("Done.\n")
		return gene_annotation
	
	def getSimilarGWResultsGivenResultsByGene(self, phenotype_method_id, call_method_id, results_directory=None, \
											analysis_method_id_ls=[1,5,6,7]):
		"""
		2008-09-30
			use phenotype_method_id and call_method_id to find out all genome-wide results
		2008-09-24
		"""
		sys.stderr.write("Getting results with phenotype=%s and call_method=%s ..."%(phenotype_method_id, call_method_id))
		if self.debug:
			analysis_method_id_set = Set([1,7])
		else:
			analysis_method_id_set = Set(analysis_method_id_ls)
		rows = Stock_250kDB.ResultsMethod.query.filter_by(phenotype_method_id=phenotype_method_id, call_method_id=call_method_id)
		analysis_method_id2gwr = {}
		for rm in rows:
			if rm.analysis_method_id in analysis_method_id_set:
				genome_wide_result = self.getResultMethodContent(rm, results_directory, min_MAF=0, construct_chr_pos2index=True)
				analysis_method_id2gwr[rm.analysis_method_id] = genome_wide_result
		sys.stderr.write("Done.\n")
		return analysis_method_id2gwr
	
	def getSNPInfo(self, db):
		"""
		2008-09-24
			in order
		"""
		sys.stderr.write("Getting info of all SNPs in chromosomal order ...")
		chr_pos_ls = []
		chr_pos2index = {}
		snps_id2index = {}
		i = 0
		block_size = 50000
		rows = db.metadata.bind.execute("select id, chromosome, position from %s where end_position is null"%Stock_250kDB.Snps.table.name)
		#.query.offset(i).limit(block_size)
		#while rows.count()!=0:
		for row in rows:
			chr_pos = (row.chromosome, row.position)
			chr_pos_ls.append(chr_pos)
			chr_pos2index[chr_pos] = len(chr_pos2index)
			snps_id2index[row.id] = len(snps_id2index)
			i += 1
		#	if self.debug and i>40000:
		#		break
		#	rows = Stock_250kDB.Snps.query.offset(i).limit(block_size)
		snp_info = PassingData()
		snp_info.chr_pos_ls = chr_pos_ls
		snp_info.chr_pos2index = chr_pos2index
		snp_info.snps_id2index = snps_id2index
		sys.stderr.write("Done.\n")
		return snp_info
	
	def add_mid_point(self, chr_pos_ls, chr_pos2adjacent_window):
		"""
		2008-09-24
			called by getSNPsAroundThisSNP()
		"""
		new_chr_pos = chr_pos_ls[-1]
		old_chr_pos = chr_pos_ls[-2]
		if old_chr_pos not in chr_pos2adjacent_window:
			chr_pos2adjacent_window[old_chr_pos] = []
		if new_chr_pos not in chr_pos2adjacent_window:
			chr_pos2adjacent_window[new_chr_pos] = []
		mid_point = (new_chr_pos[1]+old_chr_pos[1])/2.
		chr_pos2adjacent_window[old_chr_pos].append(mid_point)
		chr_pos2adjacent_window[new_chr_pos].append(mid_point)
	
	def getSNPsAroundThisSNP(self, this_snp, snp_info, min_distance=20000):
		"""
		2008-09-24
		"""
		sys.stderr.write("\t Get SNPs around this snp ...")
		#chr_pos = snp_info.chr_pos_ls[snp_info.snps_id2index[this_snp.snps_id]]
		chromosome, position = this_snp.chromosome, this_snp.position
		chr_pos_ls = []
		chr_pos2adjacent_window = {}
		j = 0
		for i in range(min_distance*2):
			new_pos = position - min_distance + i
			new_chr_pos = (chromosome, new_pos)
			if new_chr_pos in snp_info.chr_pos2index:
				chr_pos_ls.append(new_chr_pos)
				if j!=0:
					self.add_mid_point(chr_pos_ls, chr_pos2adjacent_window)
				j += 1
		#deal with the leftest point of the 1st chr_pos
		chr_pos = chr_pos_ls[0]
		window_size = chr_pos2adjacent_window[chr_pos][0]-chr_pos[1]
		chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		
		#deal with the rightest point of the 1st chr_pos
		chr_pos = chr_pos_ls[-1]
		window_size = chr_pos[1] - chr_pos2adjacent_window[chr_pos][0]
		chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		snp_region = PassingData(chr_pos_ls=chr_pos_ls, chr_pos2adjacent_window=chr_pos2adjacent_window, center_snp=this_snp)
		sys.stderr.write("Done.\n")
		return snp_region
	
	analysis_method_id2color = {1:'b',
							5:'r',
							6:'g',
							7:'k'}
	
	analysis_method_id2name = {1:'KW',
							5:'Margarita',
							6:'RF',
							7:'Emma'}
	
	def getXY(self, snps_within_this_region, analysis_method_id2gwr=None, analysis_method_id=None, LD_info=None, which_LD_statistic=2):
		"""
		2008-10-1
			return LD value (with regard to the center SNP)
		2008-09-24
			of GW results, get values for each SNP position, adjust value for analysis_method_id=5,6
		"""
		x_ls = []
		y_ls = []
		if analysis_method_id is not None:
			gwr = analysis_method_id2gwr[analysis_method_id]
		elif LD_info is None:	#2008-10-01
			return x_ls, y_ls
		
		if analysis_method_id2gwr and 1 in analysis_method_id2gwr:
			ref_gwr = analysis_method_id2gwr[1]
		elif analysis_method_id2gwr and 7 in analysis_method_id2gwr:
			ref_gwr = analysis_method_id2gwr[7]
		else:
			ref_gwr = None
		for chr_pos in snps_within_this_region.chr_pos_ls:
			if analysis_method_id is not None:
				data_obj = gwr.get_data_obj_by_chr_pos(chr_pos[0], chr_pos[1])
			else:	#2008-10-1 fake a data_obj for LD
				chr_pos1 = (snps_within_this_region.center_snp.chromosome, snps_within_this_region.center_snp.position)
				chr_pos2 = chr_pos
				if chr_pos1<chr_pos2:
					snp_pair = (chr_pos1[0], chr_pos1[1], chr_pos2[0], chr_pos2[1])
				else:
					snp_pair = (chr_pos2[0], chr_pos2[1], chr_pos1[0], chr_pos1[1])
				if snp_pair in LD_info.snp_pair2r2:
					LD_stat = getattr(LD_info.snp_pair2r2[snp_pair], LD_statistic.get_name(which_LD_statistic), None)
					LD_stat = abs(LD_stat)	#D_prime, D need abs()
				elif chr_pos1==chr_pos2:	#fake a perfect LD point
					LD_stat = 1.
				else:	#skip this point
					continue
				data_obj = PassingData(value=LD_stat)
				
			if data_obj is not None:
				x_ls.append(chr_pos[1])
				if (analysis_method_id==5 or analysis_method_id==6) and ref_gwr:
					value = (data_obj.value-gwr.min_value)/(gwr.max_value-gwr.min_value)*(ref_gwr.max_value-ref_gwr.min_value)
				else:
					value = data_obj.value
				y_ls.append(value)
		return x_ls, y_ls
	
	def drawPvalue(self, ax1, axe_LD, axe_LD_center_SNP, snps_within_this_region, analysis_method_id2gwr, LD_info=None, which_LD_statistic=2):
		"""
		2008-10-1
			refine the LD line on axe_LD_center_SNP, with markersize=2 and etc.
		2008-10-1
			draw LD of other SNPs w.r.t the center SNP in a different axe
		2008-09-24
		"""
		sys.stderr.write("\t Drawing pvalues  ...")
		analysis_method_id_ls = analysis_method_id2gwr.keys()
		analysis_method_id_ls.sort()
		pscatter_ls = []
		legend_ls = []
		for analysis_method_id in analysis_method_id_ls:
			gwr = analysis_method_id2gwr[analysis_method_id]
			x_ls, y_ls = self.getXY(snps_within_this_region, analysis_method_id2gwr, analysis_method_id)
			if x_ls and y_ls:
				pscatter = ax1.scatter(x_ls, y_ls, s=10, linewidth=0.6, edgecolor=self.analysis_method_id2color[analysis_method_id], facecolor='w')
				legend_ls.append(self.analysis_method_id2name[analysis_method_id])
				pscatter_ls.append(pscatter)
		if LD_info:	#draw LD with regard to the center SNP
			x_ls, y_ls = self.getXY(snps_within_this_region, LD_info=LD_info, which_LD_statistic=which_LD_statistic)
			if x_ls and y_ls:
				apl = axe_LD_center_SNP.plot(x_ls, y_ls, 'c-o', linewidth=0.5, markersize=2, alpha=0.6, markerfacecolor='w')
				legend_ls.append(r'$%s$ with center SNP'%LD_statistic.get_label(which_LD_statistic))
				pscatter_ls.append(apl[0])
		ax1.set_ylabel(r'-log(pvalue)/normalized score')
		axe_LD_center_SNP.set_ylabel(r'$%s$ with center SNP'%LD_statistic.get_label(which_LD_statistic))
		axe_LD.legend(pscatter_ls, legend_ls, loc='lower left', handlelen=0.02)	#cut the legend length to 0.02, default 0.05 (5% of x-axis).
		sys.stderr.write("Done.\n")
		#return legend_ls
	
	gene_desc_names = ['gene_id', 'gene_symbol', 'type_of_gene', 'chr', 'start', 'stop', 'protein_label', 'protein_comment', 'protein_text']
	
	def returnGeneDescLs(self, gene_desc_names, gene_model, gene_commentary=None, cutoff_length=200):
		"""
		2008-10-02
		"""
		#2008-10-01	get the gene descriptions
		gene_desc_ls = []
		if gene_commentary is None:	#take gene_model
			gene_commentary = gene_model
		for gene_desc_name in gene_desc_names:
			if gene_desc_name=='protein_label':	#if not available, get the rna label
				element = getattr(gene_commentary, gene_desc_name, None)
				if not element:
					element = getattr(gene_commentary, 'label', '')
					if cutoff_length is not None and type(element) ==str:
						element = element[:cutoff_length/2]
				if cutoff_length is not None and type(element) ==str:	#cut label length half
						element = element[:cutoff_length/2]
			elif gene_desc_name=='protein_comment':
				element = getattr(gene_commentary, gene_desc_name, None)
				if not element:
					element = getattr(gene_commentary, 'comment', '')
					if cutoff_length is not None and type(element) ==str:
						element = element[:cutoff_length/2]
			elif gene_desc_name=='gene_symbol' or gene_desc_name=='chromosome':	#have to get it from the gene_model
				element = getattr(gene_model, gene_desc_name, '')
			elif gene_desc_name=='type_of_gene':
				if getattr(gene_commentary, 'protein_label', None) is not None:
					element = getattr(gene_model, gene_desc_name, '')
				else:	#it doesn't have protein, get gene_commentary_type
					element = getattr(gene_commentary, 'gene_commentary_type', '')
			elif gene_desc_name=='chr':
				element = getattr(gene_model, 'chromosome', '')
			else:
				element = getattr(gene_commentary, gene_desc_name, '')
			if cutoff_length is not None and type(element)==str and len(element)>cutoff_length:	#only 20 characters
				element = element[:cutoff_length]
			gene_desc_ls.append(element)
		return gene_desc_ls
	
	def plot_one_gene(self, ax, gene_id, param_data, base_y_value=1):
		"""
		2008-10-02
			plot all different forms of a gene out.
		2008-10-01
			gene_model is changed. now it's from GenomeDB.py
		2008-09-24
			draw a single gene on the canvas, adapted from GenomeBrowser.plot_one_gene()
		"""
		gene_model = param_data.gene_id2model.get(gene_id)
		if gene_model:
			if len(gene_model.gene_commentaries)==0:
				gene_commentaries = [gene_model]	#fake one here
			else:
				gene_commentaries = gene_model.gene_commentaries
			for gene_commentary in gene_model.gene_commentaries:	#multiple commentary
				gene_desc_ls = self.returnGeneDescLs(self.gene_desc_names, gene_model, gene_commentary)
				param_data.matrix_of_gene_descriptions.append(gene_desc_ls)
				
				if getattr(gene_commentary, 'box_ls', None):
					box_ls = gene_commentary.box_ls
				elif gene_commentary.start and gene_commentary.stop:	#no box_ls, just use start, stop
					box_ls = [(gene_commentary.start, gene_commentary.stop, 'exon')]
				elif gene_model.start and gene_model.stop:	#use gene_model's coordinate
					box_ls = [(gene_model.start, gene_model.stop, 'exon')]
				else:	#no coordinate. skip
					sys.stderr.write("Warning: gene %s has either no start or no stop.\n"%(gene_commentary.gene_id))
					continue
				if gene_id in param_data.candidate_gene_set:
					gene_symbol_color = 'b'
				else:
					gene_symbol_color = 'k'
				y_value = base_y_value+param_data.no_of_genes_drawn%param_data.gene_position_cycle	#cycling through the y position to avoid clogging
				g_artist = ExonIntronCollection(box_ls, y=y_value, strand=gene_model.strand, width=param_data.gene_width, alpha=0.3, \
											picker=True, linewidths=0.7, box_line_widths=0.3)
				ax.add_artist(g_artist)
				text_start_pos = box_ls[-1][1]
				#mid_point = (c_start_ls[0]+c_end_ls[-1])/2.
				ax.text(text_start_pos+param_data.gene_box_text_gap, y_value, gene_model.gene_symbol, size=4, \
					color=gene_symbol_color, alpha=0.8, verticalalignment='center')
				param_data.no_of_genes_drawn += 1
				
					
	def drawGeneModel(self, ax, snps_within_this_region, gene_annotation, candidate_gene_set, gene_width=1.0, \
					gene_position_cycle=4, base_y_value=1, gene_box_text_gap=100):
		"""
		2008-10-01
		2008-09-24
		"""
		sys.stderr.write("\t Drawing gene model  ...")
		left_chr, left_pos = snps_within_this_region.chr_pos_ls[0]
		right_chr, right_pos = snps_within_this_region.chr_pos_ls[-1]
		left_chr = str(left_chr)
		right_chr = str(right_chr)
		param_data = PassingData(gene_id2model=gene_annotation.gene_id2model, candidate_gene_set=candidate_gene_set, \
								gene_width=gene_width, gene_position_cycle=gene_position_cycle, no_of_genes_drawn=0, \
								gene_box_text_gap=gene_box_text_gap, matrix_of_gene_descriptions = [])
		for gene_id in gene_annotation.chr_id2gene_id_ls[left_chr]:
			gene_model = gene_annotation.gene_id2model[gene_id]
			if gene_model.start!=None and gene_model.stop!=None and gene_model.stop>left_pos:
				if left_chr==right_chr:	#same chromosome
					if gene_model.start>right_pos:	#totally out of range, skip it
						continue
				self.plot_one_gene(ax, gene_id, param_data, base_y_value=base_y_value)
		if left_chr!=right_chr:
			for gene_id in gene_annotation.chr_id2gene_id_ls[right_chr]:
				gene_model = gene_annotation.gene_id2model[gene_id]
				if gene_model.start!=None and gene_model.stop!=None and gene_model.start<right_pos:
					self.plot_one_gene(ax, gene_id, param_data, base_y_value=base_y_value)
		sys.stderr.write("Done.\n")
		return param_data.matrix_of_gene_descriptions
	
	def drawLD(self, axe_gene_model, ax2, snps_within_this_region, LD_info, gene_model_min_y,\
				gene_model_max_y, which_LD_statistic=1):
		"""
		2008-09-28
			represent r2 by HSL color, rather than a grayscale intensity. matplotlib doesn't support hsl notation (i'm not aware).
			Use ImageColor.getrgb to convert hsl representation to RGB format.
		2008-09-24
			draw LD in the bottom axe
		"""
		sys.stderr.write("\t Drawing LD info  ...")
		no_of_snps = len(snps_within_this_region.chr_pos_ls)
		left_chr, left_pos = snps_within_this_region.chr_pos_ls[0]
		right_chr, right_pos = snps_within_this_region.chr_pos_ls[-1]
		#ax1.hlines(y_value, left_pos, right_pos, linewidth=0.3)
		for i in range(no_of_snps):
			chr_pos1 = snps_within_this_region.chr_pos_ls[i]
			axe_gene_model.vlines(chr_pos1[1], gene_model_min_y, gene_model_max_y, linestyle='dashed', alpha=0.3, linewidth=0.3)
			for j in range(i+1, no_of_snps):
				chr_pos2 = snps_within_this_region.chr_pos_ls[j]
				if chr_pos1<chr_pos2:
					snp_pair = (chr_pos1[0], chr_pos1[1], chr_pos2[0], chr_pos2[1])
				else:
					snp_pair = (chr_pos2[0], chr_pos2[1], chr_pos1[0], chr_pos1[1])
				if snp_pair in LD_info.snp_pair2r2:
					LD_stat = getattr(LD_info.snp_pair2r2[snp_pair], LD_statistic.get_name(which_LD_statistic), None)
					LD_stat = abs(LD_stat)	#D_prime, D need abs()
					s11, s12 = snps_within_this_region.chr_pos2adjacent_window[chr_pos1]
					s21, s22 = snps_within_this_region.chr_pos2adjacent_window[chr_pos2]
					#draw a plot to understand this. a parallegram below SNPs, it's the intersection of two adjacent windows 
					x1 = (s12+s21)/2.
					y1 = (s12-s21)/2.
					x2 = (s12+s22)/2.
					y2 = (s12-s22)/2.
					x3 = (s11+s22)/2.
					y3 = (s11-s22)/2.
					x4 = (s11+s21)/2.
					y4 = (s11-s21)/2.
					xs = [x1, x2, x3, x4]
					ys = [y1, y2, y3, y4]
					fc = self.r2ToRGBColor(LD_stat)
					poly = Polygon(zip(xs, ys), facecolor=fc, linewidth=0)
					ax2.add_patch(poly)
		sys.stderr.write("Done.\n")
	
	def r2ToRGBColor(self, r2):
		"""
		2008-09-29
			convert r2 to matplotlib RGB color
		"""
		hue_value = int(round((1-r2)*255))	#r2 is [0,1], 255 is the maximum hue value. lower r2, higher hue. r2=0 -> blue, r2=1 -> red
		fc = ImageColor.getrgb('hsl(%s'%hue_value+',100%,50%)')	#"hsl(hue, saturation%, lightness%)" where hue is the colour given as an
		# angle between 0 and 360 (red=0, green=120, blue=240),
		#saturation is a value between 0% and 100% (gray=0%, full color=100%), and lightness is a value between 0% and 100% (black=0%, normal=50%, white=100%).
		fc = [color_value/255. for color_value in fc]	#matplotlib accepts rgb in [0-1] range
		return fc
	
	def drawLDLegend(self, ax, which_LD_statistic):
		"""
		2008-09-29
			left half is 10 color patchs, right half is r2 from 0 to 1
		"""
		sys.stderr.write("\t Drawing LD legend  ...")
		xs = [0,0,0.5,0.5]	#X-axis value for the 4 points of the rectangle starting from upper left corner. 
		ys = numpy.array([0.1,0,0,0.1])
		for i in range(10):
			r2 = ys[0]	#upper bound of ys correspond to r2
			fc = self.r2ToRGBColor(r2)
			poly = Polygon(zip(xs, ys), facecolor=fc, linewidth=0)
			ax.add_patch(poly)
			if i%2==0:	#every other band
				ax.text(0.6, ys[0], '%.1f'%r2, horizontalalignment ='left', verticalalignment='top', size=4)
			ys += 0.1	#increase y-axis
		ax.text(0.6, 1.1, r"$%s=1.0$"%LD_statistic.get_label(which_LD_statistic), horizontalalignment ='left', verticalalignment='top', size=4)
		sys.stderr.write("Done.\n")
	
	def drawRegionAroundThisSNP(self, phenotype_method_id, this_snp, candidate_gene_set, gene_annotation, snp_info, analysis_method_id2gwr, \
							LD_info, output_dir, which_LD_statistic, snp_region=None, min_distance=40000, list_type_id=None):
		"""
		2008-10-01
			remove the frame of ax1 and add a grid to ax1
			leave axe_gene_model's xticks there as otherwise ax1's xticks will go with it as they share xticks.
		2008-10-01
			draw gene models on a separate axe,
			add a twinx axe to draw LD w.r.t the center SNP
			if output_dir is not a directory, it's treated as a filename.
		2008-09-24
		"""
		sys.stderr.write("Drawing region ... \n")
		phenotype = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
		if not os.path.isdir(output_dir):
			output_fname_prefix = output_dir
		else:
			#list_type = Stock_250kDB.GeneListType.get(list_type_id)
			fname_basename = 'snp_%s_%s_id_%s_phenotype_%s_%s'%\
										(this_snp.chromosome, this_snp.position, this_snp.snps_id, phenotype.id, phenotype.short_name)
			fname_basename = fname_basename.replace('/', '_')
			output_fname_prefix = os.path.join(output_dir, fname_basename)
		if snp_region:
			snps_within_this_region = snp_region
		elif getattr(this_snp, 'stop', None):
			snps_within_this_region = self.findSNPsInRegion(snp_info, this_snp.chromosome, this_snp.position, this_snp.stop)
		else:
			snps_within_this_region = self.getSNPsAroundThisSNP(this_snp, snp_info, min_distance)
		pylab.clf()
		#fig = pylab.figure()
		ax1 = pylab.axes([0.1, 0.6, 0.8, 0.35], frameon=False)	#left gap, bottom gap, width, height, axes for pvalue, gene models
		ax1.grid(True, alpha=0.3)
		ax1.set_xticklabels([])	#remove xtick labels on ax1 because axe_LD's xtick labels cover this.
		axe_LD_center_SNP = pylab.twinx()	#axes for LD with center SNP, copy ax1's
		axe_LD_center_SNP.set_xticklabels([])
		axe_gene_model = pylab.axes([0.1, 0.5, 0.8, 0.1], frameon=False, sharex=ax1)
		#axe_gene_model.set_xticks([])	#this will set ax1's xticks off as well because the x-axis is shared.
		axe_gene_model.set_yticks([])
		axe_LD = pylab.axes([0.1, 0.05, 0.8, 0.45], frameon=False)	#axes for LD
		axe_LD_legend = pylab.axes([0.8, 0.08, 0.1, 0.13], frameon=False)	#axes for the legend of LD
		axe_LD_legend.set_xticks([])
		axe_LD_legend.set_yticks([])
		fig_title = 'SNP chr %s. pos %s.'%(this_snp.chromosome, this_snp.position)
		if getattr(this_snp, 'stop', None):
				fig_title += ' - %s. '%this_snp.stop
		fig_title += "Phenotype %s (id=%s)."%(phenotype.short_name, phenotype.id)		
		ax1.title.set_text(fig_title)	#main title using this snp.
		self.drawPvalue(ax1, axe_LD, axe_LD_center_SNP, snps_within_this_region, analysis_method_id2gwr, LD_info, which_LD_statistic)
		gene_position_cycle = 5
		base_y_value = 1
		gene_width=0.8
		gene_box_text_gap = min_distance*2*0.005
		matrix_of_gene_descriptions = self.drawGeneModel(axe_gene_model, snps_within_this_region, gene_annotation, candidate_gene_set, gene_width=gene_width, \
						gene_position_cycle=gene_position_cycle, base_y_value=base_y_value, gene_box_text_gap=gene_box_text_gap)
		gene_model_min_y = base_y_value-gene_width
		gene_model_max_y = gene_position_cycle + base_y_value -1 + gene_width	#"-1" because genes never sit on y=gene_position_cycle + base_y_value
		self.drawLD(axe_gene_model, axe_LD, snps_within_this_region, LD_info, gene_model_min_y=gene_model_min_y,\
				gene_model_max_y=gene_model_max_y, which_LD_statistic=which_LD_statistic)
		self.drawLDLegend(axe_LD_legend, which_LD_statistic)
		#adjust x, y limits and etc
		ax1_ylim = ax1.get_ylim()
		ax1.set_ylim((0, ax1_ylim[1]))	#set ax1 to 0 to sit right above axe_gene_model
		axe_gene_model.set_ylim((gene_model_min_y, gene_model_max_y))	#LD panel right under gene models
		
		axe_LD.set_xlim(ax1.get_xlim())	#make the axe_LD and ax1 within the same X range
		axe_LD_x_span = (axe_LD.get_xlim()[1]-axe_LD.get_xlim()[0])
		axe_LD.set_ylim((-axe_LD_x_span/2., 0))	#has to force here, don't know why. otherwise it's (0,1)
		axe_LD.set_yticks([])	#remove all Y ticks on LD plot
		
		png_output_fname = '%s.png'%output_fname_prefix
		pylab.savefig(png_output_fname, dpi=400)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		if self.debug:
			pylab.show()
		del ax1, axe_LD_center_SNP, axe_gene_model, axe_LD, axe_LD_legend
		sys.stderr.write("Done.\n")
		after_plot_data = PassingData(png_output_fname=png_output_fname, matrix_of_gene_descriptions=matrix_of_gene_descriptions)
		return after_plot_data
	
	def generate_params(self, param_obj):
		"""
		2008-09-24
			copied from a version of MpiGeneListRankTest.py
		"""
		sys.stderr.write("Generating parameters ...")
		i = 0
		block_size = 5000
		query = Stock_250kDB.ResultsByGene.query
		if param_obj.call_method_id!=0:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(call_method_id=param_obj.call_method_id))
		if param_obj.analysis_method_id!=0 and param_obj.analysis_method_id is not None:
			query = query.filter(Stock_250kDB.ResultsByGene.results_method.has(analysis_method_id=param_obj.analysis_method_id))
		query = query.filter_by(min_distance=param_obj.min_distance).filter_by(get_closest=param_obj.get_closest)
		rows = query.offset(i).limit(block_size)
		results_id_ls = []
		while rows.count()!=0:
			for row in rows:
				results_id_ls.append(row.id)
				i += 1
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("%s results. "%(len(results_id_ls)))
		return results_id_ls
	
	def getSNPsFromInputFile(self, input_fname, snp_info):
		"""
		2008-10-04
			input could contain one more column, 'stop', indicating it is a region, not a single SNP. 'start' becomes alternative name for 'position'.
		2008-10-02
			use GeneSNPPhenotypeAssoData to hold input data
		2008-09-30
			read input from a file. check module doc for format.
		"""
		sys.stderr.write("Getting input from %s ..."%input_fname)
		input_data = GeneSNPPhenotypeAssoData()
		delimiter = figureOutDelimiter(input_fname)
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		for row in reader:
			phenotype_id_index = col_name2index.get('phenotype_id')
			if phenotype_id_index == None:
				phenotype_id_index = col_name2index.get('phenot_id')
			if phenotype_id_index is not None:
				phenotype_id = int(row[phenotype_id_index])
			else:
				continue
			
			chromosome_index = col_name2index.get('chromosome')
			if chromosome_index ==None:
				chromosome_index = col_name2index.get('chr')
			
			position_index = col_name2index.get('position')
			if position_index==None:
				position_index = col_name2index.get('pos')
				if position_index is None:
					position_index = col_name2index.get('start')
			
			stop_index = col_name2index.get('stop')
			if stop_index is None:
				stop = None
			else:
				stop = int(row[stop_index])
			
			if chromosome_index is not None and position_index is not None:
				chromosome = int(row[chromosome_index])
				position = int(row[position_index])
				if stop:
					snps_id = '%s_%s_%s'%(chromosome, position, stop)
				else:
					snps_id = '%s_%s'%(chromosome, position)
			else:
				try:
					snps_id = int(row[col_name2index['snps_id']])
					chr_pos = snp_info.chr_pos_ls[snp_info.snps_id2index[snps_id]]
					chromosome, position = chr_pos
				except:	#forget it, skip
					continue
			
			gene_id_index = col_name2index.get('gene_id')
			if gene_id_index:
				gene_id = int(row[gene_id_index])
			else:
				gene_id = None
			counter += 1
			#snps_id = int(row[col_name2index['snps_id']])
			#disp_pos = int(row[col_name2index['disp_pos']])
			this_snp = SNPPassingData(chromosome=chromosome, position=position, stop=stop, snps_id=snps_id)
			input_data.addGeneSNPPhenotype(gene_id, this_snp, phenotype_id)
		sys.stderr.write("%s genes. %s total snps. Done.\n"%(input_data.no_of_genes, input_data.no_of_snps))
		return input_data
	
	def findSNPsInRegion(self, snp_info, chromosome, start, stop, center_snp_position=None):
		"""
		2008-10-1
			called by plotSNPRegion()
			find SNPs in this region, if center_snp_position is not given, find one.
			similar to getSNPsAroundThisSNP()
		"""
		sys.stderr.write("Get SNPs in this region ...")
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
					self.add_mid_point(chr_pos_ls, chr_pos2adjacent_window)
				j += 1
		#deal with the leftest point of the 1st chr_pos
		chr_pos = chr_pos_ls[0]
		if chr_pos in chr_pos2adjacent_window:	#if chr_pos_ls has only one element, then nothing would be in chr_pos2adjacent_window
			window_size = chr_pos2adjacent_window[chr_pos][0]-chr_pos[1]
			chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		
		#deal with the rightest point of the 1st chr_pos
		chr_pos = chr_pos_ls[-1]
		if chr_pos in chr_pos2adjacent_window:	#
			window_size = chr_pos[1] - chr_pos2adjacent_window[chr_pos][0]
			chr_pos2adjacent_window[chr_pos] = [chr_pos[1]-window_size, chr_pos[1]+window_size]
		
		center_snp.snps_id = '%s_%s'%(center_snp.chromosome, center_snp.position)
		snp_region = PassingData(chr_pos_ls=chr_pos_ls, chr_pos2adjacent_window=chr_pos2adjacent_window, center_snp=center_snp)
		sys.stderr.write("Done.\n")
		return snp_region
	
	def loadDataStructure(self, gene_annotation_picklef, LD_info_picklef, LD_fname=None, min_MAF=0.1, min_distance=20000, list_type_id=None):
		"""
		2008-10-01
			wrap a few functions up, convenient for both run() and drawSNPRegion()
		"""
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		self.db = db
		snp_info = self.getSNPInfo(db)
		gene_annotation = self.dealWithGeneAnnotation(gene_annotation_picklef)
		LD_info = self.dealLD_info(LD_info_picklef, LD_fname, min_MAF, min_distance)
		if list_type_id:
			candidate_gene_list = self.getGeneList(list_type_id)
			candidate_gene_set = Set(candidate_gene_list)
		else:
			candidate_gene_set = Set()
		return_data = PassingData(gene_annotation=gene_annotation, snp_info=snp_info, LD_info=LD_info, \
								candidate_gene_set=candidate_gene_set)
		return return_data
	
	def drawSNPRegion(self, gene_annotation_picklef, LD_info_picklef, phenotype_method_id, chromosome, start, stop, \
					output_fname, center_snp_position=None, LD_fname=None, which_LD_statistic=2, call_method_id=17, min_MAF=0.1):
		"""
		2008-10-01
			for bjarni to call
		"""
		if not hasattr(self, 'grand_dataStructure'):
			self.grand_dataStructure = self.loadDataStructure(gene_annotation_picklef, LD_info_picklef, \
													LD_fname=LD_fname, min_MAF=min_MAF)
		
		gene_annotation = self.grand_dataStructure.gene_annotation
		snp_info = self.grand_dataStructure.snp_info
		LD_info = self.grand_dataStructure.LD_info
		analysis_method_id2gwr = self.getSimilarGWResultsGivenResultsByGene(phenotype_method_id, call_method_id)
		candidate_gene_set = Set()
		
		snp_region = self.findSNPsInRegion(snp_info, chromosome, start, stop, center_snp_position)
		this_snp = snp_region.center_snp
		this_snp.phenotype_id = phenotype_method_id
		self.drawRegionAroundThisSNP(phenotype_method_id, this_snp, candidate_gene_set, gene_annotation, snp_info, \
									analysis_method_id2gwr, LD_info, output_fname, which_LD_statistic, snp_region)
	
	def handleOneGeneOnePhenotype(self, phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, input_data, \
								param_data, latex_f=None, delete_gwr=False):
		"""
		2008-10-02
			wrap up towards latex output
		"""
		analysis_method_id2gwr = phenotype_id2analysis_method_id2gwr.get(phenotype_id)
		if not analysis_method_id2gwr:
			analysis_method_id2gwr = self.getSimilarGWResultsGivenResultsByGene(phenotype_id, param_data.call_method_id, param_data.results_directory)
			if not delete_gwr:
				phenotype_id2analysis_method_id2gwr[phenotype_id] = analysis_method_id2gwr
		if not analysis_method_id2gwr:	#still nothing, skip
			sys.stderr.write("No gwr resutls available. skip.\n")
			return
		after_plot_data = self.drawRegionAroundThisSNP(phenotype_id, this_snp, param_data.candidate_gene_set, param_data.gene_annotation, \
													param_data.snp_info, \
								analysis_method_id2gwr, param_data.LD_info, param_data.output_dir, param_data.which_LD_statistic, \
								min_distance=param_data.min_distance, list_type_id=param_data.list_type_id)
		param_data.no_of_snps_drawn += 1
		if delete_gwr:
			del analysis_method_id2gwr
		
		if latex_f:
			annot_table_label = 'annottable%s'%param_data.no_of_snps_drawn
			snp_figure_label = 'snpfig%s'%param_data.no_of_snps_drawn
			no_of_genes_in_annot = len(after_plot_data.matrix_of_gene_descriptions)
			snp_label_region = 'SNP at %s, %s'%(this_snp.chromosome, this_snp.position)
			if getattr(this_snp, 'stop', None):
				snp_label_region += ' to %s'%this_snp.stop
			
			caption = 'Annotation of %s genes for Figure~\\ref{%s}. %s. Phenotype is %s(id=%s).'%\
				(no_of_genes_in_annot, snp_figure_label, snp_label_region, \
				input_data.return_phenotype_short_name(phenotype_id), phenotype_id)
			
			if after_plot_data.matrix_of_gene_descriptions:	#if matrix is empty, don't output
				latex_f.write(outputMatrixInLatexTable(after_plot_data.matrix_of_gene_descriptions, caption, annot_table_label, self.gene_desc_names))
				annot_table_ref = 'Annotation of %s genes is in Table~\\ref{%s}.'%\
					(no_of_genes_in_annot, annot_table_label)
			else:	#no mentioning this table if it doesn't exist
				annot_table_ref = ''
			caption = '%s. Phenotype is %s(id=%s). %s'%(snp_label_region, \
						input_data.return_phenotype_short_name(phenotype_id), phenotype_id, annot_table_ref)
			latex_f.write(outputFigureInLatex(after_plot_data.png_output_fname, caption, snp_figure_label))
			latex_f.flush()
	
	def run(self):
		"""
		2008-10-02
			modify to enable latex output of tables and figures
		2008-09-24
		"""
		if self.debug==1:
			import pdb
			pdb.set_trace()
		import gc
		gc.set_threshold(250,10,10)
		grand_dataStructure = self.loadDataStructure(self.gene_annotation_picklef, self.LD_info_picklef, self.LD_fname, \
							self.min_MAF, self.min_distance, self.list_type_id)
		gene_annotation, snp_info, LD_info = grand_dataStructure.gene_annotation, grand_dataStructure.snp_info, grand_dataStructure.LD_info
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		"""
		if not self.results_id_ls:
			param_obj = PassingData(call_method_id=self.call_method_id, analysis_method_id=self.analysis_method_id,\
								min_distance=self.min_distance, get_closest=self.get_closest)
			self.results_id_ls = self.generate_params(param_obj)
		"""
		which_LD_statistic = self.which_LD_statistic
		input_fname = self.input_fname
		latex_fname = self.latex_output_fname
		
		phenotype_id2analysis_method_id2gwr = {}
		
		param_data = grand_dataStructure
		param_data.call_method_id = self.call_method_id
		param_data.results_directory = self.results_directory
		param_data.output_dir = self.output_dir
		param_data.which_LD_statistic = which_LD_statistic
		param_data.min_distance = self.min_distance
		param_data.list_type_id = self.list_type_id
		param_data.no_of_snps_drawn = 0
		no_of_genes_handled = 0
		while 1:
			try:
				input_data = self.getSNPsFromInputFile(input_fname, snp_info)
				param_data.which_LD_statistic = which_LD_statistic
				
				for gene_id in input_data.gene_id_ls:
					sys.stderr.write("Handling Gene %s ...\n"%gene_id)
					snp2phenotype_id_ls = input_data.return_snp_data_given_gene_id(gene_id)
					gene_model = gene_annotation.gene_id2model.get(gene_id)
					if not gene_model:
						continue
					no_of_genes_handled += 1
					snp_ls = snp2phenotype_id_ls.keys()
					snp_ls.sort()
					if latex_fname:
						latex_fname_appx = 'gene_%s_%s'%(gene_id, getattr(gene_model, 'gene_symbol', ''))
						latex_fname_appx = latex_fname_appx.replace('/', '_')
						_latex_fname = '%s_%s'%(latex_fname, latex_fname_appx)
						latex_f = open(_latex_fname, 'w')
						matrix_of_gene_descriptions = []
						if len(gene_model.gene_commentaries)>0:
							for gene_commentary in gene_model.gene_commentaries:
								gene_desc_ls = self.returnGeneDescLs(self.gene_desc_names, gene_model, gene_commentary)
								matrix_of_gene_descriptions.append(gene_desc_ls)
						else:
							gene_desc_ls = self.returnGeneDescLs(self.gene_desc_names, gene_model)
							matrix_of_gene_descriptions.append(gene_desc_ls)
						gene_table_label = 'gene%s'%no_of_genes_handled
						snps_table_label = 'snps%s'%no_of_genes_handled
						caption = 'Gene %s (id=%s) has %s SNPs (check Tabel~\\ref{%s}) and %s phenotypes associated with.'%\
							(gene_model.gene_symbol, gene_id, len(snp2phenotype_id_ls), snps_table_label, input_data.get_no_of_phenotypes_given_gene_id(gene_id))
						latex_f.write(outputMatrixInLatexTable(matrix_of_gene_descriptions, caption, gene_table_label, self.gene_desc_names))
						
						matrix_of_snp_descriptions = input_data.return_matrix_of_snp_descriptions(gene_id)
						caption = '%s SNPs are related to Gene %s, id=%s, (check Tabel~\\ref{%s}).'%\
							(len(snp_ls), gene_model.gene_symbol, gene_id, gene_table_label)
						latex_f.write(outputMatrixInLatexTable(matrix_of_snp_descriptions, caption, snps_table_label, input_data.snp_desc_names))
					else:
						latex_f = None
					for snp in snp_ls:
						phenotype_id_ls = snp2phenotype_id_ls[snp]
						this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2])
						for phenotype_id in phenotype_id_ls:
							self.handleOneGeneOnePhenotype(phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
														input_data, param_data, latex_f)
					if latex_f:
						del latex_f
				
				snps_id2latex_f = {}
				if input_data.no_of_genes==0:	#no gene in input_fname, iterate over snps
					phenotype_id2snp_ls = input_data.return_phenotype_id2snp_ls()
					for phenotype_id, snp_ls in phenotype_id2snp_ls.iteritems():
						snp_ls.sort()
						for snp in snp_ls:
							
							this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2], stop=snp[3])
							latex_f = snps_id2latex_f.get(this_snp.snps_id)
							if this_snp.snps_id not in snps_id2latex_f and latex_fname:
								latex_f = open('%s_%s'%(latex_fname, this_snp.snps_id), 'w')
								snps_id2latex_f[this_snp.snps_id] = latex_f
							self.handleOneGeneOnePhenotype(phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
														input_data, param_data, latex_f)
					del phenotype_id2analysis_method_id2gwr	#clean up the memory
					phenotype_id2analysis_method_id2gwr = {}
					"""
					
					for snp, phenotype_id_ls in input_data.snp2phenotype_id_ls.iteritems():
						this_snp = SNPPassingData(chromosome=snp[0], position=snp[1], snps_id=snp[2])
						if latex_fname:
							latex_f = open('%s_%s'%(latex_fname, this_snp.snps_id), 'w')
						else:
							latex_f = None
						for phenotype_id in phenotype_id_ls:
							self.handleOneGeneOnePhenotype(phenotype_id2analysis_method_id2gwr, phenotype_id, this_snp, \
														input_data, param_data, latex_f)
					"""
				del snps_id2latex_f
			except:
				sys.stderr.write('Except: %s\n'%repr(sys.exc_info()))
				traceback.print_exc()
				#raise
			
			to_continue = raw_input("Continue running?([y]/n): ")
			if to_continue.upper()=='N':
				break
			_input_fname = raw_input("File containing phenotype_id, chromosome, position, default is %s: "%input_fname)
			_which_LD_statistic = raw_input("which LD statistic (1=r2, 2=|D_prime|), default is %s: "%(which_LD_statistic))
			_latex_fname = raw_input("Latex Filename (or prefix) to store gene description tables and figures, default is %s: "%(latex_fname))
			if _input_fname:
				input_fname = _input_fname
			if _which_LD_statistic:
				which_LD_statistic = _which_LD_statistic
			if _latex_fname:
				latex_fname = _latex_fname

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawSNPRegion
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
