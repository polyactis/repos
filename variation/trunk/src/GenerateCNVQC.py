#!/usr/bin/env python
"""

Examples:
	#generate CNV data from 2010 SNP data
	GenerateCNVQC.py -i data/2010/data_2010_ecotype_id_y0002_version3_no_offset_only_one_ref_col.tsv -o /tmp/2010_CNV_probe_qc
	
Description:
	Program to generate QC data from a SNP data for CNV data. 3 StrainXProbe files generated:
		mismatch matrix file. containing the no_of_mismatches for each probe. -2 is NA.
		insertion matrix file. containing the no_of_insertions for each probe. -2 is NA.
		deletion matrix file. containing the no_of_deletions for each probe. -2 is NA.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt
from sets import Set
from pymodule import ProcessOptions, SNPData, PassingData
from DB_250k2Array import DB_250k2Array

class GenerateCNVQC(object):
	__doc__ = __doc__
	"""
	2009-2-12
	"""
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('probes_table', 1, ): ['probes', 'e', 1, 'table where the CNV probes are'],\
							("input_fname", 0, ): [None, 'i', 1, 'StrainXSNP matrix file to be converted into CNV probe QC matrix'],\
							('output_fname_prefix',1,): ['', 'o', 1, 'filename prefix for mismatch_matrix and deletion_matrix'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-12
		"""
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def get_chr2CNV_probe_ls(self, curs, probes_table):
		probes, xy_ls, chr_pos_ls, probes_id_ls = DB_250k2Array.get_probes(curs, probes_table, snps=None, \
																		output_CNV_intensity=True)
		sys.stderr.write("Getting chr2CNV_probe_ls ...")
		chr2CNV_probe_ls = {}
		for chr_pos in chr_pos_ls:
			chromosome, position = chr_pos
			if chromosome not in chr2CNV_probe_ls:
				chr2CNV_probe_ls[chromosome] = []
			chr2CNV_probe_ls[chromosome].append(chr_pos)
		sys.stderr.write("Done.\n")
		return chr2CNV_probe_ls
	
	def findCNVprobe(self, chr2CNV_probe_ls, snp_id_tup):
		"""
		2009-2-12
			given a snp, find the corresponding CNV probe that covers it
			position for a CNV probe is the center position of the 25mer.
		"""
		probe_ls = chr2CNV_probe_ls[snp_id_tup[0]]
		left_probe_index = 0
		right_probe_index = len(probe_ls)-1
		snp_pos = snp_id_tup[1]
		while left_probe_index !=right_probe_index:
			if right_probe_index==(left_probe_index+1):	#two adjacent probes
				left_probe = probe_ls[left_probe_index]
				right_probe = probe_ls[right_probe_index]
				if left_probe[1]+12<snp_pos<right_probe[1]-12:	#falls in between, no CNV probe for it
					return None
				elif left_probe[1]-12<=snp_pos<=left_probe[1]+12:
					return left_probe
				elif right_probe[1]-12<=snp_pos<=right_probe[1]+12:
					return right_probe
			middle_index = (left_probe_index+right_probe_index)/2
			probe = probe_ls[middle_index]
			if probe[1]-12<=snp_pos<=probe[1]+12:
				return probe
			elif snp_pos<probe[1]-12:
				right_probe_index = middle_index
			elif snp_pos>probe[1]+12:
				left_probe_index = middle_index
		
		
	def get_probe_id2snp_id_ls(self, chr2CNV_probe_ls, snp_id_ls):
		sys.stderr.write("Getting probe_id2snp_id_ls ... ")
		probe_id2snp_id_ls = {}
		snp_id2tup = {}
		for snp_id in snp_id_ls:
			snp_id_tup = snp_id.split('_')
			snp_id_tup = map(int, snp_id_tup)
			snp_id2tup[snp_id] = snp_id_tup
			CNV_probe_id = self.findCNVprobe(chr2CNV_probe_ls, snp_id_tup)
			if CNV_probe_id is not None:
				if CNV_probe_id not in probe_id2snp_id_ls:
					probe_id2snp_id_ls[CNV_probe_id] = []
				probe_id2snp_id_ls[CNV_probe_id].append(snp_id)
		sys.stderr.write("Done.\n")
		return PassingData(probe_id2snp_id_ls=probe_id2snp_id_ls, snp_id2tup=snp_id2tup)
	
	def get_SNP2Col_allele(self, snpData):
		"""
		2009-2-12
		"""
		sys.stderr.write("Getting SNP2Col_allele ...")
		SNP2Col_allele = {}
		row_index = snpData.row_id2row_index["6909"]	#6909 is col-0, the reference genome
		for snp_id in snpData.col_id_ls:
			col_index = snpData.col_id2col_index[snp_id]
			SNP2Col_allele[snp_id] = snpData.data_matrix[row_index][col_index]
		sys.stderr.write("Done.\n")
		return SNP2Col_allele
	
	def getCNVQCMatrix(self, probe_id2snp_id_ls, snp_id2tup, snpData, SNP2Col_allele):
		"""
		2009-2-12
		"""
		sys.stderr.write("Getting CNV QC matricies ...")
		mismatch_matrix = numpy.zeros([len(snpData.row_id_ls), len(probe_id2snp_id_ls)], numpy.int)
		mismatch_matrix[:] = -2
		insertion_matrix = numpy.zeros(mismatch_matrix.shape, numpy.int)
		insertion_matrix[:] = -2
		deletion_matrix = numpy.zeros(mismatch_matrix.shape, numpy.int)
		deletion_matrix[:] = -2
		cnv_probe_ls = probe_id2snp_id_ls.keys()
		cnv_probe_ls.sort()
		cnv_probe2index = dict(zip(cnv_probe_ls, range(len(cnv_probe_ls))))
		
		for i in range(mismatch_matrix.shape[0]):
			for probe_id, snp_id_ls in probe_id2snp_id_ls.iteritems():
				col_index = cnv_probe2index[probe_id]
				no_of_mismatches = 0
				no_of_deletions = 0
				no_of_insertions = 0
				is_this_probe_NA = 1
				for snp_id in snp_id_ls:
					snp_id_tup = snp_id2tup[snp_id]
					snp_col_index = snpData.col_id2col_index[snp_id]
					allele = snpData.data_matrix[i][snp_col_index]
					col_allele = SNP2Col_allele[snp_id]
					if allele==-2 or allele==0:
						continue
					else:
						is_this_probe_NA = 0
						if snp_id_tup[2]!=0:	#the offset is not 0
							if allele!=-1:	#if it's deleted, then it's nothing
								no_of_insertions += 1								
						elif allele==-1:
							no_of_deletions+=1
						elif col_allele==-2 or col_allele==0:
							sys.stderr.write("allele for this accession %s at snp %s is %s while reference allele is NA: %s.\n"%\
											(snpData.row_id_ls[i], snp_id, allele, col_allele))
						elif allele!=col_allele:
							no_of_mismatches += 1
							
				if not is_this_probe_NA:
					mismatch_matrix[i][col_index] = no_of_mismatches
					insertion_matrix[i][col_index] = no_of_insertions
					deletion_matrix[i][col_index] = no_of_deletions
		mismatchData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=cnv_probe_ls, data_matrix=mismatch_matrix)
		insertionData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=cnv_probe_ls, data_matrix=insertion_matrix)
		deletionData = SNPData(row_id_ls=snpData.row_id_ls, col_id_ls=cnv_probe_ls, data_matrix=deletion_matrix)
		sys.stderr.write("Done.\n")
		return PassingData(mismatchData=mismatchData, insertionData=insertionData, deletionData=deletionData)
	
	def run(self):
		"""
		2009-2-12
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		curs = conn.cursor()
		
		chr2CNV_probe_ls = self.get_chr2CNV_probe_ls(curs, self.probes_table)
		snpData = SNPData(input_fname=self.input_fname, turn_into_array=1, ignore_2nd_column=1)
		
		probeData = self.get_probe_id2snp_id_ls(chr2CNV_probe_ls, snpData.col_id_ls)
		SNP2Col_allele = self.get_SNP2Col_allele(snpData)
		cnvQCData = self.getCNVQCMatrix(probeData.probe_id2snp_id_ls, probeData.snp_id2tup, snpData, SNP2Col_allele)
		cnvQCData.mismatchData.tofile('%s_mismatch.tsv'%self.output_fname_prefix)
		cnvQCData.insertionData.tofile('%s_insertion.tsv'%self.output_fname_prefix)
		cnvQCData.deletionData.tofile('%s_deletion.tsv'%self.output_fname_prefix)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = GenerateCNVQC
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()