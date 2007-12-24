#!/usr/bin/env python
"""
Usage: CmpAccession2Ecotype.py [OPTIONS] -p XXX -o OUTPUT_FNAME -j XXX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete'/'at.ecotype2accession'
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-l ...,	calls table, 'stock20071008.calls'(default)
	-o ...,	output_fname
	-j ...,	sub_justin_output_fname, (output the sub_data_matrix extracted from input_fname)
	-f ...,	latex_output_fname which stores both the summary and detailed difference tables
		specify this latex_output_fname if you want output
	-c,	comparison_only, output_fname and sub_justin_output_fname are already there
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/CmpAccession2Ecotype.py -p at.ecotype2accession -o 2010pcr_with_sequenom_149snps_ecotype2accession.csv -j sequenom_with_strains_matched_to_2010pcr_ecotype2accession.csv -r
	
	./src/CmpAccession2Ecotype.py -p at.accession2ecotype_complete -o 2010pcr_with_sequenom_149snps_accession2ecotype_complete.csv -j  sequenom_with_strains_matched_to_2010pcr_accession2ecotype_complete.csv -r -f 2010pcr_vs_sequenom_diff.tex

Description:
	program to compare the 149snp calls based on the common strains inn 2010 pcr data and Justin's sequenom data.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from variation.src.QualityControl import QualityControl

class Cmp250kVs2010(QualityControl):
	"""
	2007-12-19
		QC between 250k and 2010
	"""
	def __init__(self, curs, input_fname1, input_fname2, snp_locus_table1, snp_locus_table2, ecotype_duplicate2tg_ecotypeid_table, debug=0):
		self.curs = curs
		self.input_fname1 = input_fname1
		self.input_fname2 = input_fname2
		self.snp_locus_table1 = snp_locus_table1
		self.snp_locus_table2 = snp_locus_table2
		self.ecotype_duplicate2tg_ecotypeid_table = ecotype_duplicate2tg_ecotypeid_table
		
		self.debug = int(debug)
	
	def get_col_matching_dstruc(self, header1, header2):
		"""
		2007-12-19
		"""
		sys.stderr.write("Getting col matching dstruc ...\n")
		snpid2col_index1 = {}
		for i in range(2, len(header1)):
			snpid = header1[i]
			snpid2col_index1[snpid] = i-2
		snpid2col_index2 = {}
		for i in range(2, len(header2)):
			snpid = header2[i]
			snpid2col_index2[snpid] = i-2
		
		col_id12col_id2 = {}
		for snpid in snpid2col_index1:
			if snpid in snpid2col_index2:
				col_id12col_id2[snpid] = snpid
		sys.stderr.write("Done.\n")
		return snpid2col_index1, snpid2col_index2, col_id12col_id2
	
	def get_row_matching_dstruc(self, strain_acc_list1, category_list1, strain_acc_list2):
		"""
		2007-12-19
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			ecotypeid = int(strain_acc_list1[i])
			duplicate = int(category_list1[i])
			strain_acc = (ecotypeid, duplicate)
			strain_acc2row_index1[strain_acc] = i
		
		strain_acc2row_index2 = {}
		for i in range(len(strain_acc_list2)):
			strain_acc = strain_acc_list2[i]
			ecotypeid = int(strain_acc)
			strain_acc2row_index2[ecotypeid] = i
		
		row_id12row_id2 = {}
		for strain_acc in strain_acc2row_index1:
			ecotypeid, duplicate = strain_acc
			if ecotypeid in strain_acc2row_index2:
				row_id12row_id2[strain_acc] = ecotypeid
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2
	
	def load_dstruc(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header1, strain_acc_list1, category_list1, self.data_matrix1 = FilterStrainSNPMatrix_instance.read_data(self.input_fname1)
		header2, strain_acc_list2, category_list2, self.data_matrix2 = FilterStrainSNPMatrix_instance.read_data(self.input_fname2)
	 	
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(header1, header2)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(strain_acc_list1, category_list1, strain_acc_list2)


#if __name__ == '__main__':
hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

input_fname = './script/variation/genotyping/250ksnp/data/data_250k.tsv'

input_fname1 = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
input_fname2 = os.path.expanduser('~/script/variation/data/2010/data_2010_x_250k.tsv')
snp_locus_table1 = 'snps_250k'
snp_locus_table2 = 'snps'
ecotype_duplicate2tg_ecotypeid_table = 'ecotype_duplicate2tg_ecotypeid'
Cmp250kVs2010_ins= Cmp250kVs2010(curs, input_fname1, input_fname2, snp_locus_table1, snp_locus_table2, ecotype_duplicate2tg_ecotypeid_table)
Cmp250kVs2010_ins.load_dstruc()
import pylab
pylab.plot([1],[1])
pylab.show()

Cmp250kVs2010_ins.plot_row_NA_mismatch_rate('250k vs 2010 strain-wise')
#Cmp250kVs2010_ins.plot_col_NA_mismatch_rate('250k vs 2010 snp-wise')
