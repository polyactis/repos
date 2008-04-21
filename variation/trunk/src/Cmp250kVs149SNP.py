#!/usr/bin/env python
"""
Usage: Cmp250kVs149SNP.py [OPTIONS] -p XXX -o OUTPUT_FNAME -j XXX

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
	./src/Cmp250kVs149SNP.py 
	

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

class Cmp250kVs149SNP(QualityControl):
	"""
	2007-12-14
		QC between 250k and 149snp
	"""
	def __init__(self, curs='', input_fname1='', input_fname2='', snp_locus_table_250k='snps_250k', snp_locus_table_149snp='snps', ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid', debug=0):
		"""
		2008-01-21
			add default values to all parameters.
		"""
		self.curs = curs
		self.input_fname1 = input_fname1
		self.input_fname2 = input_fname2
		self.snp_locus_table_250k = snp_locus_table_250k
		self.snp_locus_table_149snp = snp_locus_table_149snp
		self.ecotype_duplicate2tg_ecotypeid_table = ecotype_duplicate2tg_ecotypeid_table

		self.debug = int(debug)
	
	def get_col_matching_dstruc(cls, header_250k, header_149, curs, snp_locus_table_250k, snp_locus_table_149snp):
		"""
		2008-04-20
			change the one snpid by one sql matching to a full sql matching first, then link among the sql results
		2007-12-18
		2007-12-21
			snpid_250k2snpid_149 replaces col_index_250k2col_index_149
		"""
		sys.stderr.write("Getting col matching dstruc ...\n")
		snpid2col_index_250k = {}
		for i in range(2, len(header_250k)):
			snpid = header_250k[i]
			snpid2col_index_250k[snpid] = i-2
		
		snpid2col_index_149 = {}
		for i in range(2, len(header_149)):
			snpid = header_149[i]
			snpid2col_index_149[snpid] = i-2
		
		snpid_250k2snpid_149 = {}
		curs.execute("select s1.snpid, s2.snpid from %s s1, %s s2 where s1.chromosome=s2.chromosome and s1.position=s2.position"%\
					(snp_locus_table_250k, snp_locus_table_149snp))
		rows = curs.fetchall()
		for row in rows:
			snpid_250k, snpid_149 = row
			if snpid_250k in snpid2col_index_250k and snpid_149 in snpid2col_index_149:
				snpid_250k2snpid_149[snpid_250k] = snpid_149
		sys.stderr.write("Done.\n")
		return snpid2col_index_250k, snpid2col_index_149, snpid_250k2snpid_149
	
	get_col_matching_dstruc = classmethod(get_col_matching_dstruc)
	def get_row_matching_dstruc(self, strain_acc_list_250k, category_list_250k, strain_acc_list_149, curs, ecotype_duplicate2tg_ecotypeid_table):
		"""
		2007-12-18
		2007-12-21
			strain_acc_250k2strain_acc_149 replaces row_index_250k2row_index_149
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index_250k = {}
		for i in range(len(strain_acc_list_250k)):
			ecotypeid = int(strain_acc_list_250k[i])
			duplicate = int(category_list_250k[i])
			strain_acc = (ecotypeid, duplicate)
			strain_acc2row_index_250k[strain_acc] = i
		
		strain_acc2row_index_149 = {}
		for i in range(len(strain_acc_list_149)):
			strain_acc = strain_acc_list_149[i]
			ecotypeid = int(strain_acc)
			strain_acc2row_index_149[ecotypeid] = i
		
		strain_acc_250k2strain_acc_149 = {}
		for strain_acc in strain_acc2row_index_250k:
			ecotypeid, duplicate = strain_acc
			curs.execute("select tg_ecotypeid from %s where ecotypeid=%s"%(ecotype_duplicate2tg_ecotypeid_table, ecotypeid))
			rows = curs.fetchall()
			if rows:
				tg_ecotypeid = rows[0][0]
				if tg_ecotypeid in strain_acc2row_index_149:	#could be no match
					strain_acc_250k2strain_acc_149[strain_acc] = tg_ecotypeid
				else:
					print 'Failure:', strain_acc
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index_250k, strain_acc2row_index_149, strain_acc_250k2strain_acc_149
	
	def load_dstruc(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		QualityControl.load_dstruc(self)
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		self.header1, self.strain_acc_list1, self.category_list1, self.data_matrix1 = FilterStrainSNPMatrix_instance.read_data(self.input_fname1)
		self.header2, self.strain_acc_list2, self.category_list2, self.data_matrix2 = FilterStrainSNPMatrix_instance.read_data(self.input_fname2)
	 	
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.header1, self.header2, self.curs, self.snp_locus_table_250k, self.snp_locus_table_149snp)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.strain_acc_list1, self.category_list1, self.strain_acc_list2, self.curs, self.ecotype_duplicate2tg_ecotypeid_table)

if __name__ == '__main__':
	hostname='localhost'
	dbname='stock20071008'
	import MySQLdb
	conn = MySQLdb.connect(db=dbname,host=hostname)
	curs = conn.cursor()
	
	input_fname = './script/variation/genotyping/250ksnp/data/data_250k.tsv'
	
	input_fname_250k = './script/variation/genotyping/250ksnp/data/data_250k_x_149SNP.tsv'
	input_fname_149snp = './script/variation/stock20071008/data.tsv'
	snp_locus_table_250k = 'snps_250k'
	snp_locus_table_149snp = 'snps'
	ecotype_duplicate2tg_ecotypeid_table = 'ecotype_duplicate2tg_ecotypeid'
	Cmp250kVs149SNP_ins = Cmp250kVs149SNP(curs, input_fname_250k, input_fname_149snp, snp_locus_table_250k, snp_locus_table_149snp, ecotype_duplicate2tg_ecotypeid_table)
	Cmp250kVs149SNP_ins.load_dstruc()
	
	Cmp250kVs149SNP_ins.plot_row_NA_mismatch_rate('250k vs 149SNP strain-wise')
	#Cmp250kVs149SNP_ins.plot_col_NA_mismatch_rate('250k vs 149SNP snp-wise')
	
