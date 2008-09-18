#!/usr/bin/env python
"""

Examples:
	Cmp250kVs2010.py -i ~/script/variation/genotyping/250ksnp/data/data_250k.tsv -j ~/script/variation/data/2010/data_2010_x_250k.tsv
	
Description:
	This is a generic program (not just 250k vs 2010) to do QC. Row is matched by 1st column from input_fname1 and 1st column from input_fname2.
	Column is matched by the header.
	
	It is like QC.py but only handles Strain X SNP format. It will pop up a GUI window to allow better viewing of QC.
	For better usage, run QCVisualizeII.py and choose this class and run.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import read_data, QualityControl
#from variation.src.QualityControl import QualityControl

class Cmp250kVs2010(QualityControl):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['stock', 'd', 1, '',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 0, ):['yh', 'u', 1, 'database username',],\
							('db_passwd', 0, ):['', 'p', 1, 'database password', ],\
							('curs', 0, ):[None, '', 1, 'database cursor, for other program to call. no need to provide it in commandline.', ],\
							('input_fname1',1, ): [None, 'i', 1, 'file of SNPData1'],\
							('input_fname2',1, ): [None, 'j', 1, 'file of SNPData2'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
							#('input_fname1_format',1,int): [1, 'k', 1, 'Format of input_fname1. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId. 3=snp X strain with arrayId.'],\
							#('input_fname2_format',1,int): [1, 'l', 1, 'Format of input_fname2. 1=strain X snp (Yu). 2=snp X strain (Bjarni) without arrayId'],\
	"""
	2008-09-18
		become a generic program to do QC, like QC.py but only handles Strain X SNP format.
	2007-12-19
		QC between 250k and 2010
	"""
	def __init__(self, **keywords):
		"""
		2008-09-18
			use pymodule.ProcessOptions
		2008-01-21
			ecotype_duplicate2tg_ecotypeid_table, snp_locus_table1 and snp_locus_table2 are useless
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		QualityControl.__init__(self, debug=self.debug)
	
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
		2008-02-12
			2010 accession_id=75/ecotype_id=8315 (Kas-2) should be mapped to 250k/149SNP ecotype_id=8424(Kas-1)
			so add it in strain_acc2row_index2
		2007-12-19
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			ecotypeid = int(strain_acc_list1[i])
			duplicate = category_list1[i]
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
			elif ecotypeid==8424:	#2010 accession_id=75/ecotype_id=8315 (Kas-2) should be mapped to 250k/149SNP ecotype_id=8424(Kas-1)
				if 8315 in strain_acc2row_index2:
					row_id12row_id2[strain_acc] = 8315
			else:
				print 'Failure:', strain_acc
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2
	
	def load_dstruc(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		QualityControl.load_dstruc(self)
		self.header1, self.strain_acc_list1, self.category_list1, self.data_matrix1 = read_data(self.input_fname1)
		self.header2, self.strain_acc_list2, self.category_list2, self.data_matrix2 = read_data(self.input_fname2)
	 	
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.header1, self.header2)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.strain_acc_list1, self.category_list1, self.strain_acc_list2)
	
	def run(self):
		"""
		2008-09-18
			wrap up to run standalone.
			copied from old lines below __name__=='__main__'
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
		if not self.curs:
			import MySQLdb
			conn = MySQLdb.connect(db=self.database, host=self.hostname, user=self.db_user, passwd=self.db_passwd)
			self.curs = conn.cursor()
		
		self.load_dstruc()
		self.plot_row_NA_mismatch_rate('%s vs %s strain-wise'%(os.path.basename(self.input_fname1), os.path.basename(self.input_fname2)))
		
if __name__ == '__main__':
	from __init__ import ProcessOptions
	main_class = Cmp250kVs2010
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()