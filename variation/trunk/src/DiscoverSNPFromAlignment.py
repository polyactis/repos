#!/usr/bin/env python
"""
Examples:
	DiscoverSNPFromAlignment.py -a 1843 -o /tmp/alignment_1843_matrix.tsv
	
Description:
	2008-8-1 program to discover polymorphic positions (SNPs) from 2010 raw alignment data.
	The definition of polymorphic position is any position having >=2 outcomes, including deletion, excluding missing.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import PassingData, dict_map, importNumericArray, write_data_matrix
num = importNumericArray()
from variation.src.AtDB import AtDB, Sequence, Alignment
from variation.src.common import nt2number
from sets import Set

class DiscoverSNPFromAlignment(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['at', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("alignment_id", 1, ): [None, 'a', 1, 'id in table at.alignment'],\
							("output_fname", 1, ): [None, 'o', 1, 'to store the SNP matrix'],\
							("strain_id_type", 1, int): [1, 'y', 1, 'which type of id used to identify strain. 1(ecotype_id), 2(accession_id)'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-8-1
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getAlignmentMatrix(self, alignment_id):
		sys.stderr.write("Getting alignment matrix for alignment=%s ..."%(alignment_id))
		snp_pos_ls = []
		accession_id_ls = []
		name_ls = []
		data_matrix = []
		rows = Sequence.query.filter_by(alignment=alignment_id).order_by(Sequence.accession).all()
		counter = 0
		for row in rows:
			if counter == 0:
				for i in range(len(row.alignment_obj.target)):
					base_number = nt2number[row.alignment_obj.target[i]]
					if base_number!=-1:
						if i==0:
							snp_pos_ls.append((row.alignment_obj.chromosome, row.alignment_obj.start, 0))	#the 3rd position is insertion offset relative to Column position
						else:
							snp_pos_ls.append((row.alignment_obj.chromosome, snp_pos_ls[i-1][1]+1, 0))
					else:	#base is deletion
						if i==0:
							snp_pos_ls.append((row.alignment_obj.chromosome, row.alignment_obj.start-1, 1))	#this probably doesn't exist in db. it's controversal whether this insertion should be assigned to the previous or alignment's start base
						else:
							snp_pos_ls.append((row.alignment_obj.chromosome, snp_pos_ls[i-1][1], snp_pos_ls[i-1][2]+1))	#position doesn't change. offset++
			accession_id_ls.append(row.accession)
			name_ls.append(row.accession_obj.name)
			data_row = dict_map(nt2number, row.bases)
			data_matrix.append(data_row)
			counter += 1
		data_matrix = num.array(data_matrix, num.int8)
		passingdata = PassingData(snp_pos_ls=snp_pos_ls, accession_id_ls=accession_id_ls, name_ls=name_ls, data_matrix=data_matrix)
		sys.stderr.write(' %s accessions, %s bases. Done.\n'%(len(accession_id_ls), len(snp_pos_ls)))
		return passingdata
	
	def pickPolymorphicColumns(self, passingdata):
		sys.stderr.write("Picking polymorphic SNPs ...")
		no_of_rows, no_of_cols = passingdata.data_matrix.shape
		polymorphic_column_index_ls = []
		for j in range(no_of_cols):
			allele_set = Set(passingdata.data_matrix[:,j])
			if 0 in allele_set:	#remove NA if it's there
				allele_set.remove(0)
			if len(allele_set)>1:	#polymorphic
				polymorphic_column_index_ls.append(j)
		
		new_data_matrix = num.zeros([no_of_rows, len(polymorphic_column_index_ls)], num.int8)
		new_snp_pos_ls = []
		for i in range(len(polymorphic_column_index_ls)):
			index = polymorphic_column_index_ls[i]
			new_snp_pos_ls.append(passingdata.snp_pos_ls[index])
			new_data_matrix[:,i] = passingdata.data_matrix[:,index]
		del passingdata.data_matrix
		passingdata.data_matrix = new_data_matrix
		passingdata.snp_pos_ls = new_snp_pos_ls
		sys.stderr.write(' %s SNPs chosen. Done.\n'%(len(passingdata.snp_pos_ls)))
	
	def run(self):
		db = AtDB(drivername=self.drivername, username=self.db_user,
				password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		
		passingdata = self.getAlignmentMatrix(self.alignment_id)
		self.pickPolymorphicColumns(passingdata)
		
		header = ['id', 'name']
		for snp_pos in passingdata.snp_pos_ls:
			header.append('%s_%s_%s'%snp_pos)
		
		if self.strain_id_type==1:
			ecotype_id_ls = []
			for accession_id in passingdata.accession_id_ls:
				rows = db.metadata.bind.execute("select * from %s where accession_id=%s"%('accession2tg_ecotypeid', accession_id))
				row = rows.fetchone()
				ecotype_id_ls.append(row.ecotype_id)
			strain_acc_list = ecotype_id_ls
		elif self.strain_id_type==2:
			strain_acc_list = passingdata.accession_id_ls
		else:
			sys.stderr.write("strain_id_type %s not supported.\n"%(self.strain_id_type))
			sys.exit(2)
		write_data_matrix(passingdata.data_matrix, self.output_fname, header, \
						strain_acc_list, passingdata.name_ls)


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DiscoverSNPFromAlignment
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()