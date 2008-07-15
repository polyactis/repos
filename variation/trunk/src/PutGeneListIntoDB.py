#!/usr/bin/env python
"""

Examples:
	PutGeneListIntoDB.py -i /Network/Data/250k/candidate_genes_GWA/anthocyanin.csv -l 4 -u yh -c

Description:
	2008-07-03 put gene list into db, stock_250k
	
	If list_type_name is given, list_type_id would be ignored although it's still required option (put a random number there).
	
	Program will avoid redundant genes in the input files. 

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

from Stock_250kDB import Stock_250kDB, GeneList, GeneListType

class PutGeneListIntoDb(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_fname", 1, ): [None, 'i', 1, ''],\
							("delimiter", 1, ): [',', '', 1, ''],\
							("list_type_id", 1, int): [None, 'l', 1, 'Gene list type. must be in table gene_list_type beforehand.'],\
							("list_type_name", 0, ): [None, '', 1, 'Gene list type short name. if given, list_type_id would be ignored.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
		"""
		2008-07-03 put gene list into db, stock_250k
		"""
		
	def putGeneListIntoDb(self, input_fname, list_type_id, list_type_name, mysql_curs, db):
		"""
		2008-07-15
			if the list_type_name is given, forget about list_type_id. program will first search db for the given list_type_name, if search failed, create a new entry.
		2008-07-15
			use gene_id2original_name to avoid redundancy in gene list
		"""
		import csv, sys, os
		from pymodule import get_gene_symbol2gene_id_set
		gene_symbol2gene_id_set = get_gene_symbol2gene_id_set(mysql_curs, 3702, table='genome.gene_symbol2id', upper_case_gene_symbol=1)	#3702 is At's tax id
		sys.stderr.write("%s entries in gene_symbol2gene_id.\n"%len(gene_symbol2gene_id_set))
		
		session = db.session
		
		reader = csv.reader(open(input_fname), delimiter=',')
		reader.next()	#skips the 1st line
		counter = 0
		success_counter = 0
		import re
		p_acc_ver = re.compile(r'(\w+)\.(\d+)')
		
		if list_type_name:	#if the short name is given, forget about list_type_id
			glt = GeneListType.query.filter_by(short_name=list_type_name).first()	#try search the db first.
			if not glt:
				glt = GeneListType(short_name=list_type_name)
				session.save(glt)
				session.flush()
		
		gene_id2original_name = {}	#to avoid redundancy in gene list
		for row in reader:
			original_name = row[0]
			gene_symbol = row[0].upper()
			if p_acc_ver.search(gene_symbol):
				gene_symbol, version = p_acc_ver.search(gene_symbol).groups()
			gene_id_set = gene_symbol2gene_id_set.get(gene_symbol)
			if gene_id_set==None:
				sys.stderr.write("Linking to gene id failed for %s. No such gene_symbol, %s, in gene_symbol2gene_id_set.\n"%(original_name, gene_symbol))
			elif len(gene_id_set)==1:
				gene_id = list(gene_id_set)[0]
				if gene_id not in gene_id2original_name:
					gene_id2original_name[gene_id] = original_name
				success_counter += 1
			elif len(gene_id_set)>1:
				sys.stderr.write("Too many gene_ids: %s, %s.\n"%(gene_symbol, gene_id_set))
			elif len(gene_id_set)==0:
				sys.stderr.write("Linking to gene id failed for %s. There is gene_symbol, %s, in gene_symbol2gene_id_set but it's empty.\n"%(original_name, gene_symbol))
			else:
				sys.stderr.write("not supposed to happen: original_name=%s, gene_symbol=%s, gene_id_set=%s\n."%(original_name, gene_symbol, gene_id_set))
			counter += 1
		del reader
		
		for gene_id, original_name in gene_id2original_name.iteritems():
			if list_type_name:
				gl = GeneList(gene_id=gene_id, list_type=glt, original_name=original_name)
			else:
				gl = GeneList(gene_id=gene_id, list_type_id=list_type_id, original_name=original_name)
			session.save(gl)
		sys.stderr.write("%s/%s linked successfully.\n"%(success_counter, counter))
		

	def run(self):
		"""
		"""
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		
		import MySQLdb
		mysql_conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		mysql_curs = mysql_conn.cursor()
		
		session = db.session
		session.begin()
		self.putGeneListIntoDb(self.input_fname, self.list_type_id, self.list_type_name, mysql_curs, db)
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutGeneListIntoDb
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()