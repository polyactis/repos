#!/usr/bin/env python
"""

Examples:
		
	PutFTPathwayGenesIntoDB.py -i ~/doc/retreat_11_14_2008/FT_pathway_RouxEtAl2006.csv -c
Description:
	2008-11-14 
	program to put flowering time pathway genes (from RouxEtAl2006) into database.

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
from sets import Set
from Stock_250kDB import Stock_250kDB, FTGene

class PutFTPathwayGenesIntoDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("input_fname", 1, ): [None, 'i', 1, ''],\
							("delimiter", 1, ): [',', '', 1, ''],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self,  **keywords):
		"""
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		2008-08-19
		"""
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		session.begin()
		
		import MySQLdb
		mysql_conn = MySQLdb.connect(db=self.dbname, host='banyan.usc.edu', user = self.db_user, passwd = self.db_passwd)
		mysql_curs = mysql_conn.cursor()
		from pymodule import get_gene_symbol2gene_id_set
		gene_symbol2gene_id_set = get_gene_symbol2gene_id_set(mysql_curs, 3702, table='genome.gene_symbol2id', upper_case_gene_symbol=1)	#3702 is At's tax id
		
		
		reader = csv.reader(open(self.input_fname), delimiter=self.delimiter)
		reader.next()
		counter = 0
		success_counter = 0
		gene_id_pathway_id_set = Set()
		for row in reader:
			gene_symbol = row[0].upper()
			pathway_id = int(row[1])
			gene_id_set = gene_symbol2gene_id_set.get(gene_symbol)
			if gene_id_set==None:
				sys.stderr.write("Linking to gene id failed for %s. No such gene_symbol, %s, in gene_symbol2gene_id_set.\n"%(gene_symbol,gene_symbol))
			elif len(gene_id_set)==1:
				gene_id = list(gene_id_set)[0]
				gene_id_pathway_id_set.add((gene_id, pathway_id))
				success_counter += 1
			elif len(gene_id_set)>1:
				sys.stderr.write("Too many gene_ids: %s, %s.\n"%(gene_symbol, gene_id_set))
			elif len(gene_id_set)==0:
				sys.stderr.write("Linking to gene id failed for %s. There is gene_symbol, %s, in gene_symbol2gene_id_set but it's empty.\n"%(gene_symbol, gene_symbol))
			else:
				sys.stderr.write("not supposed to happen: original_name=%s, gene_symbol=%s, gene_id_set=%s\n."%(gene_symbol, gene_symbol, gene_id_set))
			counter += 1
		
		for gene_id, pathway_id in gene_id_pathway_id_set:
			rows = FTGene.query.filter_by(gene_id=gene_id).filter_by(pathway_id=pathway_id)
			if rows.count()>0:
				row = rows.first()
				sys.stderr.write('gene_id %s, pathway_id %s already in db with id=%s.\n'%(gene_id, pathway_id, row.id))
				continue
			ft_gene = FTGene(gene_id=gene_id, pathway_id=pathway_id)
			session.save_or_update(ft_gene)
			session.flush()
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutFTPathwayGenesIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
