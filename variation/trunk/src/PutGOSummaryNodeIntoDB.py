#!/usr/bin/env python
"""

Examples:
	#Summary node minimum size=150
	PutGOSummaryNodeIntoDB.py -m 150 -u yh -c
	
	#minimum size=150, exclude association evidenced by 'IEA'
	PutGOSummaryNodeIntoDB.py -m 150 -I -u yh -c

Description:
	2008-08-23 generate GO summary nodes and put them as candidate genes into db

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

from Stock_250kDB import Stock_250kDB, GeneList, GeneListType, GeneListSuperType
from annot.bin.GO.go_informative_node import go_informative_node_bfs
from sets import Set

class PutGOSummaryNodeIntoDB(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("min_size", 1, int): [100, '', 1, 'minimum GO node size'],\
							('exclude_IEA', 0): [0, 'I', 0, 'Exclude gene2go entries whose evidence is IEA'],\
							("super_type_name", 0, ): [None, '', 1, 'short name to summarize all nodes. if not given, generate on its own.'],\
							('tax_id', 0, int): [3702, 'x', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-07-10
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def saveSummaryNodes(self, session, gin_bfs_ins, super_type):
		sys.stderr.write("Saving summary nodes ... ")
		no_of_go_nodes = 0
		gene_id_in_informative_nodes_set = Set()
		for go_id,value in gin_bfs_ins.informative_node_dict.iteritems():
			if value == 1:
				no_of_go_nodes += 1
				list_type = GeneListType(short_name=gin_bfs_ins.go_id2go_name[go_id], description=go_id, type=super_type)
				session.save(list_type)
				gene_id_list = gin_bfs_ins.go_id_descendent2gene_id_dict[go_id].items()
				for gene_id in gene_id_list:
					gl = GeneList(gene_id=gene_id, list_type=list_type)
					gene_id_in_informative_nodes_set.add(gene_id)
					session.save(gl)
		sys.stderr.write(" %s GO nodes retained. covering %s genes. Done.\n"%(no_of_go_nodes, len(gene_id_in_informative_nodes_set)))
	
	def run(self):
		"""
		2008-08-19
		"""
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		
		session = db.session
		session.begin()
		
		if not self.super_type_name:
			self.super_type_name = 'GOSumNodeMinSize%s'%self.min_size
			if self.exclude_IEA:
				self.super_type_name += 'NoIEA'
		super_type = GeneListSuperType(short_name=self.super_type_name)
		session.save(super_type)
				
		gin_bfs_ins = go_informative_node_bfs(drivername=self.drivername, hostname=self.hostname, dbname=self.dbname,\
															schema=self.schema, db_user=self.db_user, db_passwd=self.db_passwd, \
															size=self.min_size, exclude_IEA=self.exclude_IEA, node_type=2, tax_id=self.tax_id)
		gin_bfs_ins.run()
		self.saveSummaryNodes(session, gin_bfs_ins, super_type)
		
		if self.commit:
			session.commit()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutGOSummaryNodeIntoDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()