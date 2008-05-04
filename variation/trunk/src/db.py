#!/usr/bin/env python
"""
2008-04-28
A wrapper on top of sqlalchemy around a database. Mostly copied from collective.lead.Database. Can't directly use it because
of trouble in understanding how to use adapter involved in TreadlocalDatabaseTransactions.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import sqlalchemy, threading
from sqlalchemy.engine.url import URL
from sqlalchemy import Table, mapper, relation
from sqlalchemy.orm.session import Session

from pymodule.db import Database, TableClass

class Results(TableClass):
	pass

class ResultsMethod(TableClass):
	pass

class PhenotypeAvg(TableClass):
	pass

class PhenotypeMethod(TableClass):
	pass

class QCMethod(TableClass):
	pass

class Call_QC(TableClass):
	pass

class Call_Info(TableClass):
	pass

class Call_Method(TableClass):
	pass

class Array_Info(TableClass):
	pass


class Stock_250kDatabase(Database):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['stock_250k', 'd', 1, '',],\
							('username', 1, ):[None, 'u', 1, 'database username',],\
							('password',1, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		
		self._threadlocal = threading.local()
		self.tables = {}
		self.mappers = {}
		self._engine = None
		
	@property
	def _url(self):
		return URL(drivername=self.drivername, username=self.username,
				   password=self.password, host=self.hostname,
				   port=self.port, database=self.database)
	
	def _setup_tables(self, metadata, tables):
		"""Map the database structure to SQLAlchemy Table objects
		"""
			
		tables['phenotype_avg'] = Table('phenotype_avg', metadata, autoload=True)
		tables['phenotype_method'] = Table('phenotype_method', metadata, autoload=True)
		tables['QC_method'] = Table('QC_method', metadata, autoload=True)
		tables['results'] = Table('results', metadata, autoload=True)
		tables['results_method'] = Table('results_method', metadata, autoload=True)
		tables['call_QC'] = Table('call_QC', metadata, autoload=True)
		tables['call_info'] = Table('call_info', metadata, autoload=True)
		tables['call_method'] = Table('call_method', metadata, autoload=True)
		tables['array_info'] = Table('array_info', metadata, autoload=True)
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		mappers['phenotype_method'] = mapper(PhenotypeMethod, tables['phenotype_method'])
		mappers['phenotype_avg'] = mapper(PhenotypeAvg, tables['phenotype_avg'],
										properties={'phenotype_method_obj': relation(PhenotypeMethod),})
		mappers['QC_method'] = mapper(QCMethod, tables['QC_method'])
		mappers['results'] = mapper(Results, tables['results'], properties={'results_method_obj': relation(ResultsMethod), 'phenotype_method_obj': relation(PhenotypeMethod)})
		mappers['results_method'] = mapper(ResultsMethod, tables['results_method'])
		mappers['call_QC'] = mapper(Call_QC, tables['call_QC'], properties={'call_info_obj': relation(Call_Info),\
																		'qc_method_obj':relation(QCMethod)})
		mappers['call_info'] = mapper(Call_Info, tables['call_info'], properties={'array_info_obj': relation(Array_Info), 'call_method_obj': relation(Call_Method)})
		mappers['call_method'] = mapper(Call_Method, tables['call_method'])
		mappers['array_info'] = mapper(Array_Info, tables['array_info'])


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Stock_250kDatabase
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	
	"""
	
	#from /usr/lib/python2.5/site-packages/zope/interface/adapter.txt
	from zope.interface.adapter import AdapterRegistry
	registry = AdapterRegistry()
	registry.register([Database], ITransactionAware, '', 12)
	
	from zope.component import registry
	from zope.component import tests
	components = registry.Components('comps')
	
	#from /usr/lib/zope2.10/lib/python/zope/, interface/ component/registry.txt
	print components.registerAdapter(TreadlocalDatabaseTransactions)
	print components.getAdapter(ITransactionAware, Database)
	"""
	if instance.debug:
		import pdb
		pdb.set_trace()
	session = instance.session
	print dir(session.query(Results))
	
	for row in session.query(ResultsMethod).list():
		print row.id
		print row.short_name
	
	i = 0
	while i <10:
		row = session.query(Results).offset(i).limit(1).list()	#all() = list() returns a list of objects. first() returns the 1st object. one() woud raise error because 'Multiple rows returned for one()'
		print len(row)
		row = row[0]
		i += 1
		print row.id
		print row.chr
		print row.start_pos
		print row.score
		print row.method_id
		print row.results_method_obj.short_name
		print row.phenotype_method_id
		print row.phenotype_method_obj.short_name
