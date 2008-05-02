#!/usr/bin/env python
"""
2008-05-02
a ORM around database dbsnp
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import sqlalchemy, threading
from sqlalchemy.engine.url import URL
from sqlalchemy import Table, mapper, relation
from sqlalchemy.orm.session import Session
from pymodule.db import Database

class SNPset(object):
	pass

class SNPs(object):
	pass

class SNPs2SNPset(object):
	pass

class CallMethod(object):
	pass

class Calls(object):
	pass

class README(object):
	pass


class DBSNP(Database):
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
		table_ls = ['snpset', 'snps', 'snps2snpset', 'call_method', 'calls', 'readme']
		for table_name in table_ls:
			tables[table_name] = Table(table_name, metadata, autoload=True)
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		standalone_table_tuple_ls = [('snpset', SNPset), ('snps', SNPs), ('call_method', CallMethod), ('readme', README)]
		for table_name, table_class in standalone_table_tuple_ls:
			mappers[table_name] = mapper(table_class, tables[table_name])
		
		mappers['snps2snpset'] = mapper(SNPs2SNPset, tables['snps2snpset'],
										properties={'snps_obj': relation(SNPs), 'snpset_obj':relation(SNPset)})
		mappers['calls'] = mapper(Calls, tables['calls'], properties={'snps_obj': relation(SNPs), 'call_method_obj': relation(CallMethod)})
