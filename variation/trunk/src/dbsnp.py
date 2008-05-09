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
from pymodule.db import Database, TableClass

class SNPset(TableClass):
	pass

class SNPs(TableClass):
	pass

class SNPs2SNPset(TableClass):
	pass

class CallMethod(TableClass):
	pass

class Calls(TableClass):
	pass

class README(TableClass):
	pass

class Accession(TableClass):
	pass

class AtEcotype2Accession(TableClass):
	pass

class AtGenotype(TableClass):
	pass

class AtAllele(TableClass):
	pass

class AtLocus(TableClass):
	pass

class SNPsABAlleleMapping(TableClass):
	pass

class DBSNP(Database):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['dbsnp', 'd', 1, '',],\
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
	
	def _setup_tables(self, metadata, tables):
		"""Map the database structure to SQLAlchemy Table objects
		"""
		table_ls = ['snpset', 'snps', 'snps2snpset', 'call_method', 'calls', 'readme', 'accession', 'snps_ab_allele_mapping']
		for table_name in table_ls:
			tables[table_name] = Table(table_name, metadata, autoload=True)
		
		"""
		schema2table_ls = {'at':['ecotype2accession', 'genotype', 'allele', 'locus']}
		for schema, table_ls in schema2table_ls.iteritems():
			for table_name in table_ls:
				tables[table_name] = Table(table_name, metadata, autoload=True, schema=schema)
		"""
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		standalone_table_tuple_ls = [('snps', SNPs), ('call_method', CallMethod), ('readme', README), ('accession', Accession)]
		for table_name, table_class in standalone_table_tuple_ls:
			mappers[table_name] = mapper(table_class, tables[table_name])
		
		mappers['snpset'] = mapper(SNPset, tables['snpset'],
										properties={'snps': relation(SNPs, secondary=tables['snps2snpset'], backref='snpset')})
		mappers['snps2snpset'] = mapper(SNPs2SNPset, tables['snps2snpset'],
										properties={'snps': relation(SNPs), 'snpset':relation(SNPset)})
		mappers['calls'] = mapper(Calls, tables['calls'], properties={'snps': relation(SNPs), 'call_method': relation(CallMethod), 'accession': relation(Accession)})
		mappers['snps_ab_allele_mapping'] = mapper(SNPsABAlleleMapping, tables['snps_ab_allele_mapping'],
										properties={'snps': relation(SNPs, backref='snps_ab_allele_mapping'),\
											'readme':relation(README)})
