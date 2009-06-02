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
from sqlalchemy import Table
from sqlalchemy.orm.session import Session
from pymodule.db import Database, TableClass
from sqlalchemy.orm import mapper, relation

from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker

from datetime import datetime
from pymodule.db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, transactional=False))
__metadata__ = MetaData()

class README(Entity):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SNPset(Entity):
	name = Field(String(200), nullable=False, unique=True)
	description = Field(String(4000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snpset',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SNPs(Entity):
	name = Field(String(200), nullable=False, unique=True)
	chromosome = Field(Integer)
	position = Field(Integer)
	offset = Field(Integer)
	probe_sequence = Field(String(4000))
	allele1 = Field(String(2))
	allele2 = Field(String(2))
	allele3 = Field(String(2))
	allele4 = Field(String(2))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snps',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class SNPs2SNPset(Entity):
	snp = ManyToOne('SNPs', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	snpset = ManyToOne('SNPset', colname='snpset_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='snps2snpset',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CallMethod(Entity):
	"""
	2009-5-19
		add ManyToMany relationship towards Accession
	"""
	short_name = Field(String(100), nullable=False, unique=True)
	method_description = Field(String(8000))
	data_description = Field(String(8000))
	comment = Field(String(8000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	accession_ls = ManyToMany("Accession", tablename='accession2call_method', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='call_method',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Calls(Entity):
	accession = ManyToOne('Accession', colname='accession_id', ondelete='CASCADE', onupdate='CASCADE')
	snp = ManyToOne('SNPs', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	genotype = Field(String(2))
	call_method = ManyToOne('CallMethod', colname='call_method_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='calls',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Accession(Entity):
	"""
	2009-5-19
		add ManyToMany relationship towards CallMethod
	"""
	name = Field(String(100), nullable=False, unique=True)
	ecotype_id = Field(Integer)
	duplicate = Field(Integer)
	accession_2010_id = Field(Integer)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	call_method_ls = ManyToMany("CallMethod", tablename='accession2call_method', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='accession',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class AtEcotype2Accession(TableClass):
	pass

class AtGenotype(TableClass):
	pass

class AtAllele(TableClass):
	pass

class AtLocus(TableClass):
	pass

class SNPsABAlleleMapping(Entity):
	snp = ManyToOne('SNPs', colname='snps_id', ondelete='CASCADE', onupdate='CASCADE')
	allele_A_nt = Field(String(2))
	allele_B_nt = Field(String(2))
	tg_snps_name = Field(String(200))
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	relative_NA_rate = Field(Float)
	relative_no_of_NAs = Field(Integer)
	relative_no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='snps_ab_allele_mapping',  metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class DBSNP_sqlalchemy(Database):
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

class DBSNP(ElixirDB):
	__doc__ = __doc__
	option_default_dict = ElixirDB.option_default_dict.copy()
	option_default_dict[('drivername', 1,)][0] = 'mysql'
	option_default_dict[('database', 1,)][0] = 'dbsnp'
	def __init__(self, **keywords):
		"""
		2009-4-10
			simplified further by moving db-common lines to ElixirDB
		2008-08-11
		"""
		from pymodule.ProcessOptions import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
		
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DBSNP
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)