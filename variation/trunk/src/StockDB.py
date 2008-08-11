#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	StockDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	StockDB.py -u yh
	
Description:
	2008-08-11
	This is a wrapper for the stock database, build on top of elixir. Several tables need to be manually created because
	exlixir can't represent 'unsigned' integer very well.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text, Boolean, Float
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint

from datetime import datetime

from pymodule.db import ElixirDB
from pymodule.db import README

class BatchEcotype(Entity):
	batch = ManyToOne('Batch', colname='batchid', ondelete='CASCADE', onupdate='CASCADE')
	ecotype = ManyToOne('Ecotype', colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	lineorder = Field(Integer, nullable=False)
	using_options(tablename='batch_ecotype')
	using_table_options(mysql_engine='InnoDB')
	
class Batch(Entity):
	batchname = Field(String(30), nullable=False)
	using_options(tablename='batch')
	using_table_options(mysql_engine='InnoDB')

class Person(Entity):
	title = Field(String(4))
	surname = Field(String(40))
	firstname = Field(String(40))
	email = Field(String(100))
	username = Field(String(10))
	password = Field(String(15))
	usertype = Field(Integer)
	donor = Field(Integer)
	using_options(tablename='person')
	using_table_options(mysql_engine='InnoDB')

class Country(Entity):
	name = Field(String(100))
	abbr = Field(String(10))
	using_options(tablename='country')
	using_table_options(mysql_engine='InnoDB')

class Address(Entity):
	person = ManyToOne("Person", colname='personid', ondelete='CASCADE', onupdate='CASCADE')
	address1 = Field(String(100))
	address2 = Field(String(100))
	city = Field(String(100))
	stateprovince = Field(String(100))
	region = Field(String(100))
	zippostal = Field(String(20))
	country = ManyToOne("Country", colname='countryid', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='address')
	using_table_options(mysql_engine='InnoDB')
	
class Site(Entity):
	address = ManyToOne("Address", colname='addressid', ondelete='CASCADE', onupdate='CASCADE')
	name = Field(String(100))
	description = Field(String(500))
	using_options(tablename='site')
	using_table_options(mysql_engine='InnoDB')
	
class Organism(Entity):
	genus = Field(String(30))
	species = Field(String(30))
	using_options(tablename='organism')
	using_table_options(mysql_engine='InnoDB')
	
class SNPs(Entity):
	organism = ManyToOne("Organism", colname='organismid', ondelete='CASCADE', onupdate='CASCADE')
	snpid = Field(String(20), unique=True)
	dir = Field(String(1))
	chromosome = Field(Integer)
	position = Field(Integer)
	refcall = Field(String(1))
	using_options(tablename='snps')
	using_table_options(mysql_engine='InnoDB')
	
class SeqInfo(Entity):
	platename = Field(String(40))
	experiment = Field(String(40))
	chip = Field(String(40))
	experiment_type = Field(String(40))
	ms_name = Field(String(40))
	creation_date = Field(DateTime, default=datetime.now)
	using_options(tablename='seqinfo')
	using_table_options(mysql_engine='InnoDB')
	
class Strain(Entity):
	ecotype = ManyToOne("Ecotype", colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	extraction = ManyToOne("Extraction", colname='extractionid', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo = ManyToOne("SeqInfo", colname='seqinfoid', ondelete='CASCADE', onupdate='CASCADE')
	plateid = Field(String(25))
	wellid = Field(String(3))
	replicate = Field(Boolean)
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='strain')
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('ecotypeid', 'plateid'))
	
class Extraction(Entity):
	namelabel = Field(String(25))
	datesubmitted = Field(DateTime, default=datetime.now)
	notes = Field(String(255))
	using_options(tablename='extraction')
	using_table_options(mysql_engine='InnoDB')

class Trip(Entity):
	collectiondate = Field(DateTime)
	notes = Field(String(500))
	using_options(tablename='trip')
	using_table_options(mysql_engine='InnoDB')

class Cross(Entity):
	maternalid = Field(Integer)
	paternalid = Field(Integer)
	name = Field(String(40))
	barcode = Field(String(25))
	generation = Field(Integer)
	using_options(tablename='crosses')
	using_table_options(mysql_engine='InnoDB')

class Ecotype(Entity):
	donor = ManyToOne("Person", colname='donorid', ondelete='CASCADE', onupdate='CASCADE')
	collector = ManyToOne("Person", colname='collectorid', ondelete='CASCADE', onupdate='CASCADE')
	site = ManyToOne("Site", colname='siteid', ondelete='CASCADE', onupdate='CASCADE')
	trip = ManyToOne("Trip", colname='tripid', ondelete='CASCADE', onupdate='CASCADE')
	organism = ManyToOne("Organism", colname='organismid', ondelete='CASCADE', onupdate='CASCADE')
	cross = ManyToOne("Cross", colname='crossid', ondelete='CASCADE', onupdate='CASCADE')
	name = Field(String(50))
	alias = Field(String(50))
	description = Field(String(270))
	barcode = Field(String(14))
	stockparent = Field(String(10))
	nativename = Field(String(50))
	collectiondate = Field(DateTime)
	latitude = Field(Float)
	longitude = Field(Float)
	locationquality = Field(Integer)
	elevation = Field(String(5))
	dnastatus = Field(Integer)
	bulkstatus = Field(Integer)
	bulkdate = Field(DateTime)
	labderived = Field(Boolean)
	incompleteplex = Field(Integer)
	using_options(tablename='ecotype')
	using_table_options(mysql_engine='InnoDB')
	
class Calls(Entity):
	strain = ManyToOne("Strain", colname='strainid', ondelete='CASCADE', onupdate='CASCADE')
	snp = ManyToOne("SNPs", colname='snpid', ondelete='CASCADE', onupdate='CASCADE')
	allele = Field(String(5))	#'call' is mysql reserved keyword
	using_options(tablename='calls')
	using_table_options(mysql_engine='InnoDB')
	
class Calls_BySeq(Entity):
	snp = ManyToOne("SNPs", colname='snpid', ondelete='CASCADE', onupdate='CASCADE')
	ecotype = ManyToOne("Ecotype", colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	extraction = ManyToOne("Extraction", colname='extractionid', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo = ManyToOne("SeqInfo", colname='seqinfoid', ondelete='CASCADE', onupdate='CASCADE')
	plateid = Field(String(25))
	wellid = Field(String(3))
	call1 = Field(String(1))
	call2 = Field(String(1))
	ext1area = Field(Float)
	ext2area = Field(Float)
	probearea = Field(Float)
	ext1height = Field(Float)
	ext2height = Field(Float)
	probeheight = Field(Float)
	ext1snr = Field(Float)
	ext2snr = Field(Float)
	probesnr = Field(Float)
	replicate = Field(Boolean)
	description = Field(String(40))
	using_options(tablename='calls_byseq')
	using_table_options(mysql_engine='InnoDB')

class StockDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['localhost', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['stock', 'd', 1, '',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('username', 1, ):[None, 'u', 1, 'database username',],\
							('password', 1, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-07-09
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		if getattr(self, 'schema', None):	#for postgres
			for entity in entities:
				using_table_options_handler(entity, schema=self.schema)
		
		metadata.bind = self._url
		setup_all(create_tables=True)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = StockDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	