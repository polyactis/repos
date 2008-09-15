#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	StockDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	StockDB.py -z papaya -u yh
	
Description:
	2008-08-11
	This is a wrapper for the stock database, build on top of elixir. Several tables need to be manually created because
	elixir can't represent 'unsigned' integer very well.
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
from elixir import OneToMany, ManyToOne, ManyToMany, OneToOne
from elixir import setup_all, entities, create_all
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from sqlalchemy import UniqueConstraint
from sqlalchemy.schema import ThreadLocalMetaData
from sqlalchemy.orm import scoped_session, sessionmaker
from datetime import datetime

from pymodule.db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, transactional=False))
__metadata__ = ThreadLocalMetaData()

class README(Entity):
	#2008-08-07
	title = Field(String(2000))
	description = Field(String(60000))
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class BatchEcotype(Entity):
	batch = ManyToOne('Batch', colname='batchid', ondelete='CASCADE', onupdate='CASCADE')
	ecotype = ManyToOne('Ecotype', colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	lineorder = Field(Integer, nullable=False)
	using_options(tablename='batch_ecotype', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Batch(Entity):
	batchname = Field(String(30), nullable=False)
	using_options(tablename='batch', metadata=__metadata__, session=__session__)
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
	using_options(tablename='person', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Country(Entity):
	name = Field(String(100))
	abbr = Field(String(10))
	capital = Field(Text)
	latitude = Field(Float)
	longitude = Field(Float)
	using_options(tablename='country', metadata=__metadata__, session=__session__)
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
	using_options(tablename='address', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Site(Entity):
	address = ManyToOne("Address", colname='addressid', ondelete='CASCADE', onupdate='CASCADE')
	name = Field(String(100))
	description = Field(String(500))
	using_options(tablename='site', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Organism(Entity):
	genus = Field(String(30))
	species = Field(String(30))
	using_options(tablename='organism', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class SNPs(Entity):
	organism = ManyToOne("Organism", colname='organismid', ondelete='CASCADE', onupdate='CASCADE')
	snpid = Field(String(20), unique=True)
	dir = Field(String(1))
	chromosome = Field(Integer)
	position = Field(Integer)
	refcall = Field(String(1))
	using_options(tablename='snps', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class SeqInfo(Entity):
	platename = Field(String(40))
	experiment = Field(String(40))
	chip = Field(String(40))
	experiment_type = Field(String(40))
	ms_name = Field(String(40))
	creation_date = Field(DateTime, default=datetime.now)
	using_options(tablename='seqinfo', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class ContaminantType(Entity):
	"""
	2008-08-11
	"""
	short_name = Field(String(100), unique=True)
	description = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='contaminant_type', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Strain(Entity):
	ecotype = ManyToOne("Ecotype", colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	extraction = ManyToOne("Extraction", colname='extractionid', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo1 = ManyToOne("SeqInfo", colname='seqinfoid1', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo2 = ManyToOne("SeqInfo", colname='seqinfoid2', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo3 = ManyToOne("SeqInfo", colname='seqinfoid3', ondelete='CASCADE', onupdate='CASCADE')
	seqinfo4 = ManyToOne("SeqInfo", colname='seqinfoid4', ondelete='CASCADE', onupdate='CASCADE')
	plateid = Field(String(25))
	wellid = Field(String(3))
	replicate = Field(Boolean)
	contaminant_type = ManyToOne("ContaminantType", colname='contaminant_type_id', ondelete='CASCADE', onupdate='CASCADE')
	call_qc_ls = OneToMany("CallQC")
	ecotypeid_strainid2tg_ecotypeid = OneToOne("EcotypeIDStrainID2TGEcotypeID", inverse="strain")
	created_by = Field(String(128))
	updated_by = Field(String(128))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='strain', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('ecotypeid', 'plateid', 'wellid'))
	
class Extraction(Entity):
	namelabel = Field(String(25))
	datesubmitted = Field(DateTime, default=datetime.now)
	notes = Field(String(255))
	using_options(tablename='extraction', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Trip(Entity):
	collectiondate = Field(DateTime)
	notes = Field(String(500))
	using_options(tablename='trip', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class Cross(Entity):
	maternalid = Field(Integer)
	paternalid = Field(Integer)
	name = Field(String(40))
	barcode = Field(String(25))
	generation = Field(Integer)
	using_options(tablename='crosses', metadata=__metadata__, session=__session__)
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
	using_options(tablename='ecotype', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	
class Calls(Entity):
	strain = ManyToOne("Strain", colname='strainid', ondelete='CASCADE', onupdate='CASCADE')
	snp = ManyToOne("SNPs", colname='snpid', ondelete='CASCADE', onupdate='CASCADE')
	allele = Field(String(5))	#'call' is mysql reserved keyword
	using_options(tablename='calls', metadata=__metadata__, session=__session__)
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
	using_options(tablename='calls_byseq', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class QCMethod(Entity):
	"""
	2008-08-17
	"""
	short_name = Field(String(30), unique=True)
	data1_type = Field(String(30), nullable=False)
	data2_type = Field(String(30), nullable=False)
	method_description = Field(Text)
	data_description = Field(Text)
	comment = Field(Text)
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='qc_method', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class CallQC(Entity):
	"""
	2008-08-17
	"""
	strain = ManyToOne("Strain", colname='strainid', ondelete='CASCADE', onupdate='CASCADE')
	ecotype = ManyToOne("Ecotype", colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	target_id = Field(Integer)
	NA_rate = Field(Float)
	no_of_NAs = Field(Integer)
	no_of_totals = Field(Integer)
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	qc_method = ManyToOne("QCMethod", colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(200))
	updated_by = Field(String(200))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='call_qc', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')
	using_table_options(UniqueConstraint('strainid', 'qc_method_id', 'target_id'))

class EcotypeIDStrainID2TGEcotypeID(Entity):
	"""
	2008-08-18
	"""
	ecotype = ManyToOne("Ecotype", colname='ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	strain = ManyToOne("Strain", colname='strainid', ondelete='CASCADE', onupdate='CASCADE')
	tg_ecotype = ManyToOne("Ecotype", colname='tg_ecotypeid', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='ecotypeid_strainid2tg_ecotypeid', metadata=__metadata__, session=__session__)
	using_table_options(mysql_engine='InnoDB')

class QCCrossMatch(Entity):
	"""
	2008-08-26
	"""
	strain = ManyToOne("Strain", colname='strainid', ondelete='CASCADE', onupdate='CASCADE')
	target_id = Field(Integer)
	qc_method = ManyToOne("QCMethod", colname='qc_method_id', ondelete='CASCADE', onupdate='CASCADE')
	mismatch_rate = Field(Float)
	no_of_mismatches = Field(Integer)
	no_of_non_NA_pairs = Field(Integer)
	readme = ManyToOne("README", colname='readme_id', ondelete='CASCADE', onupdate='CASCADE')
	using_options(tablename='qc_cross_match')
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
				if entity.__module__==self.__module__:	#entity in the same module
					using_table_options_handler(entity, schema=self.schema)
		
		__metadata__.bind = self._url
		self.metadata = __metadata__
		self.session = __session__
	
	def setup(self, create_tables=True):
		"""
		2008-08-27
			add option create_tables to disable create_tables in the case that 
				another database metadata gets added into a program because it's imported but no binding at all.
		2008-08-26
		"""
		setup_all(create_tables=create_tables)	#create_tables=True causes setup_all to call elixir.create_all(), which in turn calls metadata.create_all()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = StockDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.setup()
	import pdb
	pdb.set_trace()
