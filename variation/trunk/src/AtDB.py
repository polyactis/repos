#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	AtDB.py -v postgres -u crocea -z localhost -d graphdb -k public
	
	#setup database in mysql
	AtDB.py -u yh
	
Description:
	2008-07-31
	This is a wrapper for the at database, build on top of elixir.
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

class Population(Entity):
	name = Field(String(50))
	region_obj = ManyToOne('Region', colname='region')
	description = Field(Text)
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	latitude = Field(Float)
	longitude = Field(Float)
	using_options(tablename='population')
	using_table_options(mysql_engine='InnoDB')
	
class Region(Entity):
	name = Field(String(50))
	description = Field(Text)
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	using_options(tablename='region')
	using_table_options(mysql_engine='InnoDB')
	
class Accession(Entity):
	population_obj = ManyToOne('Population', colname='population')
	region_obj = ManyToOne('Region', colname='region')
	description = Field(Text)
	name = Field(String(20))
	origin = Field(String(50))
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	number = Field(String(10))
	using_options(tablename='accession')
	using_table_options(mysql_engine='InnoDB')

class Submitter(Entity):
	first = Field(String(50))
	last = Field(String(50))
	organization = Field(String(50))
	email = Field(String(50))
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	using_options(tablename='submitter')
	using_table_options(mysql_engine='InnoDB')

class Alignment(Entity):
	chromosome = Field(Integer)
	start = Field(Integer)
	end = Field(Integer)
	target = Field(Text)
	version = Field(Integer)
	directory = Field(String(255))
	code = Field(Integer)
	note = Field(Text)
	submitter_obj = ManyToOne('Submitter', colname='submitter')
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	plus_encode = Field(Text)
	minus_encode = Field(Text)
	using_options(tablename='alignment')
	using_table_options(mysql_engine='InnoDB')

class Sequence(Entity):
	accession_obj = ManyToOne('Accession', colname='accession')
	alignment_obj = ManyToOne('Alignment', colname='alignment')
	start = Field(Integer)
	bases = Field(Text)
	quality = Field(Text)
	_quality = Field(Text)
	created = Field(DateTime, default=datetime.now)
	modified = Field(DateTime)
	using_options(tablename='sequence')
	using_table_options(mysql_engine='InnoDB')

class AtDB(ElixirDB):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ):['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('database', 1, ):['at', 'd', 1, '',],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('username', 0, ):[None, 'u', 1, 'database username',],\
							('password', 0, ):[None, 'p', 1, 'database password', ],\
							('port', 0, ):[None, 'o', 1, 'database port number'],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-07-31
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
	main_class = AtDB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	
	#test how to do raw sql on a view
	accession_id=1
	rows = instance.metadata.bind.execute("select * from %s where accession_id=%s"%('accession2tg_ecotypeid', accession_id))
	row = rows.fetchone()
	print "ecotype_id for accession_id %s is %s.\n"%(accession_id, row.ecotype_id)
	"""
	rows = Accession.query.all()
	for row in rows:
		print row.region_obj.name
	"""