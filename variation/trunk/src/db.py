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


class a(object):
	def __init__(self, **keywords):
		self.hostname = keywords['hostname']
	
	@property
	def _initvalue(self):
		return 'abc'

class Results(object):
	pass

class ResultsMethod(object):
	pass

class PhenotypeAvg(object):
	pass

class PhenotypeMethod(object):
	pass

class QCMethod(object):
	pass

class Stock_250kDatabase(object):
	__doc__ = __doc__
	option_default_dict = {('v', 'drivername', 1, '', 1, ):'mysql',\
							('z', 'hostname', 1, '', 1, ):'papaya.usc.edu',\
							('d', 'database',1, '', 1, ):'stock_250k',\
							('u', 'username',1, '', 1, ):None,\
							('p', 'password',1, '', 1, ):None,\
							('o', 'port', 1, '', 0, ):None,\
							('c', 'commit', 0, '', 0, int):0,\
							('b', 'debug', 0, '', 0, int):0,\
							('r', 'report', 0, '', 0, int):0}
	"""
	2008-02-28
		argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
		argument_type is optional
	"""
	def __init__(self, **keywords):
		from pymodule import process_function_arguments, turn_option_default_dict2argument_default_dict
		argument_default_dict = turn_option_default_dict2argument_default_dict(self.option_default_dict)
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
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
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		mappers['phenotype_method'] = mapper(PhenotypeMethod, tables['phenotype_method'])
		mappers['phenotype_avg'] = mapper(PhenotypeAvg, tables['phenotype_avg'],
										properties={'method_id': relation(PhenotypeMethod),}, allow_column_override=True)
		mappers['QC_method'] = mapper(QCMethod, tables['QC_method'])
		mappers['results'] = mapper(Results, tables['results'], properties={'method_id': relation(ResultsMethod), 'phenotype_method_id': relation(PhenotypeMethod)}, allow_column_override=True)
		mappers['results_method'] = mapper(ResultsMethod, tables['results_method'])
	
	@property
	def _engine_properties(self):
		return {}
	
	def invalidate(self):
		self._initialize_engine()
		
	# IDatabase implementation - code using (not setting up) the database
	# uses this
	
	@property
	def session(self):
		if getattr(self._threadlocal, 'session', None) is None:
			# Without this, we may not have mapped things properly, nor
			# will we necessarily start a transaction when the client
			# code begins to use the session.
			ignore = self.engine
			self._threadlocal.session = Session()
		return self._threadlocal.session
	
	@property
	def connection(self):
		return self.engine.contextual_connect()
	
	@property
	def engine(self):
		if self._engine is None:
			self._initialize_engine()
		
		return self._engine
	
	# Helper methods
	
	def _initialize_engine(self):
		kwargs = dict(self._engine_properties).copy()
		if 'strategy' not in kwargs:
			kwargs['strategy'] = 'threadlocal'
		if 'convert_unicode' not in kwargs:
			kwargs['convert_unicode'] = True
		
		engine = sqlalchemy.create_engine(self._url, **kwargs)
		metadata = sqlalchemy.MetaData(engine)
		
		# We will only initialize once, but we may rebind metadata if
		# necessary

		if not self.tables:
			self._setup_tables(metadata, self.tables)
			self._setup_mappers(self.tables, self.mappers)
		else:
			for name, table in self.tables.items():
				self.tables[name] = table.tometadata(self._metadata)
		
		self._engine = engine
		self._metadata = metadata

if __name__ == '__main__':
	from pymodule import process_options, generate_program_doc
	main_class = Stock_250kDatabase
	opts_dict = process_options(sys.argv, main_class.option_default_dict, error_doc=generate_program_doc(sys.argv[0], main_class.option_default_dict)+main_class.__doc__)
	
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
	instance = main_class(**opts_dict)
	if instance.debug:
		import pdb
		pdb.set_trace()
	print dir(instance)
	session = instance.session
	for row in session.query(ResultsMethod).list():
		print row.id
		print row.short_name
		print row.method_description
