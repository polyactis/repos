from persistent import Persistent

from zope.interface import implements
from zope.component import getUtility

from collective.lead import Database
from Variation.SNP250k.interfaces import IStockDatabaseSettings

from sqlalchemy.engine.url import URL
from sqlalchemy import Table, mapper, relation

#from Variation.SNP250k.dbphenotype import PhenotypeAvg, PhenotypeMethod, QCMethod
#from optilux.cinemacontent.reservation import Results
from variation.src.db import Stock_250kDatabase

class StockDatabaseSettings(Persistent):
	"""Database connection settings
	
	We use raw fields here so that we can more easily use a zope.formlib
	form in the control panel to configure it. This is registered as a
	persistent local utility, with name 'optilux.reservations', which is
	then used by collective.lead.interfaces.IDatabase to find connection settings.
	"""
	
	implements(IStockDatabaseSettings)
	
	drivername = 'mysql'
	hostname = 'papaya.usc.edu'
	port = None
	username = 'nordborglab'
	password = 'papaya'
	database = 'stock_250k'

class StockDatabase(Database):
	"""The reservations database - registered as a utility providing
	collective.lead.interfaces.IDatabase and named 'optilux.reservations'
	"""
	
	@property
	def _url(self):
		settings = getUtility(IStockDatabaseSettings, name='variation.stockdatabasesettings')
		#import sys
		#sys.stderr.write(settings.username + '\n')
		#sys.stderr.write(settings.password + '\n')
		return URL(drivername=settings.drivername, username=settings.username,
				   password=settings.password, host=settings.hostname,
				   port=settings.port, database=settings.database)
	
	def _setup_tables(self, metadata, tables):
		"""Map the database structure to SQLAlchemy Table objects
		"""
		Stock_250kDatabase._setup_tables(metadata, tables)
		"""
		#2008-05-21
		tables['phenotype_avg'] = Table('phenotype_avg', metadata, autoload=True)
		tables['phenotype_method'] = Table('phenotype_method', metadata, autoload=True)
		tables['QC_method'] = Table('QC_method', metadata, autoload=True)
		#tables['reservation'] = Table('reservation', metadata, autoload=True)
		"""
	
	def _setup_mappers(self, tables, mappers):
		"""Map the database Tables to SQLAlchemy Mapper objects
		"""
		Stock_250kDatabase._setup_mappers(tables, mappers)
		"""
		#2008-05-21
		mappers['phenotype_method'] = mapper(PhenotypeMethod, tables['phenotype_method'])
		mappers['phenotype_avg'] = mapper(PhenotypeAvg, tables['phenotype_avg'],
										properties={'phenotype_method': relation(PhenotypeMethod),})
		mappers['QC_method'] = mapper(QCMethod, tables['QC_method'])
		#mappers['reservation'] = mapper(Reservation, tables['reservation'],
		#							   properties = {
		#									'screening' : relation(Screening),
		#									})
		"""