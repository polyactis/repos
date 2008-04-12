from persistent import Persistent

from zope.interface import implements
from zope.component import getUtility

from collective.lead import Database
from optilux.cinemacontent.interfaces import IDatabaseSettings

from sqlalchemy.engine.url import URL
from sqlalchemy import Table, mapper, relation

from optilux.cinemacontent.screening import Screening
from optilux.cinemacontent.reservation import Reservation

class ReservationsDatabaseSettings(Persistent):
    """Database connection settings
    
    We use raw fields here so that we can more easily use a zope.formlib
    form in the control panel to configure it. This is registered as a
    persistent local utility, with name 'optilux.reservations', which is
    then used by collective.lead.interfaces.IDatabase to find connection settings.
    """
    
    implements(IDatabaseSettings)
    
    drivername = 'mysql'
    hostname = 'localhost'
    port = None
    username = ''
    password = None
    database = ''

class ReservationsDatabase(Database):
    """The reservations database - registered as a utility providing
    collective.lead.interfaces.IDatabase and named 'optilux.reservations'
    """
    
    @property
    def _url(self):
        settings = getUtility(IDatabaseSettings)
        return URL(drivername=settings.drivername, username=settings.username,
                   password=settings.password, host=settings.hostname,
                   port=settings.port, database=settings.database)
    
    def _setup_tables(self, metadata, tables):
        """Map the database structure to SQLAlchemy Table objects
        """
            
        tables['screening'] = Table('screening', metadata, autoload=True)
        tables['reservation'] = Table('reservation', metadata, autoload=True)
    
    def _setup_mappers(self, tables, mappers):
        """Map the database Tables to SQLAlchemy Mapper objects
        """
        
        mappers['screening'] = mapper(Screening, tables['screening'])
        mappers['reservation'] = mapper(Reservation, tables['reservation'],
                                        properties = {
                                            'screening' : relation(Screening),
                                            })
        
