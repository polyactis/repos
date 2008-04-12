from zope.interface import implements
from zope.component import getUtility
from zope.app.component.hooks import getSite

from Products.CMFCore.interfaces import ISiteRoot
from Products.CMFCore.utils import getToolByName

from optilux.cinemacontent.interfaces import IFilm
from optilux.cinemacontent.interfaces import ICinema
from optilux.cinemacontent.interfaces import IScreening
from optilux.cinemacontent.interfaces import IScreeningLocator

import sqlalchemy as sql
from collective.lead.interfaces import IDatabase

class Screening(object):
    """A screening of a film at a particular cinema
    """
    
    implements(IScreening)
    
    screening_id = None
    cinema = None
    film = None
    show_time = None
    remaining_tickets = 0
    
    def __init__(self, cinema, film, show_time, remaining_tickets):
        self.cinema = cinema
        self.film = film
        self.show_time = show_time
        self.remaining_tickets = remaining_tickets
        
class ScreeningLocator(object):
    """Find screenings of films at cinemas
    """
    
    implements(IScreeningLocator)
    
    def films_at_cinema(self, cinema, from_date, to_date):
        """Return a list of all films showing at the particular ICinema
        between the specified dates.
        
        Returns a list of dictionaries with keys 'film_code', 'url', 'title' 
        and 'summary'.
        """
        
        # Set up and issue the query, making sure we use the same transaction
        # context as SQLAlchemy's session. 
        
        db = getUtility(IDatabase, name='optilux.reservations')
        connection = db.connection
                
        statement = sql.select([Screening.c.film_code],
                               sql.and_(
                                    Screening.c.cinema_code == cinema.cinema_code,
                                    Screening.c.show_time.between(from_date, to_date)
                               ),
                               distinct=True)
        
        results = connection.execute(statement).fetchall()
        
        # Now use the catalog to find films for the returned film codes
        film_codes = [row['film_code'] for row in results]
        
        # We don't have a context to acquire from here, but we can get
        # it from zope.app.component.hooks.getSite(), which returns the
        # 'site' (in the Component Architecture sense) set during traversal.
        # Since we know that the Plone site root is a "local site" we can
        # be sure that either this, or a sub-site, will be returned when
        # this function is invoked from somewhere within the Plone site.

        site = getSite()
        catalog = getToolByName(site, 'portal_catalog')
        
        return [ dict(film_code=film.film_code,
                      url=film.getURL(),
                      title=film.Title,
                      summary=film.Description,)
                 for film in 
                    catalog(object_provides=IFilm.__identifier__,
                            film_code=film_codes,
                            sort_on='sortable_title')
               ]
        
    def cinemas_for_film(self, film, from_date, to_date):
        """Return a list of all cinemas showing the given film between the
        specified dates.
        
        Returns a list of dictionaries with keys 'cinema_code', 'url', 'name' 
        and 'address'.
        """
        
        db = getUtility(IDatabase, name='optilux.reservations')
        connection = db.connection
        
        statement = sql.select([Screening.c.cinema_code],
                               sql.and_(
                                    Screening.c.film_code == film.film_code,
                                    Screening.c.show_time.between(from_date, to_date)
                               ),
                               distinct=True)
        
        results = connection.execute(statement).fetchall()
        
        cinema_codes = [row['cinema_code'] for row in results]
        
        site = getSite()
        catalog = getToolByName(site, 'portal_catalog')
        
        return [ dict(cinema_code=cinema.cinema_code,
                      url=cinema.getURL(),
                      name=cinema.Title,
                      address=cinema.Description,)
                 for cinema in 
                    catalog(object_provides=ICinema.__identifier__,
                            cinema_code=cinema_codes,
                            sort_on='sortable_title')
               ]
        
    def screenings(self, film, cinema, from_date, to_date):
        """Return all screenings of the given film, at the given cinema,
        between the given dates.
        
        Returns a list of IScreening objects.
        """
        
        db = getUtility(IDatabase, name='optilux.reservations')
        session = db.session
        
        screenings = session.query(Screening).select(sql.and_(
                                                        Screening.c.film_code==film.film_code,
                                                        Screening.c.cinema_code==cinema.cinema_code,
                                                        Screening.c.show_time.between(from_date, to_date)
                                                     ), order_by=[Screening.c.show_time])

        # Now set the 'film' and 'cinema' properties
        
        for screening in screenings:
            screening.film = film
            screening.cinema = cinema
        
        return screenings
        
    def screening_by_id(self, screening_id):
        """Get an IScreening from a screening id
        """
        
        db = getUtility(IDatabase, name='optilux.reservations')
        session = db.session
        
        screening = session.query(Screening).get(screening_id)
        
        site = getSite()
        catalog = getToolByName(site, 'portal_catalog')
        
        # Set cinema and film. We assume these are in the portal and that 
        # there is a one-to-one mapping between cinema/film codes and 
        # cinemas/films.

        screening.cinema = catalog(cinema_code=screening.cinema_code)[0].getObject()
        screening.film = catalog(film_code=screening.film_code)[0].getObject()
        
        return screening