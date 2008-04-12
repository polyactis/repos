from zope.interface import Interface
from zope.app.container.constraints import contains

from zope import schema

from optilux.cinemacontent import CinemaMessageFactory as _

# Exceptions

class ReservationError(Exception):
    """Exception raised if there is an error making a reservation
    """

    def __init__(self, message):
        Exception.__init__(self, message)
        self.error_message = message

# Basic content types

class IFilmFolder(Interface):
    """A folder containing films
    """
    contains('optilux.cinemacontent.interfaces.IFilm')
    
    title = schema.TextLine(title=_(u"Title"),
                            required=True)
                            
    description = schema.TextLine(title=_(u"Description"),
                                  description=_(u"A short summary of this folder"))
    
class IFilm(Interface):
    """A film
    """
    
    film_code = schema.ASCIILine(title=_(u"Film Code"),
                                 description=_(u"This should match the film code used by the booking system"),
                                 required=True)
    
    title = schema.TextLine(title=_(u"Film title"),
                            required=True)
    
    summary = schema.TextLine(title=_(u"Short summary"),
                              description=_(u"Plain-text blurb about the film"))
    
    teaser = schema.SourceText(title=_(u"Teaser"),
                               description=_(u"A teaser/description of the film"),
                               required=True)
    
    shown_from = schema.Date(title=_(u"Visible from"),
                             description=_(u"Date when film first appears on the website"))
    
    shown_until = schema.Date(title=_(u"Visible until"),
                             description=_(u"Date when film last appears on the website"))

class ICinemaFolder(Interface):
    """A folder containing cinemas
    """
    contains('optilux.cinemacontent.interfaces.ICinema',
             'optilux.cinemacontent.interfaces.IPromotion',)
    
    title = schema.TextLine(title=_(u"Title"),
                            required=True)
                            
    description = schema.TextLine(title=_(u"Description"),
                                  description=_(u"A short summary of this folder"))
                                  
    text = schema.SourceText(title=_(u"Descriptive text"),
                             description=_(u"Descriptive text about this cinema"),
                             required=True)
    
class ICinema(Interface):
    """A cinema
    """
    
    cinema_code = schema.ASCIILine(title=_(u"Cinema Code"),
                                   description=_(u"This should match the cinema code used by "
                                                 "the booking system"),
                                   required=True)
    
    name = schema.TextLine(title=_(u"Cinema name"),
                           required=True)
                            
    phone = schema.TextLine(title=_(u"Telephone number"),
                            description=_(u"Main contact number for this cinema"),
                            required=True)
                            
    address = schema.Text(title=_(u"Address"),
                          description=_(u"Description of this cinema"),
                          required=True)
                            
    text = schema.SourceText(title=_(u"Descriptive text"),
                             description=_(u"Descriptive text about this cinema"),
                             required=True)
                             
    highlighted_films = schema.List(title=_(u"Highlighted films"),
                                     description=_(u"Selected films to highlight"),
                                     value_type=schema.Object(title=_(u"Film"),
                                                              schema=IFilm),
                                     unique=True)
                                
class IPromotion(Interface):
    """A promotion running for one or more cinemas
    """
    
    title = schema.TextLine(title=_(u"Promotion title"),
                            required=True)
                            
    summary = schema.TextLine(title=_(u"Short summary"),
                              description=_(u"Plain-text summary of the promotion"))
    
    details = schema.SourceText(title=_(u"Details"),
                               description=_(u"Details about the promotion"),
                               required=True)
                               
    shown_from = schema.Date(title=_(u"Visible from"),
                             description=_(u"Date when promotion first appears on the website"))
    
    shown_until = schema.Date(title=_(u"Visible until"),
                             description=_(u"Date when promotion last appears on the website"))
    
# Adapters providing additional functionality for content types
    
class IBannerProvider(Interface):
    """A component which can provide an HTML tag for a banner image
    """
    
    tag = schema.TextLine(title=_(u"A HTML tag to render to show the banner image"))
    
class IRatings(Interface):
    """An object which can be rated
    """
    
    score = schema.Int(title=_(u"A score from 1-100"),
                       readonly=True)
         
    def available(user_token):
        """Whether or not rating is available for the given user
        """
                       
    def rate(user_token, positive):
        """Give a positive (True) or negative (False) vote.
        """

# Entities found in the database

class IScreening(Interface):
    """A screening of a film at a particular cinema
    """
    
    screening_id = schema.Int(title=_(u"Screening identifier"),
                              description=_(u"A unique id for this screening"),
                              required=True,
                              readonly=True)
    
    cinema = schema.Object(title=_(u"Cinema"),
                           schema=ICinema,
                           required=True,
                           readonly=True)
                           
    film = schema.Object(title=_(u"Film"),
                         schema=IFilm,
                         required=True,
                         readonly=True)
                           
    show_time = schema.Date(title=_(u"Date/time"),
                            required=True,
                            readonly=True)
    
    remaining_tickets = schema.Int(title=_(u"Remaining tickets"),
                                   description=_(u"Number of tickets available for this screening"))
    
class IReservation(Interface):
    """A ticket reservation for a particular screening
    """
    
    customer_name = schema.TextLine(title=_(u"Customer name"),
                                    description=_(u"The name of the customer making the reservation"),
                                    required=True)
                                    
    num_tickets = schema.Int(title=_(u"Number of tickets"),
                             description=_(u"Number of tickets to reserve"),
                             required=True,
                             min=1)
                             
    screening = schema.Object(title=_(u"Screening"),
                              description=_(u"Film screening to book for"),
                              schema=IScreening,
                              required=True)
                            
# Database services

class IScreeningLocator(Interface):
    """A utility used to locate appropriate screenings based on search criteria
    """
    
    def films_at_cinema(cinema, from_date, to_date):
        """Return a list of all films screening at the particular ICinema
        between the specified dates.
        
        Returns a list of dictionaries with keys 'film_code', 'url', 'title' 
        and 'summary'.
        """
        
    def cinemas_for_film(film, from_date, to_date):
        """Return a list of all cinemas screening the given film between the
        specified dates.
        
        Returns a list of dictionaries with keys 'cinema_code', 'url', 'name' 
        and 'address'.
        """
        
    def screenings(film, cinema, from_date, to_date):
        """Return all screenings of the given film, at the given cinema,
        between the given dates
        
        Returns a list of IScreening objects.
        """
        
    def screening_by_id(screening_id):
        """Get an IScreening from a screening id
        """
        
        
class ITicketReservations(Interface):
    """A utility capable of making reservations
    """
    
    def __call__(reservation):
        """Make a reservation
        """
        
# Database connectivity
        
class IDatabaseSettings(Interface):
    """Database connection settings.
    """
    
    drivername = schema.ASCIILine(title=_(u"Driver name"),
                                  description=_(u"The database driver name"),
                                  default='mysql',
                                  required=True)

    hostname = schema.ASCIILine(title=_(u"Host name"),
                                description=_(u"The database host name"),
                                default='localhost',
                                required=True)
                                
    port = schema.Int(title=_(u"Port number"),
                      description=_(u"The database port number. Leave blank to use the default."),
                      required=False)
                                
    username = schema.ASCIILine(title=_(u"User name"),
                                description=_(u"The database user name"),
                                required=True)

    password = schema.Password(title=_(u"Password"),
                                description=_(u"The database password"),
                                required=False)
                                
    database = schema.ASCIILine(title=_(u"Database name"),
                                description=_(u"The name of the database on this server"),
                                required=True)