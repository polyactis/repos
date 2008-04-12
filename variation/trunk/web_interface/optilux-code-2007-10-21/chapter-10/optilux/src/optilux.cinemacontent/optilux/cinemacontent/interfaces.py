from zope.interface import Interface
from zope import schema

from zope.app.container.constraints import contains

from optilux.cinemacontent import CinemaMessageFactory as _

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