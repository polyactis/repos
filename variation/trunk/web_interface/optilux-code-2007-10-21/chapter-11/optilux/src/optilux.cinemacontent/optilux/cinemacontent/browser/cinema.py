"""Define a browser view for the Cinema content type. In the FTI 
configured in profiles/default/types/*.xml, this is being set as the default
view of that content type.
"""

from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from plone.memoize.instance import memoize
from optilux.cinemacontent.interfaces import IBannerProvider

class CinemaView(BrowserView):
    """Default view of a cinema
    """
    
    __call__ = ViewPageTemplateFile('cinema.pt')
    
    def have_highlighted_films(self):
        return len(self.highlighted_films()) > 0
        
    @memoize
    def highlighted_films(self):
        context = aq_inner(self.context)
        return [dict(url=film.absolute_url(),
                     title=film.title,
                     summary=film.summary,
                     banner_tag=IBannerProvider(film).tag,)
                for film in context.highlighted_films]