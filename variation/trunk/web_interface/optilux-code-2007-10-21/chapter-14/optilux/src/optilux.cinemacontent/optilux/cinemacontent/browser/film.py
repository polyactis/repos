"""Define a browser view for the Film content type. In the FTI 
configured in profiles/default/types/*.xml, this is being set as the default
view of that content type.
"""

from datetime import datetime, timedelta

from zope.component import getUtility

from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from plone.memoize.instance import memoize

from optilux.cinemacontent.interfaces import IBannerProvider
from optilux.cinemacontent.interfaces import IScreeningLocator

class FilmView(BrowserView):
    """Default view of a film
    """
    
    __call__ = ViewPageTemplateFile('film.pt')
    
    @memoize
    def cinemas(self, days=14):
        context = aq_inner(self.context)
        locator = getUtility(IScreeningLocator)
        
        from_date = datetime.now()
        to_date = from_date + timedelta(days)
        return locator.cinemas_for_film(context, from_date, to_date)    
    
    def banner_tag(self):
        context = aq_inner(self.context)
        banner_provider = IBannerProvider(context)
        return banner_provider.tag