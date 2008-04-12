"""Define a browser view for the Film content type. In the FTI 
configured in profiles/default/types/*.xml, this is being set as the default
view of that content type.
"""

from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from optilux.cinemacontent.interfaces import IBannerProvider

class FilmView(BrowserView):
    """Default view of a film
    """
    
    __call__ = ViewPageTemplateFile('film.pt')
    
    def banner_tag(self):
        context = aq_inner(self.context)
        banner_provider = IBannerProvider(context)
        return banner_provider.tag