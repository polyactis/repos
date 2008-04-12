"""Define a browser view for the Promotion content type. In the FTI 
configured in profiles/default/types/*.xml, this is being set as the default
view of that content type.
"""

from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from optilux.cinemacontent.interfaces import IBannerProvider

class PromotionView(BrowserView):
    """Default view of a promotion
    """
        
    __call__ = ViewPageTemplateFile('promotion.pt')
    
    def banner_tag(self):
        context = aq_inner(self.context)
        banner_provider = IBannerProvider(context)
        return banner_provider.tag