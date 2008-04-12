"""A report of recently modified cinemas and films
"""

from DateTime import DateTime

from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from Products.CMFCore.utils import getToolByName

from optilux.cinemacontent.interfaces import IFilm
from optilux.cinemacontent.interfaces import ICinema

from plone.memoize.instance import memoize

class ActivityReportView(BrowserView):
    """View for showing recent cinema and film modifications
    """
        
    template = ViewPageTemplateFile('activityreport.pt')
        
    def __call__(self):
        # Hide the edtiable-object border
        self.request.set('disable_border', True)
        
        # Ensure we have a sensible number for days; since this is a non-
        # critical field, we fall silent back on the default if the input is
        # invalid
        try:
            self.days = int(self.request.get('days', 7))
        except ValueError:
            self.days = 7

        return self.template()
        
    def recently_modified_films(self):
        context = aq_inner(self.context)
        catalog = getToolByName(context, 'portal_catalog')
        results = []
        for r in catalog(object_provides=IFilm.__identifier__,
                         modified=dict(query=self.modified_after(), range='min'),
                         sort_on='modified',
                         sort_order='reverse',):
            results.append(dict(url=r.getURL(),
                                title=r.Title,
                                description=r.Description,
                                modified=self.localize(r.modified)))
        return results
        
    def recently_modified_cinemas(self):
        context = aq_inner(self.context)
        catalog = getToolByName(context, 'portal_catalog')
        results = []
        for r in catalog(object_provides=ICinema.__identifier__,
                         modified=dict(query=self.modified_after(), range='min'),
                         sort_on='modified',
                         sort_order='reverse',):
            results.append(dict(url=r.getURL(),
                                title=r.Title,
                                description=r.Description,
                                modified=self.localize(r.modified)))
        return results
        
    def localize(self, time):
        return self._time_localizer()(time, None, aq_inner(self.context), domain='plonelocales')
        
    def modified_after(self):
        return DateTime() - self.days
        
    @memoize
    def _time_localizer(self):
        context = aq_inner(self.context)
        translation_service = getToolByName(context, 'translation_service')
        return translation_service.ulocalized_time