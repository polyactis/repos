"""A standalone page showing screenings of a particular film at a particular
cinema. This view is registered for ICinema, and takes the film as a request
parameter.
"""

from datetime import datetime, timedelta

from zope.component import getUtility

from Acquisition import aq_inner
from AccessControl import getSecurityManager

from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from Products.CMFCore.utils import getToolByName

from plone.memoize.instance import memoize

from optilux.cinemacontent.interfaces import IScreeningLocator

from optilux.cinemacontent import CinemaMessageFactory as _
from optilux.cinemacontent import config

class CinemaScreeningsView(BrowserView):
    """List screenings of a film at a cinema
    """
        
    __call__ = ViewPageTemplateFile('screenings.pt')
    
    @memoize
    def upcoming_screenings(self, days=14):
        cinema = aq_inner(self.context)
        locator = getUtility(IScreeningLocator)
        
        from_date = datetime.now()
        to_date = from_date + timedelta(days)
        
        can_reserve = getSecurityManager().checkPermission(config.MAKE_RESERVATION_PERMISSION, cinema)
        
        film = self.film()
        return [dict(screening_id=screening.screening_id,
                     show_time=self.localize(screening.show_time),
                     remaining_tickets=screening.remaining_tickets,
                     can_reserve=(can_reserve and screening.remaining_tickets > 0))
                    for screening in 
                        locator.screenings(film, cinema, from_date, to_date)]
        
    @memoize
    def film(self):
        context = aq_inner(self.context)
        film_code = self.request['film_code']
        catalog = getToolByName(context, 'portal_catalog')
        return catalog(film_code=film_code)[0].getObject()
        
    def localize(self, time):
        return self._time_localizer()(time.isoformat(), 
                                        long_format=True, 
                                        context=aq_inner(self.context),
                                        domain='plonelocales')

    @memoize
    def _time_localizer(self):
        context = aq_inner(self.context)
        translation_service = getToolByName(context, 'translation_service')
        return translation_service.ulocalized_time