from zope.interface import implements
from zope.viewlet.interfaces import IViewlet

from Acquisition import aq_inner
from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile
from Products.CMFCore.utils import getToolByName

from optilux.cinemacontent.interfaces import IRatings

class MyCinemaViewlet(BrowserView):
    """Viewlet for allowing users to set a "home" cinema
    """
    implements(IViewlet)
    
    template = ViewPageTemplateFile('mycinema.pt')

    def render(self):
        if self.available():
            return self.template()
        else:
            return ''

    def __init__(self, context, request, view, manager):
        super(MyCinemaViewlet, self).__init__(context, request)
        self.__parent__ = view
        self.view = view
        self.manager = manager
        self.membership = getToolByName(context, 'portal_membership')

    def update(self):
        
        if not self.available():
            self.is_home = False
        else:
            cinema = aq_inner(self.context)        
            member = self.membership.getAuthenticatedMember()
    
            cinema_code = cinema.cinema_code
            home_cinemas = list(member.getProperty('home_cinemas', []))
        
            if self.request.has_key('optilux.cinemacontent.mycinema.Toggle'):
                if cinema_code in home_cinemas:
                    home_cinemas.remove(cinema_code)
                else:
                    home_cinemas.append(cinema_code)
                member.setProperties(home_cinemas=home_cinemas)
        
            self.is_home = (cinema_code in home_cinemas)
        
    def available(self):
        return not self.membership.isAnonymousUser()