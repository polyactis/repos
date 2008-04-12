from zope.interface import implements
from zope.component import getMultiAdapter
from zope.viewlet.interfaces import IViewlet

from Acquisition import aq_inner
from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from optilux.cinemacontent.interfaces import IRatings

class RatingsViewlet(BrowserView):
    """Viewlet for allowing users to rate a film
    """
    implements(IViewlet)

    render = ViewPageTemplateFile('ratings.pt')

    def __init__(self, context, request, view, manager):
        super(RatingsViewlet, self).__init__(context, request)
        self.__parent__ = view
        self.view = view
        self.manager = manager
        self.ratings = IRatings(self.context)
        self.portal_state = getMultiAdapter((context, self.request), name=u"plone_portal_state")

    def update(self):

        vote = None
        if self.request.has_key('optilux.cinemacontent.ratings.VotePositive'):
            vote = True
        elif self.request.has_key('optilux.cinemacontent.ratings.VoteNegative'):
            vote = False
            
        if vote is None or self.portal_state.anonymous():
            return

        user_token = self.portal_state.member().getId()    
        if user_token is not None and self.ratings.available(user_token):
            self.ratings.rate(user_token, vote)
    
    def have_score(self):
        return self.score() is not None
    
    def available(self):
        if self.portal_state.anonymous():
            return False
        return self.ratings.available(self.portal_state.member().getId())
    
    def score(self):
        return self.ratings.score