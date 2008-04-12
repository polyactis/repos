from zope.interface import implements, alsoProvides
from zope.component import getMultiAdapter

from zope.viewlet.interfaces import IViewlet

from kss.core import kssaction
from plone.app.kss.plonekssview import PloneKSSView

from plone.app.layout.globals.interfaces import IViewView

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
        
class DynamicRatings(PloneKSSView):
    """When the user clicks to rate a film, perform the rating and show
    the results inline
    """
    
    # By annotating the method with @kssaction are able to issue KSS
    # commands
    
    @kssaction
    def rate_film(self, vote):
        
        # The vote parameter is passed in the cinemacontent.kss file as
        # "yes" or "no"

        vote = vote.lower()
        if vote not in ("yes", "no"):
            return

        # Just to be safe, make sure we can actually rate
        portal_state = getMultiAdapter((self.context, self.request), name=u"plone_portal_state")
        if portal_state.anonymous():
            return
        
        ratings = IRatings(self.context)
        user_token = portal_state.member().getId()    
        
        if user_token is None or not ratings.available(user_token):
            return
            
        # Apply the rating
        ratings.rate(user_token, (vote == "yes"))
        
        # Now send the command back that we should refresh this viewlet
        # Because the viewlet is registered for IViewView, we need to mark
        # the KSS view (self) wit this first.
        
        alsoProvides(self, IViewView)

        ksscore = self.getCommandSet('core')
        zopecommands = self.getCommandSet('zope')

        selector = ksscore.getHtmlIdSelector('film-rating-box')
        zopecommands.refreshViewlet(selector, manager='plone.belowcontenttitle', 
                                    name='optilux.cinemacontent.ratings')