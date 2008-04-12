from zope.interface import implements, alsoProvides
from zope.component import getMultiAdapter

from zope.viewlet.interfaces import IViewlet
from zope.viewlet.interfaces import IViewletManager

from kss.core import kssaction
from plone.app.kss.plonekssview import PloneKSSView

from plone.app.layout.globals.interfaces import IViewView

from Acquisition import aq_inner
from Products.Five.browser import BrowserView
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile
from Products.CMFCore.utils import getToolByName

from optilux.cinemacontent.interfaces import ICinema
from optilux.cinemacontent.browser.interfaces import IMyCinemasBox
from optilux.cinemacontent import CinemaMessageFactory as _

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
        
class DynamicMyCinema(PloneKSSView):
    """When the user clicks to assign or unassign a "home" cinema, refresh
    the viewlet and the portlet.
    """
    
    implements(IMyCinemasBox)
    
    # By annotating the method with @kssaction are able to issue KSS
    # commands
    
    @kssaction
    def toggle(self):
        
        portal_state = getMultiAdapter((self.context, self.request), name=u"plone_portal_state")
        if portal_state.anonymous():
            return
        
        cinema = aq_inner(self.context)
        member = portal_state.member()
        
        cinema_code = cinema.cinema_code
        home_cinemas = list(member.getProperty('home_cinemas', []))
        
        # Add or remove this cinema from the list of home cinemas
        
        enabled = True
        if cinema_code in home_cinemas:
            home_cinemas.remove(cinema_code)
            enabled = False
        else:
            home_cinemas.append(cinema_code)
        
        member.setProperties(home_cinemas=home_cinemas)
    
        # Refresh the viewlet with the toggle button
        # See ratings.py for more on this set of commands

        alsoProvides(self, IViewView)

        ksscore = self.getCommandSet('core')
        zopecommands = self.getCommandSet('zope')
        
        selector = ksscore.getHtmlIdSelector('my-cinema-toggle')
        zopecommands.refreshViewlet(selector, manager='plone.belowcontentbody', 
                                    name='optilux.cinemacontent.mycinema')
        
        plonecommands = self.getCommandSet('plone')
        
        # Issue a status message that things have changed
        
        if enabled:
            plonecommands.issuePortalMessage(_(u"This cinema is now a 'home' cinema."))
        else:
            plonecommands.issuePortalMessage(_(u"This cinema is no longer a 'home' cinema."))
            
        # Refresh any instances of the "my cinema" portlet. Here, we cheat 
        # and simply re-use the template from the mycinema portlet.

        # This approach isn't entirely perfect - if the portlet wasn't 
        # displayed before, it won't be shown until the page is reloaded,
        # and we may end up with an empty (non-hidden) portlet column if
        # the portlet ends up not being rendered.
        
        if not home_cinemas:
            # If there are no home cinemas, replace the portlet <dl /> wiht
            # a blank <div /> with the same class. We do this so that we can
            # find the node again later if the user toggles the cinema back on

            ksscore.replaceHTML('.portletMyCinema', '<div class="portletMyCinema" />')
        else:
            # There are cinemas to display - render the portlet template. This
            # view will be the 'view' variable in the template. This is okay,
            # because we provide the IMyCinemasBox interface, which as the
            # template is expecting.
            
            self.home_cinemas = home_cinemas
            ksscore.replaceHTML('.portletMyCinema', self.portlet_template())
 
    portlet_template = ViewPageTemplateFile('../portlets/mycinema.pt') 
    def cinemas(self):
        """List cinemas to display
        """
        context = aq_inner(self.context)
        cinema_codes = self.home_cinemas
        catalog = getToolByName(context, 'portal_catalog')
        for brain in catalog(cinema_code = cinema_codes,
                             object_provides = ICinema.__identifier__):
            yield dict(title=brain.Title,
                       address=brain.Description,
                       url=brain.getURL())