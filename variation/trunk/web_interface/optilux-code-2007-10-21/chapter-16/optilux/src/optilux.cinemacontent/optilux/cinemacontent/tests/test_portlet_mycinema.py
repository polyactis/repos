from zope.component import getUtility, getMultiAdapter

from plone.portlets.interfaces import IPortletType
from plone.portlets.interfaces import IPortletManager
from plone.portlets.interfaces import IPortletAssignment
from plone.portlets.interfaces import IPortletDataProvider
from plone.portlets.interfaces import IPortletRenderer

from Products.CMFCore.utils import getToolByName

from plone.app.portlets.storage import PortletAssignmentMapping

from optilux.cinemacontent.portlets import mycinema

from optilux.cinemacontent.tests.base import CinemaContentTestCase

class TestPortlet(CinemaContentTestCase):

    def afterSetUp(self):
        self.setRoles(('Manager',))

    def testPortletTypeRegistered(self):
        portlet = getUtility(IPortletType, name='optilux.MyCinema')
        self.assertEquals(portlet.addview, 'optilux.MyCinema')

    def testInterfaces(self):
        portlet = mycinema.Assignment()
        self.failUnless(IPortletAssignment.providedBy(portlet))
        self.failUnless(IPortletDataProvider.providedBy(portlet.data))

    def testInvokeAddview(self):
        portlet = getUtility(IPortletType, name='optilux.MyCinema')
        mapping = self.portal.restrictedTraverse('++contextportlets++plone.leftcolumn')
        for m in mapping.keys():
            del mapping[m]
        addview = mapping.restrictedTraverse('+/' + portlet.addview)

        addview()

        self.assertEquals(len(mapping), 1)
        self.failUnless(isinstance(mapping.values()[0], mycinema.Assignment))

    def testRenderer(self):
        context = self.folder
        request = self.folder.REQUEST
        view = self.folder.restrictedTraverse('@@plone')
        manager = getUtility(IPortletManager, name='plone.rightcolumn', context=self.portal)
        assignment = mycinema.Assignment()

        renderer = getMultiAdapter((context, request, view, manager, assignment), IPortletRenderer)
        self.failUnless(isinstance(renderer, mycinema.Renderer))

class TestRenderer(CinemaContentTestCase):
    
    def afterSetUp(self):
        self.setRoles(('Manager',))
        
        self.portal.invokeFactory('Cinema Folder', 'cinemas')
        self.portal.cinemas.invokeFactory('Cinema', 'c1')
        self.portal.cinemas.c1.cinema_code = 'C1'
        self.portal.cinemas.c1.reindexObject()
        
        self.portal.cinemas.invokeFactory('Cinema', 'c2')
        self.portal.cinemas.c2.cinema_code = 'C2'
        self.portal.cinemas.c2.reindexObject()
        
        self.membership = getToolByName(self.portal, 'portal_membership')
        
        self.setRoles(('Member',))

    def renderer(self, context=None, request=None, view=None, manager=None, assignment=None):
        context = context or self.folder
        request = request or self.folder.REQUEST
        view = view or self.folder.restrictedTraverse('@@plone')
        manager = manager or getUtility(IPortletManager, name='plone.rightcolumn', context=self.portal)
        assignment = assignment or mycinema.Assignment()

        return getMultiAdapter((context, request, view, manager, assignment), IPortletRenderer)

    def test_anonymous(self):
        member = self.membership.getAuthenticatedMember()
        member.setProperties(home_cinemas=['C1'])
        self.logout()
        r = self.renderer(context=self.portal, assignment=mycinema.Assignment())
        self.failIf(r.available)

    def test_no_cinemas(self):
        member = self.membership.getAuthenticatedMember()
        member.setProperties(home_cinemas=[])
        r = self.renderer(context=self.portal, assignment=mycinema.Assignment())
        self.failIf(r.available)
        self.assertEquals(0, len(list(r.cinemas())))
        
    def test_single(self):
        member = self.membership.getAuthenticatedMember()
        member.setProperties(home_cinemas=['C1'])
        r = self.renderer(context=self.portal, assignment=mycinema.Assignment())
        self.failUnless(r.available)
        cinema_urls = [c['url'] for c in r.cinemas()]
        self.assertEquals(1, len(cinema_urls))
        self.assertEquals(self.portal.cinemas.c1.absolute_url(), cinema_urls[0])

    def test_multiple(self):
        member = self.membership.getAuthenticatedMember()
        member.setProperties(home_cinemas=['C1', 'C2'])
        r = self.renderer(context=self.portal, assignment=mycinema.Assignment())
        self.failUnless(r.available)
        cinema_urls = [c['url'] for c in r.cinemas()]
        self.assertEquals(2, len(cinema_urls))
        self.assertEquals(self.portal.cinemas.c1.absolute_url(), cinema_urls[0])
        self.assertEquals(self.portal.cinemas.c2.absolute_url(), cinema_urls[1])

def test_suite():
    from unittest import TestSuite, makeSuite
    suite = TestSuite()
    suite.addTest(makeSuite(TestPortlet))
    suite.addTest(makeSuite(TestRenderer))
    return suite