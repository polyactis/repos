from zope.component import getUtility, getMultiAdapter

from plone.portlets.interfaces import IPortletType
from plone.portlets.interfaces import IPortletManager
from plone.portlets.interfaces import IPortletAssignment
from plone.portlets.interfaces import IPortletDataProvider
from plone.portlets.interfaces import IPortletRenderer

from plone.app.portlets.storage import PortletAssignmentMapping

from optilux.cinemacontent.portlets import promotions

from optilux.cinemacontent.tests.base import CinemaContentTestCase

class TestPortlet(CinemaContentTestCase):

    def afterSetUp(self):
        self.setRoles(('Manager',))

    def testPortletTypeRegistered(self):
        portlet = getUtility(IPortletType, name='optilux.Promotions')
        self.assertEquals(portlet.addview, 'optilux.Promotions')

    def testInterfaces(self):
        portlet = promotions.Assignment()
        self.failUnless(IPortletAssignment.providedBy(portlet))
        self.failUnless(IPortletDataProvider.providedBy(portlet.data))

    def testInvokeAddview(self):
        portlet = getUtility(IPortletType, name='optilux.Promotions')
        mapping = self.portal.restrictedTraverse('++contextportlets++plone.leftcolumn')
        for m in mapping.keys():
            del mapping[m]
        addview = mapping.restrictedTraverse('+/' + portlet.addview)

        addview.createAndAdd(data={})

        self.assertEquals(len(mapping), 1)
        self.failUnless(isinstance(mapping.values()[0], promotions.Assignment))

    def testInvokeEditView(self):
        mapping = PortletAssignmentMapping()
        request = self.folder.REQUEST

        mapping['foo'] = promotions.Assignment()
        editview = getMultiAdapter((mapping['foo'], request), name='edit')
        self.failUnless(isinstance(editview, promotions.EditForm))

    def testRenderer(self):
        context = self.folder
        request = self.folder.REQUEST
        view = self.folder.restrictedTraverse('@@plone')
        manager = getUtility(IPortletManager, name='plone.rightcolumn', context=self.portal)
        assignment = promotions.Assignment()

        renderer = getMultiAdapter((context, request, view, manager, assignment), IPortletRenderer)
        self.failUnless(isinstance(renderer, promotions.Renderer))

class TestRenderer(CinemaContentTestCase):
    
    def afterSetUp(self):
        self.setRoles(('Manager',))
        self.portal.invokeFactory('Cinema Folder', 'cf1')
        self.portal.invokeFactory('Cinema Folder', 'cf2')
        self.portal.cf1.invokeFactory('Promotion', 'p1')
        self.portal.cf1.invokeFactory('Promotion', 'p2')
        self.portal.cf1.invokeFactory('Promotion', 'p3')
        self.portal.cf1.invokeFactory('Promotion', 'p4')
        self.portal.cf1.invokeFactory('Promotion', 'p5')
        self.portal.cf2.invokeFactory('Promotion', 'p6')
        self.portal.cf2.invokeFactory('Promotion', 'p7')

    def renderer(self, context=None, request=None, view=None, manager=None, assignment=None):
        context = context or self.folder
        request = request or self.folder.REQUEST
        view = view or self.folder.restrictedTraverse('@@plone')
        manager = manager or getUtility(IPortletManager, name='plone.rightcolumn', context=self.portal)
        assignment = assignment or promotions.Assignment()

        return getMultiAdapter((context, request, view, manager, assignment), IPortletRenderer)

    def test_count(self):
        r = self.renderer(context=self.portal.cf1, assignment=promotions.Assignment(count=5))
        self.assertEquals(5, len([p for p in r.promotions()]))

    def test_randomize(self):
        r = self.renderer(context=self.portal.cf1, assignment=promotions.Assignment(count=5, randomize=True))
        self.assertEquals(5, len([p for p in r.promotions()]))
        # Mmmm, hard to test for random things :)
        
    def test_sitewide(self):
        r = self.renderer(context=self.portal.cf1, assignment=promotions.Assignment(count=10, sitewide=True))
        p6_url = self.portal.cf2.p6.absolute_url()
        self.failUnless(p6_url in [p['url'] for p in r.promotions()])

def test_suite():
    from unittest import TestSuite, makeSuite
    suite = TestSuite()
    suite.addTest(makeSuite(TestPortlet))
    suite.addTest(makeSuite(TestRenderer))
    return suite
