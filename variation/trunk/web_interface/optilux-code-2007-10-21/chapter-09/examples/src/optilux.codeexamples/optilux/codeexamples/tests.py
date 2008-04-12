import unittest

from zope.testing import doctestunit
from zope.component import testing, eventtesting

from Testing import ZopeTestCase as ztc

from Products.PloneTestCase import PloneTestCase as ptc
ptc.setupPloneSite()

def test_suite():
    return unittest.TestSuite([

        # Demonstrate containment concepts
        ztc.ZopeDocFileSuite(
            'containment.txt', package='optilux.codeexamples',
            test_class=ptc.PloneTestCase),
             
        # Demonstrate acquisition concepts
        ztc.ZopeDocFileSuite(
            'acquisition.txt', package='optilux.codeexamples',
            test_class=ptc.PloneTestCase),
             
        # Demonstrate use of (un)restrictedTraverse
        ztc.ZopeDocFileSuite(
            'path_traversal.txt', package='optilux.codeexamples',
            test_class=ptc.PloneTestCase),
             
        # Demonstrate use of portal_catalog
        ztc.ZopeDocFileSuite(
            'catalog.txt', package='optilux.codeexamples',
            test_class=ptc.PloneTestCase),
            
        # Demonstrate interfaces
        doctestunit.DocFileSuite(
            'interfaces.txt', package='optilux.codeexamples'),

        # Demonstrate utilities
        doctestunit.DocFileSuite(
            'utilities.txt', package='optilux.codeexamples',
            setUp=testing.setUp, tearDown=testing.tearDown),
        
        # Demonstrate adapters
        doctestunit.DocFileSuite(
            'adapters.txt', package='optilux.codeexamples',
            setUp=testing.setUp, tearDown=testing.tearDown),
            
        # Demonstrate events (notice use of eventesting for setup!)
        doctestunit.DocFileSuite(
            'events.txt', package='optilux.codeexamples',
            setUp=eventtesting.setUp, tearDown=testing.tearDown),

        ])

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')
