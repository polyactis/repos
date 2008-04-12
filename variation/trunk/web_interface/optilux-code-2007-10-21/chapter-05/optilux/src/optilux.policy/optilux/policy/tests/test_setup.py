import unittest
from optilux.policy.tests.base import OptiluxPolicyTestCase

class TestSetup(OptiluxPolicyTestCase):
    
    def test_portal_title(self):
        self.assertEquals("Optilux Cinemas", self.portal.getProperty('title'))
        
    def test_portal_description(self):
        self.assertEquals("Welcome to Optilux Cinemas", self.portal.getProperty('description'))

def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSetup))
    return suite
