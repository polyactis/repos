import unittest

from zope.testing.doctestunit import DocTestSuite
from zope.app.testing import placelesssetup

def test_suite():
    return unittest.TestSuite((
        DocTestSuite('optilux.cinemacontent.ratings',
                     setUp=placelesssetup.setUp,
                     tearDown=placelesssetup.tearDown),
        ))

if __name__=='__main__':
    unittest.main(defaultTest='test_suite')
