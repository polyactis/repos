import unittest
import doctest

from zope.testing import doctestunit
from zope.component import testing, eventtesting

from Testing import ZopeTestCase as ztc

from optilux.cinemacontent.tests import base

def test_suite():
    return unittest.TestSuite([

        # Demonstrate the main content types
        ztc.ZopeDocFileSuite(
            'README.txt', package='optilux.cinemacontent',
            test_class=base.CinemaContentFunctionalTestCase,
            optionflags=doctest.REPORT_ONLY_FIRST_FAILURE | doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),

        # Test database connectivity
        ztc.ZopeDocFileSuite(
            'tests/database.txt', package='optilux.cinemacontent',
            test_class=base.CinemaContentFunctionalTestCase,
            optionflags=doctest.REPORT_ONLY_FIRST_FAILURE | doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),
            
        # Test the database control panel
        ztc.ZopeDocFileSuite(
            'tests/db_controlpanel.txt', package='optilux.cinemacontent',
            test_class=base.CinemaControlPanelTestCase,
            optionflags=doctest.REPORT_ONLY_FIRST_FAILURE | doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),

        # Test the KSS functionality for ratings
        ztc.ZopeDocFileSuite(
            'tests/dynamic_ratings.txt', package='optilux.cinemacontent',
            test_class=base.CinemaContentKSSTestCase,
            optionflags=doctest.REPORT_ONLY_FIRST_FAILURE | doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),

        # Test the KSS functionality for "my cinema"
        ztc.ZopeDocFileSuite(
            'tests/dynamic_mycinema.txt', package='optilux.cinemacontent',
            test_class=base.CinemaContentKSSTestCase,
            optionflags=doctest.REPORT_ONLY_FIRST_FAILURE | doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),

        ])

if __name__ == '__main__':
    unittest.main(defaultTest='test_suite')
