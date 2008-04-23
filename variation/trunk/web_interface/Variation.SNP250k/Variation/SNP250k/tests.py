import unittest
import doctest

from zope.testing import doctestunit
from zope.component import testing
from Testing import ZopeTestCase as ztc

from Products.Five import zcml
from Products.Five import fiveconfigure
from Products.PloneTestCase import PloneTestCase as ptc
from Products.PloneTestCase.layer import onsetup, PloneSite
ptc.setupPloneSite()

import Variation.SNP250k

class TestCase(ptc.FunctionalTestCase):
	class layer(PloneSite):
		@classmethod
		def setUp(cls):
			fiveconfigure.debug_mode = True
			zcml.load_config('configure.zcml',
							 Variation.SNP250k)
			fiveconfigure.debug_mode = False
			
			ztc.installPackage('Variation.SNP250k')

		@classmethod
		def tearDown(cls):
			pass
@onsetup
def setup_variation_snp250k():
    """Set up the additional products required for the Optilux Cinema Content.
    
    The @onsetup decorator causes the execution of this body to be deferred
    until the setup of the Plone site testing layer.
    """
    
    # Load the ZCML configuration for the optilux.policy package.
    # This includes the other products below as well.
    
    fiveconfigure.debug_mode = True
    zcml.load_config('configure.zcml', Variation.SNP250k)
    fiveconfigure.debug_mode = False
    
    # We need to tell the testing framework that these products
    # should be available. This can't happen until after we have loaded
    # the ZCML.
    
    ztc.installPackage('Variation.SNP250k')
    
# The order here is important: We first call the (deferred) function which
# installs the products we need for the Optilux package. Then, we let 
# PloneTestCase set up this product on installation.

setup_variation_snp250k()
ptc.setupPloneSite(products=['Variation.SNP250k'])

@onsetup
def setup_optilux_cinemacontent():
    """Set up the additional products required for the Optilux Cinema Content.
    
    The @onsetup decorator causes the execution of this body to be deferred
    until the setup of the Plone site testing layer.
    """
    
    # Load the ZCML configuration for the optilux.policy package.
    # This includes the other products below as well.
    
    fiveconfigure.debug_mode = True
    import optilux.cinemacontent
    zcml.load_config('configure.zcml', optilux.cinemacontent)
    fiveconfigure.debug_mode = False
    
    # We need to tell the testing framework that these products
    # should be available. This can't happen until after we have loaded
    # the ZCML.
    
    ztc.installPackage('optilux.cinemacontent')
    
# The order here is important: We first call the (deferred) function which
# installs the products we need for the Optilux package. Then, we let 
# PloneTestCase set up this product on installation.

#04/18/08 no longer need optilux.cinemacontent to test
#setup_optilux_cinemacontent()
#ptc.setupPloneSite(products=['optilux.cinemacontent'])

def test_suite():
	return unittest.TestSuite([

		# Unit tests
		#doctestunit.DocFileSuite(
		#	'README.txt', package='Variation.SNP250k',
		#	setUp=testing.setUp, tearDown=testing.tearDown),

		#doctestunit.DocTestSuite(
		#	module='Variation.SNP250k.mymodule',
		#	setUp=testing.setUp, tearDown=testing.tearDown),


		# Integration tests that use PloneTestCase
		ztc.FunctionalDocFileSuite(
			'README.txt', package='Variation.SNP250k',
			test_class=TestCase,
			optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS),
		ztc.FunctionalDocFileSuite(
			'database.txt', package='Variation.SNP250k',
			test_class=TestCase),
		#ztc.FunctionalDocFileSuite(
		#	'browser.txt', package='Variation.SNP250k',
		#	test_class=TestCase),

		])

if __name__ == '__main__':
	unittest.main(defaultTest='test_suite')
