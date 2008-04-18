import unittest

from zope.testing import doctestunit
from zope.component import testing
from Testing import ZopeTestCase as ztc

from Products.Five import zcml
from Products.Five import fiveconfigure
from Products.PloneTestCase import PloneTestCase as ptc
from Products.PloneTestCase.layer import PloneSite
ptc.setupPloneSite()

import Variation.SNP250k

class TestCase(ptc.PloneTestCase):
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
			test_class=TestCase),
		ztc.ZopeDocFileSuite(
			'database.txt', package='Variation.SNP250k',
			test_class=TestCase),
		#ztc.FunctionalDocFileSuite(
		#	'browser.txt', package='Variation.SNP250k',
		#	test_class=TestCase),

		])

if __name__ == '__main__':
	unittest.main(defaultTest='test_suite')
