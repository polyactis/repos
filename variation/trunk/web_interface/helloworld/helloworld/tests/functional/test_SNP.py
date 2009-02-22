from helloworld.tests import *

class TestSnpController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='SNP', action='index'))
        # Test response...
