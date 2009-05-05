from helloworld.tests import *

class TestPhenotypeController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='Phenotype', action='index'))
        # Test response...
