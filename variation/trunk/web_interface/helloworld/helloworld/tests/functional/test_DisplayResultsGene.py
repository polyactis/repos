from helloworld.tests import *

class TestDisplayresultsgeneController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='DisplayResultsGene'))
        # Test response...
