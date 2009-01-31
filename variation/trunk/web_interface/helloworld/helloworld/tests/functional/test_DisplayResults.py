from helloworld.tests import *

class TestDisplayresultsController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='DisplayResults', action='index'))
        # Test response...
