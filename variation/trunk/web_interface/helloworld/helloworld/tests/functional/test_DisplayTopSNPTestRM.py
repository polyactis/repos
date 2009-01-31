from helloworld.tests import *

class TestDisplaytopsnptestrmController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='DisplayTopSNPTestRM'))
        # Test response...
