from helloworld.tests import *

class TestGreetingController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='greeting'))
        # Test response...
