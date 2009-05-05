from helloworld.tests import *

class TestWebhelpController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='WebHelp', action='index'))
        # Test response...
