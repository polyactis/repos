from helloworld.tests import *

class TestAccessionserviceController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='AccessionService', action='index'))
        # Test response...
