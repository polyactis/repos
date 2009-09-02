from helloworld.tests import *

class TestUtilsController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='utils', action='index'))
        # Test response...
