from helloworld.tests import *

class TestHaplotypeviewController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='HaplotypeView', action='index'))
        # Test response...
