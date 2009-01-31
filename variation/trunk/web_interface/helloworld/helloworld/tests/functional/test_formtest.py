from helloworld.tests import *

class TestFormtestController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='formtest'))
        # Test response...
