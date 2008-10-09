from helloworld.tests import *

class TestCallinfoController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='CallInfo'))
        # Test response...
