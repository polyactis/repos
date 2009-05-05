from helloworld.tests import *

class TestHelpothercontrollersController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='HelpOtherControllers', action='index'))
        # Test response...
