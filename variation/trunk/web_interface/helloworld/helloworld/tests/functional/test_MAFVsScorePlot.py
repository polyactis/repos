from helloworld.tests import *

class TestMafvsscoreplotController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='MAFVsScorePlot'))
        # Test response...
