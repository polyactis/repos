from helloworld.tests import *

class TestSnpregionplotController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='SNPRegionPlot'))
        # Test response...
