from helloworld.tests import *

class TestScorerankhistogramController(TestController):

    def test_index(self):
        response = self.app.get(url_for(controller='ScoreRankHistogram'))
        # Test response...
