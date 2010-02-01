from helloworld.tests import *

class TestAssociationoverlapController(TestController):

    def test_index(self):
        response = self.app.get(url(controller='associationoverlap', action='index'))
        # Test response...
