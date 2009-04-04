from pyjamas.ui.RootPanel import RootPanel
from pyjamas.ui.TextArea import TextArea
from pyjamas.ui.Label import Label
from pyjamas.ui.Button import Button
from pyjamas.ui.HTML import HTML
from pyjamas.ui.VerticalPanel import VerticalPanel
from pyjamas.ui.HorizontalPanel import HorizontalPanel
from pyjamas.ui.ListBox import ListBox
from pyjamas.ui.TabPanel import TabPanel
from pyjamas.JSONService import JSONProxy

from pyjamas.ui.decoratorpanel import DecoratedTabPanel, DecoratorPanel
#from pyjamas.ui.decoratorpanel import DecoratorTitledPanel

from pyjamas.History import History
from AccessionByName import AccessionByName
from AccessionRawSQL import AccessionRawSQL

class Accession:
	def onModuleLoad(self):
		self.fTabs = DecoratedTabPanel()
		
		self.fTabs.add(VerticalPanel(), 'Accessions', True)
		self.fTabs.add(AccessionByName(), 'Search By Name', True)
		self.fTabs.add(VerticalPanel(), "Search By Country")
		self.fTabs.add(VerticalPanel(), "Search By Genetic Distance")
		self.fTabs.add(AccessionRawSQL(), "Raw SQL Query")
		
		self.fTabs.selectTab(1)
		#dp = DecoratorTitledPanel("Tabs", "bluetitle", "bluetitleicon",
		#						["bluetop", "bluetop2", "bluemiddle", "bluebottom"])
		#dp.add(self.fTabs)
		#RootPanel().add(panel)
		History().addHistoryListener(self)
		RootPanel().add(self.fTabs)


class EchoServicePython(JSONProxy):
	def __init__(self):
		JSONProxy.__init__(self, "/Accession/json", ["echo", "reverse", "uppercase", "lowercase"])

if __name__ == '__main__':
	Accession().onModuleLoad()
