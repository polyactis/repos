# Autocomplete component for Pyjamas
# Ported by Willie Gollino from Autocomplete component for GWT - Originally by Oliver Albers http://gwt.components.googlepages.com/

from pyjamas.ui.TextBox import TextBox
from pyjamas.ui.PopupPanel import PopupPanel
from pyjamas.ui.ListBox import ListBox
from pyjamas.ui.KeyboardListener import KeyboardListener
from pyjamas.ui.RootPanel import RootPanel

from AutoComplete import AutoCompleteTextBox
from pyjamas.JSONService import JSONProxy

class AutoCompleteByURLTextBox(AutoCompleteTextBox):
	def __init__(self, url):
		self.choicesPopup = PopupPanel(True)
		self.choices = ListBox()
		self.items = AutoCompleteURL(url)
		self.popupAdded = False
		self.visible = False
		
		TextBox.__init__(self)
		self.addKeyboardListener(self)
		self.choices.addChangeListener(self)
		self.setStyleName("AutoCompleteTextBox")
		
		self.choicesPopup.add(self.choices)
		self.choicesPopup.addStyleName("AutoCompleteChoices")
		
		self.choices.setStyleName("list")
		self.url = url
	
	"""
	def setCompletionItems(self, items):
		if not items.getCompletionItems:
			items = AutoCompletionURL(self.url)
		
		self.items = items

	def getCompletionItems(self):
		return self.items
	"""
	
	def onKeyUp(self, arg0, arg1, arg2):
		if arg1 == KeyboardListener.KEY_DOWN:
			selectedIndex = self.choices.getSelectedIndex()
			selectedIndex += 1
			if selectedIndex > self.choices.getItemCount():
				selectedIndex = 0
		
			self.choices.setSelectedIndex(selectedIndex)		   
			return

		if arg1 == KeyboardListener.KEY_UP:
			selectedIndex = self.choices.getSelectedIndex()
			selectedIndex -= 1
			if selectedIndex < 0:
				selectedIndex = self.choices.getItemCount()
			self.choices.setSelectedIndex(selectedIndex)
			return

		if arg1 == KeyboardListener.KEY_ENTER:
			if self.visible:
				self.complete()	  
			return

		if arg1 == KeyboardListener.KEY_ESCAPE:
			self.choices.clear()
			self.choicesPopup.hide()
			self.visible = False
			return

		text = self.getText()
		matches = []
		if len(text) > 0:
			matches = self.items.getCompletionItems(text)

		if len(matches) > 0:
			self.choices.clear()

			for i in range(len(matches)):
				self.choices.addItem(matches[i])
				
			if len(matches) == 1 and matches[0] == text:
				self.choicesPopup.hide()
			else:
				self.choices.setSelectedIndex(0)
				self.choices.setVisibleItemCount(len(matches) + 1)
					
				if not self.popupAdded:
					RootPanel().add(self.choicesPopup)
					self.popupAdded = True

				self.choicesPopup.show()
				self.visible = True
				self.choicesPopup.setPopupPosition(self.getAbsoluteLeft(), self.getAbsoluteTop() + self.getOffsetHeight())
				self.choices.setWidth(self.getOffsetWidth() + "px")
		else:
			self.visible = False
			self.choicesPopup.hide()
	
class CompletionService(JSONProxy):
	"""
	2009-4-1
		JSON service module
	"""
	def __init__(self, url):
		JSONProxy.__init__(self, url, ["fetch"])
		
class AutoCompleteURL:
	"""
	2009-4-1
		modelled after class SimpleAutoCompletionItems from /usr/share/pyjamas/addons/AutoComplete.py and /usr/share/pyjamas/examples/jsonrpc/JSONRPCExample.py
		
	"""
	def __init__(self, url):
		self.url = url
		self.matches = []
		self.remote_py = CompletionService(url)
	
	def onRemoteResponse(self, response, request_info):
		self.matches = response

	def onRemoteError(self, code, message, request_info):
		self.matches = []
		#self.status.setText("Server Error or Invalid Response: ERROR " + code + " - " + message)
	
	def getCompletionItems(self, match):
		id = self.remote_py.fetch(match, self)
		if id<0:
			self.matches = []
			#self.status.setText(self.TEXT_ERROR)
		return self.matches
