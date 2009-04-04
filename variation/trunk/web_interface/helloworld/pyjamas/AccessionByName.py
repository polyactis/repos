from pyjamas.ui.HTML import HTML
from pyjamas.ui.VerticalPanel import VerticalPanel
from pyjamas.ui.HorizontalPanel import HorizontalPanel
from pyjamas.ui.FormPanel import FormPanel
from pyjamas.ui.Composite import Composite
from pyjamas.ui.Button import Button
from AutoCompleteByURL import AutoCompleteByURLTextBox

class AccessionByName(Composite):
	def __init__(self):
		self.form = FormPanel()
		self.form.setAction("/Accession/searchByName/")
		
		# Because we're going to add a FileUpload widget, we'll need to set the
		# form to use the POST method, and multipart MIME encoding.
		self.form.setEncoding(FormPanel.ENCODING_MULTIPART)
		self.form.setMethod(FormPanel.METHOD_POST)
		
		# Create a panel to hold all of the form widgets.
		vpanel = VerticalPanel()
		self.form.setWidget(vpanel)
		
		self.colour_input = AutoCompleteByURLTextBox('/Accession/autoComplete')
		
		hpanel = HorizontalPanel()
		hpanel.add(HTML("Enter an ecotype name: "))
		hpanel.add(self.colour_input)

		hpanel.setSpacing(8)
		# Add a 'submit' button.
		hpanel.add(Button("Submit", self))

		vpanel.add(hpanel)
		
		# Add an event handler to the form.
		self.form.addFormHandler(self)
		
		self.setWidget(self.form)
		

	def onShow(self):
		#self.colour_input.setFocus(True)
		return