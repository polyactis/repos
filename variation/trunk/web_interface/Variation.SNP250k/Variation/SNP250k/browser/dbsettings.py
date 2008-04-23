from zope.formlib import form
from Products.Five.formlib import formbase

from plone.fieldsets.form import FieldsetsEditForm
from plone.app.form.validators import null_validator
from zope.component import getMultiAdapter
from zope.component import getUtility
from zope.component import adapter
from zope.interface import implementer
from Products.statusmessages.interfaces import IStatusMessage

from collective.lead.interfaces import IDatabase

from Variation.SNP250k.interfaces import IStockDatabaseSettings, IVariationFolder
from Variation.SNP250k import VariationMessageFactory as _

#from plone.app.controlpanel.form import ControlPanelForm
"""
2008-04-22 ControlPanelForm has 'save' and 'cancel' action/button setup. only need _on_save() (part of its interface)
but to keep things more clear, use FieldsetsEditForm or formbase.EditForm
"""

@implementer(IStockDatabaseSettings)
@adapter(IVariationFolder)
def stock_database_settings(context):
	return getUtility(IStockDatabaseSettings)

class StockDatabaseSettingsForm(FieldsetsEditForm):
	"""
	FieldsetsEditForm is same as formbase.EditForm except it could identify multiple sets by passing multiple interfaces to form.FormFields()
	"""
	form_fields = form.FormFields(IStockDatabaseSettings)
	form_name = _(u"Database Settings")
	label = _(u"Database Settings")
	description = _(u"Please enter the appropriate connection settings for the database")
	
	# This trick hides the editable border and tabs in Plone
	def __call__(self):
		self.request.set('disable_border', True)
		return super(StockDatabaseSettingsForm, self).__call__()
	
	@form.action(_(u'label_save', default=u'Save'), name=u'save')
	def handle_edit_action(self, action, data):
		if form.applyChanges(self.context, self.form_fields, data,
							 self.adapters):
			self.status = _("Changes saved.")
			self._on_save(data)
		else:
			self.status = _("No changes made.")

	@form.action(_(u'label_cancel', default=u'Cancel'),
				 validator=null_validator,
				 name=u'cancel')
	def handle_cancel_action(self, action, data):
		IStatusMessage(self.request).addStatusMessage(_("Changes canceled."),
													  type="info")
		url = getMultiAdapter((self.context, self.request),
							  name='absolute_url')()
		self.request.response.redirect(url)
		return ''
		
	def _on_save(self, data=None):
		"""This method is called when the form is successfully saved. We use
		it to inform the database engine that the settings have changed.
		"""
		db = getUtility(IDatabase, name='variation.stockdatabase')
		db.invalidate()