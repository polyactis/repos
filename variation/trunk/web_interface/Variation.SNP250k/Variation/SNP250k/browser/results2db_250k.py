from zope.formlib import form
from zope.app.form.browser import FileWidget, ChoiceDisplayWidget
from zope.app.form import CustomWidgetFactory

from Products.Five.formlib import formbase

from plone.fieldsets.form import FieldsetsEditForm
from plone.app.form.validators import null_validator
from zope.component import getMultiAdapter
from zope.component import getUtility
from zope.component import adapter
from zope.interface import implementer
from Products.statusmessages.interfaces import IStatusMessage

from collective.lead.interfaces import IDatabase

from Variation.SNP250k.interfaces import IResults2DB_250k, IVariationFolder
from Variation.SNP250k import VariationMessageFactory as _
from variation.src.db import PhenotypeMethod, ResultsMethod
from variation.src.Results2DB_250k import Results2DB_250k

from persistent import Persistent
from zope.interface import implements

@implementer(IResults2DB_250k)
@adapter(IVariationFolder)
def results2db_250k_adapter(context):
	return getUtility(IResults2DB_250k)

class Results2DB_250kFormContent(Persistent):
	implements(IResults2DB_250k)
	title = None
	description = None
	input_fname = None
	phenotype_method_id = None
	short_name = None
	method_description = None
	data_description = None

class Results2DB_250kForm(FieldsetsEditForm):
	"""
	FieldsetsEditForm is same as formbase.EditForm except it could identify multiple sets by passing multiple interfaces to form.FormFields()
	"""
	form_fields = form.FormFields(IResults2DB_250k)
	form_fields['input_fname'].custom_widget = FileWidget
	form_name = _(u"Results2DB_250k")
	label = _(u"Results2DB_250k")
	description = _(u"Submit genome-wide results to database")
	# This trick hides the editable border and tabs in Plone
	#def __call__(self):
	#	#2008-05-22 stop hiding the borders and tabs
	#	self.request.set('disable_border', True)
	#	return super(StockDatabaseSettingsForm, self).__call__()
	
	@form.action(_(u'label_save', default=u'Save'), name=u'save')
	def handle_edit_action(self, action, data):
		if form.applyChanges(self.context, self.form_fields, data,
							 self.adapters):
			self.status = _("Data submitted.")
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
		connection = db.connection
		session = db.session
		transaction = session.create_transaction()
		pm = session.query(PhenotypeMethod).get_by(id=data['phenotype_method_id'])
		rm = ResultsMethod(short_name=data['short_name'], method_description=data['method_description'], data_description=data['data_description'])
		
		#results_method_id = self.submit_results_method(curs, self.results_method_table, self.short_name, self.method_description, self.data_description)
		#Results2DB_250k.submit_results(session, data['input_fname'], rm, pm)
		#session.flush()	#not necessary as no immediate query on the new results after this and commit() would execute this.
		#transaction.commit()
		transaction.rollback()