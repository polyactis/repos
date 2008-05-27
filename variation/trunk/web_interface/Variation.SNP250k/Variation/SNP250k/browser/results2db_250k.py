from zope.formlib import form
from zope.app.form.browser import FileWidget, ChoiceDisplayWidget, RadioWidget, ChoiceCollectionInputWidget, SelectWidget
from zope.app.form import CustomWidgetFactory

from Products.Five.formlib import formbase

from plone.fieldsets.form import FieldsetsEditForm
from plone.app.form.validators import null_validator
from Products.CMFCore.utils import getToolByName
from zope.component import getMultiAdapter
from zope.component import getUtility
from zope.component import adapter
from zope.interface import implementer
from Products.statusmessages.interfaces import IStatusMessage

from collective.lead.interfaces import IDatabase

from Variation.SNP250k.interfaces import IResults2DB_250k, IVariationFolder, IStockDatabaseSettings, IDBLocator
from Variation.SNP250k import VariationMessageFactory as _
from variation.src.Results2DB_250k import Results2DB_250k

from persistent import Persistent
from zope.interface import implements

import gc

@implementer(IResults2DB_250k)
@adapter(IVariationFolder)
def results2db_250k_adapter(context):
	return getUtility(IResults2DB_250k, name='variation.results2db_250k')

class Results2DB_250kFormContent(Persistent):
	implements(IResults2DB_250k)
	
	short_name = None
	username = None
	password = None
	phenotype_method_id = None
	call_method_id = None
	data_description = None
	results_method_type_id = None
	results_method_type_short_name = None
	method_description = None
	comment = None
	input_fname = None
	commit_type = None

class MyRadioWidget(RadioWidget):
	def __init__(self, field, request):
		super(MyRadioWidget, self).__init__(field, field.vocabulary, request)


class Results2DB_250kForm(FieldsetsEditForm):
	"""
	FieldsetsEditForm is same as formbase.EditForm except it could identify multiple sets by passing multiple interfaces to form.FormFields()
	"""
	form_fields = form.FormFields(IResults2DB_250k)
	form_fields['input_fname'].custom_widget = FileWidget
	form_fields['commit_type'].custom_widget = MyRadioWidget
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
		"""
		2008-05-25
			add garbage collection
		2008-05-25
			check commit_type
		"""
		if data['results_method_type_id'] is None and not data['results_method_type_short_name']:
			IStatusMessage(self.request).addStatusMessage(_("Either results_method_type_id or results_method_type_short_name has to be something."), type='error')
		elif data['commit_type']==2:	#don't care whether contents change or not
			start_time = self.context.ZopeTime().timeTime()
			db = getUtility(IDatabase, name='variation.stockdatabase')
			if getattr(db, 'transaction', None) is not None:
				db.transaction.commit()	#commit all previous transactions
				db.session.close()
				db.transaction = None
				Results2DB_250k.reset_marker_pos2snp_id()
				self.status = _("Took %f seconds to submit all previous data."%(self.context.ZopeTime().timeTime()-start_time))
			else:
				IStatusMessage(self.request).addStatusMessage(_("Warning: No data submission as no previous transaction exists."), type='warning')
		elif data['commit_type']==3:	#don't care whether contents change or not
			start_time = self.context.ZopeTime().timeTime()
			db = getUtility(IDatabase, name='variation.stockdatabase')
			if getattr(db, 'transaction', None) is not None:
				db.transaction.rollback()	#commit all previous transactions
				db.session.close()
				db.transaction = None
				Results2DB_250k.reset_marker_pos2snp_id()
				self.status = _("Took %f seconds to rollback all previous transactions."%(self.context.ZopeTime().timeTime()-start_time))
			else:
				IStatusMessage(self.request).addStatusMessage(_("Warning: No rollback as no previous transaction exists."), type='warning')
		elif form.applyChanges(self.context, self.form_fields, data,
							 self.adapters):
			locator = getUtility(IDBLocator)
			if locator.checkIfResultsMethodExist(data['short_name']):
				IStatusMessage(self.request).addStatusMessage(_("Error: Short Name, %s already exists in database"%data['short_name']), type='error')
			else:
				start_time = self.context.ZopeTime().timeTime()
				self._on_save(data)
				del data['input_fname']	#release it from memory
				self.status = _("Took %f seconds to submit the data."%(self.context.ZopeTime().timeTime()-start_time))
		else:
			IStatusMessage(self.request).addStatusMessage(_("Warning: No data submission as nothing in the form changed (Hint: add a space somewhere)."), type='warning')
		
		gc.collect()    #run a full collection, clean up memory
		
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
		#import pdb
		#pdb.set_trace()
		db = getUtility(IDatabase, name='variation.stockdatabase')
		if data['username']:	#not None
			#update username/password first in the settings object
			settings = getUtility(IStockDatabaseSettings, name='variation.stockdatabasesettings')
			settings.username = data['username']
			settings.password = data['password']
			db.invalidate()	#need this to update the database connection object with new settings
			user = data['username']	#if user specified the username, take it as user_id
		else:
			#get the logged in user id from plone membership
			mtool = getToolByName(self.context, 'portal_membership')
			member = mtool.getAuthenticatedMember()
			user = member.getId()
		file_object = self.request.form['form.input_fname']
		#data['input_fname'] is the actual data in the file (plone read it out already). self.request.form['form.input_fname'] takes the original file object.
		comment = ''
		if data['comment']:
			comment = data['comment']
		comment += '. Original Filename = %s'%(file_object.filename)
		Results2DB_250k.plone_run(db, data['short_name'], data['phenotype_method_id'], data['call_method_id'], \
								data['data_description'], data['method_description'], comment, \
								file_object, user, data['results_method_type_id'], \
								results_method_type_short_name=data['results_method_type_short_name'], commit=data['commit_type'])
