"""Define a browser view for the Film content type. In the FTI 
configured in profiles/default/types/*.xml, this is being set as the default
view of that content type.
"""

from datetime import datetime, timedelta

from zope.component import getUtility
from zope.formlib import form
from zope.schema import vocabulary
import zope
#from zope.app.form import browser


from Acquisition import aq_inner, aq_parent

from Products.Five.browser import BrowserView
from Products.Five.browser import pagetemplatefile
from Products.Five.formlib import formbase

from Products.statusmessages.interfaces import IStatusMessage

from plone.memoize.instance import memoize

from Variation.SNP250k.interfaces import IQCOnDirectory, IStockDatabaseSettings
from Variation.SNP250k import VariationMessageFactory as _
from zope.lifecycleevent import ObjectCreatedEvent, ObjectModifiedEvent
#zope.app.event.objectevent.ObjectModifiedEvent is the old one

class QC_unit(object):
	def __init__(self, **keywords):
		argument_default_dict = {('array_id', 1, ): None,\
								('ecotypeid', 1, ): None,\
								('NA_rate', 1, ): None,\
								('mismatch_rate', 1, ): None,\
								('no_of_NAs', 1, ): None,\
								('no_of_totals', 1, ): None,\
								('no_of_mismatches', 1, ): None,\
								('no_of_non_NA_pairs', 1, ): None}
		
		from pymodule import process_function_arguments					
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self, howto_deal_with_required_none=2)


class QCOnDirectoryView(BrowserView):
	"""
	Default view of a film
	
	#@memoize
	def cinemas(self, days=14):
		context = aq_inner(self.context)
		locator = getUtility(IScreeningLocator)
		
		from_date = datetime.now()
		to_date = from_date + timedelta(days)
		return locator.cinemas_for_film(context, from_date, to_date)	
	
	def banner_tag(self):
		context = aq_inner(self.context)
		banner_provider = IBannerProvider(context)
		return banner_provider.tag
	"""
	
	__call__ = pagetemplatefile.ViewPageTemplateFile('qcondirectory.pt')
	
	@memoize
	def iter_NA_mismatch_rate(self):
		context = aq_inner(self.context)
		row_id_ls = context.row_id2NA_mismatch_rate.keys()
		row_id_ls.sort()	#try to keep them in call_info_id order
		ls_to_return = []
		for row_id in row_id_ls:
			NA_mismatch_ls = row_id2NA_mismatch_rate[row_id]
			qc_unit = QC_unit(array_id=row_id[0], ecotypeid=row_id[1], NA_rate= NA_mismatch_ls[0], mismatch_rate=NA_mismatch_ls[1],\
				no_of_NAs=NA_mismatch_ls[2], no_of_totals=NA_mismatch_ls[3], no_of_mismatches=NA_mismatch_ls[4], no_of_non_NA_pairs=NA_mismatch_ls[5])
			ls_to_return.append(qc_unit)
		return ls_to_return
		

class EditQCOnDirectoryForm(formbase.EditForm):
	form_fields = form.FormFields(IQCOnDirectory).omit('short_name', 'row_id2NA_mismatch_rate')
	
	@form.action(_(u"Save"), name='save')
	def action_save(self, action, data):
		"""
		2008-04-22 one question remains. how to use the portal factory tool. need to move the object out of the portal_factory
		"""
		if form.applyChanges(
			self.context, self.form_fields, data, self.adapters):
			context = aq_inner(self.context)
			newId = data['title'].replace(' ','-')
			newId = newId.replace('/', '-')
			#newId = context.generateNewId()
			#context.invokeFactory(id=newId, type_name='Phenotype')
			#new_context = getattr(context, newId)
			
			#import sys
			#sys.stderr.write("id: %s, title: %s, newId: %s.\n"%(context.id, context.title, newId))
			#sys.stderr.write("type(newId): %s, dir(newId): %s.\n"%(repr(type(newId)), repr(dir(newId)) ))
			
			#use below to decide if this is first time being created by portal_factory or followup editing
			grandpa = aq_parent(aq_parent(context))
			if grandpa.id=='portal_factory':
				#phenotype.generateNewId()
				
				#2008-04-23. below is probably a simpler line. don't know where i got it.
				#new_context = context.portal_factory.doCreate(context, id)
				
				#2008-04-23. got the clue from FactoryTool.doCreate() of Products/CMFPlone/FactoryTool.py
				type_name = aq_parent(context).id  # get the ID of the TempFolder
				folder = aq_parent(aq_parent(aq_parent(context)))
				folder.invokeFactory(id=newId, type_name=type_name)
				
				new_context = getattr(folder, newId)
				new_context.setId(str(newId))	#type of newId is unicode. and str is required for catalog

				#manually set all these properties
				new_context.title = data['title']
				new_context.description = data['description']
				new_context.QC_method_id = data['QC_method_id']
				new_context.input_dir = data['input_dir']
			else:
				new_context = context
				new_context.setId(str(newId))
			
			new_context.short_name = str(new_context.QC_method_id)
			
			settings = getUtility(IStockDatabaseSettings, name='variation.stockdatabasesettings')
			from variation.src.QC_250k import QC_250k
			QC_method_id2cmp_data_filename = {1: '/usr/local/home_ubuntu/crocea/script/variation/data/2010/data_2010_x_250k_y0001.tsv',
				2: '/usr/local/home_ubuntu/crocea/script/variation/data/perlegen/data_perlegen_ecotype_id_x_250k_y0101.tsv',
				3: '/usr/local/home_ubuntu/crocea/script/variation/stock20080403/data_y10001101.tsv'}
			from pymodule import process_options, generate_program_doc
			argv_list = ['QC_250k.py', '-i', QC_method_id2cmp_data_filename[int(new_context.QC_method_id)], \
				'-n', new_context.input_dir, '-m', new_context.QC_method_id, '-u', settings.username, '-p', settings.password]
			opts_dict = process_options(argv_list, QC_250k.option_default_dict, error_doc=generate_program_doc(argv_list[0], QC_250k.option_default_dict)+QC_250k.__doc__)
			
			instance = QC_250k(**opts_dict)
			_row_id2NA_mismatch_rate = instance.plone_run()
			for row_id, NA_mismatch_ls in _row_id2NA_mismatch_rate.iteritems():
				new_context.row_id2NA_mismatch_rate[row_id] = NA_mismatch_ls
			
			zope.event.notify(
				ObjectModifiedEvent(new_context)
				)
			# TODO: Needs locale support. See also Five.form.EditView.
			self.status = _(
				"Updated on ${date_time}", 
				mapping={'date_time': str(datetime.utcnow())}
				)
			#confirm = u"Phenotype Checked Out."
			#IStatusMessage(self.request).addStatusMessage(confirm, type='info')
			
			#04/18/08 old way of handling from the (The Definitive Guide to Plone)
			#self.state.set(context=new_context, portal_status_message="Phenotype Checked out.")
			
			self.request.response.redirect(new_context.absolute_url())
		else:
			#use below to decide if this is first time being created by portal_factory or followup editing
			context = aq_inner(self.context)
			grandpa = aq_parent(aq_parent(context))
			if grandpa.id=='portal_factory':
				context = aq_parent(aq_parent(aq_parent(aq_inner(self.context))))
			
			confirm = u"Cancelled."
			IStatusMessage(self.request).addStatusMessage(confirm, type='info')
			self.request.response.redirect(context.absolute_url())
			
			#self.status = _('No changes')
	
	@form.action(_(u"Cancel"), name='cancel')
	def action_cancel(self, action, data):
		"""
		Cancel the phenotype checkout
		"""
		#use below to decide if this is first time being created by portal_factory or followup editing
		context = aq_inner(self.context)
		grandpa = aq_parent(aq_parent(context))
		if grandpa.id=='portal_factory':
			context = aq_parent(aq_parent(aq_parent(aq_inner(self.context))))
			
		confirm = u"Cancellled."
		IStatusMessage(self.request).addStatusMessage(confirm, type='info')
		self.request.response.redirect(context.absolute_url())
		return ''
