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


from Acquisition import aq_inner

from Products.Five.browser import BrowserView
from Products.Five.browser import pagetemplatefile
from Products.Five.formlib import formbase

from Products.statusmessages.interfaces import IStatusMessage

from plone.memoize.instance import memoize

from Variation.SNP250k.interfaces import IPhenotype, IPhenotypeLocator
from Variation.SNP250k import VariationMessageFactory as _

class PhenotypeView(BrowserView):
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
	
	__call__ = pagetemplatefile.ViewPageTemplateFile('phenotype.pt')
	
	@memoize
	def short_name_ls(self):
		context = aq_inner(self.context)
		ls_to_return = [short_name for short_name in context.short_name_ls]
		return ls_to_return
		

class CheckoutPhenotypeForm(formbase.EditForm):
	form_fields = form.FormFields(IPhenotype).omit('short_name_ls', 'method_description_ls', 'data_matrix')
	#result_template = pagetemplatefile.ZopeTwoPageTemplateFile('search-results.pt')
	def __init__(self, *args, **kwargs):
		formbase.EditForm.__init__(self, *args, **kwargs)
		
		# a hack to make the content tab work
		#self.template.getId = lambda: 'edit'

	@form.action(_(u"save"))
	def action_save(self, action, data):
		if form.applyChanges(
			self.context, self.form_fields, data, self.adapters):
			method_id_ls = data['method_id_ls']
			#method_id = method_id_short_name.split(' ')[0]
			#method_id = int(method_id)
			context = aq_inner(self.context)
			newId = data['title'].replace(' ','-')
			newId = newId.replace('/', '-')
			#newId = context.generateNewId()
			#context.invokeFactory(id=newId, type_name='Phenotype')
			#new_context = getattr(context, newId)
			
			#import sys
			#sys.stderr.write("id: %s, title: %s, newId: %s.\n"%(context.id, context.title, newId))
			#sys.stderr.write("type(newId): %s, dir(newId): %s.\n"%(repr(type(newId)), repr(dir(newId)) ))
			
			context.setId(str(newId))	#type of newId is unicode. and str is required for catalog
			#phenotype.generateNewId()
			#new_context = context.portal_factory.doCreate(context, id)
			locator = getUtility(IPhenotypeLocator)
			context.short_name_ls, context.method_description_ls = locator.get_short_name_description_ls(method_id_ls)
			
			zope.event.notify(
				zope.app.event.objectevent.ObjectModifiedEvent(self.context)
				)
			# TODO: Needs locale support. See also Five.form.EditView.
			self.status = _(
				"Phenotype Updated on ${date_time}", 
				mapping={'date_time': str(datetime.utcnow())}
				)
			#confirm = u"Phenotype Checked Out."
			#IStatusMessage(self.request).addStatusMessage(confirm, type='info')
			
			#04/18/08 old way of handling from the (The Definitive Guide to Plone)
			#self.state.set(context=new_context, portal_status_message="Phenotype Checked out.")
			
			self.request.response.redirect(context.absolute_url())
		else:
			context = aq_inner(self.context)
			confirm = u"Phenotype Unchanged."
			IStatusMessage(self.request).addStatusMessage(confirm, type='info')
			self.request.response.redirect(context.absolute_url())
			
			#self.status = _('No changes')

	"""
	def update(self):
		if form.applyChanges(
			self.context, self.form_fields, data, self.adapters):
			method_id_ls = data['method_id_ls']
			#method_id = method_id_short_name.split(' ')[0]
			#method_id = int(method_id)
			context = aq_inner(self.context)
			newId = data['title'].replace(' ','-')
			newId = id.replace('/', '-')
			#newId = context.generateNewId()
			context.invokeFactory(newId, type_name='Phenotype')
			new_context = getattr(context, newId)
			
			#phenotype.generateNewId()
			#new_context = context.portal_factory.doCreate(context, id)
			locator = getUtility(IPhenotypeLocator)
			new_context.short_name_ls, new_context.method_description_ls = locator.get_short_name_description_ls(method_id_ls)
			
			zope.event.notify(
				zope.app.event.objectevent.ObjectModifiedEvent(self.context)
				)
			# TODO: Needs locale support. See also Five.form.EditView.
			self.status = _(
				"Updated on ${date_time}", 
				mapping={'date_time': str(datetime.utcnow())}
				)
			
		else:
			self.status = _('No changes')
	"""
	
	@form.action(_(u"Cancel"), name='cancel')
	def action_cancel(self, action, data):
		"""
		Cancel the phenotype checkout
		"""
		context = aq_inner(self.context)
		confirm = u"Phenotype Checkout cancelled."
		IStatusMessage(self.request).addStatusMessage(confirm, type='info')
		self.request.response.redirect(context.absolute_url())
		return ''

def fill_phenotype_content(obj, event):
	locator = getUtility(IPhenotypeLocator)
	obj.short_name_ls, obj.method_description_ls = locator.get_short_name_description_ls(obj.method_id_ls)
