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

from Variation.SNP250k.interfaces import IPhenotype, IPhenotypeLocator
from Variation.SNP250k import VariationMessageFactory as _
from zope.lifecycleevent import ObjectCreatedEvent, ObjectModifiedEvent
#zope.app.event.objectevent.ObjectModifiedEvent is the old one

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
		

class EditPhenotypeForm(formbase.EditForm):
	form_fields = form.FormFields(IPhenotype).omit('short_name_ls', 'method_description_ls', 'data_matrix')
	#result_template = pagetemplatefile.ZopeTwoPageTemplateFile('search-results.pt')
		
		# a hack to make the content tab work
		#self.template.getId = lambda: 'edit'

	@form.action(_(u"Save"), name='save')
	def action_save(self, action, data):
		"""
		2008-04-22 one question remains. how to use the portal factory tool. need to move the object out of the portal_factory
		"""
		if form.applyChanges(
			self.context, self.form_fields, data, self.adapters):
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
				new_context.method_id_ls = data['method_id_ls']
			else:
				new_context = context
				new_context.setId(str(newId))
			
			locator = getUtility(IPhenotypeLocator)
			new_context.short_name_ls, new_context.method_description_ls = locator.get_short_name_description_ls(new_context.method_id_ls)
			
			zope.event.notify(
				ObjectModifiedEvent(new_context)
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
			
			self.request.response.redirect(new_context.absolute_url())
		else:
			#use below to decide if this is first time being created by portal_factory or followup editing
			context = aq_inner(self.context)
			grandpa = aq_parent(aq_parent(context))
			if grandpa.id=='portal_factory':
				context = aq_parent(aq_parent(aq_parent(aq_inner(self.context))))
			
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
		#use below to decide if this is first time being created by portal_factory or followup editing
		context = aq_inner(self.context)
		grandpa = aq_parent(aq_parent(context))
		if grandpa.id=='portal_factory':
			context = aq_parent(aq_parent(aq_parent(aq_inner(self.context))))
			
		confirm = u"Phenotype Unchanged."
		IStatusMessage(self.request).addStatusMessage(confirm, type='info')
		self.request.response.redirect(context.absolute_url())
		return ''

def fill_phenotype_content(obj, event):
	"""
	2008-04-23 EditPhenotypeForm already does this. Omits them.
	"""
	pass
	#locator = getUtility(IPhenotypeLocator)
	#obj.short_name_ls, obj.method_description_ls = locator.get_short_name_description_ls(obj.method_id_ls)

class AddPhenotypeForm(formbase.AddForm):
	form_fields = form.FormFields(IPhenotype).omit('short_name_ls', 'method_description_ls', 'data_matrix')
	
	
	@form.action(_("Save"), name='save', condition=form.haveInputWidgets)
	def handle_add(self, action, data):
		if form.applyChanges(
			self.context, self.form_fields, data, self.adapters):
			method_id_ls = data['method_id_ls']
			context = aq_inner(self.context)
			newId = data['title'].replace(' ','-')
			newId = newId.replace('/', '-')
			#newId = context.generateNewId()
			obj = Phenotype(newId)
			context.add(obj)
			#context.invokeFactory(id=newId, type_name='Phenotype')
			new_context = getattr(context, newId)
			
			#import sys
			#sys.stderr.write("id: %s, title: %s, newId: %s.\n"%(context.id, context.title, newId))
			#sys.stderr.write("type(newId): %s, dir(newId): %s.\n"%(repr(type(newId)), repr(dir(newId)) ))
			
			#new_context.setId(str(newId))	#type of newId is unicode. and str is required for catalog
			
			#phenotype.generateNewId()
			#new_context = context.portal_factory.doCreate(context, id)
			locator = getUtility(IPhenotypeLocator)
			new_context.short_name_ls, new_context.method_description_ls = locator.get_short_name_description_ls(method_id_ls)
			zope.event.notify(ObjectCreatedEvent(obj))
			self._finished_add = True	#a flag to decide whether the object is added or not
			return obj	#AddForm.render() will handle the redirect
		else:
			confirm = u"No Phenotype Added."
			IStatusMessage(self.request).addStatusMessage(confirm, type='info')
			return None
	
	 
	@form.action(_(u"Cancel"), name='cancel')
	def action_cancel(self, action, data):
		"""
		Cancel the phenotype checkout
		"""
		parent = aq_parent(aq_inner(self.context))
		confirm = u"Phenotype Checkout cancelled."
		IStatusMessage(self.request).addStatusMessage(confirm, type='info')
		self.request.response.redirect(parent.absolute_url())
		return ''
