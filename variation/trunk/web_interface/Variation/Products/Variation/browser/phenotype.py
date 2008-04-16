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

#from plone.memoize.instance import memoize

from Products.Variation.interfaces import IPhenotype, IPhenotypeLocator
#from optilux.cinemacontent.interfaces import IScreeningLocator
from Products.Variation import VariationMessageFactory as _

class PhenotypeView(BrowserView):
	"""Default view of a film
	"""
	
	__call__ = pagetemplatefile.ViewPageTemplateFile('phenotype.pt')
	"""
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

def YesNoWidget(field, request):
	locator = getUtility(IPhenotypeLocator)
	v = locator.get_phenotype_method_id_ls()
	vocabulary = vocabulary.SimpleVocabulary.fromItems(v)
	#vocabulary = vocabulary.SimpleVocabulary.fromItems(((true, True), (false, False)))
	return zope.app.form.browser.MultiSelectWidget(field, vocabulary, request)
   
   
class CheckoutPhenotypeForm(formbase.PageForm):
	form_fields = form.FormFields(IPhenotype).omit('short_name_ls', 'method_description_ls', 'data_matrix')
	#result_template = pagetemplatefile.ZopeTwoPageTemplateFile('search-results.pt')
	
	@form.action(_(u"save"))
	def action_save(self, action, data):
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
		phenotype.short_name_ls, phenotype.method_description_ls = locator.get_short_name_description_ls(method_id_ls)
		"""
		catalog = cmfutils.getToolByName(self.context, 'portal_catalog')
		
		kwargs = {}
		if data['text']:
			kwargs['SearchableText'] = data['text']
		if data['description']:
			kwargs['description'] = data['description']
		
		self.search_results = catalog(**kwargs)
		self.search_results_count = len(self.search_results)
		return self.result_template()
		"""
		return self.state.set(context=new_context, portal_status_message="Phenotype Checked out.")
	
	@form.action(_(u"cancel"))
	def action_cancel(self, action, data):
		"""
		Cancel the phenotype checkout
		"""
		context = aq_inner(self.context)
		confirm = u"Phenotype Checkout cancelled."
		IStatusMessage(self.request).addStatusMessage(confirm, type='info')
		self.request.response.redirect(context.absolute_url())
		return self.state.set(context=context, portal_status_message="Phenotype Checkout Cancelled.")
