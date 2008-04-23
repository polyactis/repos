"""Definition of the Film content type. See cinemafolder.py for more
explanation on the statements below.
"""

from zope.interface import implements, directlyProvides
from zope.component import adapts

from zope.schema.interfaces import IVocabularyFactory
from zope.schema.vocabulary import SimpleVocabulary
from zope.schema import fieldproperty

from Acquisition import aq_inner
from DateTime import DateTime

from Products.Archetypes.interfaces import IObjectPostValidation

from Products.Archetypes import atapi
from Products.validation import V_REQUIRED

from Products.ATContentTypes.content import base
from Products.ATContentTypes.content import schemata
from Products.ATContentTypes.content.schemata import finalizeATCTSchema

from Products.CMFCore.utils import getToolByName

from Variation.SNP250k.interfaces import IPhenotype
from Variation.SNP250k.config import PROJECTNAME

#from Products.Variation import VariationMessageFactory as _

"""
PhenotypeSchema = schemata.ATContentTypeSchema.copy() + atapi.Schema((

	atapi.IntegerField('methodID',
		required=False,
		searchable=True,
		storage=atapi.AnnotationStorage(),
		widget=atapi.IntegerWidget(label=u"Phenotype Method ID",
								  description=u"", visible={'view':'visible', 'edit':'hidden'})
		),
	
	
	atapi.StringField('methodIDShortName',
		required=True,
		searchable=True,
		storage=atapi.AnnotationStorage(),
		enforceVocabulary=True,
		vocabulary='listMethodIDShortName',
		widget=atapi.SelectionWidget(label=u"Phenotype Method ID/Short name",
								  description=(u""))
		),
	
	atapi.TextField('shortName',
		required=False,
		searchable=True,
		storage=atapi.AnnotationStorage(),
		widget=atapi.StringWidget(label=u"Short Name",
								description=u"", visible={'view':'visible', 'edit':'hidden'})
		),

	))

PhenotypeSchema['title'].storage = atapi.AnnotationStorage()
PhenotypeSchema['title'].widget.label = u"Method Title"
PhenotypeSchema['title'].widget.description = u""

FilmSchema['description'].storage = atapi.AnnotationStorage()
FilmSchema['description'].widget.label = _(u"Description")
FilmSchema['description'].widget.description = _(u"")
"""

#finalizeATCTSchema(PhenotypeSchema, folderish=False)

class Phenotype(atapi.BaseContent):
	"""Describe a phenotype
	"""
	implements(IPhenotype)
	
	#04/19/08 these three attributes are used by Archetypes
	portal_type = "Phenotype"
	#_at_rename_after_creation = True
	#schema = PhenotypeSchema
	
	title = fieldproperty.FieldProperty(IPhenotype['title'])
	description = fieldproperty.FieldProperty(IPhenotype['description'])
	
	method_id_ls = fieldproperty.FieldProperty(IPhenotype['method_id_ls'])
	short_name_ls = fieldproperty.FieldProperty(IPhenotype['short_name_ls'])
	method_description_ls = fieldproperty.FieldProperty(IPhenotype['method_description_ls'])
	data_matrix = fieldproperty.FieldProperty(IPhenotype['data_matrix'])
	
	#def listMethodIDShortName(self):
	#	return ['1/LD-V', '2/LD+V', '3/ShD-V'] + map(repr, range(4,46))

atapi.registerType(Phenotype, PROJECTNAME)

#04/18/08 not in Archetypes anymore. need this to reindex catalog
def catalog_content(obj, event):
	obj.reindexObject()

# This simple adapter uses Archetypes' ImageField to extract an HTML tag
# for the banner image. This is used in the promotions portlet to avoid
# having a hard dependency on the AT ImageField implementation.

# Note that we adapt a class, not an interface. This means that we will only
# match adapter lookups for this class (or a subclass), which is correct in
# this case, because we are relying on internal implementation details.

"""
class BannerProvider(object):
	implements(IBannerProvider)
	adapts(Film)
	
	def __init__(self, context):
		self.context = context
	
	@property
	def tag(self):
		return self.context.getField('image').tag(self.context, scale='thumb')
"""	 
# This is a subscription adapter which is used to validate the film object.
# It will be called after the normal schema validation.
"""
class ValidateFilmCodeUniqueness(object):
	#
	implements(IObjectPostValidation)
	adapts(IFilm)
	
	field_name = 'filmCode'
	
	def __init__(self, context):
		self.context = context
	
	def __call__(self, request):
		value = request.form.get(self.field_name, request.get(self.field_name, None))
		if value is not None:
			catalog = getToolByName(self.context, 'portal_catalog')
			results = catalog(film_code = value)
			if len(results) == 0:
				return None
			elif len(results) == 1 and results[0].UID == self.context.UID():
				return None
			else:
				return {self.field_name : _(u"The film code is already in use")}
		
		# Returning None means no error
		return None
"""		
# A vocabulary factory to return currently valid films. The factory itself 
# is registered as a named utility in configure.zcml. This is referenced
# from cinema.py, in the highlightedFilms reference field.
