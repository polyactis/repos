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

from Variation.SNP250k.interfaces import IQCOnDirectory
from Variation.SNP250k.config import PROJECTNAME

from Variation.SNP250k import VariationMessageFactory as _

class QCOnDirectory(atapi.BaseContent):
	"""Describe a QCOnDirectory
	"""
	implements(IQCOnDirectory)
	
	#04/19/08 these three attributes are used by Archetypes
	#portal_type = "Phenotype"
	#_at_rename_after_creation = True
	#schema = PhenotypeSchema
	
	title = fieldproperty.FieldProperty(IQCOnDirectory['title'])
	description = fieldproperty.FieldProperty(IQCOnDirectory['description'])
	QC_method_id = fieldproperty.FieldProperty(IQCOnDirectory['QC_method_id'])
	input_dir = fieldproperty.FieldProperty(IQCOnDirectory['input_dir'])
	short_name = fieldproperty.FieldProperty(IQCOnDirectory['short_name'])
	row_id2NA_mismatch_rate = fieldproperty.FieldProperty(IQCOnDirectory['row_id2NA_mismatch_rate'])
	

atapi.registerType(QCOnDirectory, PROJECTNAME)
