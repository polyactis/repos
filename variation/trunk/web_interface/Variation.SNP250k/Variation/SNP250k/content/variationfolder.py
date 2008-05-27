"""Definition of the FilmFolder content type. See cinemafolder.py for more
explanation on the statements below.
"""

from zope.interface import implements

from Products.Archetypes import atapi

from Products.ATContentTypes.content import folder
from Products.ATContentTypes.content.schemata import finalizeATCTSchema

from Variation.SNP250k.interfaces import IVariationFolder
from Variation.SNP250k.config import PROJECTNAME

from Variation.SNP250k import VariationMessageFactory as _

VariationFolderSchema = folder.ATBTreeFolderSchema.copy() + atapi.Schema((
	atapi.TextField('text',
		required=False,
		searchable=True,
		storage=atapi.AnnotationStorage(),
		default_output_type='text/html',
		allowable_content_types = ('text/plain', 'text/restructured', 'text/html',),
		widget=atapi.RichWidget(label=_(u"Descriptive text"),
								description=_(u""),
								rows=25,
								allow_file_upload=True),
		),
	))

VariationFolderSchema['title'].storage = atapi.AnnotationStorage()
VariationFolderSchema['description'].storage = atapi.AnnotationStorage()

finalizeATCTSchema(VariationFolderSchema, folderish=True, moveDiscussion=False)

class VariationFolder(folder.ATBTreeFolder):
	"""Contains multiple films.
	"""
	implements(IVariationFolder)
	
	portal_type = "Variation Folder"
	_at_rename_after_creation = True
	schema = VariationFolderSchema
	
	title = atapi.ATFieldProperty('title')
	description = atapi.ATFieldProperty('description')
	text = atapi.ATFieldProperty('text')

atapi.registerType(VariationFolder, PROJECTNAME)