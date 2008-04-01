"""Definition of the FilmFolder content type. See cinemafolder.py for more
explanation on the statements below.
"""

from zope.interface import implements

from Products.Archetypes import atapi

from Products.ATContentTypes.content import folder
from Products.ATContentTypes.content.schemata import finalizeATCTSchema

from Products.Variation.interfaces import IFilmFolder
from Products.Variation.config import PROJECTNAME

#from Products.Variation import VariationMessageFactory as _

VariationFolderSchema = folder.ATBTreeFolderSchema.copy()

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

atapi.registerType(VariationFolder, PROJECTNAME)