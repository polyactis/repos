"""Definition of the FilmFolder content type. See cinemafolder.py for more
explanation on the statements below.
"""

from zope.interface import implements

from Products.Archetypes import atapi

from Products.ATContentTypes.content import folder
from Products.ATContentTypes.content.schemata import finalizeATCTSchema

from optilux.cinemacontent.interfaces import IFilmFolder
from optilux.cinemacontent.config import PROJECTNAME

from optilux.cinemacontent import CinemaMessageFactory as _

FilmFolderSchema = folder.ATFolderSchema.copy()

FilmFolderSchema['title'].storage = atapi.AnnotationStorage()
FilmFolderSchema['description'].storage = atapi.AnnotationStorage()

finalizeATCTSchema(FilmFolderSchema, folderish=True, moveDiscussion=False)

class FilmFolder(folder.ATFolder):
    """Contains multiple films.
    """
    implements(IFilmFolder)
    
    portal_type = "Film Folder"
    _at_rename_after_creation = True
    schema = FilmFolderSchema
    
    title = atapi.ATFieldProperty('title')
    description = atapi.ATFieldProperty('description')

atapi.registerType(FilmFolder, PROJECTNAME)
