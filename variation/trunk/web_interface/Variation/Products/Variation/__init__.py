# Register our skins directory - this makes it available via portal_skins.
from Products.CMFCore.utils import ContentInit
from Products.CMFCore.DirectoryView import registerDirectory

from Products.Archetypes import process_types
from Products.Archetypes.public import listTypes

from config import PROJECTNAME, GLOBALS

from Products.CMFCore import utils

from content import *
# Import the content types permissions
from permissions import ADD_CONTENT_PERMISSIONS
from zope.i18nmessageid import MessageFactory

VariationMessageFactory = MessageFactory(PROJECTNAME)

registerDirectory('skins', GLOBALS)

def initialize(context):
	"""Initializer called when used as a Zope 2 product."""
	import content

	contentTypes, constructors, ftis = process_types(listTypes(PROJECTNAME), PROJECTNAME)

	allTypes = zip(contentTypes, constructors)
	for atype, constructor in allTypes:
		kind = "%s: %s" % (PROJECTNAME, atype.portal_type)
		utils.ContentInit(kind,
						  content_types	  = (atype,),
						  permission		 = ADD_CONTENT_PERMISSIONS[atype.portal_type],
						  extra_constructors = (constructor,),
						  fti				= ftis,
						  ).initialize(context)
