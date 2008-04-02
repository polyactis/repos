# Python imports
from StringIO import StringIO

# CMF imports
from Products.CMFCore.utils import getToolByName

# Archetypes imports
from Products.Archetypes import listTypes
from Products.Archetypes.Extensions.utils import installTypes, install_subskin

# Products imports
from Products.Variation.config import GLOBALS, PROJECTNAME
from Acquisition import aq_base

def install(portal):
	setup_tool = getToolByName(portal, 'portal_setup')
	originalContext = setup_tool.getImportContextID()
	setup_tool.setImportContext('profile-Products.Variation:default')
	setup_tool.runAllImportSteps()
	setup_tool.setImportContext(originalContext)
	#return "Ran all import steps."
	
	"""Install content types, skin layer, enable the portal factory
	"""
	out = StringIO()
	
	print >> out, "Installing Variation"
	
	# Install types
	classes = listTypes(PROJECTNAME)
	installTypes(portal, out, 
				  classes, 
				  PROJECTNAME)
	print >> out, "Installed types"
	
	# Install skin
	install_subskin(portal, out, GLOBALS)
	print >> out, "Installed skin"
	
	# Register types with portal_factory
	factory = getToolByName(portal, 'portal_factory')
	#types = factory.getFactoryTypes().keys()
	#if 'Variation' not in types:
	#	types.append('Variation')
	#	factory.manage_setPortalFactoryTypes(listOfTypeIds = types)
	regTypes = factory.getFactoryTypes().keys()	# Registered types on portal_factory
	for atype in classes: 
		if atype['portal_type'] not in regTypes:	 # atype is a dict // atype['portal_type'] is a string (e.g. 'InstantMessage')
			regTypes.append(atype['portal_type'])
			print >> out, "Added %s to registered types" % atype['portal_type']
	factory.manage_setPortalFactoryTypes(listOfTypeIds = regTypes)
	print >> out, "Added Variation to portal_factory"
	return out.getvalue()
	 
def uninstall(self, reinstall=0):
	'''
	2008-03-07
		add this uninstall method based on uninstall() in CMFExtFile's Install.py
		Uninstalls the Variation portal types and tool.
	'''
	out = StringIO()
	types = getToolByName(self, 'portal_types', None)
	if types is None:
		return

	for portal_type in ('Variation',):
		if getattr(aq_base(types), portal_type, None) is not None:
			types._delObject(portal_type)
	skins = getToolByName(self, 'portal_skins', None)
	if skins is None:
		return

	if getattr(aq_base(skins), 'variation', None) is not None:
		skins._delObject('variation')
	
	selections = skins.getSkinSelections()
	for skin_id in selections:
		path = skins.getSkinPath(skin_id)
		path = [x.strip() for x in path.split(',')]
		path = [x for x in path if x]
		if 'variation' in path:
			path.remove('variation')
			skins.addSkinSelection(skin_id, ','.join(path))
	print >> out, "Successfully uninstalled %s." % PROJECTNAME
	return out.getvalue()