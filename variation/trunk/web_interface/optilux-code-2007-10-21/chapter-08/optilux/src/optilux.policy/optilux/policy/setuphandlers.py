from Products.CMFCore.utils import getToolByName

def setupGroups(portal):
    acl_users = getToolByName(portal, 'acl_users')
    if not acl_users.searchGroups(name='Staff'):
        gtool = getToolByName(portal, 'portal_groups')
        gtool.addGroup('Staff', roles=['StaffMember'])
            
def renameRichDocument(portal):
    portal_types = getToolByName(portal, 'portal_types')
    rich_document_fti = getattr(portal_types, 'RichDocument')
    rich_document_fti.title = "Web page"
    
def disableDocument(portal):
    portal_types = getToolByName(portal, 'portal_types')
    document_fti = getattr(portal_types, 'Document')
    document_fti.global_allow = False
    
def importVarious(context):
    """Miscellanous steps import handle
    """
    
    # Ordinarily, GenericSetup handlers check for the existence of XML files.
    # Here, we are not parsing an XML file, but we use this text file as a 
    # flag to check that we actually meant for this import step to be run.
    # The file is found in profiles/default.
    
    if context.readDataFile('optilux.policy_various.txt') is None:
        return
    
    portal = context.getSite()
    
    setupGroups(portal)
    renameRichDocument(portal)
    disableDocument(portal)