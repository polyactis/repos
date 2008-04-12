from Products.CMFCore.utils import getToolByName

def v1_0_to_v1_1(portal_setup):
    """This is a migration step, referenced from configure.zcml
    
    Here, we "fix" the portal front page.
    """
    
    portal_url = getToolByName(portal_setup, 'portal_url')
    portal = portal_url.getPortalObject()
    front_page = portal['front-page']
    front_page.setTitle('Welcome to the Guitar Shop')
    
def v1_1_to_v1_2a(portal_setup):
    """Here is another upgrade step, this one part of a two-step upgrade
    """
    
    portal_catalog = getToolByName(portal_setup, 'portal_catalog')
    for brain in portal_catalog(portal_type = 'Document'):
        brain.getObject().setTitle('All your base are belong to us!')
    
def v1_1_to_v1_2b(portal_setup):
    """This is the second step
    """
    
    portal_url = getToolByName(portal_setup, 'portal_url')
    portal = portal_url.getPortalObject()
    
    # typo in properties.xml - we should fix that one too!
    portal.manage_changeProperties(title="Guitar Shop")
    
def v1_2_to_v1_3(portal_setup):
    """This example invokes an extension profile which in turn performs
    a migration.
    """
    portal_setup.runAllImportStepsFromProfile('profile-guitar.shop:1.2_to_1.3', purge_old=False)