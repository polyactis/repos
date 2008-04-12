import unittest
from optilux.policy.tests.base import OptiluxPolicyTestCase

from Products.CMFCore.utils import getToolByName

class TestSetup(OptiluxPolicyTestCase):
    
    def afterSetUp(self):
        self.workflow = getToolByName(self.portal, 'portal_workflow')
        self.acl_users = getToolByName(self.portal, 'acl_users')
        self.types = getToolByName(self.portal, 'portal_types')
    
    def test_portal_title(self):
        self.assertEquals("Optilux Cinemas", self.portal.getProperty('title'))
        
    def test_portal_description(self):
        self.assertEquals("Welcome to Optilux Cinemas", self.portal.getProperty('description'))
        
    def test_role_added(self):
        self.failUnless("StaffMember", self.portal.validRoles())
        
    def test_workflow_installed(self):
        self.failUnless('optilux_sitecontent_workflow' in self.workflow.objectIds())
        
    def test_workflows_mapped(self):
        self.assertEquals(('optilux_sitecontent_workflow',), self.workflow.getDefaultChain())
        for portal_type, chain in self.workflow.listChainOverrides():
            if portal_type in ('File', 'Image',):
                self.assertEquals(('optilux_sitecontent_workflow',), chain)
        
    def test_view_permisison_for_staffmember(self):
        # The API of the permissionsOfRole() function sucks - it is bound too
        # closely up in the permission management screen's user interface
        self.failUnless('View' in [r['name'] for r in 
                                self.portal.permissionsOfRole('Reader') if r['selected']])
        self.failUnless('View' in [r['name'] for r in 
                                self.portal.permissionsOfRole('StaffMember') if r['selected']])
        
    def test_staffmember_group_added(self):
        self.assertEquals(1, len(self.acl_users.searchGroups(name='Staff')))
        
    def test_richdocument_installed(self):
        self.failUnless('RichDocument' in self.types.objectIds())
        
    def test_plain_document_disabled(self):
        # the internal name for "Page" is "Document"
        document_fti = getattr(self.types, 'Document')
        self.failIf(document_fti.global_allow)
        
    def test_richdocument_renamed_to_page(self):
        rich_document_fti = getattr(self.types, 'RichDocument')
        self.assertEquals("Web page", rich_document_fti.title)
        
    def test_theme_installed(self):
        skins = getToolByName(self.portal, 'portal_skins')
        layer = skins.getSkinPath('Optilux Theme')
        self.failUnless('optilux_theme_custom_templates' in layer)
        self.assertEquals('Optilux Theme', skins.getDefaultSkin())
        
    def test_types_versioned(self):
        repository = getToolByName(self.portal, 'portal_repository')
        versionable_types = repository.getVersionableContentTypes()
        for type_id in ('RichDocument', 'Film', 'Cinema', 'Promotion',):
            self.failUnless(type_id in versionable_types)
            
    def test_promotions_not_in_navtree(self):
        self.failIf('Promotion' not in self.portal.portal_properties.navtree_properties.metaTypesNotToList)
            
    def test_cinemacontent_installed(self):
        self.failUnless('Cinema' in self.types.objectIds())
        
    def test_cinemafolder_filmfolder_have_no_workflow(self):
        for portal_type, chain in self.workflow.listChainOverrides():
            if portal_type in ('Cinema Folder', 'Film Folder'):
                self.assertEquals((), chain)
                
    def test_add_promotion_permission_for_staffmember(self):
        self.failUnless('Optilux: Add Promotion' in [r['name'] for r in 
                                self.portal.permissionsOfRole('StaffMember') if r['selected']])

def test_suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSetup))
    return suite
