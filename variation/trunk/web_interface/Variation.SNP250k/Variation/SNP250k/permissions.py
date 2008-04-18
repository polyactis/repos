
# CMF and Archetypes imports
from Products.CMFCore.permissions import setDefaultRoles
from Products.Archetypes.atapi import listTypes

# Product imports
from config import PROJECTNAME




"""
types = listTypes(PROJECTNAME)
for atype in  types:
    permission = "%s: Add %s" % (PROJECTNAME, atype['portal_type'])
    ADD_CONTENT_PERMISSIONS[atype['portal_type']] = permission

    # Assign default roles for the permission
    setDefaultRoles(permission, ('Owner', 'Manager',))
"""