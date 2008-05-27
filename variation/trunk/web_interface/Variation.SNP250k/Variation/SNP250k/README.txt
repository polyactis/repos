About this product

  Polymorphism Data Portal for 250k db.
  
  See the Archetypes Developer Manual on Plone.org (http://plone.org/documentation/manual/archetypes-developer-manual) for explanation of the code.
    
Installation

  * Tested on Plone 3.0.5 with Zope 2.10.

  * After decompressing the archive, copy the product in your Zope instance's 'Products' directory.
  
  * Restart Zope
  
  * Go to Plone configuration menu "Add/Remove product", select the product and install it.
  
  * Have fun !

Contributors

  * Yu Huang
  
Missing bits to add later !

  * i18n
   
  * Allow a message to be sent to a member (like in a true Instant Message app)
  
Setting up and logging in
-------------------------

We use zope.testbrowser to simulate browser interaction in order to show
the main flow of pages. This is not a true functional test, because we also
inspect and modify the internal state of the ZODB, but it is a useful way of
making sure we test the full end-to-end process of creating and modifying
content.

    >>> from Products.Five.testbrowser import Browser
    >>> browser = Browser()
    >>> portal_url = self.portal.absolute_url()
    >>> print portal_url

The following is useful when writing and debugging testbrowser tests. It lets
us see error messages properly.

    >>> browser.handleErrors = False
    >>> self.portal.error_log._ignored_exceptions = ()

We then turn off the various portlets, because they sometimes duplicate links
and text (e.g. the navtree, the recent recent items listing) that we wish to
test for in our own views. Having no portlets makes things easier.

    >>> from zope.component import getMultiAdapter, getUtility
    >>> from plone.portlets.interfaces import IPortletManager
    >>> from plone.portlets.interfaces import IPortletAssignmentMapping

    >>> left_column = getUtility(IPortletManager, name=u"plone.leftcolumn")
    >>> left_assignable = getMultiAdapter((self.portal, left_column), IPortletAssignmentMapping)
    >>> for name in left_assignable.keys():
    ...     del left_assignable[name]

    >>> right_column = getUtility(IPortletManager, name=u"plone.rightcolumn")
    >>> right_assignable = getMultiAdapter((self.portal, right_column), IPortletAssignmentMapping)
    >>> for name in right_assignable.keys():
    ...     del right_assignable[name]

The database integration functionality is tested in tests/database.txt. We
don't want a database state dependency in this test, so we register a
dummy utility for finding screenings - otherwise, we will get errors on the
film and cinema view pages.
    >>> from zope.interface import implements
    >>> from Variation.SNP250k.interfaces import IDBLocator
    >>> class DummyPhenotypeLocator(object):
    ...     implements(IDBLocator)
    ...     
    ...     _method_id_ls = [1,2,3]
    ...     
    ...     def get_phenotype_method_id_ls(self):
    ...         return self._method_id_ls
    
We need to remember the existing utility so that we can set it back at the
end of the tests - otherwise, we may mess up other tests running later!
    
    >>> from zope.component import provideUtility
    >>> _old_screening_locator = getUtility(IDBLocator)
    >>> provideUtility(DummyPhenotypeLocator())
    
    >>> from zope.interface import implements
    >>> from Variation.SNP250k.interfaces import IDBLocator
    >>> locator = getUtility(IDBLocator)
    
    >>> v = locator.get_phenotype_method_id_ls()
	>>> print v

Clean Up
--------
We need to undo the damage we did so that it doesn't affect other tests!

    >>> provideUtility(_old_screening_locator)

Finally, we need to log in as the portal owner, i.e. an administrator user. We
do this from the login page.

    >>> from Products.PloneTestCase.setup import portal_owner, default_password
	>>> print portal_owner
	>>> print default_password
    >>> browser.open(portal_url + '/login_form?came_from=' + portal_url)
    >>> browser.getControl(name='__ac_name').value = portal_owner
    >>> browser.getControl(name='__ac_password').value = default_password
    >>> browser.getControl(name='submit').click()
    >>> "You are now logged in" in browser.contents
    True
    >>> print browser.url

Addable content
---------------

Verify that we have the links to create cinema and film folders, from the add
item menu:

	>>> browser.open(portal_url)
	>>> browser.getLink(id='variation-folder').url.endswith("createObject?type_name=Variation+Folder")
	True
	>>> browser.getLink(id='variation-folder').click()
	>>> browser.getControl(name='title').value = "variation3"
    >>> browser.getControl(name='form_submit').click()
    >>> 'variation3' in self.portal.objectIds()
    True
    >>> variation3 = self.portal['variation3']
    >>> variation3.title
    'variation3'
    >>> variation3_url = variation3.absolute_url()
    >>> browser.open(variation3_url)

Play around a bit according to /usr/lib/zope2.10/lib/python/zope/formlib/form.txt

    #>>> from zope.publisher.browser import TestRequest, BrowserResponse
    #>>> request = TestRequest(response=BrowserResponse())
    #>>> from Variation.SNP250k.browser.phenotype import AddPhenotypeForm
    #>>> print AddPhenotypeForm(None, request)() # doctest: +NORMALIZE_WHITESPACE

Create a phenotype in variation3 folder.
04/18/08 Watch: the custom edit form (inherited from zope.formlib.form) has prefix 'form' attached to each widget name.:

    >>> browser.open(variation3_url)
    >>> browser.getLink(id='phenotype').click()
    >>> browser.getControl(name='form.title').value = "first phenotype"
    >>> browser.getControl(name='form.description').value = "first phenotype"
    >>> browser.getControl(name='form.method_id_ls-empty-marker').value = "1"
    >>> browser.getControl(name='form.actions.save').click()
    >>> 'first-phenotype' in variation3.objectIds()
    True

Create a phenotype in variation3 folder.
04/18/08 Watch: the custom edit form (inherited from zope.formlib.form) has prefix 'form' attached to each widget name.:

    #>>> browser.open(variation3_url)
    #>>> browser.getLink(id='qcondirectory').click()
    #>>> browser.getControl(name='form.title').value = "first qc"
    #>>> browser.getControl(name='form.description').value = "first qc run test"
    #>>> browser.getControl(name='form.input_dir').value = "/Network/Data/250k/tmp-yh"
    #>>> browser.getControl(name='form.QC_method_id-empty-marker').value = "1"
    #>>> browser.getControl(name='form.actions.save').click()
    #>>> 'first-qc' in variation3.objectIds()
    #True

2008-05-23 test results2db_250k

	>>> browser.open(variation3_url)
	>>> browser.getLink('Results2DB_250k').click()
	>>> import StringIO
	>>> myResults = StringIO.StringIO('1\t2334\t343\n2\t4324\t83.2\n')
	>>> browser.getControl(name='form.username').value = "yh"
	>>> browser.getControl(name='form.password').value = "yh324"
    >>> browser.getControl(name='form.short_name').value = "96 with LD phenotype"
    >>> browser.getControl(name='form.phenotype_method_id-empty-marker').value = "1"
    >>> browser.getControl(name='form.call_method_id-empty-marker').value = "5"
    >>> browser.getControl(name='form.data_description').value = "best 96 picked by tina"
	>>> browser.getControl(name='form.results_method_type_id-empty-marker').value = "None"
	>>> browser.getControl(name='form.results_method_type_short_name').value = "new method"
    >>> browser.getControl(name='form.method_description').value = "kruskal wallis, -log(pvalue)"
	>>> browser.getControl(name='form.comment').value = "pvalue log "
	>>> control = browser.getControl(name='form.input_fname')
	>>> fileControl = control.mech_control
	>>> fileControl.add_file(myResults, 'text/plain', filename='myResults.tsv')
	>>> browser.getControl('No commit this transaction.').selected = True	#this doesn't work: browser.getControl(name='form.commit_type').value = ["1"]		#radio button
	>>> browser.getControl(name='form.actions.save').click()
	>>> browser.title
	>>> browser.headers
	>>> browser.contents
	>>> browser.url