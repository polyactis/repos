========================
 Optilux Cinema content
========================

This package contains content types that pertain to the Optilux Cinema
application. In this testbrowser doctest, we will demonstrate how the content 
types interact. See tests/test_doctest.py for how it is set up.

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

The following is useful when writing and debugging testbrowser tests. It lets
us see error messages properly.

    >>> browser.handleErrors = False
    >>> self.portal.error_log._ignored_exceptions = ()

We then turn off the various portlets, because they sometimes duplicate links
and text (e.g. the navtree, the recent recent items listing) that we wish to
test for in our own views. Having no portlets makes things easier.

    >>> from zope.component import getUtility, getMultiAdapter
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

Finally, we need to log in as the portal owner, i.e. an administrator user. We
do this from the login page.

    >>> from Products.PloneTestCase.setup import portal_owner, default_password

    >>> browser.open(portal_url + '/login_form?came_from=' + portal_url)
    >>> browser.getControl(name='__ac_name').value = portal_owner
    >>> browser.getControl(name='__ac_password').value = default_password
    >>> browser.getControl(name='submit').click()

Addable content
---------------

Cinema content is managed inside two root content types: A "Cinema Folder"
contains cinemas and information about them. A "Film Folder" contains films.

    >>> browser.open(portal_url)

Verify that we have the links to create cinema and film folders, from the add
item menu:

    >>> browser.getLink(id='cinema-folder').url.endswith("createObject?type_name=Cinema+Folder")
    True
    >>> browser.getLink(id='film-folder').url.endswith("createObject?type_name=Film+Folder")
    True
    
And likewise, we don't have add links for the other types
    
    >>> browser.getLink(id='cinema')
    Traceback (most recent call last):
    ...
    LinkNotFoundError
    >>> browser.getLink(id='film')
    Traceback (most recent call last):
    ...
    LinkNotFoundError
    >>> browser.getLink(id='promotion')
    Traceback (most recent call last):
    ...
    LinkNotFoundError
    
Adding film folders and films
-----------------------------

Film folders contain films. Unlike cinemas, films are self-contained and do 
not allow further pages of information. Like cinemas, films have a special 
code which relates to the existing relational database system. They also
specify start and end dates, specifying when the film is first shown and when 
the film will be taken off the program. Note that this is different from the 
Dublin Core effective/expiry date fields, which can be used to restrict when 
films show up in searches and so on. These are, however, part of the standard 
Plone metadata set, and site admins will be able to set them through the 
usual UI.

    >>> browser.open(portal_url)
    >>> browser.getLink(id='film-folder').click()
    >>> browser.getControl(name='title').value = "Films"
    >>> browser.getControl(name='form_submit').click()

    >>> 'films' in self.portal.objectIds()
    True
    >>> films = self.portal['films']
    >>> films_url = films.absolute_url()

To add a film, we need to fake a banner image:

    >>> import StringIO
    >>> dummy_banner1 = StringIO.StringIO('Dummy banner image contents')

Now let us add two films.

    >>> browser.open(films_url)
    >>> browser.getLink(id='film').click()
    >>> browser.getControl(name='title').value = "Attack of the Ninjas"
    >>> browser.getControl(name='description').value = "Ninjas attack. People angry."
    >>> browser.getControl(name='filmCode').value = "F1"
    >>> browser.getControl(name='teaser').value = "<b>Who doesn't like attacking ninjas?</b>"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner1, filename='dummy.png')
    >>> browser.getControl(name='startDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='startDate_month').getControl("March").selected = True
    >>> browser.getControl(name='startDate_day').getControl("15").selected = True
    >>> browser.getControl(name='endDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='endDate_month').getControl("September").selected = True
    >>> browser.getControl(name='endDate_day').getControl("15").selected = True
    >>> browser.getControl(name='form_submit').click()
    
    >>> 'attack-of-the-ninjas' in films.objectIds()
    True
    >>> film1 = films['attack-of-the-ninjas']
    >>> film1_url = film1.absolute_url()
    >>> film1.film_code
    'F1'
    >>> film1.title
    'Attack of the Ninjas'
    >>> film1.summary
    'Ninjas attack. People angry.'
    >>> film1.teaser
    "<b>Who doesn't like attacking ninjas?</b>"
    
    >>> # film1.shown_from
        # datetime.datetime(2007, 3, 15, 0, 0)
    >>> # film1.shown_until
        # datetime.datetime(2007, 9, 15, 0, 0)

NOTE - there is a bug in Archetypes presently: the above two tests will fail,
because the input widget for the date field does not work properly without
JavaScript (and testbrowser does not emulate JavaScript). Therefore, we need
to fix this manually for now.

    >>> from datetime import datetime
    >>> film1.shown_from = datetime(2007, 3, 15, 0, 0)
    >>> film1.shown_until = datetime(2007, 9, 15, 0, 0)
    >>> film1.reindexObject()

Create a new dummy image for this film. We could have used the same one as
for the first film, but then we'd need to remember to seek back to the
beginning of the "file".

    >>> dummy_banner2 = StringIO.StringIO('More dummy banner image contents')

    >>> browser.open(films_url)
    >>> browser.getLink(id='film').click()
    >>> browser.getControl(name='title').value = "Attack of the Zombies"
    >>> browser.getControl(name='description').value = "Zombies attack. People irritated."
    >>> browser.getControl(name='filmCode').value = "F1"
    >>> browser.getControl(name='teaser').value = "<b>Who doesn't like flesh-eating zombies?</b>"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner2, filename='dummy.png')
    >>> browser.getControl(name='startDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='startDate_month').getControl("April").selected = True
    >>> browser.getControl(name='startDate_day').getControl("10").selected = True
    >>> browser.getControl(name='endDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='endDate_month').getControl("June").selected = True
    >>> browser.getControl(name='endDate_day').getControl("10").selected = True
    >>> browser.getControl(name='form_submit').click()
    >>> "The film code is already in use" in browser.contents
    True
    
Oops - a duplicate film code. Try again:

    >>> dummy_banner2a = StringIO.StringIO('More dummy banner image contents')
    >>> browser.getControl(name='filmCode').value = "F2"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner2a, filename='dummy.png')
    >>> browser.getControl(name='form_submit').click()
    
    >>> 'attack-of-the-zombies' in films.objectIds()
    True
    >>> film2 = films['attack-of-the-zombies']
    >>> film2_url = film2.absolute_url()
    
    >>> film2.shown_from = datetime(2007, 4, 10, 0, 0)
    >>> film2.shown_until = datetime(2007, 6, 15, 0, 0)
    >>> film2.reindexObject()
    
We should now be on the view of the film. The banner provider adapter should
have ensured that we got a URL referring to the banner image as well.

    >>> browser.url == film2_url
    True
    >>> "Zombies attack. People irritated." in browser.contents
    True
    >>> film2_url + '/image_thumb' in browser.contents
    True
    
Let us quickly publish these films. We need this so that they show up in
the vocabulary for selecting highlighted films in the cinema edit screen.

    >>> self.setRoles(('Manager',))
    >>> self.portal.portal_workflow.doActionFor(film1, 'publish')
    >>> self.portal.portal_workflow.doActionFor(film2, 'publish')
    
Adding cinema folders and cinemas
---------------------------------
    
Let us now add a cinema folder and some cinemas. The cinema folder
can contain a rich-text description of the cinema folder (e.g. of a group
of cinemas), which will be displayed on the front page of that folder. 

    >>> browser.open(portal_url)
    >>> browser.getLink(id='cinema-folder').click()
    >>> browser.getControl(name='title').value = "Cinemas"
    >>> browser.getControl(name='text').value = "<b>About this cinema</b>"
    >>> browser.getControl(name='form_submit').click()

This should have added an object called 'cinemas' in the portal root, invoking
the title-to-id renaming.

    >>> 'cinemas' in self.portal.objectIds()
    True
    >>> cinemas = self.portal['cinemas']
    >>> cinemas.title
    'Cinemas'
    >>> cinemas.text
    '<b>About this cinema</b>'

    >>> cinemas_url = cinemas.absolute_url()

Cinemas include a cinema code. This will be used later to link cinemas and
films to the organization's existing booking system, which holds these values
as keys in a relational database. For now, we will just populate them with
strings. There is also an address, which is mapped to the Description Dublin
Core metadata field, a phone number, and some free text describing the cinema.

    >>> browser.open(cinemas_url)
    >>> browser.getLink(id='cinema').click()
    >>> browser.getControl(name='title').value = "First Cinema"
    >>> browser.getControl(name='description').value = "1 Star Walk, Hollywood"
    >>> browser.getControl(name='cinemaCode').value = "C1"
    >>> browser.getControl(name='phone').value = "1-800 555 111"
    >>> browser.getControl(name='text').value = "Welcome to the <em>First Cinema</em>"
    >>> browser.getControl(name='highlightedFilms:list').value = [film1.UID(), film2.UID()]
    >>> browser.getControl(name='form_submit').click()
    
    >>> 'first-cinema' in cinemas.objectIds()
    True
    >>> cinema1 = cinemas['first-cinema']
    >>> cinema1_url = cinema1.absolute_url()
    
    >>> cinema1.name
    'First Cinema'
    >>> cinema1.address
    '1 Star Walk, Hollywood'
    >>> cinema1.cinema_code
    'C1'
    >>> cinema1.phone
    '1-800 555 111'
    >>> cinema1.text
    'Welcome to the <em>First Cinema</em>'
    >>> cinema1.highlighted_films[0].title
    'Attack of the Ninjas'
    >>> cinema1.highlighted_films[1].title
    'Attack of the Zombies'

Cinema codes must be unique site-wide. This is validated on save.

    >>> browser.open(cinemas_url)
    >>> browser.getLink(id='cinema').click()
    >>> browser.getControl(name='title').value = "Second Cinema"
    >>> browser.getControl(name='description').value = "2 Star Walk, Hollywood"
    >>> browser.getControl(name='cinemaCode').value = "C1"
    >>> browser.getControl(name='phone').value = "1-800 555 222"
    >>> browser.getControl(name='text').value = "Welcome to the <em>Second Cinema</em>"
    >>> browser.getControl(name='form_submit').click()
    >>> "The cinema code is already in use" in browser.contents
    True
    
Let's correct that...

    >>> browser.getControl(name='cinemaCode').value = "C2"
    >>> browser.getControl(name='form_submit').click()
    
    >>> 'second-cinema' in cinemas.objectIds()
    True
    >>> cinema2 = cinemas['second-cinema']
    >>> cinema2_url = cinema2.absolute_url()

The cinema folder view should now be listing these two cinemas, provided we
have the rights to see them, which we do as a manager user. For other users,
this will depend on workflow permissions, of course.

    >>> browser.open(cinemas_url)
    >>> browser.getLink("First Cinema").url == cinemas_url + '/first-cinema'
    True
    >>> browser.getLink("Second Cinema").url == cinemas_url + '/second-cinema'
    True
    
    >>> '1 Star Walk, Hollywood' in browser.contents
    True
    >>> '2 Star Walk, Hollywood' in browser.contents
    True
    
A cinema is actually a folder which can contain pages of text, linked to
from its front page.

    >>> dummy_image = StringIO.StringIO('A dummy image')
    >>> browser.open(cinema2_url)
    >>> browser.getLink(id="image").click()
    >>> browser.getControl(name='title').value = "More information"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_image, filename='map.png')
    >>> browser.getControl(name="form_submit").click()
    
    >>> 'map.png' in cinema2.objectIds()
    True
    
Promotions and the promotions portlet
-------------------------------------

Let us add some promotions. We will make one that is expired and one that is 
not. We will need to do some date acrobatics for this.

    >>> from datetime import timedelta
    >>> yesterday = datetime.now() - timedelta(1)
    >>> day_before_yesterday = datetime.now() - timedelta(2)
    >>> tomorrow = datetime.now() + timedelta(1)
    >>> day_after_tomorrow = datetime.now() + timedelta(2)

NOTE - Because of the aforementioned Archetypes bug, we have to do a bit of
nasty surgery here to make submitting the edit form actually work. :-(

    >>> from optilux.cinemacontent.content.promotion import Promotion
    >>> Promotion.schema['effectiveDate'].required = False
    >>> Promotion.schema['expirationDate'].required = False

Moving swiftly on...

    >>> browser.open(cinema2_url)
    >>> dummy_banner3 = StringIO.StringIO('Third dummy banner')
    >>> browser.getLink(id='promotion').click()
    >>> browser.getControl(name='title').value = "First Promotion"
    >>> browser.getControl(name='description').value = "The first promotion is great"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner3, filename='dummy.png')
    >>> browser.getControl(name='details').value = "<b>More information here!</b>"
    >>> browser.getControl(name='effectiveDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='effectiveDate_month').getControl("January").selected = True
    >>> browser.getControl(name='effectiveDate_day').getControl("15").selected = True
    >>> browser.getControl(name='effectiveDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='effectiveDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='expirationDate_month').getControl("January").selected = True
    >>> browser.getControl(name='expirationDate_day').getControl("16").selected = True
    >>> browser.getControl(name='expirationDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='form_submit').click()

    >>> 'first-promotion' in cinema2.objectIds()
    True
    >>> promotion1 = cinema2['first-promotion']
    >>> promotion1.shown_from = day_before_yesterday
    >>> promotion1.shown_until = yesterday
    >>> promotion1.reindexObject()

    >>> browser.open(cinema2_url)
    >>> dummy_banner4 = StringIO.StringIO('Fourth dummy banner')
    >>> browser.getLink(id='promotion').click()
    >>> browser.getControl(name='title').value = "Second Promotion"
    >>> browser.getControl(name='description').value = "The second promotion is great"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner4, filename='dummy.png')
    >>> browser.getControl(name='details').value = "<b>More information here!</b>"
    >>> browser.getControl(name='effectiveDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='effectiveDate_month').getControl("January").selected = True
    >>> browser.getControl(name='effectiveDate_day').getControl("15").selected = True
    >>> browser.getControl(name='effectiveDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='effectiveDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='expirationDate_month').getControl("January").selected = True
    >>> browser.getControl(name='expirationDate_day').getControl("16").selected = True
    >>> browser.getControl(name='expirationDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='form_submit').click()

    >>> 'second-promotion' in cinema2.objectIds()
    True
    >>> promotion2 = cinema2['second-promotion']
    >>> promotion2.shown_from = yesterday
    >>> promotion2.shown_until = tomorrow
    >>> promotion2.reindexObject()

Promotions can also be added inside cinema folders.

    >>> browser.open(cinemas_url)
    >>> dummy_banner5 = StringIO.StringIO('Fifth dummy banner')
    >>> browser.getLink(id='promotion').click()
    >>> browser.getControl(name='title').value = "Third Promotion"
    >>> browser.getControl(name='description').value = "The third promotion is great"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner5, filename='dummy.png')
    >>> browser.getControl(name='details').value = "<b>More information here!</b>"
    >>> browser.getControl(name='effectiveDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='effectiveDate_month').getControl("January").selected = True
    >>> browser.getControl(name='effectiveDate_day').getControl("15").selected = True
    >>> browser.getControl(name='effectiveDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='effectiveDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='expirationDate_month').getControl("January").selected = True
    >>> browser.getControl(name='expirationDate_day').getControl("16").selected = True
    >>> browser.getControl(name='expirationDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='form_submit').click()

    >>> 'third-promotion' in cinemas.objectIds()
    True
    >>> promotion3 = cinemas['third-promotion']
    >>> promotion3.shown_from = tomorrow
    >>> promotion3.shown_until = day_after_tomorrow
    >>> promotion3.reindexObject()

    >>> browser.open(cinemas_url)
    >>> dummy_banner6 = StringIO.StringIO('Sixth dummy banner')
    >>> browser.getLink(id='promotion').click()
    >>> browser.getControl(name='title').value = "Fourth Promotion"
    >>> browser.getControl(name='description').value = "The fourth promotion is great"
    >>> browser.getControl(name='image_file').mech_control.add_file(dummy_banner6, filename='dummy.png')
    >>> browser.getControl(name='details').value = "<b>More information here!</b>"
    >>> browser.getControl(name='effectiveDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='effectiveDate_month').getControl("January").selected = True
    >>> browser.getControl(name='effectiveDate_day').getControl("15").selected = True
    >>> browser.getControl(name='effectiveDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='effectiveDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_year').getControl("2007").selected = True
    >>> browser.getControl(name='expirationDate_month').getControl("January").selected = True
    >>> browser.getControl(name='expirationDate_day').getControl("16").selected = True
    >>> browser.getControl(name='expirationDate_hour').getControl("10").selected = True
    >>> browser.getControl(name='expirationDate_minute').getControl("10").selected = True
    >>> browser.getControl(name='form_submit').click()

    >>> 'fourth-promotion' in cinemas.objectIds()
    True
    >>> promotion4 = cinemas['fourth-promotion']
    >>> promotion4.shown_from = yesterday
    >>> promotion4.shown_until = day_after_tomorrow
    >>> promotion4.reindexObject()

When the cinema folder was added, a promotions portlets was created for it.
This portlet contains effective, visible promotions, by default restricted
to promotions in the current cinema folder or cinema.

    >>> browser.open(cinemas_url)
    >>> "Promotions" in browser.contents
    True
    >>> "portletPromotions" in browser.contents
    True
    
    >>> "First Promotion" in browser.contents
    False
    >>> "Second Promotion" in browser.contents
    True
    >>> "Third Promotion" in browser.contents
    False
    >>> "Fourth Promotion" in browser.contents
    True
    
    >>> browser.open(cinema2_url)
    >>> "Promotions" in browser.contents
    True
    >>> "portletPromotions" in browser.contents
    True
    
    >>> "First Promotion" in browser.contents
    False
    >>> "Second Promotion" in browser.contents
    True
    >>> "Third Promotion" in browser.contents
    False
    >>> "Fourth Promotion" in browser.contents
    False