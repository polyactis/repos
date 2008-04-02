Introduction
============

These tests are all about seeing and testing what the browser sees.  We
make no assumptions on the innards of the code -- pretending we have indeed
never seen the code itself.

First we need to setup a browser instance.

    >>> from Products.Five.testbrowser import Browser
    >>> browser = Browser()

Now we can start checking things out.  Really all we can test here now
is that the output to the browser has an integer that increments each time.

    >>> browser.open(portal.absolute_url()+'/@@view')
    >>> browser.contents
    'Retrieved 1'

    >>> browser.open(portal.absolute_url()+'/@@view')
    >>> browser.contents
    'Retrieved 2'
