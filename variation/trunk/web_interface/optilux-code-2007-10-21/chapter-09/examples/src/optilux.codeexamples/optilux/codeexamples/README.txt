====================
optilux.codeexamples
====================

This package contains various doctests with example code outlining Zope
development concepts. They are all set up in tests.py. See the various .txt
files and follow along with the examples there to learn about different 
aspects of Zope programming. Feel free to make changes to the example code
or to write new tests if it helps your understanding.

Run the tests are normal, with:

    $ ./bin/instance test -s optilux.codeexamples
    
from a buildout, or, from a standard Zope instance:

    $ ./bin/instance test -s optilux.codeexamples
    
Some tests use PloneTestCase to bootstrap an entire Zope 2 and Plone 
environment in the style of an integration test. Other tests use only the
Zope 3 testing harnesses, which are much more lightweight.

See also http://plone.org/documentation/tutorial/testing and 
http://docs.python.org/lib/module-doctest.html.