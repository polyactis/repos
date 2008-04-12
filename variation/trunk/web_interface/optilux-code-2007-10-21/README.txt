============================================
Professional Plone Development - Sample code
============================================

This folder contains example code for the book _Professional Plone Development_
by Martin Aspeli, published in 2007 by Packt Publishing. 

There is one directory for each chapter, showing the evolution of the example
code in the book. Where a chapter's directory is missing, there were no code 
changes in that chapter.

In each chapter directory, there is a sub-directory called 'optilux'. This was
created as in Chapter 3, with 'paster', and contains a "buildout" of the code
and configuration making up the Optilux Cinemas example application. The code
is not distributed in compiled form, so you will need to bootstrap and build
any examples you wish to run. You can do that with:

 $ python bootstrap.py
 $ ./bin/buildout

This assumes you have a viable Internet connection, over which Zope and Plone
will be downloaded and then installed.

NOTE: If you are on Windows, make sure to read the README.txt file found inside
the 'optilux' directory. It contains important information about dependencies
which you may not have. Also be sure to use backslashes (\) in place of forward
slahes (/) as directory separators if following the examples here and in the
book.

Python Imaging Library (PIL)
----------------------------

Plone requires PIL, which must be installed separately. You can find downloads
and installation instructions at http://www.pythonware.com/products/pil.

Python versions
---------------

The code targets Plone 3.0 on Zope 2.10.4 or later. It assumes that your system
Python is version 2.4.x (not 2.5). You can check which is your system Python by
running (on a command line):

 $ python -V
 Python 2.4.4

If it outputs something other than 2.4.x, you need to install Python 2.4 from
http://python.org. It will be easiest if Python 2.4 is the main version, i.e.
the one you get by running 'python' on the command line. This can be achieved
by making sure that this installation comes first in the PATH environment
variable.

If that is not possible, you need to explicltly invoke the correct version of
Python for each command shown in the book, e.g. with:

 $ /opt/python2.4/python 

In this case, you will also need to edit the buildout.cfg file (there is one
per chapter, but you may not need to run every single example), and add an
'executable' option to the '[buildout]' section:

 [buildout]
 executable = /opt/python2.4/python
 # other options follow as in the file ...

As an alternative, you can add a file called default.cfg in a directory
called ~/.buildout, containing:

 [buildout]
 executable = /opt/python2.4/python

This will all buildouts that do not set the executable explicitly.

Adjust the path to the Python binary as appropriate.

Setuptools and easy_install
---------------------------

'setuptools' is a Python library for distributing software, using the 'egg'
mechanism. easy_install is a program which makes it easy to find, download and
install packages, by searching the Cheese Shop at http://cheesehop.python.org.

Several of the examples in the book use the 'paster' command with templates
provided by 'ZopeSkel'. To get ZopeSkel and its dependencies (including the
Paste tools), you can run:

 $ easy_install ZopeSkel

If you do not have setuptools and easy_install available, you can run:

 $ python ez_setup.py

The ez_setup script is in the 'extra' directory. This will download and install
setuptools for the python version used to run the command. You may need to run
this as root.

The easy_install binary (easy_install.exe on Windows) will be installed as well,
into a directory such as /usr/local/python2.4/bin or C:\Python2.4\Scripts on
Windows. You probably want this in your system PATH. If you look at the log
output that ensues when ez_setup.py is run, it should show you where the
binaries were installed.

MySQL
-----

From chapter 12, the example code depends on a live database connection. You
should download and install MySQL 5.0 from http://mysql.org.

In the 'extra' directory, you will find two database scripts. optilux-tables.sql
creates the 'optilux' database and its tables, as well as a copy called
'optilux_test', used for the automated tests. You can run this with:

 $ mysql -u root < optilux-tables.sql

Assuming there is a user called 'root' with no password. Pass -p if you need
to specify a password - you will be prompted for one.

Secondly, optilux-testdata.sql will create some test data. Since screening
data is time-dependent, it will use the system time and various offsets.
You can truncate the screening table and re-run the script to create newer
test data if necessary.

To create the test data, run:

 $ mysql -u root < optilux-testdata.sql

LDAP
----

From chapter 18, the example code uses an LDAP connection for authentication.
To test this, you need to install openldap. See http://openldap.org. The
Quick Start guide on that website is a useful starting point.

The examples use the 'cosine' OpenLDAP schema. This can be included using
the following lines (with an appropriately adjusted path) in slapd.conf, the
OpenLDAP deamon configuration file.

 include /etc/openldap/schema/core.schema
 include /etc/openldap/schema/cosine.schema

You will also need to set the root manager user in order to import the example
records. In slapd.conf, this would be:

 suffix      "dc=optilux-cinemas,dc=com"
 rootdn      "cn=Manager,dc=optilux-cinemas,dc=com"
 
Now, if you have 'slapd' running, you can import the file optilux.ldif from the
'extra' directory.

 $ ldapadd -xWD 'cn=Manager,dc=optilux-cinemas,dc=com' -f optilux.ldif

If you wish to inspect (and modify) the LDAP schema, you may wish to use the
LDAP Browser/Editor from http://ldapmanager.com. This is a Java application,
so you will need Java from http://java.sun.com.
