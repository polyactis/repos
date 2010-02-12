"""Helper functions

Consists of functions to typically be used within templates, but also
available to Controllers. This module is available to both as 'h'.
"""
from webhelpers import *

#2008-12-24 for the forms
from routes import redirect_to
from routes import url_for
from webhelpers.html.tags import *


#2009-3-5 store 250k snp dataset sucked in from filesystem
call_method_id2dataset = {}
#2009-3-5 load ecotype_info on demand
ecotype_info = None


