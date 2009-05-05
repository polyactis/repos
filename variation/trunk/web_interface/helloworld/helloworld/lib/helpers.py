"""Helper functions

Consists of functions to typically be used within templates, but also
available to Controllers. This module is available to both as 'h'.
"""
from webhelpers import *

#2008-12-24 for the forms
from routes import redirect_to
from routes import url_for
from webhelpers.html.tags import *

#2009-3-4 common URLs to link objects to external websites
NCBIGeneDBURL = 'http://www.ncbi.nlm.nih.gov/sites/entrez?cmd=search&db=gene&term=%s[uid]'
GBrowseURL = 'http://mahogany.usc.edu/cgi-bin/gbrowse/arabidopsis/?start=%s;stop=%s;ref=Chr%s;width=640;version=100;cache=on;drag_and_drop=on;show_tooltips=on;grid=on;label=BAC-ProteinCoding-Pseudogene-TEGenes-'
#2009-3-5 same as GBrowseURL above, but for handling in javascript
GBrowseURLJS = 'http://mahogany.usc.edu/cgi-bin/gbrowse/arabidopsis/?start={0};stop={1};ref=Chr{2};width=640;version=100;cache=on;drag_and_drop=on;show_tooltips=on;grid=on;label=BAC-ProteinCoding-Pseudogene-TEGenes-'

#2009-3-5 store 250k snp dataset sucked in from filesystem
call_method_id2dataset = {}
#2009-3-5 load ecotype_info on demand
ecotype_info = None


