#!/usr/bin/python
"""Example code for querying Entrez and saving the results.

This code is meant to go with the 'Connecting with biological databases'
section of the Tutorial demonstrates Entrez connectivity.

See http://www.ncbi.nlm.nih.gov/entrez/query/static/linking.html for
more help understanding the parameters passed.

This also requires a web browser to run -- either netscape or lynx
are supported in this example."""
# standard library
import os,sys,re

# biopython
from Bio.WWW import NCBI
	
def reademblaccs(inf):
	accs=[]
	tmp=inf.readline()
	anti_content=re.compile(r'\s')
	
	while tmp:
		tmp=tmp.split()[1]
		tmp=tmp.split('|')
		accs=accs+tmp
		tmp = inf.readline()
	return accs

search_command = 'Search'
search_database = 'Nucleotide'
return_format = 'GenBank'
search_term = 'Cypripedioideae'
search_term = 'CAB71104.1'
my_browser = 'lynx'

inf = open(sys.argv[1],'r')
result_file_name = os.path.join(os.getcwd(), 'results.html')
result_file = open(result_file_name, 'w')
accs = reademblaccs(inf)

for acc in accs:
	search_term = acc
	result_handle = NCBI.query(search_command, search_database, term = search_term,itemID=13,doptcmdl = return_format)

	result_file.write(result_handle.read())
	sys.stderr.write(acc+' finished\n' )

"""
if my_browser == 'lynx':
    os.system('lynx -force_html ' + result_file_name)
elif my_browser == 'netscape':
    os.system('netscape file:' + result_file_name)

"""
