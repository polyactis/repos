#!/usr/bin/python2.2

import sys
from Bio.expressions.blast import wublast
from xml.sax import handler
import re

class w2h_handler(handler.ContentHandler):
    
    def __init__(self,fname):
		handler.ContentHandler.__init__(self)
		self.outf=open(fname,'w')
	    	self.h_descrip=0
		self.t_descrip=0
    		self.score=0
    		self.acc=''

    def characters(self, s):
        dbhead=re.compile(r'^Os')
	if self.t_descrip==1 :
	    sp1=re.compile(r'\W+')
   	    self.acc=sp1.split(s)[1]
	    self.outf.write("<a href=http://10.100.113.147/est-bin/estacc.py?blastaccs="+self.acc+">")
	if self.h_descrip==1 and dbhead.search(s):
	    sp1=re.compile(r'\W+')
   	    tmp1=sp1.split(s)[1]
	    self.outf.write("<a name="+tmp1+"></a>")
	    self.outf.write("<a href=http://10.100.113.147/est-bin/estacc.py?blastaccs="+tmp1+">")

	if self.score==1:
	    	    self.outf.write("<a href=#"+self.acc+">")	
	self.outf.write(s)

    def startElement(self, name, attrs):
        if name == "record":
            self.outf.write("<pre>")
	if name == 'bioformat:search_table_description':
	    self.t_descrip=1
	    
	if name == 'bioformat:search_table_value' and attrs['name']=='score':
	    self.score=1

	if name == 'bioformat:hit_description':
	    self.h_descrip=1


    def endElement(self, name):
        if name == "record":
            self.outf.write("</pre>")
	if name == 'bioformat:search_table_description':
	    self.t_descrip=0
	    self.outf.write("</a>")
	if name == 'bioformat:hit_description':
	    self.h_descrip=0
	    self.outf.write("</a>")
	if name == 'bioformat:search_table_value':
	    self.score=0
	    self.outf.write("</a>")
if __name__=='__main__':

	#inf=open(sys.argv[1],'r')
	#outf=open(sys.argv[2],'w')

	parser = wublast.blastn.make_parser()
	parser.setContentHandler(w2h_handler(sys.argv[2]))
	parser.parse(sys.argv[1])

