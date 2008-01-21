#!/usr/bin/python

"""Format from EMBL Nucleotide Sequence Database Release 65, December 2000

"""
import sys, Martel
from Bio.expressions import embl
from xml.sax import saxutils, handler
import re

class my_handler(handler.ContentHandler):

	def __init__(self,fname):
		handler.ContentHandler.__init__(self)
		self.outf=open(fname,'w')
		self.m_pname=0
		self.m_accession=0
		self.m_seq=0
		self.m_seq_block=0
		
		self.seq=''
				
	def characters(self, s):
		anti_seq=re.compile(r'\s*')
		
		if self.m_accession==1:
			self.outf.write(s+'\n')
		if self.m_pname==1:
			pass
			#self.outf.write(s+'\t')
		if self.m_seq==1:
			self.seq=self.seq+anti_seq.sub('',s)
			#self.outf.write(self.seq)
		
	def startElement(self, name, attrs):
		if name=="accession":
			print "acc Start"
			self.m_accession=1
		if name=="pname":
			self.m_pname=1
		if name=="bioformat:sequence":
			self.m_seq=1
		if name=="bioformat:sequence_block":
			self.seq=''
			self.m_seq_block=1


	def endElement(self, name):
		if name == "accession":
			print "acc End"
			self.m_accession=0
		if name == "pname":
			self.m_pname=0
		if name=="bioformat:sequence":
			self.m_seq=0
		if name=="bioformat:sequence_block":
			#self.outf.write(self.seq+'\n')
			self.m_seq_block=0


if __name__=="__main__":
	exp = Martel.select_names(embl.format, ("record", "accession", "pname", "bioformat:sequence",))
	parser=embl.format_expression.make_parser()
	#outf=open(sys.argv[2],'w')
	parser.setContentHandler(my_handler(sys.argv[2]))
	parser.parse(sys.argv[1])
