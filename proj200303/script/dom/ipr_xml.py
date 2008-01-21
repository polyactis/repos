#!/usr/bin/python2.2

from xml.sax import saxexts
from xml.sax import saxlib
import sys

class ipr_DocumentHandler(saxlib.DocumentHandler):

	def __init__(self):
		self.i=0
		self.j=0
		self.ipr_prots=[]
		self.dom_prots=[]
		self.length=''
		self.pacc=''
		self.ipracc=''
		self.iprid=''
		self.method_acc=''
		self.start=''
		self.tail=''
	
	def startElement(self,name,attrs):
		if name=='protein':
			self.pacc=attrs['id']
			self.length=attrs['length']
			
		if name=='interpro':
			self.ipracc=attrs['id']
			self.iprid=attrs['name']
			
		if name=='match':
			self.method_acc=attrs['id']
		
		if name=='location':
			self.start=attrs['start']
			self.tail=attrs['end']
		
	def characters(self,data,start,length):
		pass

	def endElement(self,name):
		if name=='location':
			self.ipr_prots.append([])
			self.ipr_prots[self.i].append(self.pacc)
			self.ipr_prots[self.i].append(self.method_acc)
			self.ipr_prots[self.i].append(self.start)
			self.ipr_prots[self.i].append(self.tail)
			self.ipr_prots[self.i].append(self.ipracc)
			self.ipr_prots[self.i].append(self.iprid)
			self.i=self.i+1
	
			self.dom_prots.append([])
			self.dom_prots[self.j].append(self.pacc)
			self.dom_prots[self.j].append(self.start)
			self.dom_prots[self.j].append(self.tail)
			self.dom_prots[self.j].append(self.iprid)
			self.j=self.j+1

	
		if name=='protein':
			self.dom_prots.append([])
			self.dom_prots[self.j].append(self.pacc)
			length_1=int(self.length)-1
			self.dom_prots[self.j].append(repr(length_1))
			self.dom_prots[self.j].append(self.length)
			self.dom_prots[self.j].append('  ')
			self.j=self.j+1


parser=saxexts.make_parser()
dh=ipr_DocumentHandler()

parser.setDocumentHandler(dh)

infname=sys.argv[1]

parser.parse(infname)

ipr_outf=open(infname+'.ipr','w')
dom_outf=open(infname+'.dom','w')

for i in range(len(dh.ipr_prots)):
	for j in range(6):
		if j!=5:
			ipr_outf.write(dh.ipr_prots[i][j]+'\t')
		else:
			ipr_outf.write(dh.ipr_prots[i][j]+'\n')


for i in range(len(dh.dom_prots)):
	for j in range(4):
		if j!=3:
			dom_outf.write(dh.dom_prots[i][j]+'\t')
		else:
			dom_outf.write(dh.dom_prots[i][j]+'\n')


ipr_outf.close()
dom_outf.close()
 
