#!/usr/bin/python2.2

import sys
from Bio.expressions.blast import wublast
from xml.sax import handler
import re

def my_min(twodlist,n):
	min=twodlist[0][n]
	for i in range(len(twodlist)):
		if twodlist[i][n]<min:
			min =twodlist[i][n]
	return min
	

def my_max(twodlist,n):
	max=twodlist[0][n]
	for i in range(len(twodlist)):
		if twodlist[i][n]>max:
			max =twodlist[i][n]
	return max


class w2g_handler(handler.ContentHandler):
	m_query_size=0
	m_hit_descrip=0
	m_hit_length=0
	m_bits=0
	m_prob=0
	m_hit=0
	m_hsp_seqalign_query_start=0
	m_hsp_seqalign_query_end=0
	m_hsp_seqalign_subject_start=0
	m_hsp_seqalign_subject_end=0
	seqcount=0
	query_size=0
	hit_descrip=[]
	hit_length=[]
	bits=[]
	prob=[]
	seqalign=[[],[],[],[]]
	seqaligntmp=[[],[],[],[]]
	output=''

	def __init__(self):
			handler.ContentHandler.__init__(self)

	def characters(self, s):
		dbhead=re.compile(r'^EMBL')
		if self.m_query_size == 1:
			self.query_size=int(s)
		if self.m_hit_length ==1:
			self.hit_length.append(int(s))
		if self.m_hit_descrip == 1 and dbhead.search(s):
			sp1=re.compile(r'\W+')
			tmp1=sp1.split(s)[1]
			self.hit_descrip.append(tmp1)
		if self.m_bits == 1:
			self.bits.append(int(s))
		if self.m_prob == 1:
			self.prob[self.seqcount-1]+=s
		if self.m_hsp_seqalign_query_start == 1:
			self.seqaligntmp[0].append(int(s))
		if self.m_hsp_seqalign_query_end == 1:
			self.seqaligntmp[1].append(int(s))
		if self.m_hsp_seqalign_subject_start == 1:
			self.seqaligntmp[2].append(int(s))
		if self.m_hsp_seqalign_subject_end == 1:
			self.seqaligntmp[3].append(int(s))
			

	def startElement(self, name, attrs):
		if name == 'bioformat:query_size':
			self.m_query_size=1
		if name == 'bioformat:hit_description':
			self.m_hit_descrip=1
		if name == 'bioformat:hit_length':
			self.m_hit_length=1
		if name == 'bioformat:hsp_value' and attrs['name']=='bits':
			self.m_bits=1
		if name == 'bioformat:hsp_value' and attrs['name']=='P':
			self.m_prob=1
			self.seqcount+=1
			self.prob.append('')
		if name == 'bioformat:hit':
			self.m_hit=1
			self.seqaligntmp=[[],[],[],[]]
			
		if name == 'bioformat:hsp_seqalign_query_start':
			self.m_hsp_seqalign_query_start=1
		if name == 'bioformat:hsp_seqalign_query_end':
			self.m_hsp_seqalign_query_end=1
		if name == 'bioformat:hsp_seqalign_subject_start':
			self.m_hsp_seqalign_subject_start=1
		if name == 'bioformat:hsp_seqalign_subject_end':
			self.m_hsp_seqalign_subject_end=1

	def endElement(self, name):
		if name == 'bioformat:query_size':
			self.m_query_size=0
		if name == 'bioformat:hit_description':
			self.m_hit_descrip=0
		if name == 'bioformat:hit_length':
			self.m_hit_length=0
		if name == 'bioformat:hsp_value':
			self.m_bits=0
		if name == 'bioformat:hsp_value':
			pass
			self.m_prob=0
		if name == 'bioformat:hit':
			self.m_hit=0
			self.seqalign[0].append(min(self.seqaligntmp[0]))
			self.seqalign[1].append(max(self.seqaligntmp[1]))
			self.seqalign[2].append(min(self.seqaligntmp[2]))
			self.seqalign[3].append(max(self.seqaligntmp[3]))
			
		if name == 'bioformat:hsp_seqalign_query_start':
			self.m_hsp_seqalign_query_start=0
		if name == 'bioformat:hsp_seqalign_query_end':
			self.m_hsp_seqalign_query_end=0
		if name == 'bioformat:hsp_seqalign_subject_start':
			self.m_hsp_seqalign_subject_start=0
		if name == 'bioformat:hsp_seqalign_subject_end':
			self.m_hsp_seqalign_subject_end=0
	
	def query_size_return(self):
		return self.query_size
	def hit_descrip_return(self):
		return self.hit_descrip
	def hit_length_return(self):
		return self.hit_length
	def bits_return(self):
		return self.bits
	def prob_return(self):
		return self.prob
	def seqalign_return(self):
		return self.seqalign
	

	def adjust(self):
		self.axis=[[]]
		self.axis[0].append(1)
		self.axis[0].append(1)
		self.axis[0].append(0)
		self.axis[0].append(self.query_size)
		for i in range(len(self.seqalign[0])):
			self.axis.append([])
			query_start=self.seqalign[0][i]
			query_end=self.seqalign[1][i]
			subject_start=self.seqalign[2][i]
			subject_end=self.seqalign[3][i]
			self.axis[i+1].append(1-(subject_start-query_start))
			self.axis[i+1].append(query_start)
			self.axis[i+1].append(query_end)
			self.axis[i+1].append(self.hit_length[i]-subject_end+query_end)
		return self.axis		

	def ascplot(self):
		
		line_start=my_min(self.axis,0)
		line_end=my_max(self.axis,3)
		span=line_end-line_start
		rspan=70
		zoom=span/rspan
				
		query_char='='
		blank_char=' '
		dismatch_char='-'
		match_char='*'
		for i in range(rspan):
			if i==0 :
				self.output+=('1')
			elif i==abs(int(line_start/zoom)):
				self.output+=(repr(abs(line_start)))
			elif i==int((abs(line_start)+axis[0][3])/zoom):
				self.output+=(repr(abs(line_start)+axis[0][3]))
			else:
				self.output+=('-')
				
		self.output+=(repr(span))
		self.output+=('\n')
		line_start=int(line_start/zoom)
		line_end=int(line_end/zoom)
		
		for i in range(len(self.axis)):
			for j in range(len(self.axis[i])):
				self.axis[i][j]=int(axis[i][j]/zoom)
		
		self.output+='<form action="http://10.100.113.147/est-bin/batch_estacc.py" method=post enctype="multipart/form-data">\n'
		self.output+='<table>\n'

		for i in range(len(self.axis)):

			self.output+='<tr><td>'
			self.output+='<pre>'
			
			for j in range(line_start,self.axis[i][0]):
				self.output+=(blank_char)
			for j in range(self.axis[i][0],axis[i][1]):
				self.output+=(dismatch_char)
			for j in range(self.axis[i][1],axis[i][2]+1):
				self.output+=(match_char)
			for j in range(self.axis[i][2],axis[i][3]):
				if i==0:
					self.output+=(query_char)
				else:
					self.output+=(dismatch_char)
			for j in range(self.axis[i][3],line_end):
				self.output+=(blank_char)

			if i>0: 
				self.output+=repr(self.hit_length[i-1])+'\t'
				self.output+='<input type=checkbox name="item" value='+self.hit_descrip[i-1]+' >'
				self.output+=self.hit_descrip[i-1]+'\t'
				self.output+='<input type=text name='+self.hit_descrip[i-1]+' size=9>'
				self.output+='</td></tr>'
				
#				self.output+=('\t'+self.hit_descrip[i-1]+'\t')
#				self.output+=(repr(self.bits[i-1])+'\t')
#				self.output+=(self.prob[i-1])
			else:			
				self.output+=('\t'+'Query'+'\t')
			self.output+=('\n')
			self.output+='</pre>'
		self.output+='<tr>\n<td>\n'
		self.output+='<input type=submit><input type=reset>\n'
		self.output+='</td>\n</tr>\n'
		self.output+='</table>\n'
		self.output+='</form>'
	
	def finish(self):
		self.f_output='<form action="http://10.100.113.147/est-bin/batch_estacc.py" method=post enctype="multipart/form-data">\n'
		self.f_output+='<table>\n'
		for i in range(len(self.hit_descrip)):
			self.f_output+='<tr>\n<td>\n'
			self.f_output+='<b>'+self.hit_descrip[i]+'</b>\n'
			self.f_output+='<input type=checkbox name="item" value='+self.hit_descrip[i]+' >\n'
			self.f_output+='</td>\n</tr>\n'
		
		self.f_output+='<tr>\n<td>\n'
		self.f_output+='<input type=submit><input type=reset>\n'
		self.f_output+='</td>\n</tr>\n'
		self.f_output+='</table>\n'
		self.f_output+='</form>'

	
inf=open(sys.argv[1],'r')
outf=open(sys.argv[2],'w')

parser = wublast.blastn.make_parser()
handler=w2g_handler()
parser.setContentHandler(handler)
parser.parse(inf)

#handler=parser.getContentHandler()
#handler=w2g_handler()


axis=handler.adjust()

handler.output+='<pre>\n'
handler.ascplot()
handler.output+='</pre>\n'

#handler.finish()

outf.write(handler.output)



