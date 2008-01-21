#!/usr/bin/python

import string,sys,re
#usage: ./format_mir.py inputfile output_sequence output_segment
class mir:
	
	def __init__(self,infilename,seqfilename,segfilename):
		self.infile=open(infilename,'r')
		self.seqfile=open(seqfilename,'w')
		self.segfile=open(segfilename,'w')
		self.acc=''
		self.seq=''
		self.seg=[0,0,0,0]
		self.run()
	
	def run(self):
		rawseq=''
		line=self.infile.readline()
		i=0
		anti_content=re.compile(r'\W*')

		while line:
			if line[0]=='>':
				if i>0:
					self.seq_process(rawseq)
					self.output()
				self.acc=line[1:len(line)-2]
				rawseq=''
				i=i+1

			else:
				rawseq=rawseq+anti_content.sub('',line)
			line=self.infile.readline()
			
		self.seq_process(rawseq)
		self.output()

	def seq_process(self,rawseq):
		blocks=rawseq.split('_')
		for i in range(0,len(blocks)):
			self.seg[i+1]=self.seg[i]+len(blocks[i])
		self.seq=string.join(blocks,'')

	def output(self):
		segnames=('pre_mir','mir','post_mir')
		
		self.seqfile.write(self.acc+'\t'+self.seq.lower()+'\n')

		for i in range(0,3):
			start =self.seg[i]+1
			tail=self.seg[i+1]

			self.segfile.write(self.acc+'\t'+repr(start)+'\t'+repr(tail)+'\t'+segnames[i]+'\n')
		

if __name__=='__main__':
	Mir=mir(sys.argv[1],sys.argv[2],sys.argv[3])
