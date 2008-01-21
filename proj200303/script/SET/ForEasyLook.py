#!/usr/bin/python


import sys,re,string

class easylook:
	
	def __init__(self,infname,dictfname):
		self.inf=open(infname,'r')
		self.dictf=open(dictfname,'r')
		self.dict={}
		self.make_dict()
		
	def make_dict(self):
		
		line=self.dictf.readline()
		
		while line:
			item=string.split(line,'\t')
			self.dict[item[0]]=item[1][0:len(item[1])-1]

			line=self.dictf.readline()

	def run(self):
		
		line=self.inf.readline()
		fasta_begin=re.compile(r'^>')

		while line:
			
			if fasta_begin.match(line):
				item=line[1:len(line)-1]
				if self.dict.has_key(item):
					line='>'+item+'/'+self.dict[item]+'\n'
			sys.stdout.write(line)
				
					
#			else:
#				sys.stdout.write(line)

			
			line=self.inf.readline()


if __name__=='__main__':
	instance=easylook(sys.argv[1],sys.argv[2])
	instance.run()
