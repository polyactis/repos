#!/usr/bin/python

def smart2SETdb(inf,outf):
	a=inf.read(1)
	b=''
	while a:
		if a == '>':
			b=b+'\n'
			i=0
			j=0
			outf.write(b)
			b=''
		elif a == '|':
			if j == 0 or j == 1:
				b=b+'\t'
				j=j+1
			else:
				b=b+a
				j=j+1
		elif a == ' ':
			if i == 0:
				b=b+'\t'
				i=1
			else:
				b=b+a
		elif a == '[':
		
			b=b+'\t'
			
		elif a == ']':

			b=b+'\t'
		elif a == '\n':
			pass

		elif a == '\r':
			pass

		else:
			b=b+a
		
		a=inf.read(1)
		
	
	outf.write(b)



def chromdb2SETdb(inf,outf):
	a=inf.read(1)
	b=''
	while a:
		if a == '>':
			b=b+'\n'
			i=0
			outf.write(b)
			b=''
		elif a == '{':
		
			b=b
			
		elif a == '}':

			b=b+'\t'
		elif a == '[':
		
			b=b+'\t'
			
		elif a == ']':

			b=b+'\t'
		elif a == '\n':
			if i == 0:
				b=b+'\t'
				i=1
			else:
				b=b

		elif a == '\r':
			b=b

		else:
			b=b+a
		
		a=inf.read(1)
		
	
	outf.write(b)


def chromdbnt2SETdb(inf,outf):
	a=inf.read(1)
	b=''
	while a:
		if a == '>':
			b=b+'\n'
			i=0
			j=0
			outf.write(b)
			b=''
		elif a == ' ':
			if i == 0 or i==1 or i==3 or i==5:
				b=b+'\t'
				i=i+1
			else:
				b=b+a
				i=i+1
				
		elif a == '{':
		
			pass	
			
		elif a == '}':

			b=b+'\t'
		elif a == '[':
		
			pass
			
		elif a == ']':

			pass
		elif a == '\n':
			if j == 0:
				b=b+'\t'
				j=j+1
			else:
				pass

		elif a == '\r':
			pass

		else:
			b=b+a
		
		a=inf.read(1)
		
	
	outf.write(b)


def emblest2osest(inf,outf):
	a=inf.read(1)
	b=''
	i=0
	while a:
		if a == '>':
			i=i+1
			b=b+'\n'
			NuOfBlank=0
			NuOfReturn=0
			NuOfSep=0
			outf.write(b)
			b=repr(i)+'|'
			
		elif a == '|':
			NuOfSep=NuOfSep+1
			if NuOfSep>1:
				b=b+' '
			else:
				b=b+a

		elif a == ' ':
			if NuOfBlank == 0:
				b=b+'|'
				NuOfBlank=1
			else:
				b=b+a
		elif a == '\n':
			if NuOfReturn ==0:
				b=b+'|'
				NuOfReturn=1
			else:
				b=b

		elif a == '\r':
			b=b

		else:
			b=b+a
		
		a=inf.read(1)
		
	
	outf.write(b)

