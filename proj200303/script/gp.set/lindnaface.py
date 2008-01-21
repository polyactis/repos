#!/usr/bin/python

def infile(f):
	j=0
	a=''
	grouplist=[]
	tmp=f.readline()
	
	while tmp:
		grouplist.append([])
		for i in range(len(tmp)-1): 
			if tmp[i]=='*':
				grouplist[j].append(a)
				a=''
			else:
				a=a+tmp[i]
		
		grouplist[j].append(a)
		a=''
		tmp=f.readline()
		
		j=j+1
		i=0
	return grouplist


def mark(grouplist):
	j=0
	m=[]
	m.append(0)
	for i in range(len(grouplist)-1):
		if grouplist[i][0]==grouplist[i+1][0]:
			m[j]=m[j]+1
		else:
			j=j+1
			m.append(0)
	
	x=[]
	x.append(-1)
	
	for i in range(len(m)):
		x.append(0)
		x[i+1]=m[i]+x[i]+1
	return x


def sort(grouplist,begin,end):
	h=end
	while h>begin:
		p=begin
		for i in range(begin,h):
			if int(grouplist[i][1]) > int(grouplist[i+1][1]):
				t=grouplist[i]
				grouplist[i]=grouplist[i+1]
				grouplist[i+1]=t
				p=i
		h=p
		

	
def sort2(grouplist):
	h=len(grouplist)-1
	while h>0:
		p=0
		for i in range(0,h):
			if int(grouplist[i][1]) > int(grouplist[i+1][1]) and grouplist[i][0]==grouplist[i+1][0]:
				t=grouplist[i]
				grouplist[i]=grouplist[i+1]
				grouplist[i+1]=t
				p=i
		h=h-1

def labelwrite(f,gl,i):
	dcolor={'Nuclear protein SET':'0',
	'SET':'0',
	'Nuclear protein Zn2+-binding':'1',
	'Pre-SET':'1',
	'YDG_SRA':'3',
	'Myb DNA-binding domain':'4',
	'A+T-hook':'7',
	'AT_hook':'7',
	'Chromo domain':'5',
	'TPR':'9',
	'DUF260':'13',
	'PHD':'6',
	'PWWP':'11',
	'zf-CXXC':'8',
	'zf-MYND':'8',
	'zf-C2H2':'8',
	'AWS':'10',
	'SET-related region':'12',
	'  ':'15'}

	f.write('label\n')
	f.write('Block\t')
	f.write(gl[i][1]+'\t'+gl[i][2]+'\t')
	if dcolor.has_key(gl[i][3]):
		f.write(dcolor[gl[i][3]]+'\t')
	else:
		f.write('14\t')
	f.write('H\n')
	f.write(gl[i][3]+'\n')
	f.write('endlabel\n')
	


def outfile(f,gl):
	
	maxlen=0
	for i in range(len(gl)):
		if int(gl[i][2])>maxlen:
			maxlen=int(gl[i][2])
	
	f.write('Start\t1\n')
	f.write('End\t'+repr(maxlen)+'\n')
	
	f.write('group\n')
	f.write(gl[0][0]+'\n')
	f.write('label\n')
	f.write('Block\t1\t2\t15\n')
	f.write('endlabel\n')
	for i in range(len(gl)):
		if i == 0:
			labelwrite(f,gl,i)	
			
		elif gl[i][0]==gl[i-1][0]:
			labelwrite(f,gl,i)	
		else:
			f.write('endgroup\n')
			f.write('\n')
			f.write('group\n')
			f.write(gl[i][0]+'\n')
			f.write('label\n')
			f.write('Block\t1\t2\t15\n')
			f.write('endlabel\n')
			
			labelwrite(f,gl,i)	
			
			
	f.write('endgroup\n')	
