#!/usr/bin/python

import sys,psycopg,re,string

# Read Group information from 'group_set' 
#fetch domain information from SETdb/dom_prot 
# write into 'gp/gpX', where X denotes the group number.


def normal(string):
	blank=re.compile(r'\s')
	tmp=''
	for i in range(len(string)):
		if blank.search(string[i]):
			pass
		else:
			tmp=tmp+string[i]
	
	return tmp

def readpacc(inf):
	pacc=[]
	tmp=inf.readline()
	tmp_list=string.split(tmp,'\t')
	
	i=0
	while tmp:
		tmp_list=string.split(tmp,'\t')
		pacc.append([])
		pacc[i].append(int(tmp_list[0]))
		pacc[i].append(normal(tmp_list[1]))
		i=i+1
		tmp=inf.readline()
		
	return pacc

def output(outf,row):
	for i in range(len(row)):
		outf.write(row[i][0]+'\t')
		outf.write(repr(row[i][1])+'\t')
		outf.write(repr(row[i][2])+'\t')
		outf.write(row[i][3]+'\n')
		#outf.write('*'.join(tmp)+'\n')
	
def run(paccs,outf):
	paccs=paccs
	outf=outf
	conn=psycopg.connect("dbname=SETdb")
	curs=conn.cursor()
	
#	if len(sys.argv)>1:
#		infname=sys.argv[1]
#		if len(sys.argv)>2:
#			outfname=sys.argv[2]
#		else:
#			outfname=raw_input("Please specify the output filename:\n")
#	else:
#		infname=raw_input("please enter the files containing paccs:\n")
#		outfname=raw_input("Please specify the output filename:\n")
#	
#	inf=open(infname,'r')
#	outf=open(outfname,'w')
	
	
	for i in range(len(paccs)):
		curs.execute("select * from dom_prot where pname=%s order by start",(paccs[i],))
		
		row=curs.fetchall()
		if row:
			output(outf,row)
	
#	outf.close()
#	inf.close()

def group_return(pacc,i):
	#pacc=pacc
	#i=i
	paccs=[]
	for j in range(len(pacc)):
		if pacc[j][0]==i:
			paccs.append(pacc[j][1])
	return paccs

if __name__=='__main__':
	inf=open('group_set','r')
	pacc=readpacc(inf)
	GroupsNo=pacc[len(pacc)-1][0]
	for i in range(GroupsNo):
		paccs=group_return(pacc,i+1)
		outf=open('gp/gp'+repr(i+1),'w')
		run(paccs,outf)
