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
#	tmp_list=string.split(tmp,'\t')
	
	i=0
	while tmp:
		tmp_list=string.split(tmp,'\t')
#		pacc.append([])
#		pacc[i].append(int(tmp_list[0]))
		pacc.append(normal(tmp_list[1]))
		i=i+1
		tmp=inf.readline()
		
	return pacc

def output(outf,row):
	for i in range(len(row)):
		outf.write(repr(row[i][0])+'\t')
		outf.write(row[i][1]+'\t')
		outf.write(row[i][2]+'\t')
		outf.write(row[i][3]+'\n')
		#outf.write('*'.join(tmp)+'\n')
	
def run(pacc,outf):
	conn=psycopg.connect("dbname=SETdb")
	curs=conn.cursor()
	
	for i in range(len(pacc)):
		curs.execute("select g.gp,p.dbsrc,g.pname,p.pacc from gp g,protein_pname p where g.pname=p.pname and g.pname=%s",(pacc[i],))
		
		row=curs.fetchall()
		if row:
			output(outf,row)


if __name__=='__main__':
	inf=open('group_set','r')
	pacc=readpacc(inf)
	outf=open('gp_table','w')
	run(pacc,outf)
