#!/usr/bin/python

#Usage: mirdbfastaOut.py 1stfile 2ndfile table domain 
#Read pnames from the 1stfile (argv[1]).
#Fetch full-length or specific domain (argv[4]) sequences from a certain table (argv[3]) into the 2ndfile (argv[2]).

import sys,psycopg,re

def normal(string):
	blank=re.compile(r'\s')
	tmp=''
	for i in range(len(string)):
		if blank.search(string[i]):
			pass
		else:
			tmp=tmp+string[i]
	
	return tmp

def readpaccs(inf):
	paccs=[]
	tmp=inf.readline()
	tmp=normal(tmp)
	
	while tmp:
		tmp=normal(tmp)
		paccs.append(tmp)
		tmp=inf.readline()
		
	return paccs

def fastaOut(outf,rows):
	for i in range(len(rows)):
		start=int(rows[i][2])
		tail=int(rows[i][3])
		outf.write('>'+rows[i][0]+'\n')
		outf.write(rows[i][1][start-1:tail]+'\n')
	
def run(domain='', table='mir'):
	conn=psycopg.connect("dbname=mirdb")
	curs=conn.cursor()
	
	if len(sys.argv)>1:
		infname=sys.argv[1]
		if len(sys.argv)>2:
			outfname=sys.argv[2]
		else:
			outfname=raw_input("Please specify the output filename:\n")
	else:
		infname=raw_input("please enter the files containing paccs:\n")
		outfname=raw_input("Please specify the output filename:\n")
	
	inf=open(infname,'r')
	outf=open(outfname,'w')
	
	paccs=readpaccs(inf)
	
	rows=[]
	noofrows=0
	if domain:	
		for i in range(len(paccs)):
			curs.execute("select m.acc,m.sequence,s.start,s.tail from mir m, segment s where m.acc=s.acc and s.segname=%s and m.acc=%s",(domain,paccs[i],))
			
			row=curs.fetchall()
			if row:
				rows.append([])
				j=noofrows
				rows[j].append(row[0][0])
				rows[j].append(row[0][1])
				rows[j].append(row[0][2])
				rows[j].append(row[0][3])
				noofrows=noofrows+1
	else:
		for i in range(len(paccs)):
			if table=='mir':
				curs.execute("select acc,sequence from mir where acc=%s",(paccs[i],))

			elif table=='pr_chromdb':
				curs.execute("select pacc,sequence from "+table+" where pacc=%s",(paccs[i],))


			elif table=='temp_pr':
				curs.execute("select pacc,sequence from "+table+" where pacc=%s",(paccs[i],))
				
			row=curs.fetchall()
			if row:
				rows.append([])
				j=noofrows
				rows[j].append(row[0][0])
				rows[j].append(row[0][1])
				rows[j].append(1)
				rows[j].append(len(row[0][1]))
				noofrows=noofrows+1
		
	fastaOut(outf,rows)
	
	outf.close()
	inf.close()


if __name__=='__main__':
	if len(sys.argv)==5:
		run(domain=sys.argv[4],table=sys.argv[3])
	elif len(sys.argv)==4:
		run(table=sys.argv[3])
		
	else:
		run()
