#!/usr/bin/python2.2

import sys, psycopg

conn=psycopg.connect("dbname=SETdb")
curs=conn.cursor()

create=[]
cp=[]

for i in range(1,8):
	createtmp="create temp table gp"+repr(i)+"tmp as select d.* from dom_prot d, gp g where g.pname=d.pname and g.gp="+repr(i)
	cptmp="copy gp"+repr(i)+"tmp to '/home/gtkusr/script/gp.set/gp/gp"+repr(i)+"' delimiters '*'"
	create.append(createtmp)
	cp.append(cptmp)

for i in range(7):
	
	print "create the "+repr(i)+"th table"
	
	curs.execute(create[i])
	print "copy the "+repr(i)+"th table"
	curs.execute(cp[i])

#conn.commit()
	
