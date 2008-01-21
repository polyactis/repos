#!/usr/bin/python

import os,re,sys
from Extract import *
from stat import *


f= os.listdir(sys.argv[1])
print f
for i in range(0,len(f)):

	pathname=sys.argv[1]+'/'+f[i]
	print pathname
	gunzipL=['gunzip',pathname]
	gpara=os.spawnvp(os.P_WAIT,'gunzip',gunzipL)

	newpathname=pathname[0:len(pathname)-3]
	
	dir1='/home/awdong/script/'	
	outname=dir1+f[i][0:len(f[i])-3]
	instance=extract(newpathname,outname)
	instance.run()
	#extractL=['Extract.py',newpathname,outname]
	#epara=os.spawnvp(os.P_WAIT,'./Extract.py',extractL)
	
	gunzipL1=['gzip',newpathname]
	gpara1=os.spawnvp(os.P_WAIT,'gzip',gunzipL1)
	
	
