"""
Script to hide unwanted tracks
"""

WantedTracks = range(1,49)+range(57,83)+[158,159]+range(161,180)+range(182,187)+[272,273,274]+range(277,284)

fi = open("/etc/apache2/gbrowse.conf/arabidopsis.conf")
fo = open("/etc/apache2/gbrowse.conf/arabidopsis2.conf","w")

delete_flag = 0

for line in fi:
    if line[0]=="[": # start of new track
	u = line.split("_")
	if len(u)==1:
	    delete_flag = 0
	elif int(u[1]) not in WantedTracks:
            delete_flag = 1
        else:
            delete_flag = 0
    if delete_flag == 0:
        fo.write(line)
fi.close()
fo.close()
