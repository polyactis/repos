#!/usr/bin/env python2.5
"""
Usage: CompareSNPsData.py [OPTIONS] datafile1 datafile2

Option:
        -o ...,	--outputFile=...        Output file (R format)
        -p ..., --priority=...          1 for giving first file priority, 2 for second file (if SNPs values disagree). (1 is default)
        -d ..., --delim=...             default is ", "      
        -m ..., --missingval=...        default is "NA"
	-a ..., --withArrayId=...       0 for no array ID info (default), 1 if first file has array ID info, 2 if both have.
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	 CompareSNPsData.py datafile1 datafile2
	
Description:
	Compares two datafiles (of csv or deliminator separated format).

	Requires MySQLdb to be installed, and uses snpsdata.py, dataParsers.py, rfun.py, and util.py.
"""

import sys, getopt, traceback
import dataParsers
import snpsdata

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["outputFile=", "priority=", "delim=", "missingval=", "debug", "report", "help", "withArrayId="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:p:d:m:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	if len(args)!=2:
            raise Exception("Number of arguments isn't correct.")
        inputFile1 = args[0]
	inputFile2 = args[1]
	priority = 1
        output_fname = None
	method = None 
	delim = ", "
	missingVal = "NA"
	debug = None
	report = None
	withArrayIds = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-p", "--priority"):
			priority = int(arg)
		elif opt in ("-s", "--maxMissingSnps"):
			maxMissingSnps = int(arg)
		elif opt in ("-o", "--outputFile"):
			output_fname = arg
		elif opt in ("-t","--method"):
			version = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2
	
	if priority == 1:
		snpsds1 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
		snpsds2 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
        elif priority == 2:
		snpsds1 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
		snpsds2 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
        else:
		raise Exception("Priority option is wrong.")

        if len(snpsds1) != len(snpsds2):
		raise Exception("Unequal number of chromosomes.")
        
        res = []
        for i in range(0,len(snpsds1)):
		res.append(snpsds1[i].compareWith(snpsds2[i],withArrayIds=withArrayIds))
        
	import rfun
	totalCommonPos = 0
	totalAccessionCounts = [0]*len(res[0][2])
	accCallRate = [[0]*len(res[0][2]),[0]*len(res[0][2])]
	accErrorRate = [0]*len(res[0][2])
	


	for r in res:
		totalCommonPos += len(r[0])
		print rfun.plotOverlayingVectors(r[0],[r[1],r[7][0]],xlab="Position",ylab="Error rate",type="b")
		#print rfun.plotVectors(r[0],[r[1]],xlab="Position",ylab="Error rate",type="b")
		#print rfun.plotVectors(r[0],[r[7][0]],xlab="Position",ylab="Missing rate",type="b")
		for i in range(0,len(r[2])):
			totalAccessionCounts[i] += r[6][i]
			accCallRate[0][i]+=r[4][0][i]*float(len(r[0]))
			accCallRate[1][i]+=r[4][1][i]*float(len(r[0]))
			accErrorRate[i]+=r[3][i]*float(r[6][i])


	for i in range(0,len(r[2])):
		accCallRate[0][i]=accCallRate[0][i]/totalCommonPos
		accCallRate[1][i]=accCallRate[1][i]/totalCommonPos
		accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])


	print rfun.plotVectors(accCallRate[0],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate")
	print rfun.plotVectors(accCallRate[1],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate")


	rstr = 'z<-c("'+str(r[2][0])+"_ai"+str(r[5][0])+'"'
	for i in range(1, len(r[2])):
		rstr += ',"'+str(r[2][i])+"_ai"+str(r[5][i])+'"'			
	print rstr+")"
	print "text(x+0.0045,y-0.0004,z)"
			
			

if __name__ == '__main__':
	_run_()

