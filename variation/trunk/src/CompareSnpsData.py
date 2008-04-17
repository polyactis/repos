#!/usr/bin/env python2.5
"""
Usage: CompareSNPsData.py [OPTIONS] datafile1 datafile2

Option:
        -o ...,	--rFile=...             Output file (R format), a script that generates several graphs.
        -d ..., --delim=...             default is ", "
        -m ..., --missingval=...        default is "NA"
	-a ..., --withArrayId=...       0 for no array ID info (default), 1 if first file has array ID info, 2 if both have.
	-s ..., --statFile=...          Writes several comparison statistics to a file.
	-v, --verbose   Prints various comparison statistics to the screen.
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
	
	long_options_list = ["rFile=", "delim=", "missingval=", "statFile=", "debug", "report", "help", "withArrayId="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:s:vbrh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	if len(args)!=2:
            raise Exception("Number of arguments isn't correct.")
        inputFile1 = args[0]
	inputFile2 = args[1]
        rFile = None
        statFile = None
	verbose = False
	delim = ", "
	missingVal = "NA"
	debug = None
	report = None
	withArrayIds = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-o", "--rFile"):
			rFile = arg
		elif opt in ("-s", "--statFile"):
			statFile = arg
		elif opt in ("-t","--method"):
			version = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("-v", "--verbose"):
			verbose = True
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2
	
	snpsds1 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	snpsds2 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
        

        if len(snpsds1) != len(snpsds2):
		raise Exception("Unequal number of chromosomes in files.")
        
        res = []
        for i in range(0,len(snpsds1)):
		res.append(snpsds1[i].compareWith(snpsds2[i],withArrayIds=withArrayIds))
        
	import rfun
	totalCommonPos = 0
	totalAccessionCounts = [0]*len(res[0][2])
	accCallRate = [[0]*len(res[0][2]),[0]*len(res[0][2])]
	accErrorRate = [0]*len(res[0][2])
	
	statstr = "#Common SNPs positions:\n"
	rstr = "#Snps error rates\n"
	rstr = "par(mfrow=c(5,1));\n"
	snpsErrorRate = []
	for i in range(0,len(res)):
		r = res[i]
		snpsErrorRate +=r[1]
		totalCommonPos += len(r[0])
		statstr += "Chr. "+str(i+1)+":\n"
		statstr += str(r[0])+"\n"
		xname = "commonPos_ch"+str(i+1)
		ynames = ["errorRates_ch"+str(i+1)]
		rstr += rfun.plotOverlayingVectors(r[0],[r[1]],xlab="Position, chr. "+str(i+1),ylab="Error (red)",type="b",xname=xname,ynames=ynames)+"\n\n"
		for j in range(0,len(r[2])):
			totalAccessionCounts[j] += r[6][j]
			accCallRate[0][j]+=r[4][0][j]*float(len(r[0]))
			accCallRate[1][j]+=r[4][1][j]*float(len(r[0]))
			accErrorRate[j]+=r[3][j]*float(r[6][j])

	statstr += "#Number of common SNPs positions:\n"
	statstr += str(totalCommonPos)+"\n"
	statstr += "#SNPs errors:\n"
	for i in range(0,len(res)):
		r = res[i]
		statstr += "Chr. "+str(i+1)+":\n"
		statstr += str(r[1])+"\n"


	statstr += "#Average Snp Error:\n"
	statstr += str(sum(snpsErrorRate)/float(len(snpsErrorRate)))+"\n"
	
	statstr += "#Commmon accessions:\n"
	statstr += str(res[0][2])+'\n'
	statstr += "#Number of commmon accessions:\n"
	statstr += str(len(res[0][2]))+'\n'
	if withArrayIds:
		statstr += "#ArrayIds:\n"
		statstr += str(res[0][5])+'\n'


        if not verbose:
		print "In all",len(res[0][2]),"common accessions found"
		print "In all",totalCommonPos,"common snps found"
		print "Average Snp Error:",sum(snpsErrorRate)/float(len(snpsErrorRate))


	for i in range(0,len(res)):
		r = res[i]
		xname = "commonPos_ch"+str(i+1)
		ynames = ["missingRates1_ch"+str(i+1),"missingRates2_ch"+str(i+1)]
		rstr += rfun.plotOverlayingVectors(r[0],[r[7][0],r[7][1]],xlab="Position, chr. "+str(i+1),ylab="Missing (red,green)",type="b",xname=xname,ynames=ynames)+"\n\n"

	for i in range(0,len(r[2])):
		accCallRate[0][i]=accCallRate[0][i]/totalCommonPos
		accCallRate[1][i]=accCallRate[1][i]/totalCommonPos
		accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])

	accErrAndID = []
	accMissAndID = [[],[]]
	if withArrayIds:
		for i in range(0,len(r[2])):
			accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))
			accMissAndID[0].append((accCallRate[0][i], r[2][i], r[5][i]))
			accMissAndID[1].append((accCallRate[1][i], r[2][i], r[5][i]))
	else:
		for i in range(0,len(r[2])):
			accErrAndID.append((accErrorRate[i], r[2][i]))
			accMissAndID[0].append((accCallRate[0][i], r[2][i]))
			accMissAndID[1].append((accCallRate[1][i], r[2][i]))
	accErrAndID.sort(reverse=True)
	accMissAndID[0].sort(reverse=True)
	accMissAndID[1].sort(reverse=True)
	statstr += "#Sorted list, based on error rates (Error rate, ecotype id, array id):\n"
	statstr += str(accErrAndID)+'\n'
	statstr += "#Sorted list, based on missing rates of 1st file, (Missing rate, ecotype id, array id):\n"
	statstr += str(accMissAndID[0])+'\n'
	statstr += "#Sorted list, based on missing rates of 2nd file, (Missing rate, ecotype id, array id):\n"
	statstr += str(accMissAndID[1])+'\n'
 
	"""
	print "Sorted list, based on error rates: ",accErrAndID,'\n'
        accMissAndID[0].sort(reverse=True)
        print "Sorted list, based on missing rates (1st file): ",accMissAndID[0],'\n'
        accMissAndID[1].sort(reverse=True)
        print "Sorted list, based on missing rates (2nd file): ",accMissAndID[1],'\n'
	"""
		

	rstr += 'accessions<-c("'+str(r[2][0])+"_ai"+str(r[5][0])+'"'
	for i in range(1, len(r[2])):
		rstr += ',"'+str(r[2][i])+"_ai"+str(r[5][i])+'"'			
	rstr +=")\n"
	rstr += rfun.plotVectors(accCallRate[0],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate",xname="accMissingRate1",ynames=["accErrorRate"])
	rstr += "text(accMissingRate1+0.0045,accErrorRate-0.0004,accessions)\n\n"
	rstr += rfun.plotVectors(accCallRate[1],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate",xname="accMissingRate2",ynames=["accErrorRate"])
	rstr += "text(accMissingRate2+0.0045,accErrorRate-0.0004,accessions)\n\n"

	f = open(rFile,"w")
	f.write(rstr)
	f.close()
	if verbose:
		print statstr
	if statFile:
		f = open(statFile,"w")
		f.write(statstr)
		f.close()
	
if __name__ == '__main__':
	_run_()

