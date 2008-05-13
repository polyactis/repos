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
	delim = ","
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
        naRate1 = 0
        naRate2 = 0
	numSNPs1 = 0
	numSNPs2 = 0
	for i in range(0,len(snpsds1)):
		res.append(snpsds1[i].compareWith(snpsds2[i],withArrayIds=withArrayIds))
		naRate1 += snpsds1[i].countMissingSnps()*len(snpsds1[i].positions)
		naRate2 += snpsds2[i].countMissingSnps()*len(snpsds2[i].positions)
		numSNPs1 += len(snpsds1[i].positions)
		numSNPs2 += len(snpsds2[i].positions)
		
		
	naRate1 = naRate1/float(numSNPs1)
	naRate2 = naRate2/float(numSNPs2)

	import rfun
	totalCommonPos = 0
	totalPos = [0,0]
	commonAccessions = res[0][2]
	totalAccessionCounts = [0]*len(commonAccessions)
	accOverlappingCallRate = [[0]*len(commonAccessions),[0]*len(commonAccessions)]
	accCallRate = [[0]*len(commonAccessions),[0]*len(commonAccessions)]
	accErrorRate = [0]*len(commonAccessions)
	
	statstr = "#Common SNPs positions:\n"
	rstr = "#Snps error rates\n"
	rstr = "par(mfrow=c(5,1));\n"
	snpsErrorRate = []
	for i in range(0,len(res)): #for all chromosomes
		r = res[i]
		snpsErrorRate +=r[1]
		totalCommonPos += len(r[0])
		totalPos[0] += len(snpsds1[i].positions)
		totalPos[1] += len(snpsds2[i].positions)
		statstr += "Chr. "+str(i+1)+":\n"
		statstr += str(r[0])+"\n"
		xname = "commonPos_ch"+str(i+1)
		ynames = ["errorRates_ch"+str(i+1)]
		rstr += rfun.plotOverlayingVectors(r[0],[r[1]],xlab="Position, chr. "+str(i+1),ylab="Error (red)",type="b",xname=xname,ynames=ynames)+"\n\n"
		for j in range(0,len(commonAccessions)):
			totalAccessionCounts[j] += r[6][j]
			accOverlappingCallRate[0][j]+=r[4][0][j]*float(len(r[0]))
			accOverlappingCallRate[1][j]+=r[4][1][j]*float(len(r[0]))
			accCallRate[0][j]+=r[8][0][j]
			accCallRate[1][j]+=r[8][1][j]
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
	statstr += str(commonAccessions)+'\n'
	statstr += "#Number of commmon accessions:\n"
	statstr += str(len(commonAccessions))+'\n'
	statstr += "#Number of accessions (1):\n"
	statstr += str(len(snpsds1[0].accessions))+'\n'
	statstr += "#Number of accessions (2):\n"
	statstr += str(len(snpsds2[0].accessions))+'\n'

	if withArrayIds:
		commonArrayIds = res[0][5]
		statstr += "#ArrayIds:\n"
		statstr += str(commonArrayIds)+'\n'

        if not verbose:
		print "In all",len(commonAccessions),"common accessions found"
		print "In all",totalCommonPos,"common snps found"
		print "Average Snp Error:",sum(snpsErrorRate)/float(len(snpsErrorRate))
		print "NA rate (1) =",naRate1
		print "NA rate (2) =",naRate2


	for i in range(0,len(res)):
		r = res[i]
		xname = "commonPos_ch"+str(i+1)
		ynames = ["missingRates1_ch"+str(i+1),"missingRates2_ch"+str(i+1)]
		rstr += rfun.plotOverlayingVectors(r[0],[r[7][0],r[7][1]],xlab="Position, chr. "+str(i+1),ylab="Missing (red,green)",type="b",xname=xname,ynames=ynames)+"\n\n"

	for i in range(0,len(commonAccessions)):
		accOverlappingCallRate[0][i]=accOverlappingCallRate[0][i]/float(totalCommonPos)
		accOverlappingCallRate[1][i]=accOverlappingCallRate[1][i]/float(totalCommonPos)
		accCallRate[0][i]=accCallRate[0][i]/float(totalPos[0])
		accCallRate[1][i]=accCallRate[1][i]/float(totalPos[1])
		accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])

	accErrAndID = []
	accMissAndID = [[],[]]
	accOverlMissAndID = [[],[]]
	if withArrayIds:
		for i in range(0,len(commonAccessions)):
			accErrAndID.append((accErrorRate[i], commonAccessions[i], commonArrayIds[i]))
			accMissAndID[0].append((accCallRate[0][i], commonAccessions[i], commonArrayIds[i]))
			accMissAndID[1].append((accCallRate[1][i], commonAccessions[i], commonArrayIds[i]))
		accOverlMissAndID[0] = zip(accOverlappingCallRate[0],commonAccessions,commonArrayIds)
		accOverlMissAndID[1] = zip(accOverlappingCallRate[1],commonAccessions,commonArrayIds)
	else:
		for i in range(0,len(commonAccessions)):
			accErrAndID.append((accErrorRate[i], commonAccessions[i]))
			accMissAndID[0].append((accCallRate[0][i], commonAccessions[i]))
			accMissAndID[1].append((accCallRate[1][i], commonAccessions[i]))
		accOverlMissAndID[0] = zip(accOverlappingCallRate[0],commonAccessions)
		accOverlMissAndID[1] = zip(accOverlappingCallRate[1],commonAccessions)
	accErrAndID.sort()	#05/10/08 yh. sort(reverse=True) is not available in python 2.3
	accErrAndID.reverse()
	accMissAndID[0].sort()
	accMissAndID[0].reverse()
	accOverlMissAndID[1].sort()
	accOverlMissAndID[1].reverse()
	statstr += "#Sorted list, based on error rates (Error rate, ecotype id, array id):\n"
	for t in accErrAndID:
		statstr += str(t)+'\n'
	statstr += "#Sorted list, based on missing rates of 1st file, (Missing rate, ecotype id, array id):\n"
	for t in accMissAndID[0]:
		statstr += str(t)+'\n'
	statstr += "#Sorted list, based on missing rates of 2nd file, (Missing rate, ecotype id, array id):\n"
	for t in accMissAndID[1]:
		statstr += str(t)+'\n'
	statstr += "#Sorted list, based on (overlapping positions) missing rates of 1st file, (Missing rate, ecotype id, array id):\n"
	for t in accOverlMissAndID[0]:
		statstr += str(t)+'\n'
	statstr += "#Sorted list, based on (overlapping positions) missing rates of 2nd file, (Missing rate, ecotype id, array id):\n"
	for t in accOverlMissAndID[1]:
		statstr += str(t)+'\n'
 
	"""
	print "Sorted list, based on error rates: ",accErrAndID,'\n'
        accMissAndID[0].sort(reverse=True)
        print "Sorted list, based on missing rates (1st file): ",accMissAndID[0],'\n'
        accMissAndID[1].sort(reverse=True)
        print "Sorted list, based on missing rates (2nd file): ",accMissAndID[1],'\n'
	"""
		

	if withArrayIds:
		rstr += 'accessions<-c("'+str(r[2][0])+"_ai"+str(r[5][0])+'"'
	else:
		rstr += 'accessions<-c("'+str(r[2][0])+'"'		
	for i in range(1, len(r[2])):
		if withArrayIds:
			rstr += ',"'+str(r[2][i])+"_ai"+str(r[5][i])+'"'			
		else:
			rstr += ',"'+str(r[2][i])+'"'
	rstr +=")\n"
	rstr += rfun.plotVectors(accCallRate[0],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate",xname="accMissingRate1",ynames=["accErrorRate"])
	rstr += "text(accMissingRate1+0.0045,accErrorRate-0.0004,accessions)\n\n"
	rstr += rfun.plotVectors(accCallRate[1],[accErrorRate],xlab="Accession missing value rate",ylab="Accession error rate",xname="accMissingRate2",ynames=["accErrorRate"])
	rstr += "text(accMissingRate2+0.0045,accErrorRate-0.0004,accessions)\n\n"

	if rFile:
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

