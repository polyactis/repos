#!/usr/bin/env python2.5
"""
Usage: FilterSnps.py [OPTIONS] -o OUTPUT_FILE INPUT_FILE

Option:

        -o ...,	output file
	-d ..., --delim=...         default is \", \"      
        -m ..., --missingval=...    default is \"NA\"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info, 2 if comparison file also.
	--maxError=...              maximum allowed error percentage (requires a comparison file).
        --comparisonFile=...        a file which is used to claculate the SNPs error.
        --callProbFile=...          a csv file with call probabilities.
	--maxMissing=...            maximum allowed missing percentage.
	--monomorphic               filter monomorphic SNPs.
	--onlyBinary                filter all SNPs which are not bimorphic (monomorphic, tertiary and quaternary alleled SNPs are removed).
	--minCallProb=...           minimum allowed call probability. SNPs below the threshold are set to NA. (Requires a callProbFile file.)
	--minMAF=...                minimum allowed minor allele frequency.
	--output01Format            Output SNPs data in 01 format.  (Discards all SNPs that are not bimorphic.)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help
	
	--filterRegion=chr,start,stop	Retrieve only SNPs within region...

Examples:
	FilterSnps.py --maxMissing=0.5 -o /tmp/2010_filtered.csv 2010.csv
	
Description:
	Filter a csv formatted file with respect to various criteria.  
	Note that this script only removes SNPs not accessions. 

	Requires MySQLdb to be installed, as well as util.py, rfun.py, snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback
import snpsdata
import dataParsers

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["maxError=", "comparisonFile=", "maxMissing=", "monomorphic", "onlyBinary", "delim=", 
						"missingval=", "withArrayId=", "callProbFile=", "minMAF=", "minCallProb=", "debug", 
						"report", "help", "output01Format", "filterRegion="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	inputFile = args[0]
	output_fname = None
	delim = ","
	missingVal = "NA"
	comparisonFile = None
	maxMissing = 1.0
	maxError = 1.0
	monomorphic = False
	debug = None
	report = None
	help = 0
	withArrayIds = 0
	minCallProb=None
	minMAF=None
	callProbFile = None
	onlyBinary = False
	output01Format = False
	filterRegion = False
	startPos = None
	endPos = None
	chromosome = None
	chromosomes=[1,2,3,4,5]
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("--comparisonFile"):
			comparisonFile = arg
		elif opt in ("--maxError"):
			maxError = float(arg)
		elif opt in ("--maxMissing"):
			maxMissing = float(arg)
		elif opt in ("--minCallProb"):
			minCallProb = float(arg)
		elif opt in ("--minMAF"):
			minMAF = float(arg)
		elif opt in ("--callProbFile"):
			callProbFile = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("--monomorphic"):
			monomorphic = True
		elif opt in ("--onlyBinary"):
			onlyBinary = True
		elif opt in ("--output01Format"):
			output01Format = True
		elif opt in ("--filterRegion"):
			filterRegion = True
			region = arg.split(",")
			region = map(int,region)
			chromosome = region[0]
			startPos = region[1]
			endPos = region[2]
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)

	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2

	if callProbFile and minCallProb:
		#Read prob file into SNPsdatas.
		#snpsds = dataParsers.parseCSVDataWithCallProb(inputFile, callProbFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
		pass
	else:
		snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	
        #Filtering monomorphic
	if monomorphic:
		print "Filtering monomorphic SNPs"
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	if onlyBinary or output01Format:
		print "Filtering non-binary SNPs"
		for snpsd in snpsds:
			print "Removed", str(snpsd.onlyBinarySnps()),"Snps"

	#Filtering missing values
	if maxMissing<1.0 and maxMissing>=0.0:
		print "Filtering SNPs with missing values"
		numAccessions = len(snpsds[0].accessions)
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMissingSnps(int(maxMissing*numAccessions))),"Snps"

	#Filtering bad SNPs
	if comparisonFile and maxError<1.0:
		print "Filtering erroneous SNPs, with maxError=",maxError
		snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
		for i in range(0,len(snpsds)):
			snpsds[i].filterBadSnps(snpsds2[i],maxError)
			

	if minMAF:
		print "Removing SNPs withe MAF <",minMAF
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMinMAF(minMAF)),"Snps"

	#Output specific region..
	if filterRegion:
		chromosomes = [chromosome]
		snpsd = snpsds[chromosome-1]
		snpsd.filterRegion(startPos,endPos)
		snpsds = [snpsd]
		
		
	#Converting lousy calls to NAs
	if callProbFile and minCallProb:
		print "Converting base calls with call prob. lower than",minCallProb,"to NAs"
		#To avoid memory problems, the file/data is processed one line at a time.
		gInFile = open(inputFile,"r")
		pInFile = open(callProbFile,"r")
		outFile = open(output_fname,"w")
		if withArrayIds==2:
			gline = gInFile.readline()
			outFile.write(gline)
			pInFile.readline()
		gline = gInFile.readline()
		outFile.write(gline)
		pInFile.readline()
		i = 0
		totalCount = 0.0
		convertedCount = 0.0 
		
		while(1):
			i += 1
			gline = gInFile.readline()
			pline = pInFile.readline()
			#print gline
			if gline and pline:
				snp = gline.strip().split(delim) 
				probs = pline.strip().split(delim)
				probs = map(float,probs)
				newSNP = []
				totalCount += len(snp)
				for (nt,prob) in zip(snp,probs):
					if prob>minCallProb:
						newSNP.append(nt)
						convertedCount += 1.0
					else:
						newSNP.append('NA')
				outFile.write(delim.join(newSNP)+"\n")
			else:
				print i,gline,pline		
				break
			
			if i%10000==0:
				print i
		print i
		gInFile.close()
		pInFile.close()
		outFile.close()		
		print "Fraction converted =",convertedCount/totalCount
		
	else:
		if output01Format:
			snpsds01format = []
			for snpsd in snpsds:
				snpsds01format.append(snpsd.getSnpsData(missingVal=missingVal))
			#FINISH
			snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds01format,chromosomes=chromosomes, deliminator=delim, missingVal = missingVal, withArrayIds = waid1)
		else:
			snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=chromosomes, deliminator=delim, missingVal = missingVal, withArrayIds = waid1)


def filterByError(snpsds,comparisonSnpsds,maxError):
	#Filtering bad SNPs
	sys.stderr.write("Filtering erroneous SNPs, with maxError=%s ... \n"%maxError)
	for i in range(0,len(snpsds)):
		snpsds[i].filterBadSnps(comparisonSnpsds[i],maxError)
	return snpsds

def filterByNA(snpsds,maxMissing):
        #Filtering bad SNPs
	sys.stderr.write("Filtering SNPs with missing values ...\n")
	numAccessions = len(snpsds[0].accessions)
	for snpsd in snpsds:
		sys.stderr.write("Removed " + str(snpsd.filterMissingSnps(int(maxMissing*numAccessions))) + " Snps.\n")
	return snpsds

def filterMonomorphic(snpsds):
	#Filtering monomorphic
	sys.stderr.write("Filtering monomorphic SNPs ... \n")
	for snpsd in snpsds:
		sys.stderr.write("Removed" + str(snpsd.filterMonoMorphicSnps()) + "Snps.\n")
	return snpsds

def _test1_():
	import dataParsers
	snpsds = dataParsers.parseCSVData("2010_v3.csv")
	#snpsds = dataParsers.parseCSVData("250K_m3.csv",withArrayIds=1)
	#comparisonSnpsds = dataParsers.parseCSVData("2010_v3.csv")
	filterMonomorphic(snpsds)


if __name__ == '__main__':
	#_test1_()
	_run_()


