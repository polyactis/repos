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
	--minCallProb=...           minimum allowed call probability. SNPs below the threshold are set to NA. (Requires a callProbFile file.)
	--minAverageCallProb=...    minimum allowed average call probability. Bad SNPs are removed. (Requires a callProbFile file.)
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	FilterSnps.py --maxMissing=0.5 -o /tmp/2010_filtered.csv 2010.csv
	
Description:
	Filter a csv formatted file with respect to various criteria.  
	Note that this script only removes SNPs not accessions. 

	Requires MySQLdb to be installed, as well as util.py, rfun.py, snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["maxError=", "comparisonFile=", "maxMissing=", "monomorphic", "delim=", "missingval=", "withArrayId=", "callProbFile=", "minAverageCallProb=", "minCallProb=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	inputFile = args[0]
	output_fname = None
	delim = ", "
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
	minAverageCallProb=None
	callProbFile = None

	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("--comparisonFile"):
			comparisonFile = arg
		elif opt in ("--maxError"):
			maxError = float(arg)
		elif opt in ("--maxMissing"):
			maxMissing = float(arg)
		elif opt in ("--monomorphic"):
			monomorphic = True
		elif opt in ("--minCallProb"):
			minCallProb = float(arg)
		elif opt in ("--minAverageCallProb"):
			minAverageCallProb = float(arg)
		elif opt in ("--callProbFile"):
			callProbFile = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
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

	import dataParsers
	if callProbFile and (minCallProb or minAverageCallProb):
		#Read prob file into SNPsdatas.
		snpsds = dataParsers.parseCSVDataWithCallProb(inputFile, callProbFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	else:
		snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	
        #Filtering monomorphic
	if monomorphic:
		print "Filtering monomorphic SNPs"
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	
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
			
	#Converting lousy calls to NAs
	if callProbFile and minCallProb:
		print "Converting base calls with call prob. lower than",minCallProb,"to NAs"
		for i in range(0,len(snpsds)):
			print "Fraction converted =",snpsds[i].convertBadCallsToNA(minCallProb)

	import snpsdata
	snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)


if __name__ == '__main__':
	_run_()


