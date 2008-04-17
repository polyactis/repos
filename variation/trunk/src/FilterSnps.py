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
	--maxMissing=...            maximum allowed missing percentage.
	--monomorphic               filter monomorphic SNPs.
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
	
	long_options_list = ["maxError=", "comparisonFile=", "maxMissing=", "monomorphic", "delim=", "missingval=", "withArrayId=", "debug", "report", "help"]
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
	

	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)

	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2
	
	import dataParsers
	snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	#Filtering monomorphic
	if monomorphic:
		for snpsd in snpsds:
			snpsd.filterMonoMorphicSnps()
	
	#Filtering missing values
	if maxMissing<1.0 and maxMissing>=0.0:
		numAccessions = len(snpsds[0].accessions)
		for snpsd in snpsds:
			snpsd.filterMissingSnps(int(maxMissing*numAccessions))

	if comparisonFile and maxError<1.0:
		snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
		for i in range(0,len(snpsds)):
			snpsds[i].filterBadSnps(snpsds2[i],maxError)
		
	import snpsdata
	snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)


if __name__ == '__main__':
	_run_()


