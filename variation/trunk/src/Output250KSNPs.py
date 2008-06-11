#!/usr/bin/env python2.5
"""
Usage: Output250KSNPs.py [OPTIONS] -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, (papaya.usc.edu is default).
	-u ..., --user=...	the username, (otherwise it will ask for it).
	-p ..., --passwd=...	the password, (otherwise it will ask for it).
        -o ...,	output file
        -t ..., --method=...    (1 is default)
        -d ..., --delim=...         default is ", "      
        -m ..., --missingval=...    default is "NA"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info.
        --callProbFile=...,	output call probabilities in a file
	-h, --help	show this help
	--newBatch      retrieves only the new batch. (This option is volatile)

Examples:
	Output250KSNPs.py -o /tmp/250K.csv
	
Description:
	Output 250 data, in a csv format (or with another deliminator separated format).

	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["newBatch","hostname=", "user=", "passwd=", "method=", "delim=", "missingval=", "withArrayId=", "callProbFile=", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:u:p:o:t:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = 'papaya.usc.edu'
	user = None
	passwd = None
	output_fname = None
	method = 1 
	delim = ","
	missingVal = "NA"
	help = 0
	withArrayId = False
	callProbFile = None
	newBatch = False

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-t","--method"):
			method = int(arg)
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-a","--withArrayId"):
			withArrayId = bool(arg)
		elif opt in ("--callProbFile"):
			callProbFile =arg
		elif opt in ("--newBatch"):
			newBatch = True
	
	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)


	import dataParsers
	import snpsdata
	if callProbFile:
		snpsds = dataParsers.get250KDataFromDb(host=hostname,chromosomes=[1,2,3,4,5], methodId=method, user=user, passwd=passwd, withArrayIds=withArrayId, callProb=True, newBatch=newBatch)
		snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds=withArrayId, callProbFile=callProbFile)
	else:
		snpsds = dataParsers.get250KDataFromDb(host=hostname,chromosomes=[1,2,3,4,5], methodId=method, user=user, passwd=passwd, withArrayIds=withArrayId, newBatch=newBatch)
		snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds=withArrayId)


if __name__ == '__main__':
	_run_()

