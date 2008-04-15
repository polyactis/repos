#!/usr/bin/env python2.5
"""
Usage: Output250KSNPs.py [OPTIONS] -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, (banyan.usc.edu is default).
	-u ..., --user=...	the username, (otherwise it will ask for it).
	-p ..., --passwd=...	the password, (otherwise it will ask for it).
        -o ...,	output file
        -t ..., --method=...    Currently not implemented (since there is still only one method)
        -d ..., --delim=...         default is ", "      
        -m ..., --missingval=...    default is "NA"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info.
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Output250SNPs.py -o /tmp/250K.csv
	
Description:
	Output 250 data, in a csv format (or with another deliminator separated format).

	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "user=", "passwd=", "method=", "delim=", "missingval=", "withArrayId=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:u:p:o:t:d:m:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = 'banyan.usc.edu'
	user = None
	passwd = None
	output_fname = None
	method = None #TO BE IMPLEMENTED
	delim = ", "
	missingVal = "NA"
	debug = None
	report = None
	help = 0
	withArrayId = False

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
			version = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-a","--withArrayId"):
			withArrayId = bool(arg)
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


	import dataParsers
	snpsds = dataParsers.get250KDataFromDb(host=hostname,chromosomes=[1,2,3,4,5],  user=user, passwd=passwd, withArrayIds=withArrayId)
	
	import snpsdata
	snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds=withArrayId)


if __name__ == '__main__':
	_run_()

