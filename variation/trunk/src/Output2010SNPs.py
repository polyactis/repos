#!/usr/bin/env python2.5
"""
Usage: Output2010SNPs.py [OPTIONS] -o OUTPUT_FILE

Option:
	-z ..., --hostname=...	the hostname, (papaya.usc.edu is default).
	-u ..., --user=...	the username, (otherwise it will ask for it).
	-p ..., --passwd=...	the password, (otherwise it will ask for it).
        -o ...,	output file
        -v ..., --version=...   1, 2, or 3, (3 is default and the published version).
        -d ..., --delim=...     default is ", "      
        -m ..., --missingval=...    default is "NA"
        -a, --accname   use accession name instead of ecotype id. 
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help
	--only96	Output only the 96 accesions

Examples:
	Output2010SNPs.py -o /tmp/2010.csv
	
Description:
	Output 2010 data, in a csv format (or with another deliminator separated format).

	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "user=", "passwd=", "version=", "delim=", "missingval=", "accname", "debug", "report", "help", "only96"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:u:p:o:v:d:m:abrh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = 'papaya.usc.edu'
	user = None
	passwd = None
	output_fname = None
	version = "3"
	delim = ", "
	missingVal = "NA"
	useAccessionName = False
	debug = None
	report = None
	help = 0
	only96 = False
	
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
		elif opt in ("-v","--version"):
			version = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("-a","--accname"):
			useAccessionName = True		
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("--only96"):
			only96 = True
	

	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)


	import dataParsers
	snpsds = dataParsers.get2010DataFromDb(host=hostname,chromosomes=[1,2,3,4,5], dataVersion=version, only96accessions=only96, user=user, passwd=passwd)
	
	accDecoder=None
	if useAccessionName:
		accDecoder = dataParsers.ecotypeId2Name
	import snpsdata
	snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, accDecoder=accDecoder)


if __name__ == '__main__':
	_run_()

