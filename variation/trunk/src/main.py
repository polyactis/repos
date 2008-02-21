#!/usr/bin/env python
"""
Usage: main.py [OPTIONS] -y X

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input_fname
	-o ...,	output_fname
	-y ..., --type=...	which test case should be invoked.
	-a ...,	argument 1
	-e ...,	argument 2
	-f ...,	argument 3
	-g ...,	argument 4
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	#run 1st type
	main.py -y 1 -i data_250k.tsv -o data_250k.kw.pvalue -a phenotype.csv
	
	#get help of the 1st type
	main.py -y 1 -h

Description:
	This program is an interface for other classes/types.
	
	For help regarding each specific type, turn on both -y X and -h

For -y (type):

02-14-08 1: Kruskal_Wallis.py
02-20-08 2: ProcessPhenotype.py
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt
from variation.src import Kruskal_Wallis
from variation.src import ProcessPhenotype

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:o:y:a:e:f:g:br", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	ClassDict = {1:Kruskal_Wallis.Kruskal_Wallis,
				2:ProcessPhenotype.ProcessPhenotype}
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	input_fname = None
	output_fname = None
	type = None
	argument1 = None
	argument2 = None
	argument3 = None
	argument4 = None
	help = 0
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-a",):
			argument1 = arg
		elif opt in ("-e",):
			argument2 = arg
		elif opt in ("-f",):
			argument3 = arg
		elif opt in ("-g",):
			argument4 = arg
	
	if type!=None:
		if help:
			print
			print "\tClass Type: %s"%type
			print
			print ClassDict[type].__doc__
			sys.exit(2)
		elif type==1:
			ins = ClassDict[type](input_fname, argument1, output_fname, argument2, debug=debug, report=report)
			ins.run()
		elif type==2:
			ins = ClassDict[type](hostname, dbname, schema, output_fname, argument1, argument2, debug=debug, report=report)
			ins.run()
	else:
		print __doc__
		sys.exit(2)
		