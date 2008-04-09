#!/usr/bin/env python
"""
Usage: main.py [OPTIONS] -y X

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-i ...,	input
	-o ...,	output
	-y ..., --type=...	which test case should be invoked.
	-a ...,	argument 1
	-e ...,	argument 2
	-f ...,	argument 3
	-g ...,	argument 4
	-j ...,	argument 5
	-c,	commit db transaction
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
02-28-08 3: Array2DB_250k.py
03-11-08 4: LinkArrayId2EcotypeId.py
04-08-08 5: DB_250k2Array.py
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
from variation.src import Array2DB_250k
from variation.src import LinkArrayId2EcotypeId
from variation.src import DB_250k2Array

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:y:a:e:f:g:j:cbr", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	ClassDict = {1:Kruskal_Wallis.Kruskal_Wallis,
				2:ProcessPhenotype.ProcessPhenotype,
				3:Array2DB_250k.Array2DB_250k,
				4:LinkArrayId2EcotypeId.LinkArrayId2EcotypeId,
				5:DB_250k2Array.DB_250k2Array}
	
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
	argument5 = None
	help = 0
	commit = 0
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
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-a",):
			argument1 = arg
		elif opt in ("-e",):
			argument2 = arg
		elif opt in ("-f",):
			argument3 = arg
		elif opt in ("-g",):
			argument4 = arg
		elif opt in ("-j",):
			argument5 = arg
	
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
			ins = ClassDict[type](hostname=hostname, dbname=dbname, schema=schema, output_fname=output_fname, raw_phenotype_table=argument1, \
								experiment_table=argument2, phenotype_table=argument3, phenotype_avg_table=argument4,\
								method_table=argument5, commit=commit, debug=debug, report=report)
			ins.run()
		elif type==3:
			ins = ClassDict[type](hostname=hostname, dbname=dbname, schema=schema, input_dir=input_fname, array_data_table=argument1,\
								probes_table=argument2, strain_info_table=argument3, array_info_table=argument4,\
								commit=commit, debug=debug, report=report)
			ins.run()
		elif type==4:
			ins = ClassDict[type](hostname=hostname, dbname=dbname, schema=schema, array_info_table=input_fname, ecotype_table=argument1,\
								calls_comment_table=argument2,\
								commit=commit, debug=debug, report=report)
			ins.run()
		elif type==5:
			ins = ClassDict[type](hostname=hostname, dbname=dbname, schema=schema, output_dir=output_fname, snps_table=argument1, \
								probes_table=argument2, array_info_table=argument3,\
								debug=debug, report=report)
			ins.run()
	else:
		print __doc__
		sys.exit(2)
		