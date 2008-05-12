#!/usr/bin/env python2.5
"""
Usage: MergeSNPsData.py [OPTIONS] -o OUTPUT_FILE datafile1 datafile2

Option:
        -o ...,	output file
        -p ..., --priority=...          1 for giving first file priority, 2 for second file (if SNPs values disagree). (1 is default)
        -d ..., --delim=...             default is ", "      
        -m ..., --missingval=...        default is "NA"
	-a ..., --withArrayId=...       0 for no array ID info (default), 1 if first file has array ID info, 2 if both have.
        -u ..., --union=...             1: All SNPs are included in the dataset 
	                                2: All accessions are included in the merged dataset 
					3: Both all SNPs and accessions are included in the merged dataset

        -i ..., --intersection=...      1: Only common SNPs are included in the dataset (not implemented)
	                                2: Only common accessions are included in the merged dataset (not implemented)
					3: Only common SNPs and accessions are included in the merged dataset (not implemented)

	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	 MergeSNPsData.py -o /tmp/phenotype.csv datafile1 datafile2
	
Description:
	Merges two datafiles into a new file, in a csv format (or with another deliminator separated format).

	If neither the --intersection nor the --union flags are used then the output data has the same SNPs and 
	accessions as the first datafile.

	Note that the higher priority file should never have multiple accession. 

	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""
import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["priority=", "delim=", "missingval=", "union=", "intersection=", "debug", "report", "help", "withArrayId="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:p:d:m:u:i:a:brh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	if len(args)!=2:
            raise Exception("Number of arguments isn't correct.")
        inputFile1 = args[0]
	inputFile2 = args[1]
	priority = 1
        union = 0
        intersection = 0
        output_fname = None
	delim = ","
	missingVal = "NA"
	debug = None
	report = None
	withArrayIds = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-p", "--priority"):
			priority = int(arg)
		elif opt in ("-u", "--union"):
			union = int(arg)
		elif opt in ("-i", "--intersection"):
			intersection = int(arg)
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
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
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)
        
	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2

        import dataParsers
	snpsds1 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	snpsds2 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
	withArrayIds = waid1
	
	
        if len(snpsds1) != len(snpsds2):
		raise Exception("Unequal number of chromosomes.")
        
	import snpsdata
	if union==0 and intersection==0:
		for i in range(0,len(snpsds1)):
			snpsds1[i].mergeData(snpsds2[i],priority=priority)
		snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds1,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)
	elif union>0 and union<4 and intersection==0:
		for i in range(0,len(snpsds1)):
			snpsds1[i].mergeDataUnion(snpsds2[i], priority=priority, unionType=union)
		snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds1,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal)
	else:
		if help==0:
                        print "The union or intersection options used are wrong!!\n"
                        print __doc__
                sys.exit(2)




	
def merge(snpsds1,snpsds2,unionType=0,priority=1):
        if len(snpsds1) != len(snpsds2):
		raise Exception("Unequal number of chromosomes.")
        
	import snpsdata
	if unionType==0:
		for i in range(0,len(snpsds1)):
			snpsds1[i].mergeData(snpsds2[i],priority=priority)
	elif unionType>0 and unionType<4:
		for i in range(0,len(snpsds1)):
			snpsds1[i].mergeDataUnion(snpsds2[i], priority=priority, unionType=unionType)
	else:
		print "The union option used are wrong!!\n"
	return snpsds1



def _test_():
        import dataParsers
	snpsds1 = dataParsers.parseCSVData("149_v1.csv", deliminator=",")
        snpsds2 = dataParsers.parseCSVData("384.csv", deliminator=",")	
	merge(snpsds1,snpsds2,unionType=1,priority=1)
	print snpsds1[0].positions

if __name__ == '__main__':
	_run_()
	#_test_()
