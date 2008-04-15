#!/usr/bin/env python2.5
"""
Usage: MergeSNPsData.py [OPTIONS] -o OUTPUT_FILE datafile1 datafile2

Option:
        -o ...,	output file
        -p ..., --priority=...          1 for giving first file priority, 2 for second file (if SNPs values disagree). (1 is default)
        -s ..., --maxMissingSnps=...    Maximum number of allowed missing SNPs values. (no limit is default)
        -a ..., --maxMissingAcc=...     Maximum fraction of missing SNPs allowed in any accession (otherwise the accession is removed). 
                                        (No limit is default.)
        -d ..., --delim=...             default is ", "      
        -m ..., --missingval=...        default is "NA"
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	 MergeSNPsData.py -o /tmp/phenotype.csv datafile1 datafile2
	
Description:
	Merges two datafiles into a new file, in a csv format (or with another deliminator separated format).

	Requires MySQLdb to be installed, and uses snpsdata.py and dataParsers.py.
"""
import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["priority", "delim", "missingval", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:p:s:a:d:m:brh", long_options_list)

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
        maxMissingSnps = None
        maxMissingAcc = None
        output_fname = None
	method = None 
	delim = ", "
	missingVal = "NA"
	debug = None
	report = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-p", "--priority"):
			priority = int(arg)
		elif opt in ("-s", "--maxMissingSnps"):
			maxMissingSnps = int(arg)
		elif opt in ("-a", "--maxMissingAcc"):
			maxMissingAcc = int(arg)
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-t","--method"):
			version = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg	
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
        
        
        import dataParsers
	if priority == 1:
            snpsds1 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal)
            snpsds2 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal)
        elif priority == 2:
            snpsds1 = dataParsers.parseCSVData(inputFile2, format=1, deliminator=delim, missingVal=missingVal)
            snpsds2 = dataParsers.parseCSVData(inputFile1, format=1, deliminator=delim, missingVal=missingVal)
        else:
            raise Exception("Priority option is wrong.")

        if len(snpsds1) != len(snpsds2):
            raise Exception("Unequal number of chromosomes.")
        
        newSnpsds = []
        import snpsdata
        for i in range(0,len(snpsds1)):
            newSnpsds.append(snpsds1[i].mergeData(snpsds2[i],multIdenticalAcc=True))
	
        #Filtering
        if maxMissingSnps:
            numRemoved = 0
            for snpsd in newSnpsds:
                numRemoved += snpsd.filterMissingSnps(self,maxNumMissing=0)
            print numRemoved, "Snps were filtered out of data."
        
        #Filter accessions

	snpsdata.writeRawSnpsDatasToFile(output_fname,newSnpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal)


if __name__ == '__main__':
	_run_()

