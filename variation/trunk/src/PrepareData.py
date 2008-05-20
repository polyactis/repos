#!/usr/bin/env python2.5
"""
Usage: PrepareData.py [OPTIONS] -o OUTPUT_FILE INPUT_FILE

Option:

        -o ..., --outputSNPsFile=...    SNPs data output file
        -u ..., --outputPhenotFile=...  Phenotype output file
	-d ..., --delim=...             default is \", \"      
        -m ..., --missingval=...        default is \"NA\"
	-a ..., --withArrayId=...       0 for no array ID info (default), 1 if input file has array ID info.
	-f ..., --phenotypeFile=...     File with phenotypes (deliminator separated format)
        -p ..., --phenotype=...         Phenotype index
        --phenotypeName=...             Phenotype name
        --orderAccessions               Order phenotype accessions (to match genotype order)
        --rawDataFormat                 Output data in raw format (with nucleotides)
        --filterMonomorphic             Filter monomorphic SNPs.
        --calcKinshipMatrix=...         Calculate kinship matrix and output in the given file.
	-h, --help	                Show this help

Examples:
	PrepareData.py -o readyData.csv 250K.csv 
	
Description:

"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["outputSNPsFile=","outputPhenotFile=", "monomorphic", "rawDataFormat", "delim=", "missingval=", "withArrayId=", "phenotype=", "phenotypeFile=", "phenotypeName=", "calcKinshipMatrix=", "orderAccessions", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:u:d:m:a:f:p:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	inputFile = args[0]
        output_fname = None
	outputPhenotFile = None
	delim = ","
	missingVal = "NA"
	phenotypeFile = None
        kinshipMatrixFile = None
	phenotype = None
	phenotypeName = None
	rawDataFormat = False
	monomorphic = False
	help = 0
	withArrayIds = 0
	orderAccessions = False
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-a","--withArrayId"):
			withArrayIds = int(arg)
		elif opt in ("-f","--phenotypeFile"):
			phenotypeFile = arg
                elif opt in ("calcKinshipMatrix"):
                        kinshipMatrixFile = arg
		elif opt in ("--monomorphic"):
			monomorphic = True
		elif opt in ("--rawDataFormat"):
			rawDataFormat = True
		elif opt in ("--minCallProb"):
			minCallProb = float(arg)
		elif opt in ("-p","--phenotype"):
			phenotype = int(arg)
		elif opt in ("-o","--outputSNPsFile"):
			output_fname = arg
		elif opt in ("--orderAccessions"):
			orderAccessions = True
		elif opt in ("-u","--phenotypeFile"):
			outputPhenotFile = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not output_fname:
		print output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)

	waid1 = withArrayIds==1 or withArrayIds==2
	waid2 = withArrayIds==2

	import dataParsers
	snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	
	if phenotypeFile:
		import phenotypeData
		phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')  #Get Phenotype data 
		accIndicesToKeep = []			
		phenAccIndicesToKeep = []
		numAcc = len(snpsds[0].accessions)
		if phenotype>=0:
		        #Load phenotype file
			sys.stdout.write("Removing accessions which do not have a phenotype value for "+phed.phenotypeNames[phenotype]+".")
			sys.stdout.flush()
			for i in range(0,len(snpsds[0].accessions)):
				acc1 = snpsds[0].accessions[i]
				for j in range(0,len(phed.accessions)):
					acc2 = phed.accessions[j]
					if acc1==acc2 and phed.phenotypeValues[j][phenotype]!='NA':
						accIndicesToKeep.append(i)
						phenAccIndicesToKeep.append(j)
						break					

		elif phenotype==None:
			sys.stdout.write("Removing accessions which do not have any phenotype values.")
			sys.stdout.flush()
			for i in range(0,len(snpsds[0].accessions)):
				acc1 = snpsds[0].accessions[i]
				for j in range(0,len(phed.accessions)):
					acc2 = phed.accessions[j]
					if acc1==acc2:
						accIndicesToKeep.append(i)
						phenAccIndicesToKeep.append(j)
						break
			
					
		#Filter Accessions which do not have the phenotype value.
		for snpsd in snpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			snpsd.removeAccessionIndices(accIndicesToKeep)
		print ""
		print numAcc-len(accIndicesToKeep),"accessions removed, leaving",len(accIndicesToKeep),"accessions in all."
		
		if outputPhenotFile:
			print "Filtering phenotype data."
			phed.removeAccessions(phenAccIndicesToKeep)
			if orderAccessions:
				accessionMapping = []
				i = 0
				for acc in snpsds[0].accessions:
					if acc in phed.accessions:
						accessionMapping.append((phed.accessions.index(acc),i))
						i += 1
				phed.orderAccessions(accessionMapping)
			if phenotype>=0:
				phed.writeToFile(outputPhenotFile, [phenotype])
			else:
				phed.writeToFile(outputPhenotFile)

		
        #Filtering monomorphic
	if monomorphic:
		print "Filtering monomorphic SNPs"
		for snpsd in snpsds:
			print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	import snpsdata
	
	newSnpsds = []
	if not rawDataFormat:
		sys.stdout.write("Converting data format")
		for snpsd in snpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			newSnpsds.append(snpsd.getSnpsData())
		print ""
		waid1 = 0
		snpsDataset = snpsdata.SnpsDataSet(newSnpsds,[1,2,3,4,5])
		decoder = {1:1, 0:0, -1:'NA'}
	else:
		snpsDataset = snpsdata.SnpsDataSet(snpsds,[1,2,3,4,5])
		decoder=None
	
	snpsDataset.writeToFile(output_fname, deliminator=delim, missingVal = missingVal, withArrayIds = waid1, decoder=decoder)

	
if __name__ == '__main__':
	_run_()


