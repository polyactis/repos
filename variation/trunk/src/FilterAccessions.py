#!/usr/bin/env python2.5
"""
Usage: FilterAccessions.py [OPTIONS] -o OUTPUT_FILE INPUT_FILE

Option:

	-o ...,	output file
	-d ..., --delim=...	    default is \", \"	  
	-m ..., --missingval=...    default is \"NA\"
	-a ..., --withArrayId=...   1 for array ID info (default), 0 if file has array no ID info, 2 if comparison file also, 3 if only comparison file.
	--maxError=...		    maximum allowed error percentage (requires a comparison file).
	--comparisonFile=...	    a file which is used to claculate the SNPs error.
	--maxMissing=...	    maximum allowed missing percentage.
	--removeEcotypeId=...	    removes all accessions with the given ecotype ID. 
	--removeArrayId=...	    removes an accessions with the given array ID. 
	--removeIdentical	    removes redundant accessions picking the one with the least error (requires a comparison 
				    file).
	--onlyCommon		    removes all accessions which are not both in the input file and the comparison file,
				    (requires a comparison file).

        --first96                   Outputs only the first (old) 96 accessions.
	-b, --debug	            enable debug
	-h, --help                  show this help

Examples:
	FilterAccessions.py --maxMissing=0.5 -o /tmp/2010_filtered.csv 2010.csv
	
Description:
	Filter a csv formatted file with respect to various criteria.  
	Note that this script only removes Accessions not SNPs. 

	Requires MySQLdb to be installed, as well as util.py, rfun.py, snpsdata.py and dataParsers.py.
"""

import sys, getopt, traceback

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["maxError=", "comparisonFile=", "maxMissing=", "removeEcotypeId=", "removeArrayId=", "first96", "removeIdentical", "onlyCommon", "delim=", "missingval=", "withArrayId=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:bh", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	inputFile = args[0]
	output_fname = None
	delim = ","
	missingVal = "NA"
	comparisonFile = None
	maxMissing = 1.0
	maxError = 1.0
	removeEcotype = None
	removeArray = None
	removeIdentical = False
	onlyCommon = False
	debug = None
	report = None
	help = 0
	withArrayIds = 1
	first96 = False
	
	for opt, arg in opts:
		if opt in ('-o'):
			output_fname = arg
		elif opt in ("-h", "--help"):
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
		elif opt in ("--removeEcotypeId"):
			removeEcotype = float(arg)
		elif opt in ("--removeArrayId"):
			removeArray = float(arg)
		elif opt in ("--removeIdentical"):
			removeIdentical = True
		elif opt in ("--onlyCommon"):
			onlyCommon = True
		elif opt in ("--first96"):
			first96 = True
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
	waid2 = withArrayIds==2 or withArrayIds==3

	import dataParsers
	snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	
	accessionsToRemove = []
	arraysToRemove = None

	if first96:
		import dataParsers
		d = dataParsers.getEcotypeToAccessionDictionary(defaultValue='-1',user="bvilhjal",passwd="bamboo123")
		ecotd = dataParsers.getEcotypeToNameDictionary(defaultValue='-1',user="bvilhjal",passwd="bamboo123")
		print "Dictionaries loaded"
		names = []
		first96Names = []
		for i in range(0,len(snpsds[0].accessions)):
			ecotype = snpsds[0].accessions[i]
			arrayID = snpsds[0].arrayIds[i]
			names.append((arrayID,ecotd[ecotype],ecotype))
			if int(d[ecotype][0]) > 97 or int(d[ecotype][0]) < 0:
				accessionsToRemove.append(ecotype)
			else:
				first96Names.append((arrayID,d[ecotype][1],d[ecotype][0],ecotype))

		first96Names.sort()
		print "First 96 accessions, len:",len(first96Names),":"
		for name in first96Names:
			print name
		names.sort()
		print "All accessions:"
		for name in names:
			print name


	#Retrieve comparison list of accessions.  (Error rates for accessions)
	if (removeIdentical or maxError<1.0) and comparisonFile:
		sys.stderr.write("Loading comparison file:")
		snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
		res = []
		sys.stderr.write("Comparing accessions.")
		for i in range(0,len(snpsds)):
			res.append(snpsds[i].compareWith(snpsds2[i],withArrayIds=withArrayIds,verbose=False))
			sys.stderr.write(".")
		sys.stderr.write("\n")

		totalAccessionCounts = [0]*len(res[0][2])
		accErrorRate = [0]*len(res[0][2])
		for i in range(0,len(snpsds)):
			r = res[i]
			for j in range(0,len(r[2])):
				totalAccessionCounts[j] += r[6][j]
				accErrorRate[j]+=r[3][j]*float(r[6][j])
		
		for i in range(0,len(accErrorRate)):
			accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])

		accErrAndID = []
		if 0<withArrayIds<3:
			for i in range(0,len(r[2])):
				accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))
		else:
			for i in range(0,len(r[2])):
				accErrAndID.append((accErrorRate[i], r[2][i]))
		accErrAndID.sort()
		accErrAndID.reverse()


   
	#Figure out which accessions are too erroraneous
	if maxError<1.0 and comparisonFile:
		if withArrayIds:
			arraysToRemove = []
			for (error,ecotype,array) in accErrAndID:
				if error> maxError:
					accessionsToRemove.append(ecotype)
					arraysToRemove.append(array)

		else:
			for (error,ecotype) in accErrAndID:
				if error> maxError:
					accessionsToRemove.append(ecotype)


	if removeIdentical and comparisonFile and withArrayIds:
		print "Locating identical accessions"
		accErrAndID.sort()
		if not arraysToRemove:
			arraysToRemove = []
		for accession in set(snpsds[0].accessions):
			if snpsds[0].accessions.count(accession)>1:
				found = 0
				for (error,ecotype,array) in accErrAndID:
					if ecotype==accession:
						if found>0:
							accessionsToRemove.append(ecotype)
							arraysToRemove.append(array)
						found += 1

	if onlyCommon and comparisonFile:
		print "Locating accessions which are not shared"
		snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
		#print snpsds2[0].accessions,'\n',snpsds[0].accessions,'\n',len(set(snpsds2[0].accessions).intersection(set(snpsds[0].accessions)))
		if not arraysToRemove:
			arraysToRemove = []
		for i in range(0,len(snpsds[0].accessions)):
			acc = snpsds[0].accessions[i]
			if not acc in snpsds2[0].accessions:
				accessionsToRemove.append(acc)
				if 0<withArrayIds<3:
					arraysToRemove.append(snpsds[0].arrayIds[i])


	if maxMissing<1.0:
		missingCounts = [0]*len(snpsds[0].accessions)
		numSnps = 0
		for snpsd in snpsds:
			mc = snpsd.accessionsMissingCounts()
			numSnps += len(snpsd.positions)
			for i in range(0,len(snpsds[0].accessions)):
				missingCounts[i] += mc[i]
		
		missingRates = []		
		if withArrayIds:
			arraysToRemove = []
			for i in range(0,len(snpsds[0].accessions)):
				missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i],snpsds[0].arrayIds[i]))
			missingRates.sort()
			missingRates.reverse()
			for (mrate,ecotype,array) in missingRates:
				if mrate>maxMissing:
					accessionsToRemove.append(ecotype)
					arraysToRemove.append(array)
		else:
			for i in range(0,len(snpsds[0].accessions)):
				missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i]))
			missingRates.sort()
			missingRates.reverse()
			for (mrate,ecotype) in missingRates:
				if mrate>maxMissing:
					accessionsToRemove.append(ecotype)


	if removeEcotype:
		accessionsToRemove.append(removeEcotype)
	if removeArray:
		if not arraysToRemove:
			arraysToRemove = []
		arraysToRemove.append(removeArray)
		

	numAccessions = len(snpsds[0].accessions)
	sys.stderr.write("Removing accessions.")
	for snpsd in snpsds:
		snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
		sys.stderr.write(".")
	print "\n", (numAccessions-len(snpsds[0].accessions)), "accessions out of "+str(numAccessions)+" were removed."
		
	import snpsdata
	snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)



#--------------------------------------------------------------------------------#
def filterByError(snpsds,comparisonSnpsds,maxError,withArrayIds=1,onlyCommon=False):
	"""
	05/13/2008 add no_of_accessions_filtered_by_mismatch
	Removes the accessions in the snpsds object if they have error-rates greater than maxError. 
	Error is calculated by comparisonwith the comparisonSnpsds object.
	"""
	accessionsToRemove = []
	arraysToRemove = None

	res = []
	sys.stderr.write("Comparing accessions.")
	for i in range(0,len(snpsds)):
		res.append(snpsds[i].compareWith(comparisonSnpsds[i],withArrayIds=withArrayIds,verbose=False))
		sys.stderr.write(".")
	print ""

	totalAccessionCounts = [0]*len(res[0][2])
	accErrorRate = [0]*len(res[0][2])
	for i in range(0,len(snpsds)):
		r = res[i]
		for j in range(0,len(r[2])):
			totalAccessionCounts[j] += r[6][j]
			accErrorRate[j]+=r[3][j]*float(r[6][j])
		
	for i in range(0,len(accErrorRate)):
		if totalAccessionCounts[i]>0:
			accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])
		else:
			accErrorRate[i] = 1

	accErrAndID = []
	if withArrayIds: #Then include arrayID in list..
		for i in range(0,len(r[2])):
			accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))
	else:
		for i in range(0,len(r[2])):
			accErrAndID.append((accErrorRate[i], r[2][i]))
	accErrAndID.sort()
	accErrAndID.reverse()

	#Figure out which accessions are too erroraneous
	if withArrayIds:
		arraysToRemove = []
		for (error,ecotype,array) in accErrAndID:
			if error> maxError:
				accessionsToRemove.append(ecotype)
				arraysToRemove.append(array)

		if onlyCommon:
			for i in range(0,len(snpsds[0].accessions)): #Removing all arrays that are comparable
				acc = snpsds[0].accessions[i]
				if not acc in comparisonSnpsds[0].accessions:
					accessionsToRemove.append(acc)
					arraysToRemove.append(snpsds[0].arrayIds[i])

	else:
		for (error,ecotype) in accErrAndID:
			if error> maxError:
				accessionsToRemove.append(ecotype)
		if onlyCommon:
			for i in range(0,len(snpsds[0].accessions)): #Removing all arrays that are comparable
				acc = snpsds[0].accessions[i]
				if not acc in comparisonSnpsds[0].accessions:
					accessionsToRemove.append(acc)

	#Remove accessions
	sys.stderr.write("Removing accessions.")
	numAccessions = len(snpsds[0].accessions)
	no_of_accessions_removed = None
	for snpsd in snpsds:
		snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
		if no_of_accessions_removed ==None:
   			no_of_accessions_removed = numAccessions-len(snpsd.accessions)
		snpsd.no_of_accessions_filtered_by_mismatch = no_of_accessions_removed
		sys.stderr.write(".")
	sys.stderr.write("\n %s accessions out of %s were removed.\n"%(no_of_accessions_removed, numAccessions))
	return snpsds



def filterByNA(snpsds,maxMissing,withArrayIds=1, debug=0):
	"""
	2008-05-12
		cosmetic change a bit
		add no_of_accessions_removed to every snpsd
	Removes the accessions in the snpsds list if they have NA-rates greater than maxMissing. 
	"""
	sys.stderr.write("Calculating NA rate ...")
	accessionsToRemove = []
	arraysToRemove = None

	missingCounts = [0]*len(snpsds[0].accessions)
	numSnps = 0
	for snpsd in snpsds:
		mc = snpsd.accessionsMissingCounts()
		numSnps += len(snpsd.positions)
		for i in range(0,len(snpsds[0].accessions)):
			missingCounts[i] += mc[i]
		
	missingRates = []		
	if withArrayIds:
		arraysToRemove = []
		for i in range(0,len(snpsds[0].accessions)):
			missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i],snpsds[0].arrayIds[i]))
		missingRates.sort()
		missingRates.reverse()
		for (mrate,ecotype,array) in missingRates:
			if mrate>maxMissing:
				accessionsToRemove.append(ecotype)
				arraysToRemove.append(array)
	else:
		for i in range(0,len(snpsds[0].accessions)):
			missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i]))
		missingRates.sort()
		missingRates.reverse()
		if debug:
			sys.stderr.write(' (%s, %s) '%(repr(missingRates), len(missingRates)))
		
		for (mrate,ecotype) in missingRates:
			if mrate>maxMissing:
				accessionsToRemove.append(ecotype)
	sys.stderr.write("Done.\n")
	
	sys.stderr.write("Removing accessions ... ")
	numAccessions = len(snpsds[0].accessions)
	no_of_accessions_removed = None
	for snpsd in snpsds:
		snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
		if no_of_accessions_removed ==None:
   			no_of_accessions_removed = numAccessions-len(snpsd.accessions)
		snpsd.no_of_accessions_filtered_by_na = no_of_accessions_removed
		sys.stderr.write(".")
	sys.stderr.write("\n%s accessions out of %s were removed.\n"%(no_of_accessions_removed, numAccessions))
	return snpsds


def removeIdentical(snpsds,comparisonSnpsds,withArrayIds=1, debug=0):
	"""
	2008-05-17
		cosmetic change
	2008-05-12
		cosmetic change a bit
		add no_of_accessions_removed to every snpsd
	Removes the accessions in the snpsds list if they have NA-rates greater than maxMissing. 
	"""
	accessionsToRemove = []
	arraysToRemove = []

	res = []
	sys.stderr.write("Comparing accessions.")
	for i in range(0,len(snpsds)):
		res.append(snpsds[i].compareWith(comparisonSnpsds[i],withArrayIds=withArrayIds,verbose=False))
		sys.stderr.write(".")
	sys.stderr.write("\n")

	totalAccessionCounts = [0]*len(res[0][2])
	accErrorRate = [0]*len(res[0][2])
	for i in range(0,len(snpsds)):
		r = res[i]
		for j in range(0,len(r[2])):
			totalAccessionCounts[j] += r[6][j]
			accErrorRate[j]+=r[3][j]*float(r[6][j])
		
	for i in range(0,len(accErrorRate)):
		if totalAccessionCounts[i]>0:
			accErrorRate[i]=accErrorRate[i]/float(totalAccessionCounts[i])
		else:
			accErrorRate[i] = 1

	accErrAndID = []
	for i in range(0,len(r[2])):
		accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))

	accErrAndID.sort()

	print "Locating identical accessions"
	arraysToRemove = []
	for accession in set(snpsds[0].accessions):
		if snpsds[0].accessions.count(accession)>1:
			found = 0
			for (error,ecotype,array) in accErrAndID:
				if ecotype==accession:
					if found>0:
						accessionsToRemove.append(ecotype)
						arraysToRemove.append(array)
					found += 1
	
	sys.stderr.write("Removing accessions ... ")
	numAccessions = len(snpsds[0].accessions)
	no_of_accessions_removed = None
	for snpsd in snpsds:
		snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
		if no_of_accessions_removed ==None:
   			no_of_accessions_removed = numAccessions-len(snpsd.accessions)
		snpsd.no_of_accessions_filtered_by_na = no_of_accessions_removed
		sys.stderr.write(".")
	sys.stderr.write("\n%s accessions out of %s were removed.\n"%(no_of_accessions_removed, numAccessions))
	return snpsds


def _testRun1_():
	import dataParsers
	snpsds = dataParsers.parseCSVData("250K_m3.csv",withArrayIds=1)
	comparisonSnpsds = dataParsers.parseCSVData("2010_v3.csv")
	filterByError(snpsds,comparisonSnpsds,0.2,withArrayIds=1)

if __name__ == '__main__':
	#_testRun1_()
	_run_()

