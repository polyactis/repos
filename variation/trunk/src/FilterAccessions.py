#!/usr/bin/env python2.5
"""
Usage: FilterAccessions.py [OPTIONS] -o OUTPUT_FILE INPUT_FILE

Option:

        -o ...,	output file
	-d ..., --delim=...         default is \", \"      
        -m ..., --missingval=...    default is \"NA\"
	-a ..., --withArrayId=...   0 for no array ID info (default), 1 if file has array ID info, 2 if comparison file also.
	--maxError=...              maximum allowed error percentage (requires a comparison file).
        --comparisonFile=...        a file which is used to claculate the SNPs error.
	--maxMissing=...            maximum allowed missing percentage.
        --removeEcotypeId=...       removes all accessions with the given ecotype ID. 
        --removeArrayId=...         removes an accessions with the given array ID. 
	--removeIdentical           removes redundant accessions picking the one with the least error (requires a comparison 
                                    file).
        --onlyCommon                removes all accessions which are not both in the input file and the comparison file,
                                    (requires a comparison file).
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

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
	
    long_options_list = ["maxError=", "comparisonFile=", "maxMissing=", "removeEcotypeId=", "removeArrayId=", "removeIdentical", "onlyCommon", "delim=", "missingval=", "withArrayId=", "debug", "report", "help"]
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
    removeEcotype = None
    removeArray = None
    removeIdentical = False
    onlyCommon = False
    debug = None
    report = None
    help = 0
    withArrayIds = 0

	
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
    waid2 = withArrayIds==2 #or withArrayIds==0

    import dataParsers
    snpsds = dataParsers.parseCSVData(inputFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid1)
	
    accessionsToRemove = []
    arraysToRemove = None

    #Retrieve comparison list of accessions.  (Error rates for accessions)
    if (removeIdentical or maxError<1.0) and comparisonFile:
        print("Loading comparison file:")
        sys.stdout.flush()
        snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
        res = []
        sys.stdout.write("Comparing accessions.")
        for i in range(0,len(snpsds)):
            res.append(snpsds[i].compareWith(snpsds2[i],withArrayIds=withArrayIds,verbose=False))
            sys.stdout.write(".")
            sys.stdout.flush()
        print ""

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
	if withArrayIds:
            for i in range(0,len(r[2])):
                accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))
	else:
            for i in range(0,len(r[2])):
                accErrAndID.append((accErrorRate[i], r[2][i]))
	accErrAndID.sort(reverse=True)
        #print "(Error,'ecotype_id','array_id')"
        #for t in accErrAndID :
        #    print t    

   
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
        print accessionsToRemove, arraysToRemove


    if onlyCommon and comparisonFile:
        print "Locating accessions which are not shared"
        snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
        if not arraysToRemove:
            arraysToRemove = []
        for i in range(0,len(snpsds[0].accessions)):
            acc = snpsds[0].accessions[i]
            if not acc in snpsds2[0].accessions:
                accessionsToRemove.append(acc)
                if withArrayIds:
                    arraysToRemove.append(snpsds[0].arrayIds[i])
        print accessionsToRemove, arraysToRemove


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
            missingRates.sort(reverse=True)
            for (mrate,ecotype,array) in missingRates:
                if mrate>maxMissing:
                    accessionsToRemove.append(ecotype)
                    arraysToRemove.append(array)
        else:
            for i in range(0,len(snpsds[0].accessions)):
                missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i]))
            missingRates.sort(reverse=True)
            for (mrate,ecotype) in missingRates:
                if mrate>maxMissing:
                    accessionsToRemove.append(ecotype)
        #print "(NA-rate,'ecotype_id','array_id')"
        #for t in missingRates :
        #    print t



    if removeEcotype:
        accessionsToRemove.append(removeEcotype)
    if removeArray:
        if not arraysToRemove:
            arraysToRemove = []
        arraysToRemove.append(removeArray)
        

    numAccessions = len(snpsds[0].accessions)
    sys.stdout.write("Removing accessions.")
    sys.stdout.flush()
    for snpsd in snpsds:
        snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
        sys.stdout.write(".")
        sys.stdout.flush()
    print "\n", (numAccessions-len(snpsds[0].accessions)), "accessions out of "+str(numAccessions)+" were removed."
        

    """
    #Filtering missing values
    if maxMissing<1.0 and maxMissing>=0.0:
        print "Filtering SNPs with missing values"
        numAccessions = len(snpsds[0].accessions)
        for snpsd in snpsds:
            print "Removed", str(snpsd.filterMissingSnps(int(maxMissing*numAccessions))),"Snps"

    #Filtering bad SNPs
    if comparisonFile and maxError<1.0:
        print "Filtering bad SNPs"
        snpsds2 = dataParsers.parseCSVData(comparisonFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=waid2)
        for i in range(0,len(snpsds)):
            snpsds[i].filterBadSnps(snpsds2[i],maxError)
		
    import snpsdata
    snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)
    """
    import snpsdata
    snpsdata.writeRawSnpsDatasToFile(output_fname,snpsds,chromosomes=[1,2,3,4,5], deliminator=delim, missingVal = missingVal, withArrayIds = waid1)



#--------------------------------------------------------------------------------#
def filterByError(snpsds,comparisonSnpsds,maxError,withArrayIds=1):
    """
    Removes the accessions in the snpsds object if they have error-rates greater than maxError. 
    Error is calculated by comparisonwith the comparisonSnpsds object.
    """
    accessionsToRemove = []
    arraysToRemove = None

    res = []
    sys.stdout.write("Comparing accessions.")
    for i in range(0,len(snpsds)):
        res.append(snpsds[i].compareWith(comparisonSnpsds[i],withArrayIds=withArrayIds,verbose=False))
        sys.stdout.write(".")
        sys.stdout.flush()
    print ""

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
    if withArrayIds: #Then include arrayID in list..
        for i in range(0,len(r[2])):
            accErrAndID.append((accErrorRate[i], r[2][i], r[5][i]))
    else:
        for i in range(0,len(r[2])):
            accErrAndID.append((accErrorRate[i], r[2][i]))
    accErrAndID.sort(reverse=True)

    #Figure out which accessions are too erroraneous
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

    #Remove accessions
    numAccessions = len(snpsds[0].accessions)
    sys.stdout.write("Removing accessions.")
    sys.stdout.flush()
    for snpsd in snpsds:
        snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
        sys.stdout.write(".")
        sys.stdout.flush()
    print "\n", (numAccessions-len(snpsds[0].accessions)), "accessions out of "+str(numAccessions)+" were removed."
    return snpsds



def filterByNA(snpsds,maxMissing,withArrayIds=1):
    """
    Removes the accessions in the snpsds list if they have NA-rates greater than maxMissing. 
    """
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
        missingRates.sort(reverse=True)
        for (mrate,ecotype,array) in missingRates:
            if mrate>maxMissing:
                accessionsToRemove.append(ecotype)
                arraysToRemove.append(array)
    else:
        for i in range(0,len(snpsds[0].accessions)):
            missingRates.append((missingCounts[i]/float(numSnps),snpsds[0].accessions[i]))
        missingRates.sort(reverse=True)
        print missingRates
        print len(missingRates)
        for (mrate,ecotype) in missingRates:
            if mrate>maxMissing:
                accessionsToRemove.append(ecotype)

    numAccessions = len(snpsds[0].accessions)
    sys.stdout.write("Removing accessions.")
    sys.stdout.flush()
    for snpsd in snpsds:
        snpsd.removeAccessions(accessionsToRemove,arrayIds=arraysToRemove)
        sys.stdout.write(".")
        sys.stdout.flush()
    print "\n", (numAccessions-len(snpsds[0].accessions)), "accessions out of "+str(numAccessions)+" were removed."

    return snpsds


def _testRun1_():
    import dataParsers
    snpsds = dataParsers.parseCSVData("250K_m3.csv",withArrayIds=1)
    comparisonSnpsds = dataParsers.parseCSVData("2010_v3.csv")
    filterByError(snpsds,comparisonSnpsds,0.2,withArrayIds=1)

if __name__ == '__main__':
    #_testRun1_()
    _run_()


