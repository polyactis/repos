#!/usr/bin/env python
"""
Usage: RandomForest.py [OPTIONS] -o R_FILE SNPS_DATA_FILE PHENOTYPE_DATA_FILE [PHENOTYPE_INDEX]

Option:

        -o ..., --impFile=...       Random Forest output in importance score format.
        -d ..., --delim=...         default is ", "      
        -m ..., --missingval=...    default is "NA"
	--logTransform              Log transforms the phenotype values, before running RF.
        --phenotypeFileType=...     1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-a ..., --withArrayId=...   1 for array ID info (default), 0 if file has no array ID info.
	-h, --help	            show this help
	--parallel=...              Run randomForest on the cluster with standard parameters.  The arguement is used for runid as well as output files.  
	--parallelAll               Run randomForest on all phenotypes.
	--chunkSize=...             The chunk size. (default is full genome).
	--round2Size=...            The size of the second round run. (default is 5000).
	--nTrees=...                The number of trees. (default is 10000)
	--nodeSize=...              Default value is the R package default.
	--mem=..                    Memory requirements for the cluster.
	--secondRound               Do a second round RF.
        --minMAF=...                Remove all SNPs which have MAF smaller than the given argument.  (0.0 is set as default).
	

Examples:
	RandomForest.py -o imp_scores  250K.csv phenotypes.tsv phenotype_index 

	On the cluster for phenotype 0 (LD):
	python ../Projects/Python-snps/RandomForest.py -a 1 --parallel=RF_LD /home/cmb-01/bvilhjal/Projects/data/250K_method_5_after_imputation_noRedundant_051908.csv /home/cmb-01/bvilhjal/Projects/data/phenotypes_052208.tsv 0	

Description:

"""
#Relative to the home directory (not absolute directories).  These parameters are only used on the cluster.

resultDir="/home/cmb-01/bvilhjal/results/"
randomForestDir="/home/cmb-01/bvilhjal/Projects/randomForest/"
programDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"

#resultDir="/home/cmb-01/atarone/Projects/"
#randomForestDir="/home/cmb-01/atarone/randomForest/"
#programDir="/home/cmb-01/atarone/"

import sys, getopt, traceback
import os, env
import phenotypeData

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["chunkSize=", "nTrees=", "impFile=", "delim=", "missingval=", "withArrayId=", "logTransform", "phenotypeFileType=", "help", "parallel=", "parallelAll", "nodeSize=", "mem=", "round2Size=", "secondRound", "minMAF="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
        phenotypeFileType = 1
        impFile = None
	delim = ","
	missingVal = "NA"
	help = 0
	withArrayIds = 1
	parallel = None
	logTransform = False
	parallelAll = False
	chunkSize = 250000
	round2Size = 5000
	nTrees = 15000
	nodeSize = None
	mem = "8g"
	skipSecondRound = True
	minMAF = 0.0

	for opt, arg in opts:
            if opt in ("-h", "--help"):
                help = 1
                print __doc__
            elif opt in ("-a","--withArrayId"):
                withArrayIds = int(arg)
            elif opt in ("-o","--rFile"):
                impFile = arg
            elif opt in ("--phenotypeFileType"):
                phenotypeFileType = int(arg)
            elif opt in ("--parallel"):
                parallel = arg
            elif opt in ("--parallelAll"):
                parallelAll = True
            elif opt in ("--logTransform"):
                logTransform = True
            elif opt in ("--secondRound"):
                skipSecondRound = False
            elif opt in ("-d","--delim"):
                delim = arg
            elif opt in ("--chunkSize"):
                chunkSize = int(arg)
            elif opt in ("--round2Size"):
                round2Size = int(arg)		
            elif opt in ("--nTrees"):
                nTrees = int(arg)
            elif opt in ("--nodeSize"):
                nodeSize = int(arg)
            elif opt in ("--mem"):
                mem = arg
            elif opt in ("-m","--missingval"):
                missingVal = arg
            elif opt in ("-m","--minMAF"):
                minMAF = float(arg)
            else:
                if help==0:
                    print "Unkown option!!\n"
                    print __doc__
                sys.exit(2)

        if len(args)<3 and not parallel:
            if help==0:
                print "Arguments are missing!!\n"
                print __doc__
            sys.exit(2)

	
	def runParallel(phenotypeIndex):
		#Cluster specific parameters
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
		phenName = phed.getPhenotypeName(phenotypeIndex)
		phenName = phenName.replace("/","_div_")
		phenName = phenName.replace("*","_star_")
		impFileName = resultDir+"RF_"+parallel+"_"+phenName
		outFileName = impFileName
		shstr = """#!/bin/csh
#PBS -l walltime=120:00:00
"""
		shstr += "#PBS -l mem="+mem+"\n"
		shstr +="""
#PBS -q cmb
"""
		
		shstr += "#PBS -N RF"+phenName+"_"+parallel+"\n"
		shstr += "(python "+programDir+"RandomForest.py -o "+impFileName+" --chunkSize "+str(chunkSize)+" --nTrees "+str(nTrees)+" --mem "+str(mem)+" --round2Size "+str(round2Size)+""
		if nodeSize:
			shstr += " --nodeSize "+str(nodeSize)+" "
		if logTransform:
			shstr += " --logTransform "
		if not skipSecondRound:
			shstr += " --secondRound "
		shstr += " -a "+str(withArrayIds)+" "			
		shstr += snpsDataFile+" "+phenotypeDataFile+" "+str(phenotypeIndex)+" "
		shstr += "> "+outFileName+"_job"+".out) >& "+outFileName+"_job"+".err\n"

		f = open(parallel+".sh",'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	#Nested function ends

	snpsDataFile = args[0]
	phenotypeDataFile = args[1]
	if parallel:  #Running on the cluster..
		if parallelAll:
			phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
			for phenotypeIndex in phed.phenIds:
				runParallel(phenotypeIndex)
		else:
			phenotypeIndex = int(args[2])
			runParallel(phenotypeIndex)
		return
	else:
		phenotypeIndex = int(args[2])

	print "chunkSize:",chunkSize
	print "nTrees:",nTrees
	print "nodeSize:",nodeSize
	print "mem:",mem
	print "logTransform:",logTransform
	print "round2Size:",round2Size
	print "skipSecondRound:",skipSecondRound

	#Loading genotype data
	import dataParsers
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=withArrayIds)
	
	phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
	phenotype = phed.getPhenIndex(phenotypeIndex)
	accIndicesToKeep = []			
	phenAccIndicesToKeep = []
	numAcc = len(snpsds[0].accessions)

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


	#Filter accessions which do not have the phenotype value.
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snpsd.removeAccessionIndices(accIndicesToKeep)
	print ""
	print numAcc-len(accIndicesToKeep),"accessions removed, leaving",len(accIndicesToKeep),"accessions in all."
		
	print "Filtering phenotype data."
	phed.removeAccessions(phenAccIndicesToKeep) #Removing accessions that don't have genotypes or phenotype values
	
	#Ordering accessions according to the order of accessions in the genotype file
	accessionMapping = []
	i = 0
	for acc in snpsds[0].accessions:
		if acc in phed.accessions:
			accessionMapping.append((phed.accessions.index(acc),i))
			i += 1
	phed.orderAccessions(accessionMapping)

	#Log-transforming
	if logTransform:
		print "Log transforming phenotype"
		phed.logTransform(phenotype)

        #Filtering monomorphic
	print "Filtering monomorphic SNPs"
	for snpsd in snpsds:
		print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	#Remove minor allele frequencies
	if minMAF!=0:
		sys.stdout.write("Filterting SNPs with MAF<"+str(minMAF)+".")
		for snpsd in snpsds:
			sys.stdout.write(".")
			sys.stdout.flush()
			snpsd.filterMinMAF(minMAF)
		
	
      	#Converting format to 01
	import snpsdata
	newSnpsds = []
	sys.stdout.write("Converting data format")
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		newSnpsds.append(snpsd.getSnpsData())
	print ""
	snpsds = newSnpsds
	
	#Writing files
	import tempfile
	if env.user=="bjarni":
		tempfile.tempdir='/tmp'
	(fId, phenotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, genotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	
	phed.writeToFile(phenotypeTempFile, [phenotype])	
	sys.stdout.write( "Phenotype file written\n")
	sys.stdout.flush()
	
	#Retain only the correct runchunk of data.
	chromasomes = []
	positions = []
	snps = []
	for i in range(0,len(snpsds)):
		snpsd = snpsds[i]
		positions += snpsd.positions
		snps += snpsd.snps
		chrList = [i+1]*len(snpsd.positions)
		chromasomes += chrList

	#Is the phenotype binary?
	binary = phed.isBinary(phenotypeIndex)
	import util
	impFile = impFile+".imp"
	rDataFile = impFile+".rData"
	rFile = impFile+".r"
	outRfile = rFile+".out"
	errRfile = rFile+".err"
	topImpFile = impFile+"_top"+str(chunkSize)+".imp"
	topRDataFile = impFile+"_top.rData"
	try:
		os.remove(impFile)    #Removing file if it already exits.
	except Exception:
		print "Couldn't remove",impFile
	try:
		os.remove(topImpFile) #Removing file if it already exits.
	except Exception:
		print "Couldn't remove",topImpFile
	for startIndex in range(0,len(positions),chunkSize):
		if startIndex+chunkSize>=len(positions):
			endIndex = len(positions)
		else:
			endIndex = startIndex+chunkSize

	        #Writing genotype data to file.
		tmpFile = open(genotypeTempFile,"w")
		for i in range(startIndex,endIndex):
			outStr =""
			snp = util.valListToStrList(snps[i])
			outStr += str(chromasomes[i])+","+str(positions[i])+","
			outStr += ",".join(snp)
			outStr += "\n"
			tmpFile.write(outStr)
		tmpFile.close()
			
		rstr = _generateRScript_(genotypeTempFile, phenotypeTempFile, impFile, rDataFile, cluster=True, binary=binary, nTrees=nTrees, nodeSize=nodeSize)
		f = open(rFile,'w')
		f.write(rstr)
		f.close()
		#outRfile = rFile+"_"+str(startIndex/chunkSize)+".out"
		#errRfile = rFile+"_"+str(startIndex/chunkSize)+".err"
		print "Running model nr",startIndex/chunkSize,":"
		cmdStr = "(R --vanilla < "+rFile+" > "+outRfile+") >& "+errRfile
		sys.stdout.write(cmdStr+"\n")
		sys.stdout.flush()
		os.system(cmdStr)
	print "Random forest output saved in", impFile
	
	if not skipSecondRound:
		#Run on the top 'chunkSize' number of hits.
		#loading the R output file.
		impF = open(impFile,"r")
		lines=impF.readlines()
		impF.close()
		impList = list()
		for i in range(1,len(lines)):
			line = lines[i]
			line.strip()
			l = line.split(",")
			impList.append( (float(l[2]),l[0],l[1],snps[i]) )
		impList.sort()
		impList.reverse()

		#Writing genotype data to file.
		tmpFile = open(genotypeTempFile,"w")
		for i in range(0,round2Size):
			outStr = ""
			snp = util.valListToStrList(impList[i][3])
			outStr += str(impList[i][1])+","+str(impList[i][2])+","
			outStr += ",".join(snp)
			outStr += "\n"
			tmpFile.write(outStr)
		tmpFile.close()
		rstr = _generateRScript_(genotypeTempFile, phenotypeTempFile, topImpFile, topRDataFile, cluster=True, binary=binary, nTrees=nTrees, nodeSize=nodeSize)
		f = open(rFile,'w')
		f.write(rstr)
		f.close()
		print "Running randomForest on the top importance scores:"
		cmdStr = "(R --vanilla < "+rFile+" > "+outRfile+") >& "+errRfile
		sys.stdout.write(cmdStr+"\n")
		sys.stdout.flush()
		os.system(cmdStr)
	
	
def _generateRScript_(genotypeFile, phenotypeFile, impFile, rDataFile, cluster=False, binary=False, nTrees=500, nodeSize=False):
	
	if cluster:
		rstr = 'library(randomForest,lib.loc="'+randomForestDir+'");\n'
	else:
		rstr = "library(randomForest);\n"
	rstr += 'data250K <- read.csv("'+str(genotypeFile)+'", header=FALSE);\n'
	rstr += "mat250K <- as.matrix(data250K);\n"
	rstr += 'phenotData <- read.csv("'+str(phenotypeFile)+'",header=TRUE);\n'
	rstr += """
phenMat <- as.vector(as.matrix(phenotData)[,2]);
res <- list();
"""
	
	rstr += """
res[["pos"]] <- mat250K[,2];
res[["chr"]] <- mat250K[,1];
mat250K <- mat250K[,3:length(mat250K[1,])]

#Calculate MAF
maf <- c();
marf <- c();
for(i in (1:length(mat250K[,1]))){
  f = factor(mat250K[i,]);
  freq <- summary(f)[[1]]
  alleleCount <- length(f)
  v <- freq/alleleCount;
  maf[i] <- min(freq,alleleCount-freq);
  marf[i] <- min(v,1-v);
}
res[["marf"]] <- marf;
res[["maf"]] <- maf;

mat250K <- t(mat250K);
res[["m"]] <-list()
"""

	if binary: 
		print "Random Forest: Using classification"
		if nodeSize:
			rstr += 'res$m <- randomForest(mat250K,y=as.factor(phenMat),ntree='+str(nTrees)+',keep.forest=F,nodesize='+str(nodeSize)+');\n'
		else:
			rstr += 'res$m <- randomForest(mat250K,y=as.factor(phenMat),ntree='+str(nTrees)+',keep.forest=F);\n'
	else:
		print "Random Forest: Using regression"
		if nodeSize:
			rstr += 'res$m <- randomForest(mat250K,y=phenMat,ntree='+str(nTrees)+',keep.forest=F,nodesize='+str(nodeSize)+');\n'
			#rstr += 'res$m <- randomForest(phenMat~.,data=as.data.frame(mat250K),ntree='+str(nTrees)+',keep.forest=F,nodesize='+str(nodeSize)+');\n'
		else:
			rstr += 'res$m <- randomForest(mat250K,y=phenMat,ntree='+str(nTrees)+',keep.forest=F);\n'
			#rstr += 'res$m <- randomForest(phenMat~.,data=as.data.frame(mat250K),ntree='+str(nTrees)+',keep.forest=F);\n'
	rstr +="""


#write to a importance score
l <- list();
l[["Chromasome"]]<-res[["chr"]];
l[["Positions"]]<-res[["pos"]];
l[["Importance"]]<-res[["m"]]$importance;
l[["MARF"]]<-marf;
l[["MAF"]]<-maf;
dl <- as.data.frame(l)
"""
	rstr +=' write.table(dl,file="'+impFile+'", append=T, sep=", ", col.names=F, row.names=F);\n'		
	rstr += 'save(file="'+rDataFile+'",res);\n'
	return rstr	
	
if __name__ == '__main__':
	_run_()
