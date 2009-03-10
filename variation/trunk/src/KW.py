#!/usr/bin/env python2.5
"""
Usage: KW.py [OPTIONS] [-o OUT_FILE] SNPS_DATA_FILE PHENOTYPE_DATA_FILE [PHENOTYPE_INDEX]

Option:

	-o ..., --outputFile=...		Output files, one 'name'.rData file, and one 'name'.pvals file.
	-d ..., --delim=...			default is ", "	  
	-m ..., --missingval=...		default is "NA"
	--phenotypeFileType=...	 		1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-a ..., --withArrayId=...   		1 for array ID info (default), 0 if file has no array ID info.
	-h, --help				show this help
	--parallel=...				Run KW on the cluster with standard parameters.  The arguement is used for runid 
						as well as output files.  Note that if this option is used then no output file should be specified.
	--parallelAll				Run KW on all phenotypes.
	--addToDB				Adds the result to the DB.
	--callMethodID=...			What data set is used. (Only applicable if result is added to DB.)
	--comment=...				Comment for DB. (Only applicable if result is added to DB.)
	--onlyOriginal96			Run KW on the original 96 accessions only.
	--onlyOriginal192			Run KW on the original 192 2010 data accessions only.
	--onlyBelowLatidue=...			Only accessions below the given latitude are left in.
	--onlyAboveLatidue=...			Only accessions above the given latitude are left in.
	--complement				Choose the complement set of accessions (works only with onlyOriginal96 or onlyOriginal192).
	--subSample=...				Pick a random subsample of the given size.
	--subsampleTest=num_acc,num_samples	Outputs num_samples files, each with random num_acc accessions.
	--subSampleLikePhenotype=...		Pick subsample with accessions which are only in the given phenotype.
	
	--permTest=numPerm			Perform a permutation test
	--savePermutations			Stores the permutations in a file

	--testRobustness			Perform a robustness test

	--permutationFilter=...                 Only use a random sample of proportional size equal to the given number.
	
	--sr							Perform a second run.
	--srSkipFirstRun				Skip first run
	--srInput=pvalFile 		       	Use given results as input. (use with srSkipFirstRun)
	--srOutput=resultFile 		    SOutput new results in given file. 
	--srPar=topQuantile,windowSize	        Default topQuantile is 0.95, and windowSize = 20000 bases. 

Examples:
	KW.py -o outputFile  250K.csv phenotypes.tsv phenotype_index 
	
Description:
  Applies the Kruskal Wallis test to phenotypes, which are not binary.
  Applies a Chi-square test to the phenotypes which are binary!

"""
#from epydoc.docparser import str_prefixes
import sys, getopt, traceback
import os, env
import phenotypeData
import dataParsers
import gc
import random
import snpsdata
import gwaResults
import tempfile
import plotResults 
import SecondStageAnalysis
import util
from rpy import r
import math
import time
import random
import analyzePhenotype
import analyzeHaplotype


#import AddResults2DB

resultDir="/home/cmb-01/bvilhjal/results/"
scriptDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"


def _run_():
	if len(sys.argv)==1:
		print __doc__
		sys.exit(2)
	
	long_options_list=["outputFile=", "delim=", "missingval=", "withArrayId=", "phenotypeFileType=", 
					"help", "parallel=", "parallelAll", "addToDB", 
					"callMethodID=", "comment=", "onlyOriginal192","onlyOriginal96", "subSample=" , 
					"subSampleLikePhenotype=", "subsampleTest=", "complement", "onlyBelowLatidue=", 
					"onlyAboveLatidue=", "srInput=", "sr","srOutput=", "srPar=","srSkipFirstRun",
					"permTest=", "savePermutations", "permutationFilter=", "testRobustness"]
	try:
		opts, args=getopt.getopt(sys.argv[1:], "o:c:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
		phenotypeFileType=1
		outputFile=None
	delim=","
	missingVal="NA"
	help=0
	withArrayIds=1
	parallel=None
	parallelAll=False
	addToDB=False
	callMethodID=None
	comment=""
	subSample=None
	onlyOriginal96=False
	onlyOriginal192 = False
	subSampleLikePhenotype = None
	subsampleTest = False
	numSubSamples = None
	complement = False
	onlyBelowLatidue = None
	onlyAboveLatidue = None

	sr = False
	srOutput = False
	srInput = False
	srSkipFirstRun = False
	srTopQuantile = 0.95
	srWindowSize = 30000
	
	permTest = None
	savePermutations = False
	permutationFilter = 1.0
	
	testRobustness = False

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help=1
			print __doc__
		elif opt in ("-a", "--withArrayId"):
			withArrayIds=int(arg)
		elif opt in ("-o", "--outputFile"):
			outputFile=arg
		elif opt in ("--phenotypeFileType"):
			phenotypeFileType=int(arg)
		elif opt in ("--parallel"):
			parallel=arg
		elif opt in ("--parallelAll"):
			parallelAll=True
		elif opt in ("--addToDB"):
			addToDB=True
  		elif opt in ("--onlyOriginal96"):
			onlyOriginal96=True
  		elif opt in ("--onlyOriginal192"):
			onlyOriginal192=True
		elif opt in ("--complement"):
			complement=True
		elif opt in ("--subSample"):
			subSample=int(arg)
		elif opt in ("--subsampleTest"):
			subsampleTest = True
			l = arg.split(",")
			subSample=int(l[0])
			numSubSamples=int(l[1])
		elif opt in ("--onlyBelowLatidue"):
			onlyBelowLatidue=float(arg)
		elif opt in ("--onlyAboveLatidue"):
			onlyAboveLatidue=float(arg)
		elif opt in ("--subSampleLikePhenotype"):
			subSampleLikePhenotype=int(arg)
		elif opt in ("--callMethodID"):
			callMethodID=int(arg)
		elif opt in ("--comment"):
			comment=arg
		elif opt in ("-d", "--delim"):
			delim=arg
		elif opt in ("-m", "--missingval"):
			missingVal=arg
		elif opt in ("--sr"):
			sr = True
		elif opt in ("--testRobustness"):
			testRobustness = True
		elif opt in ("--permTest"):
			permTest = int(arg)
		elif opt in ("--savePermutations"):
			savePermutations = True
		elif opt in ("--permutationFilter"):
			permutationFilter = float(arg)
		elif opt in ("--srSkipFirstRun"):
			srSkipFirstRun = True
		elif opt in ("--srInput"):
			srInput = arg
		elif opt in ("--srOutput"):
			srOutput = arg
		elif opt in ("--srPar"):
			vals = arg.split(",")
			srTopQuantile = float(vals[0]) 
			srWindowSize = int(vals[1]) 
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

	snpsDataFile=args[0]
	phenotypeDataFile=args[1]

	print "Kruskal-Wallis is being set up with the following parameters:"
	print "phenotypeDataFile:",phenotypeDataFile
	print "snpsDataFile:",snpsDataFile
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "onlyOriginal96:",onlyOriginal96
	print "onlyOriginal192:",onlyOriginal192
	print "onlyBelowLatidue:",onlyBelowLatidue
	print "onlyAboveLatidue:",onlyAboveLatidue
	print "subSampleLikePhenotype:",subSampleLikePhenotype
	print "subsampleTest:",subsampleTest
	print "numSubSamples:",numSubSamples
	print "subSample:",subSample
	print "sr:",sr
	print "srSkipFirstRun:",srSkipFirstRun
	print "srInput:",srInput
	print "srOutput:",srOutput
	print "srTopQuantile:",srTopQuantile
	print "srWindowSize:",srWindowSize
	print "permTest:",permTest
	print "savePermutations:",savePermutations
	print "permutationFilter:",permutationFilter
	print "testRobustness:",testRobustness
	

	def runParallel(phenotypeIndex,id=""):
		#Cluster specific parameters
		phed=phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
		phenName=phed.getPhenotypeName(phenotypeIndex)
		phenName=phenName.replace("/", "_div_")
		phenName=phenName.replace("*", "_star_")
		outputFile=resultDir+"KW_"+parallel+"_"+phenName+id

		shstr="""#!/bin/csh
#PBS -l walltime=100:00:00
#PBS -l mem=4g 
#PBS -q cmb
"""
		
		shstr+="#PBS -N K"+phenName+"_"+parallel+"\n"
		shstr+="set phenotypeName="+parallel+"\n"
		shstr+="set phenotype="+str(phenotypeIndex)+"\n"
		shstr+="(python "+scriptDir+"KW.py -o "+outputFile+" "
		shstr+=" -a "+str(withArrayIds)+" "			
		if subSample:
			shstr+=" --subSample="+str(subSample)+" "			
		elif onlyOriginal96:
			shstr+=" --onlyOriginal96 "			
		elif onlyOriginal192:
			shstr+=" --onlyOriginal192 "
		if onlyBelowLatidue:
			shstr+=" --onlyBelowLatidue="+str(onlyBelowLatidue)+" "
		elif onlyAboveLatidue:
			shstr+=" --onlyAboveLatidue="+str(onlyAboveLatidue)+" "
		if complement: 			
			shstr+=" --complement "
		if permTest:
			shstr+=" --permTest="+str(permTest)+" "
			if savePermutations:
				shstr+=" --savePermutations "
		
		shstr+=" --permutationFilter="+str(permutationFilter)+" "
		if testRobustness:
			shstr+=" --testRobustness "
			
		if sr:
			shstr += " --sr "			
			if not srOutput:
				output = resultDir+"KW_"+parallel+"_"+phenName+".sr.pvals"				
			shstr += " --srOutput="+str(output)+" "
			if srSkipFirstRun:
				if not srInput:
					output = resultDir+"KW_"+parallel+"_"+phenName+".pvals"
				shstr += " --srInput="+str(output)+" "
				shstr += " --srSkipFirstRun "				
			shstr += " --srPar="+str(srTopQuantile)+","+str(srWindowSize)+" "


		shstr+=snpsDataFile+" "+phenotypeDataFile+" "+str(phenotypeIndex)+" "
		shstr+="> "+outputFile+"_job"+".out) >& "+outputFile+"_job"+".err\n"

		f=open(parallel+".sh", 'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	if parallel:  #Running on the cluster..
		if parallelAll:
			phed=phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
			for phenotypeIndex in phed.phenIds:
				runParallel(phenotypeIndex)
		elif subsampleTest:
			phenotypeIndex=int(args[2])
			for i in range(0,numSubSamples):
				runParallel(phenotypeIndex,id="_r"+str(subSample)+"_"+str(i))
		else:
			phenotypeIndex=int(args[2])
			runParallel(phenotypeIndex)
		return
	else:
		phenotypeIndex=int(args[2])


	print "phenotypeIndex:",phenotypeIndex
	print "output:",outputFile
	print "\nStarting program now!\n"


	#Load phenotype file
	phed=phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
	
	#If onlyOriginal96, then remove all other phenotypes..
	if onlyOriginal96: 
		print "Filtering for the first 96 accessions"
		original_96_ecotypes = phenotypeData._getFirst96Ecotypes_()
		original_96_ecotypes = map(str,original_96_ecotypes)
		keepEcotypes = []
		if complement:
			for acc in phed.accessions:
				if not acc in original_96_ecotypes:
					keepEcotypes.append(acc)
		else:
			keepEcotypes = original_96_ecotypes
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)
		
	if onlyOriginal192: 
		print "Filtering for the first 192 accessions"
		original_192_ecotypes = phenotypeData._getFirst192Ecotypes_()
		original_192_ecotypes = map(str,original_192_ecotypes)
		keepEcotypes = []
		if complement:
			for acc in phed.accessions:
				if not acc in original_192_ecotypes:
					keepEcotypes.append(acc)
		else:
			keepEcotypes = original_192_ecotypes
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)
	
	if onlyBelowLatidue:
		print "Filtering for the accessions which orginate below latitude",onlyBelowLatidue
		eiDict = phenotypeData._getEcotypeIdInfoDict_()
		print eiDict
		keepEcotypes = []
		for acc in phed.accessions:
			acc = int(acc)
			if eiDict.has_key(acc) and eiDict[acc][2] and eiDict[acc][2]<onlyBelowLatidue:
				keepEcotypes.append(str(acc))
			elif eiDict.has_key(acc) and eiDict[acc][2]==None:
				keepEcotypes.append(str(acc))
				
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)

	elif onlyAboveLatidue:
		print "Filtering for the accessions which orginate above latitude",onlyAboveLatidue
		eiDict = phenotypeData._getEcotypeIdInfoDict_()
		print eiDict
		keepEcotypes = []
		for acc in phed.accessions:
			acc = int(acc)
			if eiDict.has_key(acc) and eiDict[acc][2] and eiDict[acc][2]>onlyAboveLatidue:
				keepEcotypes.append(str(acc))
			elif eiDict.has_key(acc) and eiDict[acc][2]==None:
				keepEcotypes.append(str(acc))
				
		phed.filterAccessions(keepEcotypes)
		print "len(phed.accessions)", len(phed.accessions)

	
	if subSampleLikePhenotype:
		p_name = phed.getPhenotypeName(subSampleLikePhenotype)
		print "Picking sample as in",p_name
		ecotypes = phed.getNonNAEcotypes(subSampleLikePhenotype)
		print ecotypes
		phed.filterAccessions(ecotypes)
		print "len(phed.accessions)", len(phed.accessions)


	if subSample: 
		sample_ecotypes = []
		ecotypes = phed.getNonNAEcotypes(phenotypeIndex)
		sample_ecotypes = random.sample(ecotypes,subSample)			
		phed.filterAccessions(sample_ecotypes)
		print "len(phed.accessions)", len(phed.accessions)
		
	sys.stdout.write("Finished prefiltering phenotype accessions.\n")
	sys.stdout.flush()
	
	
	
	#Load genotype file
	snpsds=dataParsers.parseCSVData(snpsDataFile, format = 1, deliminator = delim, missingVal = missingVal, withArrayIds = withArrayIds)


	#Checking overlap between phenotype and genotype accessions. 
	phenotype=phed.getPhenIndex(phenotypeIndex)
	accIndicesToKeep=[]			
	phenAccIndicesToKeep=[]
	numAcc=len(snpsds[0].accessions)
	sys.stdout.write("Removing accessions which do not have a phenotype value for "+phed.phenotypeNames[phenotype]+".")
	sys.stdout.flush()
	for i in range(0, len(snpsds[0].accessions)):
		acc1=snpsds[0].accessions[i]
		for j in range(0, len(phed.accessions)):
			acc2=phed.accessions[j]
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
	print numAcc-len(accIndicesToKeep), "accessions removed, leaving", len(accIndicesToKeep), "accessions in all."
		
	print "Filtering phenotype data."
	phed.removeAccessions(phenAccIndicesToKeep) #Removing accessions that don't have genotypes or phenotype values
	
	#Ordering accessions according to the order of accessions in the genotype file
	accessionMapping=[]
	i=0
	for acc in snpsds[0].accessions:
		if acc in phed.accessions:
			accessionMapping.append((phed.accessions.index(acc), i))
			i+=1
	phed.orderAccessions(accessionMapping)

		#Filtering monomorphic
	print "Filtering monomorphic SNPs"
	for snpsd in snpsds:
		print "Removed", str(snpsd.filterMonoMorphicSnps()), "Snps"

	#Converting format to 01
	newSnpsds=[]
	sys.stdout.write("Converting data format")
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		newSnpsds.append(snpsd.getSnpsData())
	print ""
	
	#Double check genotype file:
	problems = 0
	for i in range(0,len(newSnpsds)):
		snpsd = newSnpsds[i]
		for j in range(0,len(snpsd.snps)):
			snp = snpsd.snps[j]
			sc = snp.count(0)
			if sc==0 or sc==len(snp):
				print "Problem in file found at chr,pos",(i+1),",",snpsd.positions[i]
				problems += 1
	if problems >0:
		print "Genotype file appears to have potential problems"
	else:
		print "Genotype file appears to be good"

	if permTest:
		print "Starting a permutation test"
		allSNPs = []
		for snpsd in newSnpsds:
			allSNPs += snpsd.snps
		phenVals = phed.getPhenVals(phenotypeIndex)
		test_type = "KW"
		if phed.isBinary(phenotypeIndex):
			test_type = "Fisher"
			permTest = 100	
		_perm_test_(allSNPs,phenVals,permTest,outputFile, test_type=test_type,savePermutations=savePermutations, filter=permutationFilter)
		sys.exit(0)
	
	if testRobustness:
		print "Starting a robustness test"
		allSNPs = []
		for snpsd in newSnpsds:
			allSNPs += snpsd.snps
		phenVals = phed.getPhenVals(phenotypeIndex)
		test_type = "KW"
		if phed.isBinary(phenotypeIndex):
			test_type = "Fisher"
		_robustness_test_(allSNPs,phenVals,outputFile, test_type=test_type, filter=permutationFilter)
		sys.exit(0)
		

	sys.stdout.flush()
	print "sr:",sr, ", srSkipFirstRun:",srSkipFirstRun
	if (not sr) or (sr and not srSkipFirstRun):
		#Writing files
		if env.user=="bjarni":
			tempfile.tempdir='/tmp'
		(fId, phenotypeTempFile)=tempfile.mkstemp()
		os.close(fId)
		(fId, genotypeTempFile)=tempfile.mkstemp()
		os.close(fId)
		
		phed.writeToFile(phenotypeTempFile, [phenotype])	
		sys.stdout.write("Phenotype file written\n")
		sys.stdout.flush()
		snpsDataset=snpsdata.SNPsDataSet(newSnpsds, [1, 2, 3, 4, 5])
		decoder={1:1, 0:0,-1:'NA'}	
		snpsDataset.writeToFile(genotypeTempFile, deliminator = delim, missingVal = missingVal, withArrayIds = 0, decoder = decoder)
		sys.stdout.write("Genotype file written\n")
		sys.stdout.flush()
	
		phenotypeName=phed.getPhenotypeName(phenotypeIndex)
	
		rDataFile=outputFile+".rData"
		pvalFile=outputFile+".pvals"
		#Is the phenotype binary?
		binary=phed.isBinary(phenotypeIndex)
		rstr=_generateRScript_(genotypeTempFile, phenotypeTempFile, rDataFile, pvalFile, name = phenotypeName, binary = binary)
		rFileName=outputFile+".r"
		f=open(rFileName, 'w')
		f.write(rstr)
		f.close()
		outRfile=rFileName+".out"
		errRfile=rFileName+".err"
		print "Running R file:"
		cmdStr="(R --vanilla < "+rFileName+" > "+outRfile+") >& "+errRfile
		sys.stdout.write(cmdStr+"\n")
		sys.stdout.flush()	
		gc.collect() 
		os.system(cmdStr)
		#print "Emma output saved in R format in", rDataFile
		print "Generating a GW plot."
		res = gwaResults.Result(pvalFile,name="KW_"+phenotypeName, phenotypeID=phenotypeIndex)
		res.negLogTransform()
		pngFile = pvalFile+".png"
		plotResults.plotResult(res,pngFile=pngFile,percentile=90,type="pvals",ylab="$-$log$_{10}(p)$", plotBonferroni=True,usePylab=False)	
		srInput = pvalFile
		
	else:
		print "Skipping first stage analysis."
		sys.stdout.flush()

	if sr:
		_secondRun_(srOutput,srInput,srTopQuantile,srWindowSize,newSnpsds,phed,phenotypeIndex,binary=binary)
		print "Generating second run GW plot."
		res = gwaResults.Result(srInput,name="KW_"+phenotypeName, phenotypeID=phenotypeIndex)
		res.negLogTransform()
		srRes = gwaResults.Result(srOutput,name="KW_SR_"+phenotypeName, phenotypeID=phenotypeIndex)
		srRes.negLogTransform()
		srPngFile = pvalFile+".sr.png"
		plotResults.plotResultWithSecondRun(res,srRes,pngFile=srPngFile,ylab="$-$log$_{10}(p)$", plotBonferroni=True)	
	
	
def _generateRScript_(genotypeFile, phenotypeFile, rDataFile, pvalFile, name = None, binary = False, chiSquare = False):
	
	rstr='data250K <- read.csv("'+str(genotypeFile)+'", header=TRUE);\n'
	rstr+='phenotData <- read.csv("'+str(phenotypeFile)+'",header=TRUE);\n'
	rstr+="""
phenMat <- t(as.matrix(as.matrix(phenotData)[,2]));
res <- list();
pvals <- c();
positions <- c();
chrs <- c();
maf <- c();
marf <- c();
for (chr in (1:5)){
  mat250K <- as.matrix(data250K);
  mat250K <- mat250K[mat250K[,1]==chr,];
  pos <- mat250K[,2];
  mat250K <- mat250K[,3:length(mat250K[1,])];
  res[[chr]] <- list();
  res[[chr]][["ps"]] <- c();
  res[[chr]][["stat"]] <- c();
  for (j in (1:length(mat250K[,1]))){
"""
	if binary: 
		if chiSquare:
			rstr+="	v <- chisq.test(as.vector(phenMat),as.vector(mat250K[j,]));"
			rstr+="""
	res[[chr]]$ps[j] <- as.double(v$p.value);
	res[[chr]]$stat[j] <- as.double(v[1]);		
"""
		else:
			rstr+="	v <- fisher.test(as.vector(phenMat),as.vector(mat250K[j,]));"			
			rstr+="""
	res[[chr]]$ps[j] <- as.double(v$p.value);
	res[[chr]]$stat[j] <- as.double(v$p.value);		
"""

	else:
		rstr+="	v <- kruskal.test(as.vector(phenMat),as.vector(mat250K[j,]));"
		rstr+="""
	res[[chr]]$ps[j] <- as.double(v$p.value);
	res[[chr]]$stat[j] <- as.double(v[1]);		
"""
	rstr+="""
  }  
  res[[chr]][["pos"]] <- pos;
  res[[chr]][["chr_pos"]] <- pos;

  pvals <- append(pvals,res[[chr]]$ps);
  positions <- append(positions,pos);
  chrs <- append(chrs,rep(chr,length(pos)));

  af <- c();
  arf <- c();
  for(i in (1:length(mat250K[,1]))){
	f = factor(mat250K[i,]);
	freq <- summary(f)[[1]]
	alleleCount <- length(f)
	v <- freq/alleleCount;
	af[i] <- min(freq,alleleCount-freq);
	arf[i] <- min(v,1-v);
  }
  res[[chr]][["maf"]] <- af;
  res[[chr]][["marf"]] <- arf;

  marf <- append(marf,arf);
  maf <- append(maf,af);
}

#write to a pvalue-file
l <- list();
l[["Chromasome"]]<-chrs;
l[["Positions"]]<-positions;
l[["Pvalues"]]<-pvals;
l[["MARF"]]<-marf;
l[["MAF"]]<-maf;
dl <- as.data.frame(l)
"""
	rstr+='write.table(dl,file="'+pvalFile+'", sep=", ", row.names = FALSE);\n'		
	rstr+="""
#Save data as R object.
res[[2]]$pos <- res[[2]]$pos+res[[1]]$pos[length(res[[1]]$pos)];
res[[3]]$pos <- res[[3]]$pos+res[[2]]$pos[length(res[[2]]$pos)];
res[[4]]$pos <- res[[4]]$pos+res[[3]]$pos[length(res[[3]]$pos)];
res[[5]]$pos <- res[[5]]$pos+res[[4]]$pos[length(res[[4]]$pos)];

res[["ylim"]] <- c(min(min(-log(res[[1]]$ps)),min(-log(res[[2]]$ps)),min(-log(res[[3]]$ps)),min(-log(res[[4]]$ps)),min(-log(res[[5]]$ps))), max(max(-log(res[[1]]$ps)),max(-log(res[[2]]$ps)),max(-log(res[[3]]$ps)),max(-log(res[[4]]$ps)),max(-log(res[[5]]$ps))));

res[["xlim"]] <- c(min(res[[1]]$pos),max(res[[5]]$pos));

res[["FRI"]]<-c(269026+res[[3]]$pos[length(res[[3]]$pos)],271503+res[[3]]$pos[length(res[[3]]$pos)]);
res[["FLC"]]<-c(3173498+res[[4]]$pos[length(res[[4]]$pos)],3179449+res[[4]]$pos[length(res[[4]]$pos)]);


"""
	if name:
		if binary:
			if chiSquare:
				rstr+='res[["lab"]]= "Chi-Square p-values for '+name+'";\n'
			else:
				rstr+='res[["lab"]]= "Fisher p-values for '+name+'";\n'
		else:
			rstr+='res[["lab"]]= "Kruskal Wallis p-values for '+name+'";\n'
	else:
		rstr+='res[["lab"]]="";\n'
	rstr+='save(file="'+rDataFile+'",res);\n'
	return rstr	
	

def run_kw(snps,phenotypeValues,verbose=False):
	print "Running KW on",len(snps),"snps, and",len(phenotypeValues),"phenotype values."
	pvals = []
	for snp in snps:
		res = r.kruskal_test(phenotypeValues,snp)
		#print snp,res,phenotypeValues
		#pdb.set_trace()
		pvals.append(res["p.value"])
	return pvals

	
def run_fet(snps,phenotypeValues,verbose=False):
	"""
	Fisher's exact test.
	"""
	print "Running Fisher's exact test on",len(snps),"snps, and",len(phenotypeValues),"phenotype values."
	pvals = []
	for snp in snps:
		#if snp.count(0)!=0 and snp.count(1)!=0:
		res = r.fisher_test(phenotypeValues,snp)
		#print snp,res,phenotypeValues
		#pdb.set_trace()
		pval = res["p.value"]
		#if pval == 1:
		#	print phenotypeValues,snp
		#else:
		#pval = -1.0
		#print "Found a monomorphic SNP.."
		pvals.append(pval)
	return pvals
	

def _secondRun_(srOutput,srInput,srTopQuantile,srWindowSize,snpsds,phed,p_i,binary = False):
	"""
	Sets up a second run, using and/or on SNPs.
	"""
	print "Preparing a second stage analyzis."
	sys.stdout.flush()
	#Collect top n% SNPs.
	print "Loading first stage results."
	sys.stdout.flush()
	result = gwaResults.SNPResult(resultFile=srInput,snpsds=snpsds,phenotypeID=p_i)
	# - filter results
	result.negLogTransform()
		
	res = SecondStageAnalysis.retrieveSecondRunSNPs(result,srTopQuantile,srWindowSize)
	snps = res['snps']
	snpTypes = res['snpTypes']
	chromosomes = res['chromosomes']
	positions = res['positions']
	marfs = res['marfs']
	mafs = res['mafs']


	if binary:
		#Run Fisher's Exact Test:
		sys.stdout.write( "Running Fisher's exact test (in R).\n")
		sys.stdout.flush()
		phenVals = phed.getPhenVals(p_i)
		pvals = run_fet(snps,phenVals)
		#print zip(positions,pvals,genotype_var_perc)
	else:
		#Run Emma:
		sys.stdout.write( "Running Kruskal-Wallis (in R).\n")
		sys.stdout.flush()
		phenVals = phed.getPhenVals(p_i)
		pvals = run_kw(snps,phenVals)
		#print zip(positions,pvals,genotype_var_perc)
	
	sys.stdout.write("Writing results to file.\n")
	sys.stdout.flush()
	#Write results to file!
	f = open(srOutput,"w")
	f.write("Chromosome,position,p-value,marf,maf,snpType,second_pos\n")
	for (chr,pos,pval,marf,maf,snpType) in zip(chromosomes,positions,pvals,marfs,mafs,snpTypes):
		f.write(str(chr)+","+str(pos[0])+","+str(pval)+","+str(marf)+","+str(maf)+","+str(snpType)+","+str(pos[1])+"\n")
	f.close()

	#Plot results?
	


def _kruskal_wallis_(snps,phenVals):
	from scipy import stats
	import numpy
	"""
	Uses Scipy to perform KW.
	"""
	ds = []
	ps = []
	for snp in snps:
		group1 = []
		group2 = []
		for (nt,val) in zip(snp,phenVals):
			if int(nt)==0:
				group1.append(val)
			elif int(nt)==1:
				group2.append(val)
			else:
				raise Exception, "An invalid nt: "+str(nt)
		g1 = numpy.array(group1)
		g2 = numpy.array(group2)
		(d,p) = stats.kruskal(g1, g2)
		ds.append(d)
		ps.append(p)
	return {"ps":ps,"ds":ds}
		
		
	
	
def _perm_test_(all_snps,phenVals,numPerm,outputFile,filter=0.1,test_type = "KW",savePermutations=False,useSameSnps=False):

	def _calc_statistics_(pvals,exp_quantiles,exp_median=0.5,exp_pvals=None):
		m = analyzePhenotype._calcMedian_(pvals,exp_median)
		ks_res = analyzePhenotype._calcKS_(pvals,exp_pvals)
		s = analyzePhenotype._estLogSlope_(pvals,exp_pvals)-1.0
		ks_stat = ks_res["D"]
		ks_pvalue = ks_res["p.value"]
		quantiles = analyzePhenotype._getQuantiles_(pvals, 1000)
		#exp_quantiles = analyzePhenotype.__getExpectedPvalueQuantiles__(1000)
		a = analyzePhenotype._estAreaBetweenCurves_(quantiles,exp_quantiles)
		
		return (m,a,ks_stat,ks_pvalue,s)
	
	if filter <1.0:
		snps = random.sample(all_snps,int(len(all_snps)*filter))
		print "Number of SNPs:",len(snps)
	else:
		snps = all_snps 

	#Calc norm stats, and est. p-value 
#	print "running old KW"
#	t1 = time.time()
#	pvals = analyzeHaplotype._run_kw_(snps,phenVals)
#	t2 = time.time()
#	print "Took",t2-t1,"seconds."
	if test_type=="KW":
		print "running KW"
		t1 = time.time()
		true_pvals = util.kruskal_wallis(snps,phenVals)["ps"]
		t2 = time.time()
		print "Took",t2-t1,"seconds."
	elif test_type=="Fisher":
		print "running Fisher's exact test"
		t1 = time.time()
		true_pvals = run_fet(snps,phenVals)
		t2 = time.time()
		print "Took",t2-t1,"seconds."
		

	
	
	perm_pvalues_list = []
	for i in range(0,numPerm):#For every perm
		if filter <1.0:
			snps = random.sample(all_snps,int(len(all_snps)*filter))
			print "Number of SNPs:",len(snps)	
		print i
		random.shuffle(phenVals) #Permute phenotype
		#pvals = analyzeHaplotype._run_kw_(snps,phenVals)	#Run KW
		if test_type=="KW":
			print "running KW"
			t1 = time.time()
			pvals = util.kruskal_wallis(snps,phenVals)["ps"]
			t2 = time.time()
			print "Took",t2-t1,"seconds."
		elif test_type=="Fisher":
			print "running Fisher's exact test"
			t1 = time.time()
			pvals = run_fet(snps,phenVals)
			t2 = time.time()
			print "Took",t2-t1,"seconds."
		perm_pvalues_list.append(pvals)

	
	print "Combining p-values"
	quantiles = []
	all_pvals = []
	for pvals in perm_pvalues_list:
		for pval in pvals:
			all_pvals.append(pval)
	print len(all_pvals),"permuted pvals in all"
	quantiles = analyzePhenotype._getQuantiles_(all_pvals, 1000)
	print "len(quantiles):", len(quantiles)
	exp_median = (quantiles[499]+quantiles[500])/2.0

	(true_m,true_a,true_ks_stat,true_ks_pvalue,true_s) = _calc_statistics_(true_pvals,quantiles,exp_median,all_pvals)

	m_list = []
	a_list = []
	ks_stat_list = []
	ks_pvalue_list = []
	s_list = []
	for i in range(0,numPerm):
		pvals = perm_pvalues_list[i]
		(m,a,ks_stat,ks_pvalue,s) = _calc_statistics_(pvals,quantiles,exp_median,all_pvals) #Calc. statistic
		m_list.append(m)
		a_list.append(a)
		s_list.append(s)
		ks_stat_list.append(ks_stat)
		ks_pvalue_list.append(ks_pvalue)
	
	del all_pvals,quantiles

		
	if savePermutations:
		permOutputFile = outputFile+".perm.pvals"
		print "Writing to",permOutputFile
		f = open(permOutputFile,"w")
		i = 0
		for pvals in perm_pvalues_list:
			pvals_str = map(str,pvals)
			f.write(",".join(pvals_str)+"\n")
		print "Done writing to",permOutputFile
	
	f.close()


	#Output results
	outputFile = outputFile+".perm.stat.txt"
	f = open(outputFile,"w")
	f.write("Perm_nr, median, area, ks_stat, s_stat \n")
	for i in range(0,numPerm):
		str_l = map(str,[i, m_list[i],a_list[i],ks_stat_list[i],s_list[i]])
		f.write(", ".join(str_l)+"\n")
	
	f.write("\n"+"Observed values: "+str((true_m,true_a,true_ks_stat,true_s))+"\n")

	pvals = [0.0,0.0,0.0,0.0]
	
	#M stat p-value (two sided)
	#Assuming symm. dist.
	for i in range(0,numPerm):
		if abs(true_m) <= abs(m_list[i]):
			pvals[0]+=1.0/numPerm
	
	#A stat p-value (one tailed)
	for i in range(0,numPerm):
		if true_a <= a_list[i]:
			pvals[1]+=1.0/numPerm
		

	#KS stat p-value (one tailed)
	for i in range(0,numPerm):
		if true_ks_stat <= ks_stat_list[i]:
			pvals[2]+=1.0/numPerm
		
	#S stat p-value (one tailed)
	for i in range(0,numPerm):
		if abs(math.log(true_s+1.0)) <= abs(math.log(s_list[i]+1.0)):
			pvals[3]+=1.0/numPerm
		


	for i in range(0,len(pvals)):
		if pvals[i] == 0.0:
			pvals[i] = 0.5*(1.0/numPerm)
			
	str_pvals = map(str,pvals)
	f.write("\n"+"Estimated p-values: "+",".join(str_pvals)+"\n")
	f.close()
	
	#Plot results
	pngFile_median = outputFile+".perm.m.png"
	pngFile_area = outputFile+".perm.a.png"
	pngFile_ks = outputFile+".perm.ks.png"
	pngFile_s = outputFile+".perm.s.png"
	
	def _getBinning_(n_bins,min_val,max_val):
		bins = []
		delta = (max_val-min_val)/n_bins
		start_val = min_val-delta*0.5
		for i in range(0,n_bins+2):
			bins.append(start_val+delta*i)
		return (bins,delta)
		
	n_bins = 20+int(4*(math.log(numPerm)))
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	min_val = min(min(m_list),true_m)
	max_val = max(max(m_list),true_m)
	(bins,delta) = _getBinning_(n_bins,min_val,max_val)
	print (bins,delta)

	plt.figure(figsize=(10,7))
	plt.hist(m_list+[true_m], bins = bins)#, range=[start_val,end_val])
	plt.hist([true_m], bins = bins)#, range=[start_val,end_val])
	plt.savefig(pngFile_median, format = "png")
	plt.legend()
	plt.clf()

	min_val = min(min(a_list),true_a)
	max_val = max(max(a_list),true_a)
	(bins,delta) = _getBinning_(n_bins,min_val,max_val)
	print (bins,delta)
	plt.figure(figsize=(10,7))
	plt.hist(a_list+[true_a], bins = bins)
	plt.hist([true_a], bins = bins)
	plt.savefig(pngFile_area, format = "png")
	plt.clf()

	min_val = min(min(ks_stat_list),true_ks_stat)
	max_val = max(max(ks_stat_list),true_ks_stat)
	(bins,delta) = _getBinning_(n_bins,min_val,max_val)
	print (bins,delta)
	plt.figure(figsize=(10,7))
	plt.hist(ks_stat_list+[true_ks_stat], bins = bins)
	plt.hist([true_ks_stat], bins = bins)
	plt.savefig(pngFile_ks, format = "png")
	plt.clf()

	min_val = min(min(s_list),true_s)
	max_val = max(max(s_list),true_s)
	(bins,delta) = _getBinning_(n_bins,min_val,max_val)
	print (bins,delta)
	plt.figure(figsize=(10,7))
	plt.hist(s_list+[true_s], bins = bins)
	plt.hist([true_s], bins = bins)
	plt.savefig(pngFile_s, format = "png")
	plt.clf()


		
	
def _robustness_test_(all_snps,phenVals,outputFile,filter=0.1,test_type = "KW",):
	"""
	Leave one out test..
	"""
	
	new_all_snps = []
	for snp in all_snps:
		if snp.count(0)>1 and snp.count(1)>1:
			new_all_snps.append(snp)
	print "Filtered",len(all_snps)-len(new_all_snps)," with minor allele count <2."
	all_snps = new_all_snps

	if filter <1.0:
		snps = random.sample(all_snps,int(len(all_snps)*filter))
		print "Number of SNPs:",len(snps)
	else:
		snps = all_snps 

	if test_type=="KW":
		print "running KW"
		t1 = time.time()
		true_pvals = util.kruskal_wallis(snps,phenVals)["ps"]
		t2 = time.time()
		print "Took",t2-t1,"seconds."
	elif test_type=="Fisher":
		print "running Fisher's exact test"
		t1 = time.time()
		true_pvals = run_fet(snps,phenVals)
		t2 = time.time()
		print "Took",t2-t1,"seconds."
	
	log_true_pvals = []
	for pval in true_pvals:
		log_true_pvals.append(-math.log(pval,10))
	
	perm_pvalues_list = []
	for i in range(0,len(phenVals)):
		newPhenvals = phenVals[:]
		newPhenvals.pop(i)
				
		newSNPs = []
		for snp in snps:
			newSNP = snp[:]
			newSNP.pop(i)
			newSNPs.append(newSNP)
		
		print i
		if test_type=="KW":
			print "running KW"
			t1 = time.time()
			pvals = util.kruskal_wallis(newSNPs,newPhenvals)["ps"]
			t2 = time.time()
			print "Took",t2-t1,"seconds."
		elif test_type=="Fisher":
			print "running Fisher's exact test"
			t1 = time.time()
			pvals = run_fet(newSNPs,newPhenvals)
			t2 = time.time()
			print "Took",t2-t1,"seconds."
		perm_pvalues_list.append(pvals)
		
		
	
	delta_pvals_list = []
	delta_log_pvals_list = []
	for perm_pvals in perm_pvalues_list:
		log_pvals = []
		delta_pvals = []
		delta_log_pvals = []
		for i in range(0,len(true_pvals)):
			pval = perm_pvals[i]
			true_pval = true_pvals[i]
			delta_pvals.append(true_pval-pval)

			log_true_pval = log_true_pvals[i]
			if pval > 0.0:
				log_pval = -math.log(pval,10)
			else:
				print "Damn those random 0 prob. events: event #", i
				log_pval = -math.log(true_pval,10)
				
			log_pvals.append(log_pval)
			delta_log_pvals.append(log_true_pval-log_pval)
		
		delta_pvals_list.append(delta_pvals)
		delta_log_pvals_list.append(delta_log_pvals)
	
	sd_log_pvals = []
	sd_pvals = []
	t_delta_log_pvals_list = map(list,zip(*delta_log_pvals_list))
	t_delta_pvals_list = map(list,zip(*delta_pvals_list))
	for i in range(0,len(true_pvals)):
		sd_log_pvals.append(util.calcSD(t_delta_log_pvals_list[i]))
		sd_pvals.append(util.calcSD(t_delta_pvals_list[i]))
	
	
	#Write SDs out to file, to be able to replot, or plot together with other methods... etc
	import csv
	sd_log_pval_file = outputFile+".rob.log_pvals_sd"
	f = open(sd_log_pval_file,"w")
	w = csv.writer(f)
	w.writerow(["log_true_pval","sd_log_pvals"])
	l = zip(log_true_pvals,sd_log_pvals)
	w.writerows(l)
	f.close()

	
	#Plot things....
	pngFile_log_pvals = outputFile+".rob.log_pval.png"
	pngFile_pval = outputFile+".rob.pval.png"
	pngFile_sd_log_pval = outputFile+".rob.sd_log_pval.png"
	pngFile_sd_pval = outputFile+".rob.sd_pval.png"


	min_val = min(true_pvals)
	max_val = max(true_pvals)
	val_range = max_val-min_val

	min_log_val = min(log_true_pvals)
	max_log_val = max(log_true_pvals)
	log_val_range = max_val-min_val


	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plt.figure(figsize=(10,7))
	max_perm_val = 0
	min_perm_val = 0
	for i in range(0,len(perm_pvalues_list)):
		delta_log_pvals = delta_log_pvals_list[i]
		plt.plot(log_true_pvals,delta_log_pvals,"b.")
		max_perm_val = max(max_perm_val,max(delta_log_pvals))
		min_perm_val = min(min_perm_val,min(delta_log_pvals))
	perm_val_range = max_perm_val - min_perm_val
	plt.axis([min_log_val-0.02*log_val_range, max_log_val+0.02*log_val_range, min_perm_val-0.02*perm_val_range, max_perm_val+0.02*perm_val_range])
	plt.savefig(pngFile_log_pvals, format = "png")

	plt.figure(figsize=(10,7))
	max_perm_val = 0
	min_perm_val = 0
	for i in range(0,len(perm_pvalues_list)):
		delta_pvals = delta_pvals_list[i]
		plt.plot(true_pvals,delta_pvals,"b.")
		max_perm_val = max(max_perm_val,max(delta_pvals))
		min_perm_val = min(min_perm_val,min(delta_pvals))
	perm_val_range = max_perm_val - min_perm_val
	plt.axis([min_val-0.02*val_range, max_val+0.02*val_range, min_perm_val-0.02*perm_val_range, max_perm_val+0.02*perm_val_range])
	plt.savefig(pngFile_pval, format = "png")

	plt.figure(figsize=(10,7))
	max_sd_log_pval = max(sd_log_pvals)
	min_sd_log_pval = min(sd_log_pvals)
	sd_val_range = max_sd_log_pval-min_sd_log_pval
	plt.plot(log_true_pvals,sd_log_pvals,"b.")
	plt.axis([min_log_val-0.02*log_val_range, max_log_val+0.02*log_val_range, min_sd_log_pval-0.02*sd_val_range, max_sd_log_pval+0.02*sd_val_range])
	plt.savefig(pngFile_sd_log_pval, format = "png")

	plt.figure(figsize=(10,7))
	max_sd_pval = max(sd_pvals)
	min_sd_pval = min(sd_pvals)
	sd_val_range = max_sd_pval-min_sd_pval
	plt.plot(true_pvals,sd_pvals,"b.")
	plt.axis([min_val-0.02*val_range, max_val+0.02*val_range, min_sd_pval-0.02*sd_val_range, max_sd_pval+0.02*sd_val_range])
	plt.savefig(pngFile_sd_pval, format = "png")
	
	
	
	
if __name__=='__main__':
	_run_()
	print "Done!"
