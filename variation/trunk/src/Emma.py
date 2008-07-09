#!/usr/bin/env python2.5
"""
Usage: Emma.py [OPTIONS] -o R_FILE SNPS_DATA_FILE PHENOTYPE_DATA_FILE PHENOTYPE_INDEX

Option:

        -o ..., --rFile=...         Emma output in R format.
        -d ..., --delim=...         default is ", "      
        -m ..., --missingval=...    default is "NA"
        -c ..., --chr=...           default is all chromosomes
	--logTransform              Log transforms the phenotype values, before running Emma.
	--kinshipDatafile=          Datafile which is used to calculate the kinship matrix.  
	                            (default is /home/cmb-01/bvilhjal/Projects/data/2010_01format_070808.csv)
        --BoundaryStart=...         Only the region within the boundary is considered in the GWA. (Default is no boundaries)
        --minMAF=...                Remove all SNPs which have MAF smaller than the given argument.  (0.05 is set as default).
	--LRT                       Use ML and LRT instead of REML and a t-test.
        --BoundaryEnd=...           
        --phenotypeFileType=...     1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-a ..., --withArrayId=...   1 for array ID info (default), 0 if file has no array ID info.
	-h, --help	            show this help
	--parallel=...              Run Emma on the cluster with standard parameters.  The arguement is used for runid as well as output files.  
	--parallelAll               Run Emma on all phenotypes.

Examples:
	Emma.py -o emma_result.r  250K.csv phenotypes.tsv phenotype_index 
	
Description:

"""

import sys, getopt, traceback
import os, env
import phenotypeData

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["rFile=","chr=", "delim=", "missingval=", "withArrayId=", "BoundaryStart=", "logTransform", "BoundaryEnd=", "phenotypeFileType=", "help", "parallel=", "parallelAll", "LRT", "minMAF=", "kinshipDatafile="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:c:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
        phenotypeFileType = 1
        rFile = None
	delim = ","
	missingVal = "NA"
	help = 0
	minMAF=0.05
	withArrayIds = 1
        boundaries = [-1,-1]
        chr=None
	parallel = None
	logTransform = False
	parallelAll = False
	lrt = False
	kinshipDatafile = "/home/cmb-01/bvilhjal/Projects/data/2010_01format_070808.csv"  #This default path should perhaps be something else..

	for opt, arg in opts:
            if opt in ("-h", "--help"):
                help = 1
                print __doc__
            elif opt in ("-a","--withArrayId"):
                withArrayIds = int(arg)
            elif opt in ("-o","--rFile"):
                rFile = arg
            elif opt in ("--phenotypeFileType"):
                phenotypeFileType = int(arg)
            elif opt in ("--BoundaryStart"):
                boundaries[0] = int(arg)
            elif opt in ("--BoundaryEnd"):
                boundaries[1] = int(arg)
            elif opt in ("--parallel"):
                parallel = arg
            elif opt in ("--minMAF"):
                minMAF = float(arg)
            elif opt in ("--parallelAll"):
                parallelAll = True
            elif opt in ("--logTransform"):
                logTransform = True
            elif opt in ("--LRT"):
                lrt = True
            elif opt in ("-c","--chr"):
                chr = int(arg)
            elif opt in ("-d","--delim"):
                delim = arg
            elif opt in ("-m","--missingval"):
                missingVal = arg
            elif opt in ("--kinshipDatafile"):
                kinshipDatafile = arg
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
		emmadir = '/home/cmb-01/bvilhjal/Projects/Python-snps/'
		resultDir = '/home/cmb-01/bvilhjal/results/'
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
		phenName = phed.phenotypeNames[phenotypeIndex]
		phenName = phenName.replace("/","_div_")
		phenName = phenName.replace("*","_star_")
		outFileName = resultDir+"Emma_"+parallel+"_"+phenName
		rFile = outFileName 

		shstr = """#!/bin/csh
#PBS -l walltime=72:00:00
#PBS -l mem=4g 
#PBS -q cmb
"""

		shstr += "#PBS -N E"+phenName+"_"+parallel+"\n"
		shstr += "set phenotypeName="+parallel+"\n"
		shstr += "set phenotype="+str(phenotypeIndex)+"\n"
		shstr += "(python "+emmadir+"Emma.py -o "+rFile+" "
		if logTransform:
			shstr += " --logTransform "
		shstr += " -a "+str(withArrayIds)+" "			
		shstr += " --kinshipDatafile="+str(kinshipDatafile)+" "			
		shstr += snpsDataFile+" "+phenotypeDataFile+" "+str(phenotypeIndex)+" "
		shstr += "> "+outFileName+"_job"+".out) >& "+outFileName+"_job"+".err\n"

		f = open(parallel+".sh",'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	snpsDataFile = args[0]
	phenotypeDataFile = args[1]
	if parallel:  #Running on the cluster..
		if parallelAll:
			phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
			for phenotypeIndex in range(0,len(phed.phenotypeNames)):
				runParallel(phenotypeIndex)
		else:
			phenotypeIndex = int(args[2])
			runParallel(phenotypeIndex)
		return
	else:
		phenotype = int(args[2])



	print "Emma is being set up with the following parameters:"
	print "withArrayId:",withArrayIds
	print "phenotypeFileType:",phenotypeFileType
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "minMAF:",minMAF
	print "LRT:",lrt
	print "delim:",delim
	print "missingval:",missingVal
	print "kinshipDatafile:",kinshipDatafile
	print ""


	import dataParsers
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=withArrayIds)
	kinshipSnpsds = dataParsers.parseCSVData(kinshipDatafile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=withArrayIds)

	#Load phenotype file
	phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
	accIndicesToKeep = []			
	numAcc = len(snpsds[0].accessions)

	#Checking which accessions to keep and which to remove (genotype data).
	sys.stdout.write("Removing accessions which do not have a phenotype value for "+phed.phenotypeNames[phenotype]+".")
	sys.stdout.flush()
	for i in range(0,len(kinshipSnpsds[0].accessions)):
		acc1 = kinshipSnpsds[0].accessions[i]
		for j in range(0,len(phed.accessions)):
			acc2 = phed.accessions[j]
			if acc1==acc2 and phed.phenotypeValues[j][phenotype]!='NA':
				accIndicesToKeep.append(i)
				break	

	for snpsd in kinshipSnpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snpsd.removeAccessionIndices(accIndicesToKeep)
	print ""
	print numAcc-len(accIndicesToKeep),"accessions removed from kinship genotype data, leaving",len(accIndicesToKeep),"accessions in all."

	accIndicesToKeep = []			
	phenAccIndicesToKeep = []
	#Checking which accessions to keep and which to remove (kinship genotype data).
	for i in range(0,len(snpsds[0].accessions)):
		acc1 = snpsds[0].accessions[i]
		for j in range(0,len(phed.accessions)):
			acc2 = phed.accessions[j]
			if acc1==acc2 and phed.phenotypeValues[j][phenotype]!='NA':
				accIndicesToKeep.append(i)
				phenAccIndicesToKeep.append(j)
				break	


	#Filter accessions which do not have the phenotype value (from the genotype data).
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		snpsd.removeAccessionIndices(accIndicesToKeep)
	print ""
	print numAcc-len(accIndicesToKeep),"accessions removed from genotype data, leaving",len(accIndicesToKeep),"accessions in all."
		
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

	accessionMapping = []
	i = 0
	for acc in snpsds[0].accessions:
		if acc in kinshipSnpsds[0].accessions:
			accessionMapping.append((kinshipSnpsds[0].accessions.index(acc),i))
			i += 1

	print "Ordering kinship data accessions."
	for snpsd in kinshipSnpsds:
		snpsd.orderAccessions(accessionMapping)
	
	
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
		newSnpsds.append(snpsd.getSnpsData(missingVal=missingVal))
	print ""

	newKinshipSnpsds = []
	sys.stdout.write("Converting data format")
	for snpsd in kinshipSnpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		newKinshipSnpsds.append(snpsd.getSnpsData(missingVal=missingVal))  #This data might have NAs
	print ""

	#Writing files
	cluster = "/home/cmb-01/bvilhjal/"==env.homedir #Am I running on the cluster?
	import tempfile
	#tempfile.tempdir = "/home/cmb-01/bvilhjal/tmp/" #(Temporary) debug hack...
	if not cluster:
		tempfile.tempdir='/tmp'
	(fId, phenotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, genotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, kinshipTempFile) = tempfile.mkstemp()
	os.close(fId)
	
	phed.writeToFile(phenotypeTempFile, [phenotype])	
	sys.stdout.write( "Phenotype file written\n")
	sys.stdout.flush()
	snpsDataset = snpsdata.SnpsDataSet(newSnpsds,[1,2,3,4,5])
	snpsDataset.writeToFile(genotypeTempFile, deliminator=delim, missingVal = missingVal, withArrayIds = 0)
	sys.stdout.write( "Genotype file written\n")
	sys.stdout.flush()
	kinshipSnpsDataset = snpsdata.SnpsDataSet(newKinshipSnpsds,[1,2,3,4,5])
	kinshipSnpsDataset.writeToFile(kinshipTempFile, deliminator=delim, missingVal = missingVal, withArrayIds = 0)
	sys.stdout.write( "Kinship genotype file written\n")
	sys.stdout.flush()

	phenotypeName = phed.phenotypeNames[phenotype].split("_")[1]
	phenotypeName = phenotypeName.replace("/","_div_")
	phenotypeName = phenotypeName.replace("*","_star_")

	rDataFile = rFile+".rData"
	pvalFile = rFile+".pvals"
	rstr = _generateRScript_(genotypeTempFile, phenotypeTempFile, rDataFile, pvalFile, kinshipTempFile, name = phenotypeName, boundaries = boundaries, chr=chr,cluster=cluster, lrt=lrt)
	f = open(rFile,'w')
	f.write(rstr)
	f.close()
	outRfile = rFile+"_R.out"
	errRfile = rFile+"_R.err"
	print "Running R file:"
        cmdStr = "(R --vanilla < "+rFile+" > "+outRfile+") >& "+errRfile
	sys.stdout.write(cmdStr+"\n")
	sys.stdout.flush()
	os.system(cmdStr)
	print "Emma output saved in R format in", rDataFile
	
	
def _generateRScript_(genotypeFile, phenotypeFile, rDataFile, pvalFile, kinshipDatafile, name=None, boundaries = [-1,-1],chr=None, cluster=False, lrt=False):
	
	if cluster:
		rstr = 'library(emma,lib.loc="/home/cmb-01/bvilhjal/Projects/emma");\n'
	else:
		rstr = "library(emma);\n"
	rstr += 'data250K <- read.csv("'+str(genotypeFile)+'", header=TRUE);\n'
	rstr += "mat250K <- as.matrix(data250K);\n"
	rstr += 'data2010 <- read.csv("'+kinshipDatafile+'", header=TRUE);\n'
	rstr += "mat2010 <- as.matrix(data2010);\n"
	rstr += 'phenotData <- read.csv("'+str(phenotypeFile)+'",header=TRUE);\n'
	rstr += """
phenMat <- t(as.matrix(as.matrix(phenotData)[,2]));
res <- list();
mat250KAll <- mat250K[,3:length(mat250K[1,])]
mat2010All <- mat2010[,3:length(mat2010[1,])]
K2010<- emma.kinship(mat2010All,use="all") #Kinship matrix
"""
	
	if chr:
		pass #Not implemented yet
	else:
				
		rstr += """
pvals <- c();
positions <- c();
chrs <- c();
maf <- c();
for (chr in (1:5)){
  mat250K <- as.matrix(data250K);
  mat250K <- mat250K[mat250K[,1]==chr,];

  pos <- mat250K[,2];
  mat250K <- mat250K[,3:length(mat250K[1,])];
"""

		if lrt: 
			rstr += "  res[[chr]] <- emma.ML.LRT(phenMat,mat250K,K2010);"
		else:
			rstr += "  res[[chr]] <- emma.REML.t(phenMat,mat250K,K2010);"
		rstr +="""
  res[[chr]][["pos"]] <- pos;
  res[[chr]][["chr_pos"]] <- pos;

  pvals <- append(pvals,res[[chr]]$ps);
  positions <- append(positions,pos);
  chrs <- append(chrs,rep(chr,length(pos)));
  af <- c();
  for(i in (1:length(mat250K[,1]))){
    f = factor(mat250K[i,]);
    v <- summary(f)[[1]]/length(f);
    af[i] <- min(v,1-v);
  }
  res[[chr]][["maf"]] <- af;

  maf <- append(maf,af);
}
res[["K"]] <- K2010;

#write to a pvalue-file
l <- list();
l[["Chromasome"]]<-chrs;
l[["Positions"]]<-positions;
l[["Pvalues"]]<-pvals;
l[["MAF"]]<-maf;
dl <- as.data.frame(l);
"""
		rstr +=' write.table(dl,file="'+pvalFile+'", sep=", ", row.names = FALSE);\n'		
		rstr += """
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
		rstr += 'res[["lab"]]= "Emma p-values for '+name+'";\n'
	else:
		rstr += 'res[["lab"]]="";\n'
	rstr += 'save(file="'+rDataFile+'",res);\n'
	return rstr	
	
if __name__ == '__main__':
	_run_()


