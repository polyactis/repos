#!/usr/bin/env python2.5
"""
Usage: KW.py [OPTIONS] [-o OUT_FILE] SNPS_DATA_FILE PHENOTYPE_DATA_FILE [PHENOTYPE_INDEX]

Option:

        -o ..., --outputFile=...         Output files, one 'name'.rData file, and one 'name'.pvals file.
        -d ..., --delim=...         default is ", "      
        -m ..., --missingval=...    default is "NA"
        --phenotypeFileType=...     1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-a ..., --withArrayId=...   1 for array ID info (default), 0 if file has no array ID info.
	-h, --help	            show this help
	--parallel=...              Run emma on the cluster with standard parameters.  The arguement is used for runid 
                                    as well as output files.  Note that if this option is used then no output file should be specified.
	--parallelAll               Run Emma on all phenotypes.
	--addToDB                   Adds the result to the DB.
	--callMethodID=...          What data set is used. (Only applicable if result is added to DB.)
	--comment=...               Comment for DB. (Only applicable if result is added to DB.)

Examples:
	KW.py -o outputFile  250K.csv phenotypes.tsv phenotype_index 

Cluster examples:
	KW.py -o --parallel=runIdentifier --parallelAll /global_path_to_data/250K.csv /global_path_to_data/phenotypes.tsv
	

Description:
  Applies the Kruskal Wallis test to phenotypes, which are not binary.
  Applies a Chi-square test to the phenotypes which are binary!

"""

import sys, getopt, traceback
import os, env
import phenotypeData


resultDir=env.resultDir
scriptDir=env.scriptDir


def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["outputFile=", "delim=", "missingval=", "withArrayId=", "phenotypeFileType=", "help", "parallel=", "parallelAll", "addToDB","callMethodID=", "comment="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:c:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
        phenotypeFileType = 1
        outputFile = None
	delim = ","
	missingVal = "NA"
	help = 0
	withArrayIds = 1
	parallel = None
	parallelAll = False
	addToDB=False
	callMethodID=None
	comment=""

	for opt, arg in opts:
            if opt in ("-h", "--help"):
                help = 1
                print __doc__
            elif opt in ("-a","--withArrayId"):
                withArrayIds = int(arg)
            elif opt in ("-o","--outputFile"):
                outputFile = arg
            elif opt in ("--phenotypeFileType"):
                phenotypeFileType = int(arg)
            elif opt in ("--parallel"):
                parallel = arg
            elif opt in ("--parallelAll"):
                parallelAll = True
            elif opt in ("--addToDB"):
                addToDB = True
            elif opt in ("--callMethodID"):
                callMethodID = int(arg)
            elif opt in ("--comment"):
                comment = arg
            elif opt in ("-d","--delim"):
                delim = arg
            elif opt in ("-m","--missingval"):
                missingVal = arg
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
		phenName = phed.phenotypeNames[phenotypeIndex]
		phenName = phenName.replace("/","_div_")
		phenName = phenName.replace("*","_star_")
		outputFile = resultDir+"KW_"+parallel+"_"+phenName

		shstr = """#!/bin/csh
#PBS -l walltime=72:00:00
#PBS -l mem=4g 
#PBS -q cmb
"""

		shstr += "#PBS -N K"+phenName+"_"+parallel+"\n"
		shstr += "set phenotypeName="+parallel+"\n"
		shstr += "set phenotype="+str(phenotypeIndex)+"\n"
		shstr += "(python "+scriptDir+"KW.py -o "+outputFile+" "
		shstr += " -a "+str(withArrayIds)+" "			
		shstr += snpsDataFile+" "+phenotypeDataFile+" "+str(phenotypeIndex)+" "
		shstr += "> "+outputFile+"_job"+".out) >& "+outputFile+"_job"+".err\n"

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





	import dataParsers
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=delim, missingVal=missingVal, withArrayIds=withArrayIds)
	
	phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
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

        #Filtering monomorphic
	print "Filtering monomorphic SNPs"
	for snpsd in snpsds:
		print "Removed", str(snpsd.filterMonoMorphicSnps()),"Snps"

	#Converting format to 01
	import snpsdata
	newSnpsds = []
	sys.stdout.write("Converting data format")
	for snpsd in snpsds:
		sys.stdout.write(".")
		sys.stdout.flush()
		newSnpsds.append(snpsd.getSnpsData())
	print ""

	#Writing files
	import tempfile
	if env.user == "bjarni":
		tempfile.tempdir='/tmp'
	(fId, phenotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	(fId, genotypeTempFile) = tempfile.mkstemp()
	os.close(fId)
	
	phed.writeToFile(phenotypeTempFile, [phenotype])	
	sys.stdout.write( "Phenotype file written\n")
	sys.stdout.flush()
	snpsDataset = snpsdata.SnpsDataSet(newSnpsds,[1,2,3,4,5])
	decoder = {1:1, 0:0, -1:'NA'}	
	snpsDataset.writeToFile(genotypeTempFile, deliminator=delim, missingVal = missingVal, withArrayIds = 0, decoder=decoder)
	sys.stdout.write( "Genotype file written\n")
	sys.stdout.flush()

	phenotypeName = phed.phenotypeNames[phenotype].split("_")[1]
	phenotypeName = phenotypeName.replace("/","_div_")
	phenotypeName = phenotypeName.replace("*","_star_")

	rDataFile = outputFile+".rData"
	pvalFile = outputFile+".pvals"
	binary = phed.isBinary(phenotype)
	rstr = _generateRScript_(genotypeTempFile, phenotypeTempFile, rDataFile, pvalFile, name = phenotypeName, binary=binary)
	rFileName = outputFile+".r"
	f = open(rFileName,'w')
	f.write(rstr)
	f.close()
	outRfile = rFileName+".out"
	errRfile = rFileName+".err"
	print "Running R file:"
        cmdStr = "(R --vanilla < "+rFileName+" > "+outRfile+") >& "+errRfile
	sys.stdout.write(cmdStr+"\n")
	sys.stdout.flush()
	os.system(cmdStr)
	print "Emma output saved in R format in", rDataFile
	
	
def _generateRScript_(genotypeFile, phenotypeFile, rDataFile, pvalFile, name=None, binary=False, chiSquare=False):
	
	rstr = 'data250K <- read.csv("'+str(genotypeFile)+'", header=TRUE);\n'
	rstr += 'phenotData <- read.csv("'+str(phenotypeFile)+'",header=TRUE);\n'
	rstr += """
phenMat <- t(as.matrix(as.matrix(phenotData)[,2]));
res <- list();
pvals <- c();
positions <- c();
chrs <- c();
maf <- c();
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
			rstr += "    v <- chisq.test(as.vector(phenMat),as.vector(mat250K[j,]));"
			rstr += """
    res[[chr]]$ps[j] <- as.double(v$p.value);
    res[[chr]]$stat[j] <- as.double(v[1]);		
"""
		else:
			rstr += "    v <- fisher.test(as.vector(phenMat),as.vector(mat250K[j,]));"			
			rstr += """
    res[[chr]]$ps[j] <- as.double(v$p.value);
    res[[chr]]$stat[j] <- as.double(v$p.value);		
"""

	else:
		rstr += "    v <- kruskal.test(as.vector(phenMat),as.vector(mat250K[j,]));"
		rstr += """
    res[[chr]]$ps[j] <- as.double(v$p.value);
    res[[chr]]$stat[j] <- as.double(v[1]);		
"""
	rstr +="""
  }  
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

#write to a pvalue-file
l <- list();
l[["Chromasome"]]<-chrs;
l[["Positions"]]<-positions;
l[["Pvalues"]]<-pvals;
l[["MAF"]]<-maf;
dl <- as.data.frame(l)
"""
	rstr +='write.table(dl,file="'+pvalFile+'", sep=", ", row.names = FALSE);\n'		
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
		if binary:
			if chiSquare:
				rstr += 'res[["lab"]]= "Chi-Square p-values for '+name+'";\n'
			else:
				rstr += 'res[["lab"]]= "Fisher p-values for '+name+'";\n'
		else:
			rstr += 'res[["lab"]]= "Kruskal Wallis p-values for '+name+'";\n'
	else:
		rstr += 'res[["lab"]]="";\n'
	rstr += 'save(file="'+rDataFile+'",res);\n'
	return rstr	
	
if __name__ == '__main__':
	_run_()


