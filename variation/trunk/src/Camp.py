
"""
A wrapper for the camp program (Kimmel et al., 2008).  

	-o ..., --outputFile=...		CAMP output.
	-d ..., --delim=...			default is ", "	  
	-m ..., --missingval=...		default is "NA"
	-n ..., --sampleNum=...		Num of samples
	--parallel=...				Run Camp on the cluster with standard parameters.  The arguement is used for runid 
	--parallelAll				Run Camp on all phenotypes.
	-h, --help				show this help
	--useFloats				Use floats in phenotype values.
"""


import sys, getopt, traceback
import os, env
import phenotypeData
import tempfile
tempfile.tempdir = "/home/cmb-01/bvilhjal/tmp/" #(Temporary) debug hack...
import dataParsers
import snpsdata

resultDir="/home/cmb-01/bvilhjal/results/"
scriptDir="/home/cmb-01/bvilhjal/Projects/Python-snps/"


def _run_():
	if len(sys.argv)==1:
		print __doc__
		sys.exit(2)
	
	long_options_list=["outputFile=", "delim=", "missingval=", "sampleNum=", "parallel=", "parallelAll", "useFloats"]
	try:
		opts, args=getopt.getopt(sys.argv[1:], "o:d:m:n:h", long_options_list)

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
	sampleNum = None
	chromosomes=[1,2,3,4,5]	
	useFloats = False

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help=1
			print __doc__
		elif opt in ("-o", "--outputFile"):
			outputFile=arg
		elif opt in ("--parallel"):
			parallel=arg
		elif opt in ("--parallelAll"):
			parallelAll=True
		elif opt in ("-d", "--delim"):
			delim=arg
		elif opt in ("-m", "--missingval"):
			missingVal=arg
		elif opt in ("n", "--sampleNum"):
			sampleNum = int(arg)
		elif opt in ("--useFloats"):
			useFloats = True
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

	print "CAMP is being set up with the following parameters:"
	print "phenotypeDataFile:",phenotypeDataFile
	if len(args)>2:
		print "Phenotype_id:",args[2]
	print "snpsDataFile:",snpsDataFile
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "sampleNum:",sampleNum


	def runParallel(phenotypeIndex,id=""):
		#Cluster specific parameters
		phed=phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 
		phenName=phed.getPhenotypeName(phenotypeIndex)
		phenName=phenName.replace("/", "_div_")
		phenName=phenName.replace("*", "_star_")
		outputFile=resultDir+"CAMP_"+parallel+"_"+phenName+id

		shstr="""#!/bin/csh
#PBS -l walltime=24:00:00
#PBS -l mem=6g 
#PBS -q cmb
"""
		
		shstr+="#PBS -N C"+phenName+"_"+parallel+"\n"
		shstr+="set phenotypeName="+parallel+"\n"
		shstr+="set phenotype="+str(phenotypeIndex)+"\n"
		shstr+="(python "+scriptDir+"Camp.py -o "+outputFile+" "
		if sampleNum:
			shstr+=" -n "+str(sampleNum)+" "			
		if useFloats:
			shstr+=" --useFloats "			

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
		else:
			phenotypeIndex=int(args[2])
			runParallel(phenotypeIndex)
		return
	else:
		phenotypeIndex=int(args[2])
		
	#Load phenotype file
	phed=phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter = '\t')  #Get Phenotype data 


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

	#Writing phenotype data to CAMP format.
	(fId, phenotypeFile) = tempfile.mkstemp()
	os.close(fId)
	phenVals = phed.getPhenVals(phenotypeIndex,asString=False)
	if not useFloats:
		phenVals = map(int,phenVals)
	phenFile = open(phenotypeFile,"w")
	for value in phenVals:
		phenFile.write(str(value)+"\n")
	phenFile.close()

	
	chromosome_list = [] 
	positions_list = [] 
	scores_list = [] 
	interaction_positions_list = []
	mafs = []
	marfs = [] 
	#Writing SNP data to CAMP format.
	for chromosome in chromosomes:
		(fId, snpsFile) = tempfile.mkstemp()
		os.close(fId)
		(fId, posFile) = tempfile.mkstemp()
		os.close(fId)
		sf = open(snpsFile,"w")
		pf = open(posFile,"w")
		snpsd = newSnpsds[chromosome-1]
		for i in range(0,len(snpsd.snps)):
			snp = snpsd.snps[i]
			(marf,maf) = snpsdata.getMAF(snp)
			marfs.append(marf)
			mafs.append(maf)
			str_snp = map(str,snp)
			double_snp = []
			for nt in str_snp:
				double_snp.append(nt)
				double_snp.append(nt)
			sf.write("".join(double_snp)+"\n")
			pf.write(str(snpsd.positions[i])+"\n")
		sf.close()
		pf.close()
		
		outFile = outputFile+"_job_"+str(chromosome)+".out"
		errFile = outputFile+"_job_"+str(chromosome)+".err"
		resFile = outputFile+"_"+str(chromosome)+".out"
		print "resFile,outFile,errFile,snpsFile,posFile,phenotypeFile:",resFile,outFile,errFile,snpsFile,posFile,phenotypeFile
		results = _runCAMP_(resFile,outFile,errFile,snpsFile,posFile,phenotypeFile,sampleNum)
		
		positions_list += results["positions"]
		scores_list += results["scores"]
		for (i,j) in results["snpIndices"]:
			if not (j<0 or i<0):
				marfs.append(0.5)  #An ugly hack!!!
				mafs.append(0.5)
			chromosome_list.append(chromosome)
	
	scoreFile = outputFile+".scores"
	f = open(scoreFile,"w")
	f.write("Chromosome,Position,Score,MARF,MAF,Second_Position\n")
	for i in range(0,len(positions_list)):
		chromosome = chromosome_list[i]
		(pos1,pos2) = positions_list[i]
		score = scores_list[i]
		marf = marfs[i]
		maf = mafs[i]
		l = map(str,[chromosome,pos1,score,marf,maf,pos2])
		f.write(",".join(l)+"\n")
	f.close()


def _runCAMP_(resFile,outFile,errFile,snpsFile,posFile,phenotypeFile,sampleNum=None,windowSize=None):
	cmdStr = "(/home/cmb-01/bvilhjal/Projects/camp/camp -i "+snpsFile+" -f -p "+posFile+" -d "+phenotypeFile+" -w 30000  -o "+resFile+" "
	if windowSize:
		cmdStr+= " -w "+str(windowSize)+" "
	if sampleNum:
		cmdStr+= " -n "+str(sampleNum)+" "
	cmdStr+=" > "+outFile+") >& "+errFile
	sys.stdout.write(cmdStr+"\n")
	sys.stdout.flush()
	os.system(cmdStr)

	#Parse results... and return
	rf = open(resFile,"r")
	#pf = open(posFile,"r")
	lines = rf.readlines()
	lines.pop(0)
	positions = []
	interactions = []
	snpIndices = []
	scores = []
	i = 1
	for line in lines:
		line_list = line.split()
		if i%10000==0:
			print i,"lines read"
			sys.stdout.flush()
		
		if len(line_list)==5:
			positions.append(int(line_list[2]))
			interaction = int(line_list[3])
			interactions.append(interaction)
			snpIndices.append([int(line_list[0])-1,int(line_list[1])-1])
			scores.append(float(line_list[4]))
		else:
			print line_list, len(line_list)
			break
	print i,"lines read out of",len(lines)
	sys.stdout.flush()
	rf.close()
	
	return {"scores":scores, "positions": positions, "interactions":interactions, "snpIndices":snpIndices}	
		
if __name__ == '__main__':
	_run_()
	print "Done!"

