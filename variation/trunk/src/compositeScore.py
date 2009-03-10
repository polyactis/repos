
"""
Usage: compositeScore.py [OPTIONS] [-o R_FILE] [PHENOTYPE_INDEX_1 PHENOTYPE_INDEX_2 PHENOTYPE_INDEX_3 .. ]

Option:

	-o ...					  compositeScore output.
	-d ..., --delim=...		 default is ", "	  
	-m ..., --missingval=...	default is "NA"
	--candGeneListID=...		The database ID of the cand. gene list. used.
	--phenotypeCategory=...	 The database biology category of the phenotypes to be tested. (Only with --parallelAll, FT (1) is default.) 
	--testDataFraction=...	  The fraction of data reserved for testing. (1/3 is default)
	--gridSize=...			  The size of the grid on which the parameters are tested. (6 is default)
	--windowSize=...			The radius of the window around cand. genes to look for SNPs. (default is 10000)
	--phenotypeFileType=...	 1 (default) if file has tsv format, 2 if file has csv format and contains accession names (instead of ecotype ID)
	-h, --help				show this help
	--parallel=...			  Run compositeScore on the cluster with standard parameters.  The arguement is used for runid as well as output files.  
	--parallelAll			   Run compositeScore on all phenotypes.

Examples:
	compositeScore.py -parallel=newDataset --parallelAll
	
Description:

Warning:
Remember to specify where to search for the result files.  And how to treat each of the results in the code!

"""

import sys, getopt, traceback
import os, env, gwaResults
import MySQLdb
import phenotypeData

#The following parameters should be set accordingly.
res_path="/home/cmb-01/bvilhjal/results/"
resultsDirs = [res_path+"kw_results/",res_path+"emma_results/",res_path+"marg_results/",res_path+"rf_results/"]
methods=["KW","Emma","Marg","RF"]
fileTypes=[".pvals",".pvals",".score",".imp"]
datasetNames=["newDataset","newDataset","newDataset","newDataset"]
logTransform = [True,True,False,False]
mafCutoffs = [0,15,0,15]
onlyQuantitative = [0,1,0,0]
perlScriptDir = "/home/cmb-01/bvilhjal/Projects/Python-snps/composite_rank3.pl"
phenotypeDataFile="/home/cmb-01/bvilhjal/Projects/data/phenotypes_transformed_publishable_v2.tsv" 


def _run_():
	import sys
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["delim=", "missingval=", "phenotypeFileType=", "help", "parallel=", "parallelAll","candGeneListID=","windowSize=","testDataFraction=","gridSize=","phenotypeCategory="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:c:d:m:a:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	phenotypeCategory = 1
	phenotypeFileType = 1
	testDataFraction = 1.0/3.0
	gridSize = 6
	outFile = None
	delim = ","
	missingVal = "NA"
	help = 0
	parallel = None
	parallelAll = False
	candGeneListID = 129
	windowSize=10000

	host = "papaya.usc.edu"
	user = "bvilhjal"
	passwd = "bamboo123"
	db = "T8_annotation_TH"
	

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-o"):
			outFile = arg
		elif opt in ("--gridSize"):
			gridSize = int(arg)
		elif opt in ("--phenotypeFileType"):
			phenotypeFileType = int(arg)
		elif opt in ("--phenotypeCategory"):
			phenotypeCategory = int(arg)
		elif opt in ("--testDataFraction"):
			testDataFraction = float(arg)
		elif opt in ("--candGeneListID"):
			candGeneListID = int(arg)
		elif opt in ("--windowSize"):
			windowSize = int(arg)
		elif opt in ("--parallel"):
			parallel = arg
		elif opt in ("--parallelAll"):
			parallelAll = True
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if len(args)<1 and not parallel:
		if help==0:
			print "Arguments are missing!!\n"
			print __doc__
		sys.exit(2)

	def runParallel(phenotypeIndex):
		#Cluster specific parameters
		scriptDir = env.scriptDir
		resultDir = env.resultDir
		phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
		phed.onlyBiologyCategory(phenotypeCategory,host=host,user=user,passwd=passwd)
		phenName = phed.getPhenotypeName(phenotypeIndex)
		phenName = phenName.replace("/","_div_")
		phenName = phenName.replace("*","_star_")
		outFile = resultDir+"CS_"+parallel+"_"+phenName
		shstr = """#!/bin/csh
#PBS -l walltime=72:00:00
#PBS -l mem=4g 
#PBS -q cmb
"""
		shstr += "#PBS -N CS"+phenName+"_"+parallel+"\n"
		shstr += "(python "+scriptDir+"compositeScore.py -o"+outFile+" "
		shstr += "--candGeneListID="+str(candGeneListID)+" --testDataFraction="+str(testDataFraction)+" --gridSize="+str(gridSize)+" --windowSize="+str(windowSize)+" --phenotypeCategory="+str(phenotypeCategory)+" "+str(phenotypeIndex)+" "

		shstr += "> "+outFile+"_job"+".out) >& "+outFile+"_job"+".err\n"

		f = open(parallel+".sh",'w')
		f.write(shstr)
		f.close()

		#Execute qsub script
		os.system("qsub "+parallel+".sh ")

	if parallel:  #Running on the cluster..
		if parallelAll:
			phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
			phed.onlyBiologyCategory(phenotypeCategory,host=host,user=user,passwd=passwd)
			for phenotypeIndex in phed.phenIds:
				runParallel(phenotypeIndex)
		else:
			for arg in args:
				runParallel(int(arg))
		return
	else:
		
		phenotype = int(args[0])
		if len(args)>1:
			print "Warning multiple phenotype_id arguments were ignored (use --parallel)."


	print "compositeScore is being set up with the following parameters:"
	print "candGeneListID:",candGeneListID
	print "phenotypeCategory:",phenotypeCategory
	print "phenotype:",phenotype
	print "gridSize:",gridSize
	print "testDataFraction:",testDataFraction
	print "phenotypeFileType:",phenotypeFileType
	print "parallel:",parallel
	print "parallelAll:",parallelAll
	print "delim:",delim
	print "missingval:",missingVal
	print ""


	#Now the algorithm!!!

	#Load phenotype file
	categoricalNames = ["158_Sil_length_16","159_Sil_length_22","161_Germ_10","163_Germ_22","173_Leaf_serr_10","174_Leaf_serr_16","175_Leaf_serr_22","179_Roset_erect_22","180_Chlor_16","181_Chlor_22"]
	phed = phenotypeData.readPhenotypeFile(phenotypeDataFile, delimiter='\t')  #Get Phenotype data 
	phed.onlyBiologyCategory(phenotypeCategory,host=host,user=user,passwd=passwd)
	#Check whether phenotype is quantitative..
	isQuantitative =  not (phed.isBinary(phenotype) or phed.getPhenotypeName(phenotype) in categoricalNames)
	if isQuantitative:
		print "Phenotype",phed.getPhenotypeName(phenotype),"is quantitaive."

	phenName = phed.getPhenotypeName(phenotype)
	phenName = phenName.replace("/","_div_")
	phenName = phenName.replace("*","_star_")

	#Load the result files:
	results = []
	resultFiles = []
	mafCutoff = max(mafCutoffs)
	for j in range(0,len(methods)):
		if isQuantitative or not onlyQuantitative[j]:
			resultFile=resultsDirs[j]+methods[j]+"_"+datasetNames[j]+"_"+phenName+fileTypes[j]
			resultFiles.append(resultFile)
			print "Loading result file",resultFile
			result = gwaResults.Result(resultFile,)
			if logTransform[j]:
				print "Log transformed the p-values"
				result.negLogTransform()
			
			result.filterMAF(minMaf=mafCutoff)
			results.append(result)
   
	#Write the results to a file.
	import tempfile
	#if os.getenv("USER")=="bjarni"
	#	tempfile.tempdir='/tmp'
	#tempfile.tempdir='/home/cmb-01/bvilhjal/tmp'
	(fId, resultsTempFile) = tempfile.mkstemp()
	os.close(fId)
	
	f = open(resultsTempFile,'w')
	for i in range(0,len(results[0].scores)):
		out_str = str(results[0].chromosomes[i])+"_"+str(results[0].positions[i])
		for result in results:
			out_str += ","+str(result.scores[i])		
		out_str += "\n"
		f.write(out_str)
	f.close()
		
	#Load cand. gene list.

	print "Connecting to db, host="+host
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = db)
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"  

	#select c.locustag, b.start, b.stop, a.comment from genome.gene_commentary a, genome.entrezgene_mapping b, genome.gene c where b.start > 25000 and b.stop < 75000 and b.chromosome=1 and b.gene_id = c.gene_id and c.gene_id = a.gene_id and a.gene_commentary_type_id = 8
	
	#select distinct t8_fd.tair_id, t8.chromosome, t8.start, t8.end, t8_fd.type, t8_fd.short_description from T8_annotation_TH.t8_063008 t8, T8_annotation_TH.t8_func_desc t8_fd, stock_250k.candidate_gene_list cgl where t8.pub_locus+'.1' = t8_fd.tair_id and cgl.list_type_id=129  and cgl.original_name=t8.pub_locus and t8.chromosome =1 order by t8.chromosome, t8.start

	#select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop

	
	numRows = int(cursor.execute("select distinct gm.chromosome, gm.start, gm.stop, g.locustag from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id=129 and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop"))
	candSnpRegions = []
	currTairID=""
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		candSnpRegion = (int(row[0]),int(row[1])-windowSize,int(row[2])+windowSize,row[3])
		candSnpRegions.append(candSnpRegion)
	cursor.close ()
	conn.close ()
	print "Candiate genelists fetched"

	chrAndSnpPosList = []
	currChr = 0
	currPos = 0
	i = -1
	for candSnpRegion in candSnpRegions:
		while currChr<=candSnpRegion[0] and currPos <candSnpRegion[1]:
			i += 1
			currChr = results[0].chromosomes[i]
			currPos = results[0].positions[i]

		while currChr==candSnpRegion[0] and currPos >candSnpRegion[1] and currPos <=candSnpRegion[2]:
			chrAndSnpPosList.append(str(currChr)+"_"+str(currPos))
			i +=1
			currChr = results[0].chromosomes[i]
			currPos = results[0].positions[i]
	
	#Write to file
	(fId, candidateTempFile) = tempfile.mkstemp()
	os.close(fId)
	f = open(candidateTempFile,'w')
	for chrAndSnpPos in chrAndSnpPosList:
		f.write(chrAndSnpPos+"\n")
	f.close()

	print "Candidate SNPs written out."
		
	#Execute perl..
	jobOutFile = outFile+"_script.out"
	jobErrFile = outFile+"_script.err"
	outFile = outFile+".score"
	cmdStr = "(perl "+perlScriptDir+" -o "+outFile+" -p "+resultsTempFile+" -c "+candidateTempFile+" -f "+str(testDataFraction)+" -i "+str(gridSize)+" > "+jobOutFile+") >& "+jobErrFile
	print "Executing perl :",cmdStr
	sys.stdout.flush()
	os.system(cmdStr)


	#Parse results?

if __name__ == '__main__':
	_run_()
