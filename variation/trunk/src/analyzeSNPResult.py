#!/usr/bin/env python2.5
"""
Usage: analyzeSNPResult.py [OPTIONS] chromosome position

Option:

		-o ...,						Output file
	--callMethodID=...		  Use the call_method_id in the DB. 
		--phenotypeMethodID=...	 The phenotype_method_id from the DB.
		--analysisMethodID=...	  The analysis_method_id from the DB.
	--snpsDataFile=...		  SNPs data file which is used for LD calculations.
		--boxPlotFile=...		   Draw a box plot in the given file.
		--calcEffect				Calculate the variance explained by this SNP and other statistics.
	--logTransform			  Negative log transform scores (necessary if using p-values).
	--minMAF=...				minimum allowed minor allele frequency.
	-h, --help				show this help

Examples:

	
Description:


For each region:

1. Plot p-values, scores (on a relative scale).
2. Plot genes.
3. Plot LD information.

4. Plot Haplotype information.
5. Calc. statistics:
	 - P-value (using mixed models) (w. and w.out the kinship matrix)
	 - Distribution of residuals, plot
	 - Variance explained (w. and w.out the kinship matrix)

(DONE) 6. Do a box-plot of most sign. SNP in region.

"""


import sys, getopt, traceback, util, pdb, gc
import dataParsers
import phenotypeData
	
import plotResults, gwaResults

#Database variables
_user_="bvilhjal"
_passwd_="bamboo123"
_host_="papaya.usc.edu"

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["callMethodID=", "phenotypeMethodID=","analysisMethodID=", "snpsDataFile=", "boxPlotFile=", "calcEffect", "minMAF=", "help", "logTransform"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:i:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	output_fname = None
	callMethodID = None
	phenotypeMethodID = None
	analysisMethodID = None
	boxPlotFile = None
	snpsDataFile = None
	help=0
	logTransform = False
	calcEffect = False
	minMAF = 0

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("--callMethodID"):
			callMethodID = int(arg)
		elif opt in ("--phenotypeMethodID"):
			phenotypeMethodID = int(arg)
		elif opt in ("--analysisMethodID"):
			analysisMethodID = int(arg)
		elif opt in ("--minMAF"):
			minMAF = float(arg)
		elif opt in ("--logTransform"):
			logTransform = True
		elif opt in ("--calcEffect"):
			calcEffect = True
		elif opt in ("--snpsDataFile"):
			snpsDataFile = arg
		elif opt in ("--boxPlotFile"):
			boxPlotFile = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not (output_fname and phenotypeMethodID and analysisMethodID and callMethodID) :
		if help==0:
			print "\nArguments missing!!\n"
			print __doc__
		sys.exit(2)

	print "Retrieving results"

	# Retrieve the correct results file.
	host = "papaya.usc.edu"
	user = "bvilhjal"
	passwd = "bamboo123"
	resultFiles = gwaResults.getResultsFilename(host,user,passwd,callMethodID,phenotypeMethodID,analysisMethodID)
	print "Results found in",resultFiles[0]
	# Load the results file.
	print "Loading results"
	if snpsDataFile: #Load the sequence data.
		result = gwaResults.SNPResult(resultFiles[0],snpsDataFile)
	else:
		result = gwaResults.SNPResult(resultFiles[0])
	
	if minMAF:
		result.filterMAF(minMAF)
		
	if logTransform:
		print "Log transforming the results"
		result.negLogTransform()



def drawBoxPlot(snp,pData,accIndices,phenotypeName="SNP",runId="",filename="/tmp/test.pdf"):
	"""
	Draws a box-plot
	
	requires pylab to be installed.
	"""
	phenoG1=[]
	phenoG2=[]
	allele1 = snp.alleles[accIndices[0][1]]
	for (i1,i2) in accIndices:
		if snp.alleles[i2] == allele1:
			phenoG1.append(float(pData[i1]))
		else:
			phenoG2.append(float(pData[i1]))				
	import pylab
	pylab.figure()
	pylab.boxplot([phenoG1,phenoG2])
	minVal = min((min(phenoG1),min(phenoG2)))
	maxVal = max((max(phenoG1),max(phenoG2)))
	rangeVal = maxVal-minVal
	pylab.axis([0.2,2.8,minVal-rangeVal*0.3,maxVal+rangeVal*0.3])
	pylab.text(0.4, minVal-0.15*rangeVal, "# of obs.: ", color='k')
	pylab.text(0.95, minVal-0.15*rangeVal, str(len(phenoG1)), color='k')
	pylab.text(1.95, minVal-0.15*rangeVal, str(len(phenoG2)), color='k')
	pylab.text(0.9, maxVal+0.15*rangeVal, "-log(p-value)/score: "+str(snp.score), color='k')
	pylab.title(phenotypeName+": chromosome "+str(snp.chromosome)+", position "+str(snp.position))
	pylab.savefig(filename,format="pdf")
	
#import yh_matplotlib_artists
#import pylab
import gwaResults

#
#def drawGenes(genes, gene_position_cycle=6):
#	"""
#	2008-09-27 More or less borrowed from Yu Huang..
#	"""
#	#print "\t Drawing gene model  ..."
#	no_of_genes_drawn = 0
#	for gene in genes:
#		y_value = no_of_genes_drawn%gene_position_cycle	#cycling through the y position to avoid clogging
#		plot_one_gene(gene, y_value=0.5+y_value)
#		no_of_genes_drawn += 1
#	#print "Done drawing genes..."
#
#
#
#def plot_one_gene(gene , y_value=1,buffer=0.1):  #ax: pylab.axis obj.
#	"""
#	2008-09-29: Code largely borrowed from Yu Huang..		  
#	"""
#	y_value = 0.5+y_value/1.5
#	pylab.text(gene.startPos, -y_value+buffer, gene.tairID, size=8)
#	pylab.plot([gene.startPos,gene.endPos],[-y_value,-y_value],"k",linewidth=3)
#
#def _testGenePlot_():
#	genes = []
#	gene = gwaResults.Gene()
#	gene.startPos=1000
#	gene.endPos=5000
#	gene.tairID="At234213412"
#	genes.append(gene)
#	gene = gwaResults.Gene()
#	gene.startPos=2000
#	gene.endPos=6000
#	gene.tairID="At2d13412"
#	genes.append(gene)
#	gene = gwaResults.Gene()
#	gene.startPos=5000
#	gene.endPos=7000
#	gene.tairID="At2d13412"
#	genes.append(gene)
#	f1 = pylab.figure(figsize=(20,6))
#	pylab.plot([0,100,1000,10000],[0,5,7,9])
#	drawGenes(genes)
#	pylab.axis([0,10000,-5,10])
#	pylab.show()
#
#
#def drawSNPPlot(signSnp,region,results,phenotypeName="SNP",runId="",genes=None,filename="/tmp/test.pdf",displayRanks=True):
#	import pylab
#	"""
#	Draws a snp-plot
#	
#	requires pylab to be installed.
#	"""
#	(snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,regionRank) = region
#
#	print region
#
#	maxVal = 12
#	minVal = 0
#	rangeVal = maxVal-minVal
#	pvals = [True,True,False,False,False]
#
#	positionsList = []
#	scoreList = []
#	snps = []
#	
#	result = results[0]
#	positions = []
#	scores = []
#	i = 0
#	currPos = result.positions[0]
#	currChr = result.chromosomes[0]
#
#	while currChr<chr: 
#		i += 1
#		currChr = result.chromosomes[i]
#	#Found chromsome..
#
#	while currChr==chr and currPos<startPos:
#		i += 1
#		currPos = result.positions[i]
#		currChr = result.chromosomes[i]
#	#Found start..
#	
#	while currChr==chr and currPos<endPos:
#		currPos = result.positions[i]
#		currChr = result.chromosomes[i]
#		snps.append(result.snps[i])
#		positions.append(currPos)			
#		if result.scores[i]>maxVal and pvals[0]:
#			scores.append(maxVal)
#		else:
#			scores.append(result.scores[i])			
#		i += 1
#
#	positionsList.append(positions)
#	scoreList.append(scores)
#
#
#	for j in range(1,len(results)):
#		result = results[j]
#		positions = []
#		scores = []
#		i = 0
#		currPos = result.positions[0]
#		currChr = result.chromosomes[0]
#		while currChr<chr: 
#			i += 1
#			currChr = result.chromosomes[i]
#		#Found chromsome..
#
#		while currChr==chr and currPos<startPos:
#			i += 1
#			currPos = result.positions[i]
#			currChr = result.chromosomes[i]
#		#Found start..
#	
#		while currChr==chr and currPos<endPos:
#			currPos = result.positions[i]
#			currChr = result.chromosomes[i]
#			positions.append(currPos)			
#			if result.scores[i]>10 and pvals[j]:
#				scores.append(10.0)
#			else:
#				scores.append(result.scores[i])			
#			i += 1
#
#		positionsList.append(positions)
#		scoreList.append(scores)
#		#Found the end
#		
#	startPos = positionsList[0][0]
#	endPos = positionsList[0][len(positionsList[0])-1]
#	for i in range(1,len(positionsList)):
#		positions = positionsList[i]
#		if positions[0]<startPos:
#			startPos = positions[0]
#		if positions[len(positions)-1]>endPos:
#			endPos = positions[len(positions)-1]
#	posRange = endPos-startPos
#
#	r2Values = []
#	for i in range(0,len(snps)):
#		snp = snps[i]
#		pos = positionsList[0][i]
#		r2Values.append(r2_ld(snp,signSnp.alleles)*rangeVal+minVal)
#	
#	pylab.figure(1,figsize=(20,4))
#	pylab.plot(positionsList[0],r2Values,"k-")
#
#	margMin = min(results[2].scores)
#	margMax = max(results[2].scores)
#	margRange = margMax-margMin
#
#	scores = []
#	for score in scoreList[2]:
#		scores.append(((score-margMin)/margRange)*rangeVal+minVal)
#	pylab.plot(positionsList[2],scores,"g.")
#
#	rfMin = min(results[3].scores)
#	rfMax = max(results[3].scores)
#	rfRange = rfMax-rfMin
#
#	scores = []
#	for score in scoreList[3]:
#		scores.append(((score-rfMin)/rfRange)*rangeVal+minVal)
#	pylab.plot(positionsList[3],scores,"c.")
#
#	csMin = min(results[4].scores)
#	csMax = max(results[4].scores)
#	csRange = csMax-csMin
#
#	scores = []
#	for score in scoreList[4]:
#		scores.append(((score-csMin)/csRange)*rangeVal+minVal)
#	pylab.plot(positionsList[4],scores,"y.")
#
#	
#	pylab.plot(positionsList[1],scoreList[1],"b.")
#	pylab.plot(positionsList[0],scoreList[0],"r.")
#
#	drawGenes(genes, gene_position_cycle=6)
#
#	pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05-4,maxVal+rangeVal*0.05])
#	if displayRanks:
#		pylab.title(phenotypeName+": chromosome "+str(signSnp.chromosome)+", position "+str(signSnp.position)+", snp rank"+str(snpRank)+", region rank"+str(regionRank)+".")
#	else:
#		pylab.title(phenotypeName+": chromosome "+str(signSnp.chromosome)+", position "+str(signSnp.position)+".")
#		
#	pylab.subplots_adjust(right=0.98)
#	pylab.subplots_adjust(left=0.03)
#	pylab.subplots_adjust(bottom=0.15)
#	pylab.subplots_adjust(top=0.9)
#	pylab.savefig(filename,format="pdf")
#
#	pylab.clf()
#
#def _to01Format_(snp):
#	all1 = snp[0]
#	tSnp = [0]*len(snp)
#	for i in range(1,len(snp)):
#		allele = snp[i]
#		if allele != all1:
#			tSnp[i]=1
#	return tSnp
#			
#		
#def r2_ld(snp1,snp2):
#	tSnp1 = _to01Format_(snp1)
#	tSnp2 = _to01Format_(snp2)
#	#print tSnp1,tSnp2
#	#Calculating the frequencies
#	delta = 1.0/float(len(snp1))
#	freqs = [0.0]*4
#	for i in xrange(0,len(snp1)):
#		val = tSnp1[i]*2+tSnp2[i]
#		freqs[val] += delta
#	
#	
#	f1 = freqs[1]+freqs[3]
#	f2 = freqs[2]+freqs[3]
#	D = freqs[3]-f1*f2
#	divisor = f1*f2*(1-f1)*(1-f2)
#	if divisor != 0:
#		r2= D*D/divisor
#	else:
#		r2 = -1
#
#	#print freqs, sum(freqs), r2
#	return r2 
#
##def getTairAnn(startPos,endPos,chr):
#	host = "papaya.usc.edu"
#	user = "bvilhjal"
#	passwd = "bamboo123"
#	db = "T8_annotation_TH"
#
#	import MySQLdb
#	print "Connecting to db, host="+host
#	if not user:
#		import sys
#		sys.stdout.write("Username: ")
#		user = sys.stdin.readline().rstrip()
#	if not passwd:
#		import getpass
#		passwd = getpass.getpass()
#	try:
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = db)
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Fetching data"  
#	numRows = int(cursor.execute("select distinct pub_locus, start, end from t8_063008 where start > "+str(startPos)+" and end < "+str(endPos)+" and chromosome="+str(chr)+" and segment_type='gene' order by start"))
#	
#	genes = []
#	currTairID=""
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		gene = gwaResults.Gene()
#		gene.startPos = int(row[1])
#		gene.endPos = int(row[2])
#		gene.tairID=row[0]
#		gene.chromosome=chr
#		genes.append(gene)
#
#	numRows = int(cursor.execute("select distinct t8.tair_id, t8_fd.short_description, t8_fd.description from t8_063008 t8, t8_func_desc t8_fd where t8.pub_locus=t8_fd.tair_id and t8.start > "+str(startPos)+" and t8.end < "+str(endPos)+" and t8.chromosome="+str(chr)+" order by t8.tair_id"))
#
#	functionDescriptions = []
#	while(1):
#		row = cursor.fetchone()
#		if not row:
#			break;
#		functionDescriptions.append(row)   
#	cursor.close ()
#	conn.close ()
#
#	for gene in genes:
#		for fdesc in functionDescriptions:
#			if gene.tairID==fdesc[0]:
#				gene.shortDescriptions.append(fdesc[1])
#				gene.functionDescriptions.append(fdesc[2])
#			
#
#	return genes
#


	
	
def getQuantitativeIndices():
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	categoricalNames = ["158_Sil_length_16","159_Sil_length_22","161_Germ_10","163_Germ_22","173_Leaf_serr_10","174_Leaf_serr_16","175_Leaf_serr_22","179_Roset_erect_22","180_Chlor_16","181_Chlor_22"]


	print "Reading phenotype file:",phenotypeFile
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
	print "Found",len(phed.phenotypeNames),"phenotypes."
	quantitativeIndices = []
	for i in phed.phenIds:
		phenName = phed.getPhenotypeName(i)
		if not phed.isBinary(i) and not phenName in categoricalNames :
			print phenName
			quantitativeIndices.append(i)
	print len(quantitativeIndices)
	return quantitativeIndices


#
#def loadAllResults(phenotypeIndices=[1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,80,81,82]):
#	"""
#	Load all the results needed.
#	"""
#
#	#ftIndices = [1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,80,81,82]
#	#ionomicsIndices = range(14,32)+range(83,101)
#
#	res_path="/Network/Data/250k/tmp-bvilhjal/"
#
#	resultsDirs = [res_path+"kw_results/",res_path+"emma_results/",res_path+"marg_results/",res_path+"rf_results/",res_path+"cs_results/"]
#	methods=["KW","Emma","Marg","RF","CS"]
#	fileTypes=[".pvals",".pvals",".score",".imp",".score"]
#	#csDataName = "newDataset_ker_v3_bc1_vs_cgl129_g8_f33_w5000"
#	#csDataName = "newDataset_ker_v3_bc3_vs_cgl43_g10_f33_w10000"
#	csDataName = "newDataset_simpleCS"
#	datasetNames=["newDataset","newDataset","newDataset","newDataset",csDataName]
#	logTransform = [True,True,False,False,False]
#	mafCutoffs = [0,15,0,15,0]
#
#	#mrIndex = 1  #The guiding (main) result
#
#	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
#	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
#	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
#	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
#
#	results_ls = []
#	for i in phenotypeIndices:
#		phenName=phed.getPhenotypeName(i)
#		phenName = phenName.replace("/","_div_")
#		phenName = phenName.replace("*","_star_")
#		phenIndex = phed.getPhenIndex(i)
#		
#		results = []
#		resultFiles = []
#		for j in range(0,len(methods)):
#			resultFile=resultsDirs[j]+methods[j]+"_"+datasetNames[j]+"_"+phenName+fileTypes[j]
#			resultFiles.append(resultFile)
#			print "Loading result file",resultFile
#			result = gwaResults.SNPResult(resultFile,snpsds,name=methods[j]+"_"+datasetNames[j]+"_"+phenName)
#			if logTransform[j]:
#				print "Log transformed the p-values"
#				result.negLogTransform()
#
#			result.filterMAF(minMaf=mafCutoffs[j])
#			results.append(result)
#			
#		results_ls.append(results)
#
#	return results_ls
#
#def loadResults(phenIndex=1,snpsds=None,phed=None,methodIndices=None, mafCutoffs = [0,15,0,15,0]):
#	"""
#	Load all the results needed.
#	"""
#
#	#ftIndices = [1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,80,81,82]
#	#ionomicsIndices = range(14,32)+range(83,101)
#
#	res_path="/Network/Data/250k/tmp-bvilhjal/"
#
#	resultsDirs = [res_path+"kw_results/",res_path+"emma_results/",res_path+"marg_results/",res_path+"rf_results/"]#,res_path+"cs_results/"]
#	methods=["KW","Emma","Marg","RF"]#,"CS"]
#	fileTypes=[".pvals",".pvals",".score",".imp"]#,".score"]
#	#csDataName = "newDataset_ker_v3_bc1_vs_cgl129_g8_f33_w5000"
#	#csDataName = "newDataset_ker_v3_bc3_vs_cgl43_g10_f33_w10000"
#	#csDataName = "newDataset_simpleCS"
#	csDataName = "newDataset_simpleCS_kerm_ws"
#	datasetNames=["newDataset","newDataset","newDataset","newDataset"]#,csDataName]
#	logTransform = [True,True,False,False]#,False]
#
#
#	#mrIndex = 1  #The guiding (main) result
#
#	#if not snpsds:
#	#	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
#	#	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
#	
#	if not phed:
#		phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
#		phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
#
#	if not methodIndices:
#		methodIndices=range(0,len(methods))
#	
#	phenName=phed.getPhenotypeName(phenIndex)
#	phenName = phenName.replace("/","_div_")
#	phenName = phenName.replace("*","_star_")
#
#		
#	results = []
#	resultFiles = []
#	for j in methodIndices:
#		try:
#			resultFile=resultsDirs[j]+methods[j]+"_"+datasetNames[j]+"_"+phenName+fileTypes[j]
#			resultFiles.append(resultFile)
#			print "Loading result file",resultFile
#			if not snpsds:
#				result = gwaResults.Result(resultFile,name=methods[j]+"_"+datasetNames[j]+"_"+phenName)
#			else:
#				result = gwaResults.SNPResult(resultFile,snpsds,name=methods[j]+"_"+datasetNames[j]+"_"+phenName)
#			if logTransform[j]:
#				print "Log transformed the p-values"
#				result.negLogTransform()
#			
#			result.filterMAF(minMaf=mafCutoffs[j])
#			results.append(result)
#		except Exception:
#			print "File not found error:",resultFile
#			
#	return results


#def drawRegionsPlots(signResult,results,snps,phenName,runId,res_path,displayRanks=True):
#	#Draw plots and write TAIR file..
#	for j in range(0,len(snps)):
#		snp = snps[j]
#		region = signResult.regions[j]   #(snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,regionRank)
#		#Print TAIR description
#		genes = getTairAnn(region[1],region[2],region[3])
#		if displayRanks:
#			snpPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_regionRank"+str(region[8])+"_"+signResult.name+"SNPplot.pdf"
#			tairFile = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_regionRank"+str(region[8])+"_"+signResult.name+"_tair.txt"
#		else:
#			snpPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_"+signResult.name+"SNPplot.pdf"
#			tairFile = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_"+signResult.name+"_tair.txt"
#			
#		f = open(tairFile,"w")
#		f.write("\n"+phenName+": "+str(region[1])+", "+str(region[2])+","+str(region[3])+"\n")
#		for gene in genes:
#			#print "Tair_id:"+str(gene.tairID)+", startPos:"+str(gene.startPos)+", endPos:"+str(gene.endPos)+", function description:"+str(gene.functionDescriptions)
#			f.write( "Tair_id:"+str(gene.tairID)+", startPos:"+str(gene.startPos)+", endPos:"+str(gene.endPos)+", function description:"+str(gene.functionDescriptions)+".\n")
#		f.close()
#		if displayRanks:
#			snpPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_regionRank"+str(region[8])+"_"+signResult.name+"SNPplot.pdf"
#		drawSNPPlot(snp,region,results,phenotypeName=signResult.name,genes=genes,filename=snpPlotFilename,displayRanks=displayRanks) #Plot SNPs
#
#		
#def plotRegion(region,results,snp,methodName,runId,res_path,genes,displayRanks=True):
#	#Draw a plot and write TAIR file..
#	#Print TAIR description
#	(chr,startPos,endPos,size) = region
#	if displayRanks:
#		snpPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(chr)+"_position"+str(startPos)+"_"+str(endPos)+"_regionRank"+str(snp.regionRank)+"_"+methodName+"_SNPplot.pdf"
#	else:
#		snpPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(chr)+"_position"+str(startPos)+"_"+str(endPos)+"_"+methodName+"SNPplot.pdf"
#	region = (None,startPos,endPos,chr,size,snp.score,snp.position,snp.rank,snp.regionRank)
#	drawSNPPlot(snp,region,results,phenotypeName=methodName,genes=genes,filename=snpPlotFilename,displayRanks=displayRanks) #Plot SNPs
#			
#
#def generateTairFile(region,runId,res_path):
#	#Draw a plot and write TAIR file..
#	#Print TAIR description
#	(chr,startPos,endPos,size) = region
#	genes = getTairAnn(startPos,endPos,chr) 
#	tairFile = res_path+"rid_"+str(runId)+"_chromosome"+str(chr)+"_position"+str(startPos)+"_"+str(endPos)+"_tair.txt"
#			
#	f = open(tairFile,"w")
#	f.write("\nRegion: "+str(chr)+", "+str(startPos)+", "+str(endPos)+"\n")
#	for gene in genes:
#		#print "Tair_id:"+str(gene.tairID)+", startPos:"+str(gene.startPos)+", endPos:"+str(gene.endPos)+", function description:"+str(gene.functionDescriptions)
#		f.write( "Tair_id:"+str(gene.tairID)+", startPos:"+str(gene.startPos)+", endPos:"+str(gene.endPos)+", function description:"+str(gene.functionDescriptions)+".\n")
#	f.close()
#
#	return genes
#
#
#
#
#
#def drawAllPlots(phenotypeIndices,n,method="rank_snps",res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/",runId="",statId="",window=[50000,50000]):
#	"""
#	1. Draws the following region-plots using top n ranked snps.
#
#	"""
#	if n:
#		runId = "top_"+str(n)+"_"+str(runId)
#	
#	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
#	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
#			
#	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
#	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
#	
#	statFilename = res_path+"rid_"+str(runId)+"_sid"+str(statId)+"_stat.txt"
#	statFile = open(statFilename,"w")
#	method_phenotype_stats = []
#	#phenotype_stats = []
#	#method_stats = []
#	#stats = []
#	totalRegionCount = 0
#	totalSNPsCount = 0
#	
#	unionRegions = []
#
#	methodUnions = []
#
#	for ind in range(0,len(phenotypeIndices)):
#		p_i = phenotypeIndices[ind]
#		print "\nNow working on phenotype id:",p_i
#		results = loadResults(p_i,snpsds,phed)
#		signResult_ls= []
#		stats_ls = []
#		for m_i in range(0,len(results)-1): #For all methods 
#			result = results[m_i]
#			print "\nInspecting result",result.name,":"
#			statFile.write("\nStatistics for "+str(result.name)+"\n")
#
#			#if method=="rank_snps":
#			#	signResult = results[m_i].getTopSnps(n)  #returns a result object.
#			#elif method=="rank_regions": 
#			#	signResult = results[m_i].getTopRegions(n,window=window)  #returns a result object.
#			
#			if m_i==1: #Emma
#				signResult = results[m_i].getTopSnps(200)
#				signResult.filterScoreCutoff(6)
#			elif m_i==0: #KW
#				signResult = results[m_i].getTopSnps(200)
#				signResult.filterScoreCutoff(8)
#			elif m_i==2: #Marg
#				signResult = results[m_i].getTopSnps(50)
#			elif m_i==3: #RF
#				signResult = results[m_i].getTopSnps(25)
#			snps = signResult.getRegionSNPs(window=window)
#			if len(signResult.scores):
#				print "Found",len(signResult.scores),"snps and",len(signResult.regions),"regions."
#			else:
#				print "Found no SNPs"
#			
#			#Recording various statistics
#			totalRegSize = 0
#			maxScore = 0
#			maxPos = None
#			maxChr = None
#			totalRegSize = 0
#			avgSize = None
#			numRegions = 0
#			if len(signResult.scores):
#				for reg in signResult.regions:
#					totalRegSize += reg[4]
#					if reg[5]>maxScore:
#						maxScore = reg[5]
#						maxPos = reg[6]
#						maxChr = reg[3]
#				avgSize = totalRegSize/float(len(signResult.regions))
#				numRegions = len(signResult.regions)
#				totalRegionCount +=len(signResult.regions)
#				totalSNPsCount += len(signResult.scores)
#			stats = (len(signResult.scores),numRegions,avgSize,maxScore,maxChr,maxPos)
#			statFile.write("Number of sign. SNPs: "+str(len(signResult.scores))+", number of sign. regions: "+str(numRegions)+", ave. region size: "+str(avgSize)+", max. score: "+str(maxScore)+", max. pos.: "+str((maxChr,maxPos))+".\n")
#			stats_ls.append(stats)
#
#			
#			if len(signResult.scores):
#				#Retrieve the accessions for this phenotype. (For the box-plot)
#				sAccessions = snps[0].accessions		
#				phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
#				phenIndex = phed.getPhenIndex(p_i)
#				pData = zip(*phed.phenotypeValues)[phenIndex]
#				pAccessions = phed.accessions
#				accIndices =[]
#				for pIndex1 in range(0,len(pAccessions)):
#					for pIndex2 in range(0,len(sAccessions)):
#						if pAccessions[pIndex1]==sAccessions[pIndex2] and pData[pIndex1]!="NA":
#							accIndices.append((pIndex1,pIndex2))
#							break			
#				for j in range(0,len(snps)):
#					snp = snps[j]
#					region = signResult.regions[j]
#					boxPlotFilename = res_path+"rid_"+str(runId)+"_chromosome"+str(region[3])+"_position"+str(region[1])+"_"+str(region[2])+"_snpRank"+str(region[7])+"_"+result.name+"Boxplot.pdf"
#					drawBoxPlot(snp,pData,accIndices,result.name,filename=boxPlotFilename) #Draw Box plot
#
#				#drawRegionsPlots(signResult,results,snps,result.name,runId,res_path)
#
#			signResult_ls.append(signResult)
#
#		
#		if ind==0:
#			for i in range(0,len(signResult_ls)):
#				methodUnions.append([signResult_ls[i].clone()])
#		else:
#			for i in range(0,len(signResult_ls)):
#				methodUnions[i].append(signResult_ls[i].clone())
#			
#		#Now union of all regions..
#		unionTopSnps = signResult_ls[0]
#		for i in range(1,len(signResult_ls)):
#			unionTopSnps.mergeWith(signResult_ls[i])  #FIXME: implement mergeWith 
#		snps = unionTopSnps.getRegionSNPs(window=window)
#		print "The union of the results contained",len(unionTopSnps.positions),"snps and",len(unionTopSnps.regions),"regions."
#
#		phenName = phed.getPhenotypeName(p_i)
#		phenName = phenName.replace("/","_div_")
#		phenName = phenName.replace("*","_star_")
#
#		unionTopSnps.name = "Union_"+phenName
#		
#		unionRegions.append(unionTopSnps.clone())
#		
#
#		if ind==0:
#			methodUnions.append([unionTopSnps.clone()])
#		else:
#			methodUnions[len(methodUnions)-1].append(unionTopSnps.clone())
#
#		#Recording various statistics
#		totalRegSize = 0
#		avgSize = None
#		numRegions = 0
#		if len(unionTopSnps.scores):
#			for reg in unionTopSnps.regions:
#				totalRegSize += reg[4]
#			avgSize = totalRegSize/float(len(unionTopSnps.regions))
#			numRegions = len(unionTopSnps.regions)
#			drawRegionsPlots(unionTopSnps,results,snps,unionTopSnps.name,runId,res_path,displayRanks=False)
#		stats = (len(unionTopSnps.scores),numRegions,avgSize,maxScore,maxChr,maxPos)
#		statFile.write("\nUnion: Number of sign. SNPs: "+str(len(unionTopSnps.scores))+", number of sign. regions: "+str(numRegions)+", ave. region size: "+str(avgSize)+".\n")		
#		statFile.flush()
#		stats_ls.append(stats) 
#
#
#
#		
#		method_phenotype_stats.append(stats_ls)
#
#
#	print len(unionRegions)
#	unionSnps = unionRegions[0]
#	for i in range(1,len(unionRegions)):
#		print i,"th union"
#		unionSnps.mergeWith(unionRegions[i])  #FIXME: implement mergeWith 
#	snps = unionSnps.getRegionSNPs(window=window)
#	totalRegSize = 0
#	if len(unionSnps.positions):
#		for reg in unionSnps.regions:
#			totalRegSize += reg[4]
#	
#	statStr = "The union of the results: Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"
#	print statStr
#	statFile.write("\n"+statStr)
#
#	##Unions accross methods.
#	j = 0
#	for results in methodUnions:
#		j += 1
#		unionSnps = unionRegions[0]
#		for i in range(1,len(unionRegions)):
#			print i,"th union"
#			unionSnps.mergeWith(unionRegions[i])   
#		snps = unionSnps.getRegionSNPs(window=window)
#		totalRegSize = 0
#		if len(unionSnps.positions):
#			for reg in unionSnps.regions:
#				totalRegSize += reg[4]
#		statStr = "The union of the results for method "+str(j)+": Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"
#		print statStr
#		statFile.write("\n"+statStr)
#	
#
#	statFile.close()




def getRegions(phenotypeIndices,runId="regions",statId="",window=[25000,25000],res_path="/Network/Data/250k/tmp-bvilhjal/snp_res_111208/",writeToFiles=True):
	"""
	1. Draws the following region-plots using top n ranked snps.

	"""
	runId = "final_"+str(runId)
	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
			
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)

	resultTypes = gwaResults._getStandardResultTypes_()

	results_map = gwaResults.loadResults(phenotypeIndices,resultTypes=resultTypes,phed=phed,snpsds=snpsds)
	
	import regionPlotter
	rp = regionPlotter.RegionPlotter(phenotypeIndices,snpsds,results_map=results_map)

	totalRegionCount = 0
	totalSNPsCount = 0
	
	result_ls = []
	
	full_results_ls = []

	if writeToFiles:
		overallStatFile=open(res_path+runId+"_over_all_stat.txt","w")
		
	
	for p_i in phenotypeIndices:
		print "\nNow working on phenotype id:",p_i
		results = results_map[p_i]
		newResults = []
		regionSet = set()  #regionSet
		for result in results:
			print "\nFiltering result",result.resultType.resultType,":"
			if result.resultType.resultType=="Emma":
				regionSet = regionSet.union(result.getTopRegions(10,window=window,minScore=5))
				
			elif result.resultType.resultType=="KW":
				regionSet = regionSet.union(result.getTopRegions(10,window=window,minScore=5))

			elif result.resultType.resultType=="Marg" or result.resultType.resultType=="RF":
				regionSet = regionSet.union(result.getTopRegions(5,window=window))
			else:
				print "Result wasn't recognized!"
				
		
		regionList = gwaResults.getRegions(regionSet,window=[25000,25000])
		totalRegSize = 0
		for region in regionList:
			totalRegSize += region.size
			
		averageSize = "NA"
		if len(regionList) and totalRegSize:
			averageSize = totalRegSize/float(len(regionList))
			
		if writeToFiles:
			overallStatFile.write("Phenotype "+phed.getPhenotypeName(p_i)+":  ")
			overallStatFile.write("Found "+str(len(regionSet))+" snps and "+str(len(regionList))+" regions with average size = "+str(averageSize)+"\n")
		else:
			print "Found",len(regionSet),"snps and",len(regionList),"regions with average size =", averageSize
		
			
		#Retrieving/updating various statistics about the regions
		totalRegSize = 0
		maxScore = 0
		for result in results:
			result.updateRegions(regionList)
		
		if writeToFiles:
			phenSpecificFile = open(res_path+runId+phed.getPhenotypeName(p_i)+"_stat.txt",'w')
		for region in regionList:
			maxPositions = {}
			if writeToFiles:
				 phenSpecificFile.write("\n"+str(region)+":\n")
			for result in results:
				maxScore = region.snpsInfo[result.resultType.name]["maxScore"]
				maxPos = region.snpsInfo[result.resultType.name]["maxPos"]
				maxRank = region.snpsInfo[result.resultType.name]["maxRank"]
				if maxPositions.has_key(maxPos):
					maxPositions[maxPos].append((result.resultType.name,maxScore,maxRank))
				else:
					maxPositions[maxPos] = [(result.resultType.name,maxScore,maxRank)]
				chr = region.chromosome
				if writeToFiles:
					phenSpecificFile.write("Result type: "+str(result.resultType.name)+", max. score: "+str(maxScore)+", max. rank.: "+str(maxRank)+", max. pos.: "+str(chr)+", "+str(maxPos)+".\n")
				else:
					print "Result type:",result.resultType.name,"Max. score:",maxScore,", max. rank.:",maxRank,", max. pos.:",(chr,maxPos),".\n"
			
			chr_pos_str = str(region.chromosome)+"_"+str(region.startPos)+"_"+str(region.endPos) 
			plotFileName = res_path+runId+"_"+phed.getPhenotypeName(p_i)+"_"+chr_pos_str+".pdf"
			tairFileName = res_path+runId+"_"+phed.getPhenotypeName(p_i)+"_"+chr_pos_str+"_stat.txt"
			rp.plotReg(region,p_i,pdfFile=plotFileName,tairFile=tairFileName)
			tairFile = open(tairFileName,"a")
			tairFile.write("\n\n"+"Max. SNPs:\n")
			for pos in maxPositions:
				plotFileName = res_path+runId+"_"+phed.getPhenotypeName(p_i)+"_"+chr_pos_str+"_wLD_"+str(pos)+".pdf"
				for (resType,score,rank) in maxPositions[pos]:
					tairFile.write(resType+": pos="+str(pos)+", rank="+str(rank)+", score="+str(score)+"\n")
				rp.plotReg(region,p_i,snpPos=pos,printTairInfo=False,pdfFile=plotFileName)
			tairFile.close()
				
		if writeToFiles:
			phenSpecificFile.close()

	if writeToFiles:
		overallStatFile.close()
	



#	elif len(phenotypeIndices)>1:
#		unionSnps = result_ls[0]
#		for i in range(1,len(result_ls)):
#			print i,"th union"
#			unionSnps.mergeWith(result_ls[i])  
#		unionSnps.getRegions(window=window)
#		totalRegSize = 0
#		if len(unionSnps.positions):
#			for reg in unionSnps.regions:
#				totalRegSize += reg[4]
#		print "The union of the results: Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"
#		del unionSnps
#
#
#		reg_met_table = gwaResults.RegionsTable(result_ls,window=window)
#		del result_ls
#		print "Garbage to be collected:",gc.garbage
#		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
#		
#		for i in range(0,len(reg_met_table.regions)):  #For all regions
#			region = reg_met_table.regions[i]
#			methods_snps_ls = reg_met_table.region_by_methods_table[i]
#			genes = generateTairFile(region,runId,res_path)
#			for j in range(0,len(methods_snps_ls)):
#				methods_snps = methods_snps_ls[j]
#				result_name = reg_met_table.resultNames[j]
#				if len(methods_snps):
#					maxSNP = methods_snps[0]
#					for s_i in range(1,len(methods_snps)):
#						if methods_snps[s_i].score>maxSNP.score:
#							maxSNP = methods_snps[s_i]
#					plotRegion(region,full_results_ls[j],maxSNP,result_name,runId,res_path,genes)		
#			reg_met_table.regions[i] = [] #An attempt to clean up the memory.
#			gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
				 
				


def getAlexMethodRegions(phenotypeIndices,runId="",statId="",window=[50000,50000],res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/"):
	"""
	1. Draws the following region-plots using top n ranked snps.

	"""
	runId = "Alex_method_"+str(runId)
	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
			
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")

	totalRegionCount = 0
	totalSNPsCount = 0
	
	result_ls = []
	
	full_results_ls = []

	for ind in range(0,len(phenotypeIndices)):
		p_i = phenotypeIndices[ind]
		print "\nNow working on phenotype id:",p_i
		results = loadResults(p_i,phed=phed,snpsds=snpsds)

		csResult = results[0].clone()  #FIXME this should be CS not KW
		emmaResult = results[1].clone()
		print "\nFiltering result",csResult.name,":"
		csResult.alexFiltering(emmaResult.scores,cutoff=6,window=window)
		del emmaResult

		csResult.getRegions(window=window)
		
		print "Found",len(csResult.scores),"snps and",len(csResult.regions),"regions."
			
		#Recording various statistics
		totalRegSize = 0
		maxScore = 0
		maxPos = None
		maxChr = None
		if len(csResult.scores):
			for reg in csResult.regions:
				totalRegSize += reg[4]
				if reg[5]>maxScore:
					maxScore = reg[5]
					maxPos = reg[6]
					maxChr = reg[3]
			
		print "Number of sign. SNPs: "+str(len(csResult.scores))+", number of sign. regions: "+str(len(csResult.regions))+", ave. region size: "+str(totalRegSize/float(len(csResult.regions)))+", max. score: "+str(maxScore)+", max. pos.: "+str((maxChr,maxPos))+".\n"
			
		result_ls.append(csResult)

		full_results_ls.append(results)
			
	unionSnps = result_ls[0]
	for i in range(1,len(result_ls)):
		print i,"th union"
		unionSnps.mergeWith(result_ls[i])  
	unionSnps.getRegions(window=window)
	totalRegSize = 0
	if len(unionSnps.positions):
		for reg in unionSnps.regions:
			totalRegSize += reg[4]
	print "The union of the results: Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"
	del unionSnps


	reg_met_table = gwaResults.RegionsTable(result_ls,window=window)
	del result_ls
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	print "Garbage:",gc.garbage

	for i in range(0,len(reg_met_table.regions)):  #For all regions
		region = reg_met_table.regions[i]
		methods_snps_ls = reg_met_table.region_by_methods_table[i]
		genes = generateTairFile(region,runId,res_path)
		for j in range(0,len(methods_snps_ls)):
			methods_snps = methods_snps_ls[j]
			result_name = reg_met_table.resultNames[j]
			if len(methods_snps):
				maxSNP = methods_snps[0]
				for s_i in range(1,len(methods_snps)):
					if methods_snps[s_i].score>maxSNP.score:
						maxSNP = methods_snps[s_i]
				plotRegion(region,full_results_ls[j],maxSNP,result_name,runId,res_path,genes)		
		reg_met_table.regions[i] = [] #An attempt to clean up the memory.
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
		 
				



def plotSignificantRegions(phenotypeIndex,runId="",statId="",window=[50000,50000],res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/"):
	"""
	Plots significant regions for the given phentope
	"""
	pass


def plotCommonSignificantRegions(phenotypeIndices,runId="",statId="",window=[50000,50000],res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/"):
	"""
	

	"""
	runId = "Common_regions_"+str(runId)
	import phenotypeData
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds

		
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f12_120508.csv"
	import dataParsers,snpsdata,phenotypeData
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
	
	results_map = gwaResults.loadResults(phenotypeIndices,phed=phed,filterPercentile=filterPercentile)


	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
			
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")

	totalRegionCount = 0
	totalSNPsCount = 0
	
	result_ls = []
	
	full_results_ls = []

	phenotype_2_regions = {}

	for ind in range(0,len(phenotypeIndices)):
		p_i = phenotypeIndices[ind]
		print "\nNow working on phenotype id:",p_i
		results = loadResults(p_i,phed=phed,snpsds=snpsds)

		kwResult = results[0].clone()  
		emmaResult = results[1].clone()
		print "\nFiltering result",kwResult.name,":"
		kwResult.filterScoreCutoff(10)
		emmaResult.filterScoreCutoff(5)
		
		kwResult.mergeWith(emmaResult)		

		kwResult.getRegions(window=window)

		phenotype_2_regions[ind] = kwResult.regions

		print "Found",len(kwResult.scores),"snps and",len(kwResult.regions),"regions."
			
		result_ls.append(kwResult)

		full_results_ls.append(results)
			
	unionSnps = result_ls[0]
	for i in range(1,len(result_ls)):
		print i,"th union"
		unionSnps.mergeWith(result_ls[i])  
	unionSnps.getRegions(window=window)
	totalRegSize = 0
	if len(unionSnps.positions):
		for reg in unionSnps.regions:
			totalRegSize += reg[4]
	print "The union of the results: Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"
	del unionSnps


	reg_met_table = gwaResults.RegionsTable(result_ls,window=window)
	del result_ls
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	for i in range(0,len(reg_met_table.regions)):  #For all regions
		region = reg_met_table.regions[i]
		methods_snps_ls = reg_met_table.region_by_methods_table[i]
		genes = generateTairFile(region,runId,res_path)
		for j in range(0,len(methods_snps_ls)):
			methods_snps = methods_snps_ls[j]
			result_name = reg_met_table.resultNames[j]
			if len(methods_snps):
				maxSNP = methods_snps[0]
				for s_i in range(1,len(methods_snps)):
					if methods_snps[s_i].score>maxSNP.score:
						maxSNP = methods_snps[s_i]
				plotRegion(region,full_results_ls[j],maxSNP,result_name,runId,res_path,genes)		
		reg_met_table.regions[i] = [] #An attempt to clean up the memory.
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
		 
				








def countingRegions(phenotypeIndices,n,method="rank_snps",res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/",runId="",statId="",window=[50000,50000]):
	"""
	"""
	runId = "top_"+str(n)+"_"+str(runId)
	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
			
	statFilename = res_path+"rid_"+str(runId)+"_sid"+str(statId)+"_stat.txt"
	statFile = open(statFilename,"w")
	method_phenotype_stats = []

	totalRegionCount = 0
	totalSNPsCount = 0
	
	unionRegions = []

	methodUnions = []

	for ind in range(0,len(phenotypeIndices)):
		p_i = phenotypeIndices[ind]
		print "\nNow working on phenotype id:",p_i
		results = loadResults(p_i,phed=phed)
		signResult_ls= []
		stats_ls = []
		for m_i in range(0,len(results)): #For all methods 
			result = results[m_i]
			print "\nInspecting result",result.name,":"
			statFile.write("\nStatistics for "+str(result.name)+"\n")

			if method=="rank_snps":
				signResult = results[m_i].getTopSnps(n)  #returns a result object.
			elif method=="rank_regions": 
				signResult = results[m_i].getTopRegions(window=window)  #returns a result object.
			signResult.getRegions()
			print "Found",len(signResult.scores),"snps and",len(signResult.regions),"regions."
			
			#Recording various statistics
			totalRegSize = 0
			maxScore = 0
			maxPos = None
			maxChr = None
			if len(signResult.scores):
				for reg in signResult.regions:
					totalRegSize += reg[4]
					if reg[5]>maxScore:
						maxScore = reg[5]
						maxPos = reg[6]
						maxChr = reg[3]
			
			statFile.write("Number of sign. SNPs: "+str(len(signResult.scores))+", number of sign. regions: "+str(len(signResult.regions))+", ave. region size: "+str(totalRegSize/float(len(signResult.regions)))+", max. score: "+str(maxScore)+", max. pos.: "+str((maxChr,maxPos))+".\n")
			
			signResult_ls.append(signResult)

		
		if ind==0:
			for i in range(0,len(signResult_ls)):
				methodUnions.append([signResult_ls[i].clone()])
		else:
			for i in range(0,len(signResult_ls)):
				methodUnions[i].append(signResult_ls[i].clone())
			
		#Now union of all regions..
		unionTopSnps = signResult_ls[0]
		for i in range(1,len(signResult_ls)):
			unionTopSnps.mergeWith(signResult_ls[i])  #FIXME: implement mergeWith 
		unionTopSnps.getRegions()
		print "The union of the results contained",len(unionTopSnps.positions),"snps and",len(unionTopSnps.regions),"regions."

		phenName = phed.getPhenotypeName(p_i)
		phenName = phenName.replace("/","_div_")
		phenName = phenName.replace("*","_star_")

		unionTopSnps.name = "Union_"+phenName
		
		unionRegions.append(unionTopSnps.clone())
		

		if ind==0:
			methodUnions.append([unionTopSnps.clone()])
		else:
			methodUnions[len(methodUnions)-1].append(unionTopSnps.clone())

		#Recording various statistics
		totalRegSize = 0
		if len(unionTopSnps.scores):
			for reg in unionTopSnps.regions:
				totalRegSize += reg[4]
		stats = (len(unionTopSnps.scores),len(unionTopSnps.regions),totalRegSize/float(len(unionTopSnps.regions)),maxScore,maxChr,maxPos)
		statFile.write("\nUnion: Number of sign. SNPs: "+str(len(unionTopSnps.scores))+", number of sign. regions: "+str(len(unionTopSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionTopSnps.regions)))+".\n")		
		statFile.flush()
		stats_ls.append(stats) 

		
		method_phenotype_stats.append(stats_ls)
		#FIXME calc. variance explained
	statFile.close()

	print len(unionRegions)
	unionSnps = unionRegions[0]
	for i in range(1,len(unionRegions)):
		print i,"th union"
		unionSnps.mergeWith(unionRegions[i])  
	unionSnps.getRegions()
	totalRegSize = 0
	if len(unionSnps.positions):
		for reg in unionSnps.regions:
			totalRegSize += reg[4]
	print "The union of the results: Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"

	j = 0
	for results in methodUnions:
		j += 1
		unionSnps = results[0]
		for i in range(1,len(unionRegions)):
			print i,"th union"
			unionSnps.mergeWith(unionRegions[i])  
		unionSnps.getRegions()
		totalRegSize = 0
		if len(unionSnps.positions):
			for reg in unionSnps.regions:
				totalRegSize += reg[4]
		print "The union of the results for method "+str(j)+": Number of sign. SNPs: "+str(len(unionSnps.scores))+", number of sign. regions: "+str(len(unionSnps.regions))+", ave. region size: "+str(totalRegSize/float(len(unionSnps.regions)))+".\n"




def drawCutoffPlots():
	"""
	1. Draws the following region-plots using top n ranked snps.
	2. Reports statistics about the regions, for each phenotype.

	Methods:
	- KW
	- Emma
	- Marg
	- Emma
	- CS
	- Union
	"""
	pass


def drawAllGWPlots(phenotypeIndices=None,res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/",runId="gwPlot"):
	"""
	Draws all the GWA plots for 6 methods.
	"""
	import plotResults

	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds
		
#	p_i1 = 5
#	p_i2 = 31
#	l = [str(p_i1)+"_o96acc",str(p_i1)+"_o192acc"]
#	for s in l:
#		rt = gwaResults.ResultType(resultType="KW",name=s)
#		resultFile = "/Users/bjarni/tmp/KW_"+s+".pvals"
#		pdfFile = "/Users/bjarni/tmp/KW_"+s+".pdf"
#		pngFile = "/Users/bjarni/tmp/KW_"+s+"acc.png"
#		result = gwaResults.Result(resultFile,name=s, resultType=rt)
#		result.negLogTransform()
#		print "\nPlotting result",result.name,":"
#		plotResults.plotResult(result,pdfFile,pngFile,ylab="- log10 pvalue",plotBonferroni=(result.resultType.resultType=="Emma" or result.resultType.resultType=="KW"))
	
	results_map = gwaResults.loadResults(phenotypeIndices,phed=phed)
			
	for p_i in phenotypeIndices:
		print "\nNow working on phenotype id:",p_i
		phenName = phed.getPhenotypeName(p_i)
		results = results_map[p_i]
		for result in results: #For all methods 
			pdfFile = res_path+"rid_"+runId+"_"+phenName+"_"+result.resultType.name+".pdf" 
			pngFile = res_path+"rid_"+runId+"_"+phenName+"_"+result.resultType.name+".png" 
			print "\nPlotting result",result.name,":"
			plotResults.plotResult(result,pdfFile,pngFile,ylab=result.resultType.name,plotBonferroni=(result.resultType.resultType=="Emma" or result.resultType.resultType=="KW"))
			

		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

def drawSecondRunGWPlots(phenotypeIndices=None,res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/",runId="gwPlot"):
	"""
	Draws all the GWA plots for 6 methods.
	"""
	import plotResults

	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds
		
#	p_i1 = 5
#	p_i2 = 31
#	l = [str(p_i1)+"_o96acc",str(p_i1)+"_o192acc"]
#	for s in l:
#		rt = gwaResults.ResultType(resultType="KW",name=s)
#		resultFile = "/Users/bjarni/tmp/KW_"+s+".pvals"
#		pdfFile = "/Users/bjarni/tmp/KW_"+s+".pdf"
#		pngFile = "/Users/bjarni/tmp/KW_"+s+"acc.png"
#		result = gwaResults.Result(resultFile,name=s, resultType=rt)
#		result.negLogTransform()
#		print "\nPlotting result",result.name,":"
#		plotResults.plotResult(result,pdfFile,pngFile,ylab="- log10 pvalue",plotBonferroni=(result.resultType.resultType=="Emma" or result.resultType.resultType=="KW"))
	
	results_map = gwaResults.loadResults(phenotypeIndices,phed=phed,secondRun=True)
			
	for p_i in phenotypeIndices:
		print "\nNow working on phenotype id:",p_i
		phenName = phed.getPhenotypeName(p_i)
		results = results_map[p_i]
		pdfFile = res_path+"rid_"+runId+"_"+phenName+"_"+results[0].resultType.name+".pdf" 
		pngFile = res_path+"rid_"+runId+"_"+phenName+"_"+results[0].resultType.name+".png" 
		print "\nPlotting result",results[0].name,":"
		plotResults.plotResultWithSecondRun(results[0],results[1],pngFile=pngFile,ylab=results[0].resultType.name,plotBonferroni=(results[0].resultType.resultType=="Emma" or results[0].resultType.resultType=="KW"),)
			

		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

def drawPowerGWPlots(phenotypeIndices=None,res_path="/Network/Data/250k/tmp-bvilhjal/power_analysis/results/",runId="gwPlot"):
	"""
	Draws all the GWA plots for 6 methods.
	"""
	import plotResults

	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds

		
	l = ["original96","original192","permTest","permTest"]
	perm_count = [1,1,10,10]
	tailString = ["","","r65","r112"]
	for p_i	in phenotypeIndices:
		for i in range(0,len(l)):
			for j in range(0,perm_count[i]):
				phenName = phed.getPhenotypeName(p_i)
				if perm_count[i]==1: 
					filename = res_path+"KW_raw_"+l[i]+"_"+phenName+".pvals"
				else:
					filename = res_path+"KW_raw_"+l[i]+"_"+phenName+"_"+tailString[i]+"_"+str(j)+".pvals"
				rt = gwaResults.ResultType(resultType="KW",name=l[i])
				pdfFile = filename+".pdf"
				pngFile = filename+".png"
				result = gwaResults.Result(filename,name=l[i], resultType=rt)
				result.negLogTransform()
				print "\nPlotting result",result.name,":"
				plotResults.plotResult(result,pdfFile,pngFile,ylab="- log10 pvalue",plotBonferroni=True)

		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..



def drawSegregatingGWAPlots(phenotypeIndices=None,ecotypePairList=None,res_path="/Network/Data/250k/tmp-bvilhjal/cross_plots/",runId="gwPlot"):
	"""
	Draws all the GWA plots for 6 methods.
	"""
	
	if not ecotypePairList:
		#ecotypePairList = [("6944","6977","NFA8_Van0"),("6046", "6962","Lov5_Sha"),("6977","6903","Van0_Bor4"),
		#			("6916","7514","Est1_RRS7"),("6916","6904","Est1_Br0"),("6904","6906","Br0_C24"),
		#			("6911","6906","Cvi0_RRS7"),("6899","6046","Bay0_Lov5"),("6944","6903","NFA8_Bor4"),
		#			("6906","7515","C24_RRS10"),("6962","8215","Sha_Fei0"),("6970","6972","Ts1_Tsu1")]
		ecotypePairList = [("6916","6904","Est1_Br0"),("6904","6906","Br0_C24"),
					("6911","6906","Cvi0_RRS7"),("6899","6046","Bay0_Lov5"),("6944","6903","NFA8_Bor4"),
					("6906","7515","C24_RRS10"),("6962","8215","Sha_Fei0"),("6970","6972","Ts1_Tsu1")]
					
	import plotResults,snpsdata,regionPlotter,pdb,gwaResults
	
			
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_transformed_publishable_v2.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
	#pdb.set_trace()

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds
			
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f11_100708.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])

	resultTypes = gwaResults._getStandardResultTypes_()

	results_map = gwaResults.loadResults(phenotypeIndices,resultTypes=resultTypes,phed=phed,snpsds=snpsds)

	regions = [gwaResults.fri_region_small, gwaResults.flc_region_small]
	rp = regionPlotter.RegionPlotter(snpsds=snpsds)
			

	for (ecotype1,ecotype2,crossName) in ecotypePairList:
		for p_i in phenotypeIndices:
			print "\nNow working on phenotype id:",p_i
			phenName = phed.getPhenotypeName(p_i)
			results = results_map[p_i]		
			filteredResults = []
			for result in results: #For all methods 
				res = result.clone()
				res.filterNonSegregatingSnps(ecotype1,ecotype2,snpsd)
				res2 = res.clone()
				filteredResults.append(res2)
				pdfFile = res_path+"rid_"+runId+"_"+phenName+"_"+crossName+"_"+result.resultType.name+".pdf" 
				print "\nPlotting result",result.name,":"
				plotResults.plotResult(res,pdfFile,ylab=result.resultType.name,plotBonferroni=(result.resultType.resultType=="Emma" or result.resultType.resultType=="KW"))
	
			for region in gwaResults.fri_flc_regions:
				chr_pos_str = region.get_chr_pos_str()
				plotFileName = res_path+runId+"_"+phenName+"_"+crossName+"_"+region.name+"_"+chr_pos_str+".pdf"
				tairFileName = res_path+runId+"_"+region.name+"_"+chr_pos_str+"_stat.txt"
				rp.plotReg(region,p_i,pdfFile=plotFileName,tairFile=tairFileName,results=filteredResults)
			gc.collect()  #Calling garbage collector, in an attempt to clean up memory..


def plotPeakWidths(phenotypeIndices,pdfFile=None,pngFile=None,ylab="",cgl_id=None, specialGenes=None, shade_scale=0.9,filterPercentile=0.9996):
	import phenotypeData
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds

		
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f12_120508.csv"
	import dataParsers,snpsdata,phenotypeData
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
	
	#results_map = gwaResults.loadResults(phenotypeIndices,phed=phed,filterPercentile=filterPercentile)
	results_map = gwaResults.loadResults(phenotypeIndices,phed=phed,filterCutoffs=[8,4])


	#Load CGL.
	if cgl_id:
		candGenes = gwaResults.getCandidateGeneList(cgl_id)

		print "Determining the dot positions"
		chr_pos_set = set()
		for p_i in phenotypeIndices:
			results = results_map[p_i]
			for result in results: #For all methods 
				print "\nRetrieving chr_pos list for ",result.name,":"
				chrPosList = result.getChrPos()
				chr_pos_set = chr_pos_set.union(set(chrPosList))
		chr_pos_list = list(chr_pos_set)
		chr_pos_list.sort()
		print "len(chr_pos_list):",len(chr_pos_list)
		chr_pos_map = {1:[],2:[],3:[],4:[],5:[]}
		for (chr,pos) in chr_pos_list:
			 chr_pos_map[chr].append(pos)

		window = 10000
		newCandGenes = []
		#Check whether cg are close to dots...
		for gene in candGenes:
			for pos in chr_pos_map[gene.chromosome]:
				if pos>gene.startPos-window and pos<gene.endPos+window:
					print gene.startPos, gene.endPos, pos
					newCandGenes.append(gene)
					break
		print "len(newCandGenes):",len(newCandGenes),", len(candGenes):",len(candGenes)
		candGenes = newCandGenes
		
		candGeneInfoFile = "/Users/bjarni/tmp/flowering_time_cand_genes.txt"
		f = open(candGeneInfoFile,"w")
		f.write("chr, pos, tairID, name, description\n")
		for gene in candGenes:
			st = str(gene)+"\n"
			f.write(st)
		f.close()
			
		
	import pylab
	for chromosome in[1,2,3,4,5]:
		y_shift = 0
		print "Working on chromosome",chromosome
		x_size=20.0*snpsds[chromosome-1].positions[-1]/float(snpsds[0].positions[-1])
		y_size=7.0*len(phenotypeIndices)/30.0+1
		pylab.figure(figsize=(x_size,y_size))
		
				


		if cgl_id:
			for gene in candGenes:
				if gene.chromosome==chromosome:
					pylab.axvspan(gene.startPos,gene.endPos, facecolor=(0.1,0.6,0.1),edgecolor=(0.1,0.6,0.1), alpha=1,linewidth=1.2)
					#pylab.plot([gene.startPos,gene.startPos],[-0.4,len(phenotypeIndices)-0.4],"--",color=(0.2,0.2,0.2),linewidth=0.8)
					#pylab.plot([gene.endPos,gene.endPos],[-0.4,len(phenotypeIndices)-0.4],"--",color=(0.2,0.2,0.2),linewidth=0.8)
		
		if specialGenes:
			for (gene_name,gene) in specialGenes:
				if gene.chromosome==chromosome:
					specialColor = "g" #(0.7,0.5,0.2)
					pylab.plot([gene.startPos,gene.startPos],[-0.4,len(phenotypeIndices)-0.4],"-",color=specialColor,linewidth=1.2)
					pylab.plot([gene.endPos,gene.endPos],[-0.4,len(phenotypeIndices)-0.4],"-",color=specialColor,linewidth=1.2)
					x_shift = 200000
					pylab.text(gene.startPos-x_shift,(len(phenotypeIndices)-1)*1.06, gene_name,color=specialColor)
				
		
		for p_i in phenotypeIndices:
			print "\nNow working on phenotype id:",p_i
			phenName = phed.getPhenotypeName(p_i)
			results = results_map[p_i]
			for result in results: #For all methods 
				print "\nPlotting result",result.name,":"
				scorePosList = result.getChrScorePos(chromosome)
				scorePosList.sort()
				print "len(scorePosList):",len(scorePosList)
				for (score,pos) in scorePosList:
					shade = shade_scale
					if result.resultType.resultType=='KW':
						shade = shade_scale*(max(min(score,10.0)-8.0,0))/2.0
						col = (shade_scale-shade,shade_scale-shade,1)
						pylab.plot([pos],[y_shift+0.2],".",color=col)
					if result.resultType.resultType=='Emma':
						shade = shade_scale*(max(min(score,6.0)-4.0,0))/2.0
						col = (1,shade_scale-shade,shade_scale-shade)
						pylab.plot([pos],[y_shift],".",color=col)
				print "\nFinished plotting result",result.name,":"
			y_shift += 1
		
		x_shift = 250000
		pylab.axis([0-x_shift,snpsds[chromosome-1].positions[-1]+x_shift,0-0.05*(y_shift),(y_shift-1)+0.05*(y_shift)])
		ticks = range(0,snpsds[chromosome-1].positions[-1]+x_shift,1000000)
		s_ticks = []
		for tick in ticks:
			if tick%5000000==0:
				s_ticks.append(str(int(tick/1000000)))
			else: 
				s_ticks.append("")
		pylab.xticks(ticks,s_ticks)
		s_ticks = []
		for p_i in phenotypeIndices:
			phenName = phed.getPhenotypeName(p_i)
			phenNameList = phenName.split("_")[1:]
			phenName = " ".join(phenNameList)
			s_ticks.append(phenName)
		#s_ticks = [""]*len(phenotypeIndices)
		ticks = [i+0.12 for i in range(0,len(phenotypeIndices))]
		pylab.yticks(ticks,s_ticks,size="medium")
		pylab.ylabel(ylab,size='x-large')
		x_shift = 1600000
		pylab.text(snpsds[chromosome-1].positions[-1]/2.0-x_shift,(len(phenotypeIndices)-1)*1.08,"Chromosome "+str(chromosome),size='x-large')
		pylab.xlabel("Position (mb)",size='x-large')
		pylab.axhspan(-2, len(phenotypeIndices)+1, facecolor='0.8', alpha=0.5)
		
		if pngFile:
			pylab.savefig(pngFile+"_"+str(chromosome)+".png",format="png")    	
		if not (pngFile or pdfFile):
			pylab.show()
		pylab.clf()
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	




def _test1_():
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f11_100708.csv"
	import dataParsers,snpsdata,phenotypeData
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
	
	
	resdir = "/Network/Data/250k/tmp-bvilhjal/"
	filename = resdir+"ecotype_acc_stockParent.csv"
	f = open(filename,"w")
	ecotype_dict = phenotypeData._getEcotypeIdToStockParentDict_()
	l = []
	for acc in snpsd.accessions:
		acc = int(acc)
		l.append(ecotype_dict[acc])
	l.sort()
	for (name,csID) in l:
		f.write(str(name)+", "+str(csID)+"\n")
	f.close()


	
ftIndices = [1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,80,81,82]  # pid 60 was removed since it's not really a flowering time
ionomicsIndices = range(14,32)
#getQuantitativeIndices()[0:10]


		
if __name__ == '__main__':
	#_test1_()
	#combinedFormat("/Network/Data/250k/tmp-bvilhjal/results_final_090809/",ftIndices)
	#getRegions([161,163,164,165,166,173,174,175]+range(179,187))
	#drawSegregatingGWAPlots(phenotypeIndices=[1,7],runId="gwPlot")
	#drawAllGWPlots([224,225,227,229,231,236,239,241],res_path="/Users/bjarni/tmp/",runId="r")
	#drawAllGWPlots([226,228,230,232,233,234,235,237,238,240,242],res_path="/Users/bjarni/tmp/",runId="r")
	drawSecondRunGWPlots([17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],res_path="/Users/bjarni/tmp/",runId="r")
	#drawPowerGWPlots([187,188,189,190],runId="r")
	#getAlexMethodRegions(ftIndices,runId="test",statId="",window=[50000,50000],res_path="/Network/Data/250k/tmp-bvilhjal/snp_res_alex/")
	#drawAllPlots(ftIndices,None,method="rank_snps",res_path="/Network/Data/250k/tmp-bvilhjal/snp_res_final/w25000/",runId="final",statId="",window=[25000,25000])
	#countingRegions(ftIndices,50,res_path="/Network/Data/250k/tmp-bvilhjal/snp_res/top10Regions/ft/",runId="ft",statId="1",method="rank_snps")
	#_run_()
	#_testGenePlot_()
	#specialGenes = [("SVP",gwaResults.Gene(2,9586954,9590973)),("FRI",gwaResults.Gene(4,269026,271503)), ("FLC",gwaResults.Gene(5,3173498,3179449))]
	#plotPeakWidths(phenotypeData.categories_2_phenotypes[1],pngFile="/Users/bjarni/tmp/FT_width_cutoff_KW10_8_E6_4_l",ylab="Phenotypes",cgl_id=28,shade_scale=0.6) #
	#plotPeakWidths([1,2],pngFile="/Users/bjarni/tmp/test",ylab="Phenotypes",cgl_id=28,shade_scale=0.6) #
	#plotPeakWidths(phenotypeData.categories_2_phenotypes[2],pngFile="/Users/bjarni/tmp/FT_related_width",ylab="Flowering time related")
	#plotPeakWidths([1,2],pngFile="/Users/bjarni/tmp/FT_related_width",ylab="Flowering time related",cgl_id=129)
	#plotPeakWidths(phenotypeData.categories_2_phenotypes[3],pngFile="/Users/bjarni/tmp/Disease_resistance_width",ylab="Disease resistance",cgl_id=130)
	#plotPeakWidths(phenotypeData.categories_2_phenotypes[4],pngFile="/Users/bjarni/tmp/Ionomics_width",ylab="Ionomics",cgl_id=20)
	#plotPeakWidths(ionomicsIndices,pngFile="/Users/bjarni/tmp/Ionomics_old_width",cgl_id=43)



#FRI deletions
#abline(v=268809)
#abline(v=269962)
