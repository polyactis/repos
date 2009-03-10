#!/usr/bin/env python2.5
"""
11/20/08

Description:
This file analyzes the haplotype structure 

"""
from numpy import *
from rpy import r
import random, pylab, gc, pdb, math
import dataParsers, snpsdata, phenotypeData

def calcKinship(snps):
	"""
	Requires EMMA to be installed.
	"""
	a = array(snps)
	r.library("emma")
	return r.emma_kinship(a)

def runEmma(phed,p_i,k,snps):
	#Assume that the accessions are ordered.
	i = phed.getPhenIndex(p_i)
	r.library("emma")
	phenValues = []
	for vals in phed.phenotypeValues:
		phenValues.append(float(vals[i]))
	phenArray = array([phenValues])
	snpsArray = array(snps)
	res = r.emma_REML_t(phenArray,snpsArray,k)
	#print res
	return res


def simulateSNPs(k,n):
	#[pos][acc]
	pass
	
def sampleSNPs(snps,n,withReplacement=True):
	l = len(snps)
	i = 0
	snpSample = []
	if withReplacement:
		while i<n:
			r = random.randint(0,l-1)
			snpSample.append(snps[r])
			i += 1
	else:
		while i<n:
			r = random.randint(0,len(snps)-1)
			snpSample.append(snps[r])
			snps.remove(snps[r])
			i += 1		
	return snpSample
	


def getKinshipDiffs(snpsd,normalGlobalKinship,windowSize=500000,binSize=50000,minSNPs=200):	
	kDiffs = []
	binPos = []
	maxPos = snpsd.positions[len(snpsd.positions)-1]
	numBins = int(maxPos/binSize-windowSize/binSize+1)
	startPos = 0
	endPos = windowSize
	i = 0
	j = 0
	localKinships = []
	for b_i in range(0,numBins):
		while i < len(snpsd.positions)-1 and snpsd.positions[i]<startPos:
			i += 1
		while j < len(snpsd.positions)-1 and snpsd.positions[j]<=endPos:
			j += 1
		if j-1-i<minSNPs:
			print "there are two few SNPs"
			kDiffs.append(0)
			localKinships.append(None)
		else:
			print j-1-i, "snps found in region ",snpsd.positions[i],"-",snpsd.positions[j-1]
			if j-1-i<minSNPs:
				kDiff = 0
				print "Too few SNPs to estimate K"
			else:
				snps = sampleSNPs(snpsd.snps[i:j],minSNPs,withReplacement=False)
				k = calcKinship(snps)
				localKinships.append(k)
				normal_k = k/mean(k)
				kDiff = mean(abs(normal_k-normalGlobalKinship))
				print "Average k matrix difference:",kDiff
			kDiffs.append(kDiff)
		binPos.append((b_i+(float(windowSize)/binSize)*0.5)*binSize)
		startPos += binSize
		endPos += binSize
	return (kDiffs,binPos,localKinships)
	

def getEmmaDiffs(snpsd,phed,p_i,globalKinship,localKinships=None,windowSize=500000,binSize=100000,nSNPs=200,simulate=False,minSNPs=200):
	emmaDiffs = []
	binPos = []
	maxPos = snpsd.positions[len(snpsd.positions)-1]
	numBins = int(maxPos/binSize-windowSize/binSize+1)
	startPos = 0
	endPos = windowSize
	i = 0
	j = 0
	for b_i in range(0,numBins):
		while i < len(snpsd.positions)-1 and snpsd.positions[i]<startPos:
			i += 1
		while j < len(snpsd.positions)-1 and snpsd.positions[j]<=endPos:
			j += 1
		if j-1-i<minSNPs:
			print "there are two few SNPs"
			emmaDiffs.append(0)
		else:
			print j-1-i, "snps found in region ",snpsd.positions[i],"-",snpsd.positions[j-1]
			if localKinships!=None and localKinships[b_i]!=None:
				k = localKinships[b_i]
			else:
				snps = sampleSNPs(snpsd.snps[i:j],minSNPs)
				k = calcKinship(snps)
			numTries = 0
			unsuccessful = True
			while numTries < 10 and unsuccessful:
				try:
					snps = sampleSNPs(snpsd.snps[i:j],nSNPs)
					emma_res_local = runEmma(phed, p_i, k, snps)
					pvals_local = list(emma_res_local["ps"])
					pvals_local = [-math.log10(pval) for pval in pvals_local]
					emma_res_global = runEmma(phed, p_i, globalKinship, snps)
					pvals_global = list(emma_res_global["ps"])
					pvals_global = [-math.log10(pval) for pval in pvals_global]
					avgPvalDiff = (sum(pvals_global)-sum(pvals_local))/len(pvals_local)
					print "Average pvalue difference:",avgPvalDiff
					unsuccessful = False
				except Exception:
					print "Emma failed on",(numTries+1),"try."
					avgPvalDiff = 0
				numTries +=1
			emmaDiffs.append(avgPvalDiff)
		binPos.append((b_i+(float(windowSize)/binSize)*0.5)*binSize)
		startPos += binSize
		endPos += binSize
	return (emmaDiffs,binPos)





def _locate_max_kinship_acc_(k):
	ai1 = 0
	ai2 = 0
	maxKinship = 0.0
	for i in range(1,len(k)):
		for j in range(0,i):
			if k[i,j]>=maxKinship:
				ai1 = j
				ai2 = i
				maxKinship = k[i,j]
	return (maxKinship,ai1,ai2)

def _merge_accessions_(ai1,ai2,acc_groups):
	acc_groups[ai1] = acc_groups[ai1]+acc_groups[ai2]
	acc_groups.remove(acc_groups[ai2])

def _mergeSNPs_(snpList):
	new_snp = []
	t_snps = zip(*snpList)
	for allele in t_snps:  #For all accessions
		if mean(allele) ==0.5:
			new_snp.append(round(random.random()))
		else:
			new_snp.append(round(mean(allele)))
	return new_snp
	
def _update_snps_(snps,acc_groups):
	t_snps = map(list,zip(*snps))
	new_snps = []
	for group in acc_groups:
		snpList = []
		for i in group:
			snpList.append(t_snps[i])
		m_snp = _mergeSNPs_(snpList)
		new_snps.append(m_snp)
	new_snps = map(list,zip(*new_snps))
	return new_snps
		
	
def _update_phenotype_(phenVals,acc_groups):
	new_phenVals = []
	for group in acc_groups:
		goodCount = 0
		totSum = 0.0
		for a_i in group:
			if phenVals[a_i] !="NA":
				goodCount += 1
				totSum += float(phenVals[a_i])
		#print totSum/goodCount
		phen_val = totSum/goodCount
		new_phenVals.append(phen_val)
	return new_phenVals
	
def _run_kw_(snps,phenotypeValues,verbose=False):
	print "Running KW on",len(snps),"snps, and",len(phenotypeValues),"phenotype values."
	pvals = []
	for snp in snps:
		res = r.kruskal_test(phenotypeValues,snp)
		#print snp,res,phenotypeValues
		#pdb.set_trace()
		pvals.append(res["p.value"])
	return pvals

def get_KW_pvals(snps,positions,phed,p_i,kinshipThreshold = 0.95, method = "EMMA"):
	"""
	Takes a set of SNPs, which based on calculates a kinship, and then starts merging strains that look alike, until a threshold is reached.
	
	Methods:
	"Emma"
	"KW"
	"""
	accessions = phed.accessions
	phenVals = phed.getPhenVals(p_i)  #FIXME: make sure that the accessions are ordered correctly.
	pre_count = len(phenVals)
	k = calcKinship(snps)
	acc_groups = []
	for i in range(0,len(accessions)):
		acc_groups.append([i])
		
	(maxKinship,ai1,ai2) = _locate_max_kinship_acc_(k)
	print "Starting max kinship is",maxKinship
	while maxKinship > kinshipThreshold and len(acc_groups)>2:
		_merge_accessions_(ai1,ai2,acc_groups) #It updates the acc_groups and map automatically
		merged_snps= _update_snps_(snps,acc_groups)
		print len(merged_snps[0])
		k = calcKinship(merged_snps)
		(maxKinship,ai1,ai2) = _locate_max_kinship_acc_(k)
		print maxKinship
		
	merged_phenVals= _update_phenotype_(phenVals,acc_groups)
	print "Grouping reduced the # of indiv. by",pre_count-len(merged_phenVals)
	
	new_snps = []
	new_positions = []
	for i in range(0,len(merged_snps)):
		snp = merged_snps[i]
		if snp.count(0)>0 and snp.count(1)>0:
			new_snps.append(snp)
			new_positions.append(positions[i])
	
	pvals = _run_kw_(new_snps,merged_phenVals)
	#print pvals
	
	return (pvals,new_positions,acc_groups)

def _plotKinshipDiffs_():
	
	filterProb = 0.2
	p_i=1
	res_dir = "/Users/bjarni/tmp/"
	runId = "full_"
	
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	snpsdata.coordinateSnpsAndPhenotypeData(phed,p_i,snpsds)
	
	for snpsd in snpsds:
		snpsd.filterMinMAF(0.1)
		snpsd.filterMonoMorphicSnps()
		
		
	totalSNPs = []
	for i in range(len(snpsds)):
		snpsds[i] = snpsds[i].getSnpsData()
		totalSNPs += snpsds[i].snps
	
	#For memory, remove random SNPs
	snps = []
	for snp in totalSNPs:
		if random.random()<filterProb:
			snps.append(snp)
	totalSNPs = snps
	
	print "Calculating the global kinship..."
	globalKinship = calcKinship(totalSNPs)
	print "done."
	normalizedGlobalKinship = globalKinship/mean(globalKinship)
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..


	for i in range(4,5):#len(snpsds)):
		chr = i+1
		snpsd = snpsds[i]
		#pylab.subplot(5,1,chr)
#		pylab.figure(figsize=(18,4))
#		(kinshipDiffs,binPos,local300Kinships) = getKinshipDiffs(snpsd,normalizedGlobalKinship,windowSize=300000)
#		pylab.plot(binPos,kinshipDiffs,"r",label='ws$=300000$')
#		(kinshipDiffs,binPos,local500Kinships) = getKinshipDiffs(snpsd,normalizedGlobalKinship,windowSize=500000)
#		pylab.plot(binPos,kinshipDiffs,"b",label='ws$=500000$')
#		pylab.legend(numpoints=2,handlelen=0.005)
#		pylab.title("Kinship diff. chr. "+str(chr))
#		pylab.savefig(res_dir+runId+"kinshipDiffs_500_300kb_chr"+str(chr)+".pdf",format="pdf")
#		pylab.clf()
		pylab.figure(figsize=(18,4))
		(emmaDiffs,binPos) = getEmmaDiffs(snpsd,phed,p_i,globalKinship,windowSize=300000)
		pylab.plot(binPos,emmaDiffs,"r",label='ws$=300000$')
		pylab.title("Emma avg. p-value diff. 500kb on chr. "+str(chr))
		(emmaDiffs,binPos) = getEmmaDiffs(snpsd,phed,p_i,globalKinship,windowSize=500000)
		pylab.plot(binPos,emmaDiffs,"b",label='ws$=500000$')
		pylab.title("Emma avg. p-value diff. on chr. "+str(chr))
		pylab.legend(numpoints=2,handlelen=0.005)
		pylab.savefig(res_dir+runId+"EmmaPvalDiffs_500_300kb_chr"+str(chr)+".pdf",format="pdf")
		pylab.clf()
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..


		
def _plotKW_():
	"""
	Analyze how population structure affects KW.
	"""
	filterProb = 0.1
	p_i=1
	res_dir = "/Users/bjarni/tmp/"
	runId = "_full_quick_"
	
	
	snpsDataFile="/Network/Data/250k/dataFreeze_080608/250K_f10_080608.csv"
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	snpsdata.coordinateSnpsAndPhenotypeData(phed,p_i,snpsds)
	
	totalSNPs = []
	for i in range(len(snpsds)):
		snpsds[i] = snpsds[i].getSnpsData()
		totalSNPs += snpsds[i].snps
	
	#For memory, remove random SNPs
	snps = []
	for snp in totalSNPs:
		if random.random()<filterProb:
			snps.append(snp)
	totalSNPs = snps
	
	#globalKinship = calcKinship(totalSNPs)
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	#chr = 1
	#for snpsd in snpsds:
	
	snpsd = snpsds[3]


	k = calcKinship(snpsd.snps[200:1400])
	res = runEmma(phed,p_i,k,snpsd.snps[200:1400]) #runEmma(phed,p_i,k,snps):
	pvals = res["ps"]
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	pylab.plot(snpsd.positions[200:1400],log_pvals,"c.",label="Emma (local)")
	
	k = calcKinship(totalSNPs)
	res = runEmma(phed,p_i,k,snpsd.snps[200:1400]) #runEmma(phed,p_i,k,snps):
	pvals = res["ps"]
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	pylab.plot(snpsd.positions[200:1400],log_pvals,"g.",label="Emma (global)")
	
	phenVals = phed.getPhenVals(p_i)
	pvals = _run_kw_(snpsd.snps[200:1400],phenVals)
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	
	pylab.plot(snpsd.positions[200:1400],log_pvals,"r.",label="KW (full data)")
	
	(pvals,new_positions,acc_groups) = get_KW_pvals(snpsd.snps[200:1400],snpsd.positions[200:1400],phed,p_i,kinshipThreshold = 0.95, method = "KW")
	ecot_map = phenotypeData._getEcotypeIdToStockParentDict_()
	
	for i in range(0,len(acc_groups)):
		acc_list = []
		for a_i in acc_groups[i]:
			e_i = snpsd.accessions[a_i]
			#print e_i
			acc_list.append(ecot_map[int(e_i)][0])
		print "group",i,":",acc_list
	
	log_pvals = []
	for pval in pvals:
		#print pval
		log_pvals.append(-math.log10(pval))
	
	pylab.plot(new_positions,log_pvals,"b.",label="KW (merged data)")
	
	pylab.legend(numpoints=2,handlelen=0.005)
	
	pylab.show()
		
	
		
		
if __name__ == "__main__":
	_plotKinshipDiffs_()
	#_plotKW_()
