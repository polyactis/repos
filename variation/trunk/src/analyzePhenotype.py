"""
This class provides functionality (in python) to evaluate which transformation to choose.
"""
import phenotypeData, gwaResults, gc, plotResults
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#import pylab as plt



def drawHistogram(phed, p_i, title = None , pdfFile = None, pngFile = None):
	plt.figure(figsize=(10,7))
	p_i = phed.getPhenIndex(p_i)
	phenValues = []
	for values in phed.phenotypeValues:
		if values[p_i] != 'NA':
			phenValues.append(float(values[p_i]))
	minVal = min(phenValues)
	maxVal = max(phenValues)
	valRange = maxVal - minVal
	 
	histRes= plt.hist(phenValues, bins = len(phenValues)/3)
	if title:
		plt.title(title)
	if pdfFile:
		plt.savefig(pdfFile, format = "pdf")
	if pngFile:
		plt.savefig(pngFile, format = "png")
	elif not pdfFile:
		plt.show()
	plt.clf()
		

def _getQuantiles_(scores, numQuantiles):
	scores.sort()
	quantiles = []
	for i in range(1, numQuantiles + 1):
		j = int(len(scores) * i / (numQuantiles + 2))
		quantiles.append(scores[j])
	return quantiles
		
def __getExpectedPvalueQuantiles__(numQuantiles):
	quantiles = []
	for i in range(1, numQuantiles + 1):
		quantiles.append(float(i) / (numQuantiles + 2))
	return quantiles
		
def _calcMedian_(scores,exp_median=0.5):
        scores.sort()
	median = scores[len(scores)/2]
	return (exp_median-median)


def drawQQPlot(results, numQuantiles, phenName = None, pdfFile = None, pngFile = None, perm_pvalues=None, **kwargs):
	plt.figure(figsize=(10,8))
	plt.plot([0, 1], [0, 1],"k",label="Expected")
	areas = []
	medians = []
	for j in range(0, len(results)):
		result = results[j]
		label = kwargs['resultTypes'][j]
		label = " ".join(label.split("_"))
		newScores = result.scores[:]
		quantiles = _getQuantiles_(newScores, numQuantiles)
		if perm_pvalues and j==0: #j==0 is KW..
			print "Getting exp. quantiles for permuted p-values"
			expQuantiles = _getQuantiles_(perm_pvalues,numQuantiles)
			q_i = numQuantiles/2
			if numQuantiles%2==0: #even
				exp_median = (expQuantiles[q_i-1]+expQuantiles[q_i])/2.0
			else: #odd
				exp_median = expQuantiles[q_i]

		else:
			exp_median = 0.5
			expQuantiles = __getExpectedPvalueQuantiles__(numQuantiles)
		area = _estAreaBetweenCurves_(quantiles,expQuantiles)
		median = _calcMedian_(newScores,exp_median)
		plt.plot(quantiles,expQuantiles, label = label+", A="+str(round(area,4))+", M="+str(round(median,4)))
		areas.append(area)
		medians.append(median)
			
	if phenName:
		plt.title(phenName)
	fontProp = matplotlib.font_manager.FontProperties(size=10)
	plt.legend(loc = 2, numpoints = 4, handlelen = 0.05, markerscale = 0.5,prop=fontProp,pad=0.02)
	plt.axis([-0.01,1.01,-0.01,1.01])
	plt.ylabel("Expected $p$-value")
	plt.xlabel("Observed $p$-value")
	if pdfFile:
		plt.savefig(pdfFile, format = "pdf")
	if pngFile:
		plt.savefig(pngFile, format = "png")
	elif not pdfFile:
		plt.show()
	plt.clf()
	return (areas,medians)


def _estAreaBetweenCurves_(quantiles,expQuantiles):
	area = 0
	for i in range(0,len(quantiles)-1):
		area += (expQuantiles[i+1]-expQuantiles[i])*(abs(quantiles[i+1]-expQuantiles[i+1]+quantiles[i]-expQuantiles[i]))/2.0
	#area = area*(expQuantiles[1]-expQuantiles[0])
	return area
	
def _calcKS_(scores,exp_scores=None):
	ret = {}
	ret["D"] = -1
	try:
		from rpy import r
		if exp_scores:
			res = r.ks_test(scores,exp_scores)
		else:
			res = r.ks_test(scores,"punif")
		ret = res["statistic"]
		ret["p.value"] = res["p.value"]
	except Exception, message:
		print "Calculating KS failed??",message
	return ret

	
def _getLogQuantilesMaxVal_(scores,maxScore=None):
	scores.sort()
	i = 0
	new_score = -math.log(scores[i],10)
	score = new_score
	#print score, maxScore
	while i < len(scores)-1 and new_score > maxScore:
		score = new_score
		i += 1
		new_score = -math.log(scores[i],10)

	maxVal = math.log((len(scores))/float(i+1),10)
	#print maxVal,i, score, maxScore
	return(maxVal)

def _getLogQuantiles_(scores, numDots, maxVal=None):
	scores.sort()
	quantiles = []
	for i in range(0,numDots):
		j = int(round(math.pow(10,-(float(i)/(numDots-1))*maxVal)*len(scores)))
		quantiles.append(-math.log10(scores[j-1]))
	#print quantiles[(numDots-100):]
	return quantiles
		


def __getExpectedLogQuantiles__(numDots,maxVal):
	quantiles = []
	for i in range(1, numDots + 1):
		quantiles.append((float(i)/(numDots+2.0))*maxVal)
	return quantiles



def _estLogSlope_(ys,xs=None):
	if xs:
		q1 = _getQuantiles_(xs,1000)
	else:
		q1 = __getExpectedPvalueQuantiles__(1000)
	q2 = _getQuantiles_(ys,1000)
	b_sum = 0.0
	num_valid = 0
	for (x,y) in zip(q1,q2):
		if x<1.0:
			b_sum += math.log(y,10)/math.log(x,10)
			num_valid += 1
	return(b_sum/num_valid)


def drawLogQQPlot(results, numDots, maxVal, phenName = None, pdfFile = None, pngFile = None, perm_pvalues=None, **kwargs):

	#FIXME: FINSIH numdots and maxVal!!!
	plt.figure(figsize=(10,8))
	maxVal = min(math.log10(len(results[0].scores)), maxVal)
	minVal = (1.0/numDots)*maxVal
	valRange = maxVal-minVal
	plt.plot([minVal, maxVal], [minVal, maxVal], "k",label="Expected")
	maxObsVals = []
	areas = []
	ds = []
	slopes = []
	for j in range(0, len(results)):
		result = results[j]
		label = kwargs['resultTypes'][j]
		label = " ".join(label.split("_"))
		if perm_pvalues and j==0:#j==0 is KW..
			exp_maxVal = _getLogQuantilesMaxVal_(perm_pvalues[:],maxVal)
			expQuantiles = _getLogQuantiles_(perm_pvalues[:],numDots,exp_maxVal)
			ks_res = _calcKS_(result.scores,perm_pvalues)
			quantiles = _getLogQuantiles_(result.scores[:], numDots, exp_maxVal)
			slope = _estLogSlope_(result.scores[:],perm_pvalues)
		else:
			quantiles = _getLogQuantiles_(result.scores[:], numDots, maxVal)
			expQuantiles = __getExpectedLogQuantiles__(numDots,maxVal)
			ks_res = _calcKS_(result.scores)
			slope = _estLogSlope_(result.scores[:])

		area = _estAreaBetweenCurves_(quantiles,expQuantiles)
		areas.append(area)
		slopes.append(slope)
		ds.append(ks_res["D"])
		plt.plot(expQuantiles, quantiles, label = label+", A="+str(round(area,2))+", D="+str(round(ks_res["D"],3))+", S="+str(round(slope,3)))
		maxObsVals.append(max(quantiles))
		
	maxObsVal = max(maxObsVals)
	obsValRange = maxObsVal-minVal
	plt.axis([minVal-0.025*valRange,maxVal+0.025*valRange,minVal-0.025*obsValRange,maxObsVal+0.025*obsValRange])
	plt.ylabel("Observed $log_{10}(p)$")
	plt.xlabel("Expected $log_{10}(p)$")
	if phenName:
		plt.title(phenName)
	fontProp = matplotlib.font_manager.FontProperties(size=10)
	plt.legend(loc = 2, numpoints = 4, handlelen = 0.05, markerscale = 0.5,prop=fontProp,pad=0.02)
	if pdfFile:
		plt.savefig(pdfFile, format = "pdf")
	if pngFile:
		plt.savefig(pngFile, format = "png")
	elif not pdfFile:
		plt.show()
	plt.clf()
	return (ds,areas,slopes)


def drawPermLogQQPlot(results, permutedResultsList, numDots=1000, maxVal=5, phenName = None, pdfFile = None, pngFile=None, **kwargs):
		
	maxVal = min(math.log10(len(results[0].scores)), maxVal)
	minVal = (1.0/numDots)*maxVal
	valRange = maxVal-minVal
	plt.plot([minVal, maxVal], [minVal, maxVal], "k",label="Expected")
	maxObsVals = []
	for (permutedResults,label,color) in permutedResultsList:
		for result in permutedResults:
			quantiles = _getLogQuantiles_(result.scores[:], numDots, maxVal)
			expQuantiles = __getExpectedLogQuantiles__(numDots,maxVal)
			plt.plot(expQuantiles, quantiles,color=color)
			maxObsVals.append(max(quantiles))
		plt.plot(expQuantiles, quantiles,label=label,color=color)

	for result in results:
		ks_res = _calcKS_(result.scores)
		quantiles = _getLogQuantiles_(result.scores[:], numDots, maxVal)
		expQuantiles = __getExpectedLogQuantiles__(numDots,maxVal)
		area = _estAreaBetweenCurves_(quantiles,expQuantiles)
		plt.plot(expQuantiles, quantiles, label = result.name+", A="+str(round(area,2))+", D="+str(round(ks_res["D"],3)))
		maxObsVals.append(max(quantiles))
		
	maxObsVal = max(maxObsVals)
	obsValRange = maxObsVal-minVal
	plt.axis([minVal-0.025*valRange,maxVal+0.025*valRange,minVal-0.025*obsValRange,maxObsVal+0.025*obsValRange])
	if phenName:
		plt.title(phenName)
	fontProp = matplotlib.font_manager.FontProperties(size=10)
	plt.legend(loc = 2, numpoints = 2, handlelen = 0.01, markerscale = 0.5,prop=fontProp,pad=0.1)
	if pdfFile:
		plt.savefig(pdfFile, format = "pdf")
	if pngFile:
		plt.savefig(pngFile, format = "png")
	elif not pdfFile:
		plt.show()
	plt.clf()


def _drawPowerQQPlots_(phenotypeIndices=None,res_path="/Network/Data/250k/tmp-bvilhjal/power_analysis/results/",runId="gwPlot"):
	"""
	Draws all the GWA plots for 6 methods.
	"""
	import plotResults

	
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')

	if not phenotypeIndices:
		phenotypeIndices = phed.phenIds

	#mainRTs = ["","_original192", "_original192_inverse","_original96","_original96_inverse"]  #FIXME: add full result 
	#mainLabels = ["Full data", "192 acc. overlap", "192 acc. complement", "96 acc. overlap", "96 acc. complement"]
	mainRTs = ["","_original192","_original96", "_latitude60","_latitude55", "_original192_latitude60", "_original192_latitude55"]  #FIXME: add full result 
	mainLabels = ["Full data", "192 acc. overlap", "96 acc. overlap", "latitude < 60", "Latitude < 55", "192 acc. overl. and lat. < 60", "192 acc. overl. and lat. < 55"]
	permRTs = []#["permTest","permTest"]
	colors = []#[[0.6,0.8,0.6],[0.6,0.6,0.8]]
	perm_counts = []#[10,10]
	perm_sample_sizes = []#[65 ,112] #[170,96] #
	permLabels = []#["random 65","random 112"]
	for p_i	in phenotypeIndices:
		mainResults = []
		phenName = phed.getPhenotypeName(p_i)
		pdfFile = res_path+phenName+"_log_QQplot.pdf"
		pngFile = res_path+phenName+"_log_QQplot.png"
		for i in range(0,len(mainRTs)):
			mainRT = mainRTs[i]
			name = mainLabels[i]
			filename = res_path+"KW_raw"+mainRT+"_"+phenName+".pvals"
			rt = gwaResults.ResultType(resultType="KW",name=name)
			print "Loading",filename
			result = gwaResults.Result(filename, name=name, resultType=rt)
			mainResults.append(result)
	
		permResultsList = []
		for i in range(0,len(permRTs)):
			permResults = []
			permRT = permRTs[i]
			for j in range(0,perm_counts[i]):
				filename = res_path+"KW_raw_"+permRT+"_"+phenName+"_r"+str(perm_sample_sizes[i])+"_"+str(j)+".pvals"
				rt = gwaResults.ResultType(resultType="KW",name=permRT)
				print "Loading",filename
				result = gwaResults.Result(filename, name=permRT, resultType=rt)
				permResults.append(result)
			permResultsList.append((permResults,permLabels[i],colors[i]))

		drawPermLogQQPlot(mainResults, permResultsList, phenName = phenName,pdfFile=pdfFile,pngFile=pngFile)
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	#plotResults.plotResult(result,pdfFile,pngFile,ylab="- log10 pvalue",plotBonferroni=True)




def _test_():
	phenotypeFile = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_102208.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	for p_i in range(215, 216):
		phenName = phed.getPhenotypeName(p_i)
		drawHistogram(phed, p_i, title = phenName)
		phed.logTransform(p_i)
		drawHistogram(phed, p_i, title = phenName)
		
def _loadData_(phed, phenotypeIndices):

	res_path = "/Network/Data/250k/tmp-bvilhjal/"

	resultsDirs = [res_path + "kw_results/", res_path + "emma_results/"]
	methods = ["KW", "Emma",]
	datasetNames = ["new_raw", "new_trans"]
	mafCutoffs = [0, 15,]

	#mrIndex = 1  #The guiding (main) result

	results_map = {}
	resultTypes_map = {}
	for i in phenotypeIndices:
		try:
			phenName = phed.getPhenotypeName(i)
			phenName = phenName.replace("/", "_div_")
			phenName = phenName.replace("*", "_star_")
			phenIndex = phed.getPhenIndex(i)
			
			results = []
			resultTypes = []
			for j in range(0, len(methods)):
				#if not (methods[j]=="Emma" and phed.isBinary(i)):
					resultFile = resultsDirs[j] + methods[j] + "_" + datasetNames[j] + "_" + phenName + ".pvals"
					try:
						print "Loading result file", resultFile
						result = gwaResults.Result(resultFile, name = methods[j] + "_" + datasetNames[j] + "_" + phenName)
						result.filterMAF(minMaf = mafCutoffs[j])
						results.append(result)
						resultTypes.append(methods[j] + "_" + datasetNames[j])
					except Exception:
						print "Couldn't load", resultFile
									
			results_map[i] = results
			resultTypes_map[i] = resultTypes
		except Exception:
			print "Couldn't load the result file"			
	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	return (results_map, resultTypes_map)


def _getPermPvalues_(phenName,perm_pval_dir="/Network/Data/250k/tmp-bvilhjal/perm_tests/"):
	perm_pvals = []
	filename = perm_pval_dir+"KW_perm_f1_n1000_"+phenName+".perm.pvals"  #FIXME finish 
	print "Getting permuted p-values:",filename
	f = open(filename,"r")
	lines = f.readlines()
	for line in lines:
		pval_str_lst = line.split(",")
		pvals = map(float,pval_str_lst)
		for pval in pvals:
			perm_pvals.append(pval)
	return perm_pvals
	
def _testQQplot_(includeEmmaInBinary=False,usePvalueFiles=True):
	resdir = "/Users/bjarni/tmp/"
	#resdir = "/Network/Data/250k/tmp-bvilhjal/phenotype_analyzis/"
	#resdir = "/Network/Data/250k/tmp-bvilhjal/qq_plots/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phed2 = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phenotypeIndices = phenotypeData.categories_2_phenotypes[4]#+phenotypeData.categories_2_phenotypes[2]+phenotypeData.categories_2_phenotypes[3]+phenotypeData.categories_2_phenotypes[4]
	#(results_map, resultTypes_map) = _loadData_(phed, phenotypeIndices)
	q_pvalues = None
	stat_dict = {}
	for p_i in phenotypeIndices:
		(results_map, resultTypes_map) = _loadData_(phed, [p_i])
		#try:
		phenName = phed.getPhenotypeName(p_i)
		phenNamePrint = " ".join(phenName.split("_")[1:])
		print "\nWorking on phenotype",phenName
		if usePvalueFiles:
			q_pvalues = _getPermPvalues_(phenName)
			print len(q_pvalues),"permuted pvalues found"

		valCount = phed.countValues(p_i)
		print valCount,"values found."
		if (not phed.isBinary(p_i)) or includeEmmaInBinary:
			histogramFile = resdir + phenName +"_hist.pdf"
			histogramFile_png = resdir + phenName +"_hist.png"
			drawHistogram(phed, p_i, title = phenNamePrint, pngFile = histogramFile_png)
			if phed.logTransform(p_i):
				histogramFile = resdir + phenName + "_hist_logTransformed.pdf"
				histogramFile_png = resdir + phenName + "_hist_logTransformed.png"
				drawHistogram(phed, p_i, title = phenNamePrint, pngFile = histogramFile_png)
			elif not phed.isBinary(p_i):
				print "adding scaled const."
				phed.addSDscaledConstant(p_i)
				if phed.logTransform(p_i):
					histogramFile = resdir + phenName + "_hist_logTransformed_const.pdf"
					histogramFile_png = resdir + phenName + "_hist_logTransformed_const.png"
					drawHistogram(phed, p_i, title = phenNamePrint, pngFile = histogramFile_png)

#				phed2.naOutliers(p_i,10)
#				histogramFile = resdir + phenName + "_hist_noOutliers.pdf"
#				histogramFile_png = resdir + phenName + "_hist_noOutliers.png"
#				drawHistogram(phed2, p_i, title = phenName, pdfFile = histogramFile, pngFile = histogramFile_png)
#				if phed2.logTransform(p_i):
#					histogramFile = resdir + phenName + "_hist_logTransformed_noOutliers.pdf"
#					histogramFile_png = resdir + phenName + "_hist_logTransformed_noOutliers.png"
#					drawHistogram(phed2, p_i, title = phenName, pdfFile = histogramFile, pngFile = histogramFile_png)
		results = results_map[p_i]
		resultTypes = resultTypes_map[p_i]
		qqplotFile = resdir + phenName + "_qqplot.pdf"
		qqplotFile_png = resdir + phenName + "_qqplot.png"
		s_dict={}
		(As,Ms)=drawQQPlot(results, 1000, phenName = phenNamePrint, resultTypes = resultTypes, pngFile=qqplotFile_png, perm_pvalues = q_pvalues)
		s_dict["A"]=As
		s_dict["M"]=Ms
		
		qqplotFile = resdir + phenName + "_qqplot_log.pdf"
		qqplotFile_png = resdir + phenName + "_qqplot_log.png"
		(ds,areas,slopes) = drawLogQQPlot(results, 1000,5, phenName = phenNamePrint, resultTypes = resultTypes, pngFile=qqplotFile_png, perm_pvalues = q_pvalues)
		s_dict["A2"]=areas
		s_dict["D"]=ds
		s_dict["S"]=slopes
		stat_dict[p_i] = s_dict
		for i in range(0,len(results)):
			result = results[i]
			result.negLogTransform()
			pngFile = resdir + phenName + "_gwplot_" +resultTypes[i]+".png"
			plotResults.plotResult(result,pngFile=pngFile,percentile=90,type="pvals", plotBonferroni=True)	
		#except Exception:
		#	print "\nPhenotype index", p_i, "failed."
		del results_map
	       	gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
		
	print stat_dict
	stat_file_name = resdir + "confounding_stat_4.txt"
	f = open(stat_file_name,"w")
	methods = ["KW","Emma"]
	f.write("phenotype_name, method_name, is_binary, D, A, B, M, S\n")
	for p_i in phenotypeIndices:
		if stat_dict.has_key(p_i):
			s_dict = stat_dict[p_i]
			phenName = phed.getPhenotypeName(p_i)
			phenName = " ".join(phenName.split("_")[1:])
			for i in range(0,len(methods)):
				st = phenName+", "+methods[i]+", "+str(phed.isBinary(p_i))+", "+str(s_dict["D"][i])+", "+str(s_dict["A"][i])+", "+str(s_dict["A2"][i])+", "+str(s_dict["M"][i])+", "+str(s_dict["S"][i])+"\n"
				f.write(st)
	f.close()
			

def _countVals_():
	resdir = "/Network/Data/250k/tmp-bvilhjal/phenotype_analyzis/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phenotypeIndices = phenotypeData.categories_2_phenotypes[1]+phenotypeData.categories_2_phenotypes[2]+phenotypeData.categories_2_phenotypes[3]+phenotypeData.categories_2_phenotypes[4]
	print "total # of phenotypes:", phed.countPhenotypes()
	print "# of phenotypes analyzed:", len(phenotypeIndices)
	
	totalCounts = []
	for p_i in phenotypeIndices:
		valCount = phed.countValues(p_i)
		totalCounts.append(valCount)

	snpsDataFile="/Network/Data/250k/dataFreeze_011209/250K_f13_012509.csv"
	import dataParsers,snpsdata
	snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")#,debug=True)
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
	phed.removeAccessionsNotInSNPsData(snpsd)
	
	overlappingCounts = []
	for p_i in phenotypeIndices:
		valCount = phed.countValues(p_i)
		overlappingCounts.append(valCount)


	#ecotypes_192 = phenotypeData._getFirst192Ecotypes_()
	ecotypes_192 = _get192Ecotypes_()
	ecotypes_192 = [str(e) for e in ecotypes_192]
	print "len(ecotypes_192):",len(ecotypes_192)
	print ecotypes_192
	phed.filterAccessions(ecotypes_192)

	filename = resdir+"phen_value_count_new_data_012509_v2.txt"
	f = open(filename,"w")
	f.write("Phenotype,  total_count, overlapping_count, 192_overlap_count\n")
	
	for i in range(0,len(phenotypeIndices)):
		p_i = phenotypeIndices[i]
		try:
			phenName = phed.getPhenotypeName(p_i)
			valCount = phed.countValues(p_i)
			f.write(str(phenName)+", "+str(totalCounts[i])+", "+str(overlappingCounts[i])+", "+str(valCount)+"\n")
		except Exception:
			print "\nPhenotype index", p_i, "failed."

	f.close()


def _get192Ecotypes_():
	resdir = "/Network/Data/250k/tmp-bvilhjal/phenotype_analyzis/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phenotypeIndices = phenotypeData.categories_2_phenotypes[1]+phenotypeData.categories_2_phenotypes[2]+phenotypeData.categories_2_phenotypes[3]+phenotypeData.categories_2_phenotypes[4]
	
	
	total_accessions = set()
	for p_i in phenotypeIndices:
		if not p_i in [5,6,7]:
			accessions = phed.getAccessionsWithValues(p_i)
			total_accessions = total_accessions.union(accessions)

	ecotypes_192 = phenotypeData._getFirst192Ecotypes_()
	ecotypes_192 = [str(e) for e in ecotypes_192]
	print "len(ecotypes_192):",len(ecotypes_192)
	#print ecotypes_192
	phed.filterAccessions(ecotypes_192)

        for p_i in [5,6,7]:
		accessions = phed.getAccessionsWithValues(p_i)
		total_accessions = total_accessions.union(accessions)
		
	total_accessions = list(total_accessions)
	print len(total_accessions)
	total_accessions.sort()
	print total_accessions
	
	ecotype_info_dict = phenotypeData._getEcotypeIdInfoDict_()
	ets = []
	
	i = 0
	for et in total_accessions:
		et = int(et)
		if ecotype_info_dict.has_key(et):
			print str(et)+", "+str(ecotype_info_dict[et][0])+", "+str(ecotype_info_dict[et][1])
			i += 1
			ets.append(et)
		else:
			print et,"is missing in genotype data."
	print i
	return ets

def _plotConfoundingStats_():
	#import pylab as plt

	resdir = "/Network/Data/250k/tmp-bvilhjal/perm_tests/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	stat_file_dir = "/Users/bjarni/tmp/"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phenotypeIndices = phenotypeData.categories_2_phenotypes[1]+phenotypeData.categories_2_phenotypes[2]+phenotypeData.categories_2_phenotypes[3]+phenotypeData.categories_2_phenotypes[4]
	
	m_pvals = {}
	a_pvals = {}
	ks_pvals = {}
	for p_i in phenotypeIndices:
		#if not phed.isBinary(p_i):
		phenName = phed.getPhenotypeName(p_i)
		print "Loading permutation stats data for",phenName
		filename = resdir+"KW_perm_f1_n1000_"+phenName+".perm.stat.txt"
		f = open(filename,"r")
		lines = f.readlines()
		pvals = (lines[-1].strip()).split(',')
		m_pvals[p_i] = float(pvals[0].split(" ")[-1])
		a_pvals[p_i] = float(pvals[1])
		ks_pvals[p_i] = float(pvals[2])
	

	x_ticks = []
	s_ticks = []
	x_pos = 1
	for cat in [1,2,3,4]:
		for p_i in phenotypeData.categories_2_phenotypes[cat]:
			s_ticks.append(phed.getPhenotypeName(p_i))
			#plt.text(x_pos+shift,min_stat-0.1*stat_range,p_i,rotation="vertical",size="xx-small")				
			x_ticks.append(x_pos-0.5)
			x_pos += 1
		x_pos = x_pos+1

	
	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.16,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0
	colors = {1:"b",2:"r",3:"g",4:"c"}
	for i in [1,2,3,4]:
		
		phenotypeIndices = phenotypeData.categories_2_phenotypes[i]
		newPhenotypeIndices = []
		for p_i in phenotypeIndices:
			#if not phed.isBinary(p_i):
				newPhenotypeIndices.append(p_i)
		phenotypeIndices = newPhenotypeIndices
			
		m_list = []
		for p_i in phenotypeIndices:
			m_list.append(m_pvals[p_i])
		plt.bar(range(x_pos,len(m_list)+x_pos),m_list,color = colors[i])
		x_pos = x_pos+len(m_list)+1
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),-0.02,1.02])
	plt.xticks(x_ticks,s_ticks,size="x-small",rotation="vertical")
	plt.ylabel("M stat. p-value")
	plt.savefig(stat_file_dir+"confounding_M_pvalues.png", format = "png")
	plt.clf()

	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.16,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0
	for i in [1,2,3,4]:
		
		phenotypeIndices = phenotypeData.categories_2_phenotypes[i]
		newPhenotypeIndices = []
		for p_i in phenotypeIndices:
			#if not phed.isBinary(p_i):
				newPhenotypeIndices.append(p_i)
		phenotypeIndices = newPhenotypeIndices
		a_list = []
		for p_i in phenotypeIndices:
			a_list.append(a_pvals[p_i])
		plt.bar(range(x_pos,len(a_list)+x_pos),a_list,color = colors[i])
		x_pos = x_pos+len(a_list)+1
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),-0.02,1.02])
	plt.xticks(x_ticks,s_ticks,size="x-small",rotation="vertical")
	plt.ylabel("A stat. p-value")
	plt.savefig(stat_file_dir+"confounding_A_pvalues.png", format = "png")
	plt.clf()

	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.16,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0
	for i in [1,2,3,4]:
		
		phenotypeIndices = phenotypeData.categories_2_phenotypes[i]
		newPhenotypeIndices = []
		for p_i in phenotypeIndices:
			#if not phed.isBinary(p_i):
				newPhenotypeIndices.append(p_i)
		phenotypeIndices = newPhenotypeIndices
		a_list = []
		for p_i in phenotypeIndices:
			a_list.append(ks_pvals[p_i])
		plt.bar(range(x_pos,len(a_list)+x_pos),a_list,color = colors[i])
		x_pos = x_pos+len(a_list)+1
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),-0.02,1.02])
	plt.xticks(x_ticks,s_ticks,size="x-small",rotation="vertical")
	plt.ylabel("KS stat. p-value")
	plt.savefig(stat_file_dir+"confounding_KS_pvalues.png", format = "png")
	plt.clf()


	print m_pvals, a_pvals, ks_pvals
	

def _plotConfoundingStats2_():
	resdir = "/Network/Data/250k/tmp-bvilhjal/perm_tests/"
	stat_file_dir = "/Users/bjarni/tmp/"
	
	stats_dict = {1:{},2:{},3:{},4:{}}
	for cat in [1,2,3,4]:
		f = open(stat_file_dir+"confounding_stat_"+str(cat)+".txt","r")
		lines = f.readlines()
		p_dict = {}
		m_dict = {}
		for line in lines[1:]:
			p_name = line[0].strip()
			line = line.strip()
			line = line.split(",")
			#"phenotype_name, method_name, is_binary, D, A, B, M, S\n"
			method = line[1].strip()
			
			m_dict[method] = {}
			m_dict[method]["is_binary"] = (line[2].strip()=="True")
			m_dict[method]["D"]= float(line[3])
			m_dict[method]["A"]= float(line[4])
			m_dict[method]["B"]= float(line[5])
			m_dict[method]["M"]= float(line[6])
			m_dict[method]["S"]= float(line[7])-1.0
			
			p_dict[line[0].strip()] = m_dict.copy()
		stats_dict[cat] = p_dict
	print stats_dict[1]
	print stats_dict[2]
	print stats_dict[3]
	print stats_dict[4]
	
			
	#import pylab as plt
	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.15,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0

	colors = {1:"b",2:"r",3:"g",4:"y"}
	max_stats = []
	min_stats = []
	for cat in [1,2,3,4]:
		
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		s_list = {0:[],1:[]}
		for p_i in phenotypes:
			if not phenotypes_dict[p_i]["KW"]["is_binary"]:
				s_list[0].append(phenotypes_dict[p_i]["KW"]["B"])
				s_list[1].append(phenotypes_dict[p_i]["Emma"]["B"])
		max_stats.append(max(max(s_list[0]),max(s_list[1])))
		min_stats.append(min(min(s_list[0]),min(s_list[1])))
		for method in [0,1]:
			plt.bar(range(method+x_pos,method+2*len(s_list[method])+x_pos,2),s_list[method],color = colors[1+2*((cat-1)%2)+(method)%2])
		x_pos = x_pos+2*len(s_list[method])+1
	max_stat = max(max_stats)
	min_stat = min(0,min(min_stats))
	print min_stat
	stat_range = max_stat-min_stat
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),min_stat-stat_range*0.02,max_stat+stat_range*0.02])

	x_ticks = []
	s_ticks = []
	x_pos = 1
	for cat in [1,2,3,4]:
		shift = 0
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		for p_i in phenotypes:
			if not phenotypes_dict[p_i]["KW"]["is_binary"]:
				s_ticks.append(p_i)
				#plt.text(x_pos+shift,min_stat-0.1*stat_range,p_i,rotation="vertical",size="xx-small")				
				x_ticks.append(x_pos+shift)
				shift += 2
		x_pos = x_pos+shift+1
	
	plt.xticks(x_ticks,s_ticks,size="x-small",rotation="vertical")
	plt.ylabel("Area between expected and observed $log(p)$ curve.")
	#plt.ylabel("Area between expected and observed p-value curve.")
	#plt.ylabel("Kolmogorov-Smirnov $D$-statistic")
	#plt.ylabel("(est. slope of $log(p)$)$- 1$")
	#plt.ylabel("(median p-value)$- 0.5$")
	
	#plt.show()
	plt.savefig(stat_file_dir+"confounding_B.png", format = "png")
	plt.clf()
		

def _plotConfoundingStats3_():
	resdir = "/Network/Data/250k/tmp-bvilhjal/perm_tests/"
	stat_file_dir = "/Users/bjarni/tmp/"
	
	stats_dict = {1:{},2:{},3:{},4:{}}
	for cat in [1,2,3,4]:
		f = open(stat_file_dir+"confounding_stat_"+str(cat)+".txt","r")
		lines = f.readlines()
		p_dict = {}
		m_dict = {}
		for line in lines[1:]:
			p_name = line[0].strip()
			line = line.strip()
			line = line.split(",")
			#"phenotype_name, method_name, is_binary, D, A, B, M, S\n"
			method = line[1].strip()
			
			m_dict[method] = {}
			m_dict[method]["is_binary"] = (line[2].strip()=="True")
			m_dict[method]["D"]= float(line[3])
			m_dict[method]["A"]= float(line[4])
			m_dict[method]["B"]= float(line[5])
			m_dict[method]["M"]= float(line[6])
			m_dict[method]["S"]= float(line[7])-1.0
			
			p_dict[line[0].strip()] = m_dict.copy()
		stats_dict[cat] = p_dict
	print stats_dict[1]
	print stats_dict[2]
	print stats_dict[3]
	print stats_dict[4]
	
			
	#import pylab as plt
	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.15,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0

	colors = {1:"b",2:"r",3:"g",4:"y"}
	max_stats = []
	min_stats = []
	for cat in [1,2,3,4]:
		
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		s_list = {0:[],1:[]}
		for p_i in phenotypes:
			if phenotypes_dict[p_i]["KW"]["is_binary"]:
				s_list[0].append(phenotypes_dict[p_i]["KW"]["M"])
				s_list[1].append(phenotypes_dict[p_i]["Emma"]["M"])
		if len(s_list[0]):
			max_stats.append(max(max(s_list[0]),max(s_list[1])))
			min_stats.append(min(min(s_list[0]),min(s_list[1])))
			for method in [0,1]:
				plt.bar(range(method+x_pos,method+2*len(s_list[method])+x_pos,2),s_list[method],color = colors[1+2*((cat-1)%2)+(method)%2])
			x_pos = x_pos+2*len(s_list[method])+1
	max_stat = max(max_stats)
	min_stat = min(0,min(min_stats))
	print min_stat
	stat_range = max_stat-min_stat
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),min_stat-stat_range*0.02,max_stat+stat_range*0.02])

	x_ticks = []
	s_ticks = []
	x_pos = 1
	for cat in [1,2,3,4]:
		shift = 0
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		for p_i in phenotypes:
			if phenotypes_dict[p_i]["KW"]["is_binary"]:
				s_ticks.append(p_i)
				#plt.text(x_pos+shift,min_stat-0.1*stat_range,p_i,rotation="vertical",size="xx-small")				
				x_ticks.append(x_pos+shift)
				shift += 2
				
		if shift:
			x_pos = x_pos+shift+1
	
	plt.xticks(x_ticks,s_ticks,size="x-small",rotation="vertical")
	#plt.ylabel("Area between expected and observed $log(p)$ curve.")
	#plt.ylabel("Area between expected and observed p-value curve.")
	#plt.ylabel("Kolmogorov-Smirnov $D$-statistic")
	#plt.ylabel("(est. slope of $log(p)$)$- 1$")
	plt.ylabel("(median p-value)$- 0.5$")
	
	#plt.show()
	plt.savefig(stat_file_dir+"confounding_binary_M.png", format = "png")
	plt.clf()
		


def _plotConfoundingStats4_():
	resdir = "/Network/Data/250k/tmp-bvilhjal/perm_tests/"
	stat_file_dir = "/Users/bjarni/tmp/"
	
	stats_dict = {1:{},2:{},3:{},4:{}}
	for cat in [1,2,3,4]:
		f = open(stat_file_dir+"confounding_stat_"+str(cat)+".txt","r")
		lines = f.readlines()
		p_dict = {}
		m_dict = {}
		for line in lines[1:]:
			p_name = line[0].strip()
			line = line.strip()
			line = line.split(",")
			#"phenotype_name, method_name, is_binary, D, A, B, M, S\n"
			method = line[1].strip()
			
			m_dict[method] = {}
			m_dict[method]["is_binary"] = (line[2].strip()=="True")
			m_dict[method]["D"]= float(line[3])
			m_dict[method]["A"]= float(line[4])
			m_dict[method]["B"]= float(line[5])
			m_dict[method]["M"]= float(line[6])
			m_dict[method]["S"]= float(line[7])-1.0
			
			p_dict[line[0].strip()] = m_dict.copy()
		stats_dict[cat] = p_dict
	print stats_dict[1]
	print stats_dict[2]
	print stats_dict[3]
	print stats_dict[4]
	
			
	#import pylab as plt
	figure = plt.figure(figsize=(14,8))
	axes = plt.Axes(figure, [.06,.15,.91,.81])
	figure.add_axes(axes) 
	x_pos = 0

	colors = {1:"b",2:"r",3:"g",4:"y"}
	max_stats = []
	min_stats = []
	for cat in [1,2,3,4]:
		
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		s_list = {0:[],1:[]}
		for p_i in phenotypes:
			s_list[0].append(phenotypes_dict[p_i]["KW"]["B"])
			s_list[1].append(phenotypes_dict[p_i]["Emma"]["B"])
		max_stats.append(max(s_list[0]))
		min_stats.append(min(s_list[0]))
		plt.bar(range(x_pos,len(s_list[0])+x_pos),s_list[0],color = colors[1+2*((cat-1)%2)])
		x_pos = x_pos+len(s_list[0])+1
	max_stat = max(max_stats)
	min_stat = min(0,min(min_stats))
	print min_stat
	stat_range = max_stat-min_stat
	plt.axis([0-0.02*(x_pos-1),1.02*(x_pos-1),min_stat-stat_range*0.02,max_stat+stat_range*0.02])

	x_ticks = []
	s_ticks = []
	x_pos = 1
	for cat in [1,2,3,4]:
		shift = 0
		phenotypes_dict = stats_dict[cat]
		phenotypes = phenotypes_dict.keys()
		phenotypes.sort()
		for p_i in phenotypes:
			s_ticks.append(p_i)
				#plt.text(x_pos+shift,min_stat-0.1*stat_range,p_i,rotation="vertical",size="xx-small")				
			x_ticks.append(x_pos+shift-0.5)
			shift += 1
		x_pos = x_pos+shift+1
	
	plt.xticks(x_ticks,s_ticks,size="xx-small",rotation="vertical")
	plt.ylabel("Area between expected and observed $log(p)$ curve.")
	#plt.ylabel("Area between expected and observed p-value curve.")
	#plt.ylabel("Kolmogorov-Smirnov $D$-statistic")
	#plt.ylabel("(est. slope of $log(p)$)$- 1$")
	#plt.ylabel("(median p-value)$- 0.5$")
	
	#plt.show()
	plt.savefig(stat_file_dir+"confounding_KWonly_B.png", format = "png")
	plt.clf()
		

def _plotRobustnessTests_():
	import csv
	resdir = "/Network/Data/250k/tmp-bvilhjal/robustness_test/"
	fig_dir = "/Users/bjarni/tmp/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	print "Loading phenotype data"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	phenotypeIndices = phenotypeData.categories_2_phenotypes[1]+phenotypeData.categories_2_phenotypes[2]+phenotypeData.categories_2_phenotypes[3]+phenotypeData.categories_2_phenotypes[4]
	
	#First Emma
	emma_sd_dict = {}
	emma_sd_list = []
	emma_log_pvalues = []
	found_phenotypes = []
	for p_i in phenotypeIndices:
		try:
			phenName = phed.getPhenotypeName(p_i)
			filename = resdir+"KW_rob_f1_"+phenName+".rob.log_pvals_sd"
			print "Loading", filename, "..."
			reader = csv.reader(open(filename, "rb"))
			reader.next()
			for row in reader:
				emma_log_pvalues.append(float(row[0]))
				emma_sd_list.append(float(row[1]))
			found_phenotypes.append(p_i)
		except Exception:
			print p_i,"failed."
	
	import numpy as np
	import matplotlib.cm as cm
	import  matplotlib.pyplot as plt

	xs = np.array(emma_log_pvalues)
	ys = np.array(emma_sd_list)
	print len(emma_sd_list),len(emma_log_pvalues)
	xmin = xs.min()
	xmax = xs.max()
	ymin = ys.min()
	ymax = ys.max()
	
	#plt.subplots_adjust(hspace=0.5)
	#plt.subplot(121)
	plt.hexbin(xs,ys,bins='log',cmap=cm.jet)
	plt.axis([xmin, xmax, ymin, ymax])
	cb = plt.colorbar()
	cb.set_label('$log_{10}(N)$')
	plt.ylabel("SD$(\Delta log(p)))$")
	plt.xlabel("$log(p)$")
	plt.savefig(fig_dir+"KW_overall_robustness.png", format = "png")
	plt.clf()

	emma_sd_dict = {}
	emma_sd_list = []
	emma_log_pvalues = []
	found_phenotypes = []
	for p_i in phenotypeIndices:
		try:
			phenName = phed.getPhenotypeName(p_i)
			filename = resdir+"Emma_rob_f1_"+phenName+".rob.log_pvals_sd"
			print "Loading", filename, "..."
			reader = csv.reader(open(filename, "rb"))
			reader.next()
			for row in reader:
				emma_log_pvalues.append(float(row[0]))
				emma_sd_list.append(float(row[1]))
			found_phenotypes.append(p_i)
		except Exception:
			print p_i,"failed."
	
	import numpy as np
	import matplotlib.cm as cm
	import  matplotlib.pyplot as plt

	xs = np.array(emma_log_pvalues)
	ys = np.array(emma_sd_list)
	print len(emma_sd_list),len(emma_log_pvalues)
	xmin = xs.min()
	xmax = xs.max()
	ymin = ys.min()
	ymax = ys.max()
	
	#plt.subplots_adjust(hspace=0.5)
	#plt.subplot(121)
	plt.hexbin(xs,ys,bins='log',cmap=cm.jet)
	plt.axis([xmin, xmax, ymin, ymax])
	cb = plt.colorbar()
	cb.set_label('$log_{10}(N)$')
	plt.ylabel("SD$(\Delta log(p)))$")
	plt.xlabel("$log(p)$")
	plt.savefig(fig_dir+"Emma_overall_robustness.png", format = "png")
	plt.clf()




if __name__ == '__main__':
	#_get192Ecotypes_()
	#_countVals_()
	#_testQQplot_(includeEmmaInBinary=True)
	#_drawPowerQQPlots_([5])
	#_plotConfoundingStats_()
	_plotRobustnessTests_()
	#_test_()
	print "Done!"
