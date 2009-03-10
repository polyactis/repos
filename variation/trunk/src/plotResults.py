"""
Functions for ploting GWA results in pylab.
"""

import gwaResults

def plotFilteredResult(result,pdfFile,minScore=0,maxScore=10, plotBonferroni=False,usePylab=True):
	if usePylab:
		import pylab as plt
	else:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt

	newScores = []
	for score in result.scores:
		if score>maxScore:
			score = maxScore
		newScores.append(score)
	newSecondaryScores = []
	for score in result.secondaryScores:
		if score>maxScore:
			score = maxScore
		newSecondaryScores.append(score)
	minScore = min([min(newScores),min(newSecondaryScores),minScore])
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	secondaryChrSplits = result.getSecondaryChromosomeSplit()
	print chromosomeSplits
	print secondaryChrSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(1,figsize=(20,4))
	for i in range(0,len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i+1][0]
		secondaryIndex1 = secondaryChrSplits[i][0]
		secondaryIndex2 = secondaryChrSplits[i+1][0]
		scoreList = newScores[index1:index2]
		posList = result.positions[index1:index2]
		secScoreList = newSecondaryScores[secondaryIndex1:secondaryIndex2]
		secPosList = result.secondaryPositions[secondaryIndex1:secondaryIndex2]
		
		newPosList = []
		for pos in posList:
			newPosList.append(offset+pos)
		
		newSecPosList = []
		for pos in secPosList:
			newSecPosList.append(offset+pos)

		plt.plot(newSecPosList,secScoreList,".",color=(0.75,0.75,0.75,0.2))
		plt.plot(newPosList,scoreList,".")
		oldOffset = offset
		textPos.append(offset+result.chromosomeEnds[i]/2-2000000)
		offset += result.chromosomeEnds[i]
		if i<4:
			plt.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
		for j in range(oldOffset,offset,1000000):
			ticksList1.append(j)
		for j in range(0,result.chromosomeEnds[i],1000000):
			if j%5000000 == 0 and j < result.chromosomeEnds[i]-1500000 :
				ticksList2.append(j/1000000)
			else:
				ticksList2.append("")
				
		
	if plotBonferroni:
		plt.plot([0,sum(result.chromosomeEnds)],[6.68,6.68],"k--")


	plt.axis([0,sum(result.chromosomeEnds),minScore-0.05*scoreRange,maxScore+0.05*scoreRange])
	plt.xticks(ticksList1,ticksList2)
	for i in range(0,len(textPos)):
		plt.text(textPos[i],minScore-scoreRange*0.2,"Chr "+str(i+1))
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(left=0.03)
	plt.subplots_adjust(bottom=0.15)
	plt.subplots_adjust(top=0.9)
	plt.text(offset/2,maxScore+scoreRange*0.1,'Position')
	plt.ylabel('-log(p-value)')
	plt.savefig(pdfFile,format="pdf")
		
	#   result.chromosomeEnds[i]

	"""
	plt.title(statistic.name)
	plt.xlabel('Observed theta')
	plt.ylabel('Estimated theta')
	plt.clf()
	"""

def plotResult(result,pdfFile=None,pngFile=None,minScore=None,maxScore=None,percentile=95,type="pvals",ylab="$-$log$_{10}(p)$", plotBonferroni=False,usePylab=False):
	if usePylab:
		import pylab as plt
	else:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
	
	"""
	type is either 'pvals' or 'scores'.
	"""

	result.filterPercentile(percentile/100.0)

	if maxScore:
		newScores = []
		for score in result.scores:
			if score>maxScore:
				score = maxScore
			newScores.append(score)
	else:
		newScores = result.scores
		maxScore = max(newScores)
	if minScore:
	   minScore = min([min(newScores),minScore])
	else: 
		if type=="pvals":
			minScore = 0
		else:
			minScore = min(newScores)
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	print chromosomeSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(figsize=(20,4))
	for i in range(0,len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i+1][0]
		scoreList = newScores[index1:index2]
		posList = result.positions[index1:index2]
		
		newPosList = []
		for pos in posList:
			newPosList.append(offset+pos)
		

		plt.plot(newPosList,scoreList,".")
		oldOffset = offset
		textPos.append(offset+result.chromosomeEnds[i]/2-2000000)
		offset += result.chromosomeEnds[i]
		if i<4:
			plt.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
		for j in range(oldOffset,offset,1000000):
			ticksList1.append(j)
		for j in range(0,result.chromosomeEnds[i],1000000):
			if j%5000000 == 0 and j < result.chromosomeEnds[i]-1500000 :
				ticksList2.append(j/1000000)
			else:
				ticksList2.append("")
				
		
	if plotBonferroni:
		plt.plot([0,sum(result.chromosomeEnds)],[6.68,6.68],"k-.")

	plt.axis([0,sum(result.chromosomeEnds),minScore-0.05*scoreRange,maxScore+0.05*scoreRange])
	plt.xticks(ticksList1,ticksList2)
	for i in range(0,len(textPos)):
		plt.text(textPos[i],minScore-scoreRange*0.2,"Chr "+str(i+1))
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(left=0.05)
	plt.subplots_adjust(bottom=0.15)
	plt.subplots_adjust(top=0.9)
	plt.text(offset/2,maxScore+scoreRange*0.1,'Position')
	if not ylab:
		if type=="pvals":
			plt.ylabel('-log(p-value)')
		else:
			plt.ylabel('score')
	else:
		plt.ylabel(ylab)

	if pdfFile:
		plt.savefig(pdfFile,format="pdf")
	if pngFile:
		plt.savefig(pngFile,format="png")		
	if not (pdfFile or pngFile):
		plt.show()
		
	plt.clf()

	"""
	plt.title(statistic.name)
	plt.xlabel('Observed theta')
	plt.ylabel('Estimated theta')
	"""
	
def plotResultWithSecondRun(result,secondRunResult,pdfFile=None,pngFile=None,minScore=None,
						maxScore=None,percentile=90,srPercentile=90,type="pvals",ylab=None, 
						plotBonferroni=False,bonferroniCutoffs=(6.64,9.02),usePylab=False):
	if usePylab:
		import pylab as plt
	else:
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt


	result.filterPercentile(percentile/100.0)
	secondRunResult.filterPercentile(srPercentile/100.0)
	#print len(secondRunResult.scores)

	if maxScore:
		newScores = []
		newSecondRunScores = []
		for score in result.scores:
			if score>maxScore:
				score = maxScore
			newScores.append(score)
		for score in secondRunResult.scores:
			if score>maxScore:
				score = maxScore
			newSecondRunScores.append(score)
	else:
		newScores = result.scores
		newSecondRunScores = secondRunResult.scores
		maxScore = max(max(newScores),max(newSecondRunScores))
	if minScore:
	   minScore = min([min(newScores),minScore,min(newSecondRunScores)])
	else: 
		if type=="pvals":
			minScore = 0
		else:
			minScore = min(newScores)
	scoreRange = maxScore - minScore
	print result.chromosomeEnds
	offset = 0
	chromosomeSplits = result.getChromosomeSplit()
	secondRunChromosomeSplits = secondRunResult.getChromosomeSplit()
	print chromosomeSplits
	print secondRunChromosomeSplits
	ticksList1 = []
	ticksList2 = []
	textPos = []
	plt.figure(figsize=(20,4))
	for i in range(0,len(result.chromosomeEnds)):
		index1 = chromosomeSplits[i][0]
		index2 = chromosomeSplits[i+1][0]
		srIndex1 = secondRunChromosomeSplits[i][0]
		srIndex2 = secondRunChromosomeSplits[i+1][0]
		scoreList = newScores[index1:index2]
		srScoreList = newSecondRunScores[srIndex1:srIndex2]
		posList = result.positions[index1:index2]
		srPosList = secondRunResult.positions[srIndex1:srIndex2]
		
		newPosList = []
		for pos in posList:
			newPosList.append(offset+pos)
		
		newSRPosList = []
		for pos in srPosList:
			newSRPosList.append(offset+pos)

		if i%2==0:
			basicColor = 'b'
			signColor = 'r'
			srBasicColor = 'y'
			srSignColor = 'm'
		else:
			basicColor = 'g'
			signColor = 'r'
			srBasicColor = 'c'
			srSignColor = 'm'
			 
		for (pos,score) in zip(newSRPosList,srScoreList):
			if score>=bonferroniCutoffs[1]:
				plt.plot([pos],[score],"+",color=srSignColor)
			else:
				plt.plot([pos],[score],"+",color=srBasicColor)

		for (pos,score) in zip(newPosList,scoreList):
			if score>=bonferroniCutoffs[0]:
				plt.plot([pos],[score],".",color=signColor)
			else:
				plt.plot([pos],[score],".",color=basicColor)
				
		oldOffset = offset
		textPos.append(offset+result.chromosomeEnds[i]/2-2000000)
		offset += result.chromosomeEnds[i]
		if i<4:
			plt.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
		for j in range(oldOffset,offset,1000000):
			ticksList1.append(j)
		for j in range(0,result.chromosomeEnds[i],1000000):
			if j%5000000 == 0 and j < result.chromosomeEnds[i]-1500000 :
				ticksList2.append(j/1000000)
			else:
				ticksList2.append("")
				
		
	if plotBonferroni:
		plt.plot([0,sum(result.chromosomeEnds)],[bonferroniCutoffs[0],bonferroniCutoffs[0]],"k-.")
		plt.plot([0,sum(result.chromosomeEnds)],[bonferroniCutoffs[1],bonferroniCutoffs[1]],"k-.")

	plt.axis([0,sum(result.chromosomeEnds),minScore-0.05*scoreRange,maxScore+0.05*scoreRange])
	plt.xticks(ticksList1,ticksList2)
	for i in range(0,len(textPos)):
		plt.text(textPos[i],minScore-scoreRange*0.2,"Chr "+str(i+1))
	plt.subplots_adjust(right=0.98)
	plt.subplots_adjust(left=0.05)
	plt.subplots_adjust(bottom=0.15)
	plt.subplots_adjust(top=0.9)
	plt.text(offset/2,maxScore+scoreRange*0.1,'Position')
	if not ylab:
		if type=="pvals":
			plt.ylabel('-log(p-value)')
		else:
			plt.ylabel('score')
	else:
		plt.ylabel(ylab)

	if pdfFile:
		plt.savefig(pdfFile,format="pdf")
	if pngFile:
		plt.savefig(pngFile,format="png")		
	if not (pdfFile or pngFile):
		plt.show()
		
	plt.clf()

	"""
	plt.title(statistic.name)
	plt.xlabel('Observed theta')
	plt.ylabel('Estimated theta')
	"""
