"""
Functions for ploting GWA results in pylab.
"""

import pylab, gwaResults
def plotFilteredResult(result,pdfFile,minScore=0,maxScore=10, plotBonferroni=False):
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
    pylab.figure(1,figsize=(20,4))
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

        pylab.plot(newSecPosList,secScoreList,".",color=(0.75,0.75,0.75,0.2))
        pylab.plot(newPosList,scoreList,".")
        oldOffset = offset
        textPos.append(offset+result.chromosomeEnds[i]/2-2000000)
        offset += result.chromosomeEnds[i]
        if i<4:
            pylab.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
        for j in range(oldOffset,offset,1000000):
            ticksList1.append(j)
        for j in range(0,result.chromosomeEnds[i],1000000):
            if j%5000000 == 0 and j < result.chromosomeEnds[i]-1500000 :
                ticksList2.append(j/1000000)
            else:
                ticksList2.append("")
                
        
    if plotBonferroni:
        pylab.plot([0,sum(result.chromosomeEnds)],[6.68,6.68],"k--")


    pylab.axis([0,sum(result.chromosomeEnds),minScore-0.05*scoreRange,maxScore+0.05*scoreRange])
    pylab.xticks(ticksList1,ticksList2)
    for i in range(0,len(textPos)):
        pylab.text(textPos[i],minScore-scoreRange*0.2,"Chr "+str(i+1))
    pylab.subplots_adjust(right=0.98)
    pylab.subplots_adjust(left=0.03)
    pylab.subplots_adjust(bottom=0.15)
    pylab.subplots_adjust(top=0.9)
    pylab.text(offset/2,maxScore+scoreRange*0.1,'Position')
    pylab.ylabel('-log(p-value)')
    pylab.savefig(pdfFile,format="pdf")
        
    #   result.chromosomeEnds[i]

    """
    pylab.title(statistic.name)
    pylab.xlabel('Observed theta')
    pylab.ylabel('Estimated theta')
    pylab.clf()
    """

def plotResult(result,pdfFile=None,minScore=None,maxScore=None,percentile=95,type="pvals",ylab=None, plotBonferroni=False):
    
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
    pylab.figure(figsize=(20,4))
    for i in range(0,len(result.chromosomeEnds)):
        index1 = chromosomeSplits[i][0]
        index2 = chromosomeSplits[i+1][0]
        scoreList = newScores[index1:index2]
        posList = result.positions[index1:index2]
        
        newPosList = []
        for pos in posList:
            newPosList.append(offset+pos)
        

        pylab.plot(newPosList,scoreList,".")
        oldOffset = offset
        textPos.append(offset+result.chromosomeEnds[i]/2-2000000)
        offset += result.chromosomeEnds[i]
        if i<4:
            pylab.plot([offset,offset],[minScore-0.05*scoreRange,maxScore+0.05*scoreRange],"k--")
        for j in range(oldOffset,offset,1000000):
            ticksList1.append(j)
        for j in range(0,result.chromosomeEnds[i],1000000):
            if j%5000000 == 0 and j < result.chromosomeEnds[i]-1500000 :
                ticksList2.append(j/1000000)
            else:
                ticksList2.append("")
                
        
    if plotBonferroni:
        pylab.plot([0,sum(result.chromosomeEnds)],[6.68,6.68],"k-.")

    pylab.axis([0,sum(result.chromosomeEnds),minScore-0.05*scoreRange,maxScore+0.05*scoreRange])
    pylab.xticks(ticksList1,ticksList2)
    for i in range(0,len(textPos)):
        pylab.text(textPos[i],minScore-scoreRange*0.2,"Chr "+str(i+1))
    pylab.subplots_adjust(right=0.98)
    pylab.subplots_adjust(left=0.05)
    pylab.subplots_adjust(bottom=0.15)
    pylab.subplots_adjust(top=0.9)
    pylab.text(offset/2,maxScore+scoreRange*0.1,'Position')
    if not ylab:
        if type=="pvals":
            pylab.ylabel('-log(p-value)')
        else:
            pylab.ylabel('score')
    else:
        pylab.ylabel(ylab)

    if pdfFile:
        pylab.savefig(pdfFile,format="pdf")
    else:
        pylab.show()
        
    #   result.chromosomeEnds[i]

    """
    pylab.title(statistic.name)
    pylab.xlabel('Observed theta')
    pylab.ylabel('Estimated theta')
    pylab.clf()
    """

