ft = [1,2,3,4,5,6,7,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,80,81,82]
ionomics = range(14,32)+range(83,101)
#resDir = "/Network/Data/250k/CEGS_meeting_plots/"
   
import sys, getopt, traceback, util, pdb, gc
import dataParsers
import phenotypeData
	
import plotResults, gwaResults,pylab

class RegionPlotter():
	
	def __init__(self,phenotypeIndices=None,snpsds=None,results=None,results_map=None):
		phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
		self.phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
		#if snpsds:
		self.snpsds = snpsds
		#else:
		#	snpsDataFile="/Network/Data/250k/dataFreeze_011209/250K_f13_012509.csv"
		#	self.snpsds = dataParsers.parseCSVData(snpsDataFile, format=1, deliminator=",")
		self.results_map = {}	
		if results_map:
			self.results_map=results_map
		elif results:
			for result in results:
				if self.results_map.has_key(result.phenotypeID):
					self.results_map[result.phenotypeID].append(result)
				else:
					self.results_map[result.phenotypeID] = [result]					
		elif phenotypeIndices:
			self.loadData(phenotypeIndices)
			
		#print len(results_map[phenotypeIndices[0]])


	def loadData(self,phenotypeIndices):
	
		res_path="/Network/Data/250k/tmp-bvilhjal/"

		resultsDirs = [res_path+"kw_results/",res_path+"kw_results/",res_path+"emma_results/",res_path+"emma_results/"]
		methods=["KW","KW","Emma","Emma"]
		fileTypes=[".pvals",".sr.pvals",".pvals",".sr.pvals",]
		datasetNames=["new_raw","new_raw","new_trans","new_trans"]
		logTransform = [True,True,True,True]
		mafCutoffs = [0,0,15,15]
		percentiles = [0.0,0.9,0.0,0.9]
		interactionResults = [False,True,False,True]

		#mrIndex = 1  #The guiding (main) result

		self.results_map = {}
		for i in phenotypeIndices:
			phenName=self.phed.getPhenotypeName(i)
			
			results = []
			resultTypes = []
			for j in range(0,len(methods)):
				#try:
					resultFile=resultsDirs[j]+methods[j]+"_"+datasetNames[j]+"_"+phenName+fileTypes[j]
					print "Loading result file",resultFile
					rt = gwaResults.ResultType(resultType=methods[j],fileType="pvals",mafCutoff=mafCutoffs[j])
					if self.snpsds:
						result = gwaResults.SNPResult(resultFile,self.snpsds,name=methods[j]+"_"+datasetNames[j]+"_"+phenName, resultType=rt, phenotypeID=i)
					else:
						result = gwaResults.Result(resultFile,name=methods[j]+"_"+datasetNames[j]+"_"+phenName, resultType=rt, phenotypeID=i,interactionResult=interactionResults[j])
					if logTransform[j]:
						print "Log transforming the p-values.."
						result.negLogTransform()

					result.filterMAF(minMaf=mafCutoffs[j])
					if percentiles[j]>0:
						result.filterPercentile(percentiles[j])
					results.append(result)
					resultTypes.append(methods[j])
				#except Exception:
				#	print "Couldn't load",resultFile
								
			self.results_map[i] = results
			gc.collect()  #Calling garbage collector, in an attempt to clean up memory..

	
	def plotReg(self,region,phenotypeID,snpPos=None,pdfFile=None,pngFile=None,tairFile=None,plotGenes=True,printTairInfo=True,binary=False,results=None):
		return self.plotRegion(region.chromosome,region.startPos,region.endPos,phenotypeID,snpPos=snpPos,pdfFile=pdfFile,pngFile=pngFile,tairFile=tairFile,plotGenes=plotGenes,printTairInfo=printTairInfo,binary=binary,results=results)

	def plotRegion(self,chr,startPos,endPos,phenotypeID,snpPos=None,pdfFile=None,pngFile=None,tairFile=None,plotGenes=True,printTairInfo=True,binary=False,results=None,maxLogPval=None,plotColors=None,plotOrder=None,bonferoniCorrections=[6.64,9.02],labels=None):
		if not results:
			results = self.results_map[phenotypeID]
		resultTypes = [result.resultType.resultType for result in results]
		phenName=self.phed.getPhenotypeName(phenotypeID)
		tairGenes = self.getTairAnn(startPos,endPos,chr)
		
		if tairFile:	
			f = open(tairFile,'w')
			for gene in tairGenes:
				f.write(str(gene)+"\n")
			f.close()
		if printTairInfo:
			print "\nGene annotation:"
			for gene in tairGenes:
				print gene
			print '\n'

		tairInfo = []
		for gene in tairGenes:
			tairInfo.append(str(gene))

		genes = None
		if plotGenes:
			genes = tairGenes
		if not snpPos:
			self.drawSNPPlot(chr,startPos,endPos,phenName,results,resultTypes,genes=genes,pdfFile=pdfFile,pngFile=pngFile,binary=binary,
					maxLogPval=maxLogPval,plotColors=plotColors,plotOrder=plotOrder,bonferoniCorrections=bonferoniCorrections,labels=labels)
		else:
			snpsd = self.snpsds[chr-1]
			pos = snpsd.positions[0]
			posDiff = abs(pos-snpPos)
			i = 0
			while i+1 < len(snpsd.positions) and abs(snpsd.positions[i+1]-snpPos)<posDiff:
				i += 1
				posDiff = abs(snpsd.positions[i]-snpPos)
				
			if posDiff == 0:
				print "SNP was found."
			else:
				print "SNP was not found.  Using the nearest SNP instead."
				print chr,snpsd.positions[i]
				
			snp = gwaResults.SNP(snpsd.positions[i],chr,alleles=snpsd.snps[i])
			self.drawSNPPlot(chr,startPos,endPos,phenName,results,resultTypes,genes=genes,ldSNP=snp,pdfFile=pdfFile,pngFile=pngFile,binary=binary,
					maxLogPval=maxLogPval,plotColors=plotColors,plotOrder=plotOrder,bonferoniCorrections=bonferoniCorrections,labels=labels)
		
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
		
		return tairInfo

	
	
	def drawSNPPlot(self,chr,startPos,endPos,phenName,results,resultTypes,genes=None,ldSNP=None,pdfFile=None,
			pngFile=None,binary=False,maxLogPval=None,plotColors=None,plotOrder=None,bonferoniCorrections=[6.64,9.02],labels=None):
		"""
		Draws a snp-plot
		
		At least on result is required.
		
		requires pylab to be installed.
		"""
  		#FIXME: Handle ldSNP....
   
   		if not plotColors:
   			plotColors = []
   			defaultColors = ['b','r','c','g','m','y']
   			for i in range(0,len(results)):
   				plotColors.append(defaultColors[i%len(defaultColors)])
   
		if not plotOrder:
			plotOrder = range(1,len(results)+1)
			
		if not labels:
   			labels = [""]*len(results)
			
   		if not maxLogPval:
   			maxList = []
   			for result in results:
   				maxList.append(max(result.scores))
   			maxLogPval = max(maxList)
   			
   		
		import pylab
		size = endPos - startPos
	 
			
		maxVal = maxLogPval
		minVal = 0
		rangeVal = maxVal-minVal

		
		positionsList = []
		scoreList = []
#		if result.interactionResult:
		
		newMaxPos=0 
		for result in results:
			positions = []
			scores = []
			i = 0
			
			currPos = result.positions[0]
			currChr = result.chromosomes[0]
		
			while currChr<chr: 
				i += 1
				currChr = result.chromosomes[i]
			#Found chromsome..
	
			while currChr==chr and currPos<startPos:
				i += 1
				currPos = result.positions[i]
				currChr = result.chromosomes[i]
			#Found start..
	
			while currChr==chr and currPos<endPos:
				currPos = result.positions[i]
				pos=currPos
				if result.interactionResult:
					pos = []
					pos.append(currPos)
					pos += result.interactionPositions[i]
					
				positions.append(pos)			
				if maxLogPval and result.scores[i]>maxVal:
					scores.append(maxVal)
				else:
					scores.append(result.scores[i])			
				i += 1
				currChr = result.chromosomes[i]
				currPos = result.positions[i]
	
			positionsList.append(positions)
			scoreList.append(scores)
			if result.interactionResult:
				for p_list in positions:
					if max(p_list)> newMaxPos:
						newMaxPos=max(p_list)
		
		endPos = max(endPos,newMaxPos)
		posRange = endPos-startPos


		#Now plotting...
		print "Now plotting.."
		if genes:
			numGeneLines= int(2+len(genes)/14)
			pylab.figure(1,figsize=(18,4+0.4*numGeneLines))
		else:
			pylab.figure(1,figsize=(18,4))


		for i in plotOrder:
			print "plotting",i," result."
			positions = positionsList[i-1]
			scores = scoreList[i-1]
			c = plotColors[i-1]
			l = labels[i-1]
			print "interactionResult:",results[i-1].interactionResult
			if results[i-1].interactionResult and len(positions)>0:
				for (pos,score) in zip(positions,scores):
					p_list = list(pos)
					s_list = [score]*len(p_list)
					pylab.plot(p_list,s_list,":",color=c,ms=1.8,lw=0.8)
					for (p,s) in zip(p_list,s_list):
						pylab.plot([p],[s],".",color=c,ms=1.8)	
				pylab.plot([p],[s],".",color=c,label=l,ms=1.8)							
			else:
				pylab.plot(positions,scores,".",color=c,label=l)

		for bfc in bonferoniCorrections:
			pylab.plot([startPos-0.05*posRange,endPos+0.05*posRange],[bfc,bfc],"k-.")

		if genes:
			numGeneLines= int(2+len(genes)/14)
			self.drawGenes(genes, gene_position_cycle=numGeneLines)
			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05-0.49*numGeneLines,maxVal+rangeVal*0.05])
		else:
			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05,maxVal+rangeVal*0.05])

		pylab.title(phenName+": chromosome "+str(chr)+".")

		pylab.subplots_adjust(right=0.98)
		pylab.subplots_adjust(left=0.03)
		pylab.subplots_adjust(bottom=0.15)
		pylab.subplots_adjust(top=0.9)
		pylab.legend(numpoints=2,handlelen=0.005)

		if pdfFile:
			pylab.savefig(pdfFile,format="pdf")
		if pngFile:
			pylab.savefig(pngFile,format="png")
		if not (pdfFile or pngFile):
			pylab.show()

		pylab.close(1)

		
	def drawSNPPlot_old(self,chr,startPos,endPos,phenName,results,resultTypes,genes=None,ldSNP=None,pdfFile=None,
			pngFile=None,binary=False,maxLogPval=12,interactionResult=[False]):
		"""
		Draws a snp-plot
		
		At least on result is required.
		
		requires pylab to be installed.
		"""
   
		import pylab
		size = endPos - startPos
	 
		maxVal = maxLogPval
		minVal = 0
		rangeVal = maxVal-minVal
		pvals = [True,True,True,True]

		positionsList = []
		scoreList = []
		snps = []
		
		result = results[0]
		positions = []
		scores = []
		i = 0
		currPos = result.positions[0]
		currChr = result.chromosomes[0]
		if interactionResults[0]:
			iPositions = []
		
		while currChr<chr: 
			i += 1
			currChr = result.chromosomes[i]
		#Found chromsome..

		while currChr==chr and currPos<startPos:
			i += 1
			currPos = result.positions[i]
			currChr = result.chromosomes[i]
		#Found start..

		while currChr==chr and currPos<endPos:
			currPos = result.positions[i]
			pos = []
			pos.append(currPos)
			if interactionResults[0]:
				pos += result.interactionPositions
			currChr = result.chromosomes[i]
			if ldSNP:
				snps.append(result.snps[i])
			positions.append(pos)			
			if result.scores[i]>maxVal and pvals[0]:
				scores.append(maxVal)
			else:
				scores.append(result.scores[i])			
			i += 1

		positionsList.append(positions)
		scoreList.append(scores)


		for j in range(1,len(results)):
			result = results[j]
			positions = []
			scores = []
			i = 0
			currPos = result.positions[0]
			currChr = result.chromosomes[0]
			while currChr<chr: 
				i += 1
				currChr = result.chromosomes[i]
			#Found chromsome..

			while currChr==chr and currPos<startPos:
				i += 1
				currPos = result.positions[i]
				currChr = result.chromosomes[i]
			#Found start..
	
			while currChr==chr and currPos<endPos:
				currPos = result.positions[i]
				pos = []
				pos.append(currPos)
				if interactionResults[0]:
					pos.append(result.interactionResults[i])
				currChr = result.chromosomes[i]
				positions.append(pos)			
				if result.scores[i]>10 and pvals[j]:
					scores.append(10.0)
				else:
					scores.append(result.scores[i])			
				i += 1

			positionsList.append(positions)
			scoreList.append(scores)
			#Found the end
		
		startPos = positionsList[0][0]
		endPos = positionsList[0][len(positionsList[0])-1]
		for i in range(1,len(positionsList)):
			positions = positionsList[i]
			if positions[0][0]<startPos:
				startPos = positions[0][0]
			if positions[len(positions)-1][0]>endPos:
				endPos = positions[len(positions)-1][0]
		posRange = endPos-startPos

		if genes:
			numGeneLines= int(2+len(genes)/14)
			pylab.figure(1,figsize=(18,4+0.4*numGeneLines))
		else:
			pylab.figure(1,figsize=(18,4))
		if ldSNP:		  
			r2Values = []
			for i in range(0,len(snps)):
				snp = snps[i]
				pos = positionsList[0][i]
				r2Values.append(self.r2_ld(snp,ldSNP.alleles)*rangeVal+minVal)
				
			pylab.plot(positionsList[0],r2Values,"y-",label='$r^2$')

		for j in range(2,len(results)):
			minScore = min(results[j].scores)
			maxScore = max(results[j].scores)
			scoreRange = maxScore-minScore
			scores = []
			for score in scoreList[j]:
				scores.append(((score-minScore)/scoreRange)*rangeVal+minVal)
			if resultTypes[j] == "Marg":
				pylab.plot(positionsList[j],scores,"g",label='Marg')
			elif resultTypes[j] == "RF":
				pos = positionsList[j][0]
				score = scores[0]
				pylab.plot([pos,pos],[0,score+0.03],"c-",label='RF',lw=1.8)
				for k in range(1,len(scores)):
					pos = positionsList[j][k]
					score = scores[k]
					pylab.plot([pos,pos],[0,score+0.03],"c-",lw=1.8)
			elif resultTypes[j] == "CS":
				pylab.plot(positionsList[j],scores,"b.",label='CS')

	
		pylab.plot(positionsList[0],scoreList[0],"r.",label='KW')
#		if not binary:
#			pylab.plot(positionsList[1],scoreList[1],"b.",label='Emma')

		pylab.plot([startPos-0.05*posRange,endPos+0.05*posRange],[6.68,6.68],"k:")

		if genes:
			numGeneLines= int(2+len(genes)/14)
			self.drawGenes(genes, gene_position_cycle=numGeneLines)
			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05-0.49*numGeneLines,maxVal+rangeVal*0.05])
		else:
			pylab.axis([startPos-0.05*posRange,endPos+0.05*posRange,minVal-rangeVal*0.05,maxVal+rangeVal*0.05])

		if ldSNP:
			pylab.title(phenName+": chromosome "+str(chr)+", position "+str(ldSNP.position)+".")
		else:
			pylab.title(phenName+": chromosome "+str(chr)+".")

		pylab.subplots_adjust(right=0.98)
		pylab.subplots_adjust(left=0.03)
		pylab.subplots_adjust(bottom=0.15)
		pylab.subplots_adjust(top=0.9)
		pylab.legend(numpoints=2,handlelen=0.005)

		if pdfFile:
			pylab.savefig(pdfFile,format="pdf")
		if pngFile:
			pylab.savefig(pngFile,format="png")
		if not (pdfFile or pngFile):
			pylab.show()

		pylab.close(1)



	def _to01Format_(self,snp):
		all1 = snp[0]
		tSnp = [0]*len(snp)
		for i in range(1,len(snp)):
			allele = snp[i]
			if allele != all1:
				tSnp[i]=1
		return tSnp
			
		
	def r2_ld(self,snp1,snp2):
		tSnp1 = self._to01Format_(snp1)
		tSnp2 = self._to01Format_(snp2)
		delta = 1.0/float(len(snp1))
		freqs = [0.0]*4
		for i in xrange(0,len(snp1)):
			val = tSnp1[i]*2+tSnp2[i]
			freqs[val] += delta
	
	
		f1 = freqs[1]+freqs[3]
		f2 = freqs[2]+freqs[3]
		D = freqs[3]-f1*f2
		divisor = f1*f2*(1-f1)*(1-f2)
		if divisor != 0:
			r2= D*D/divisor
		else:
			r2 = -1
		return r2 


	def getTairAnn(self,startPos,endPos,chr):
		host = "papaya.usc.edu"
		user = "bvilhjal"
		passwd = "bamboo123"
		db = "T8_annotation_TH"

		import MySQLdb
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
		import warnings
		
		warnings.filterwarnings("ignore")
       		cursor.execute("select distinct pub_locus, start, end from t8_063008 where start > "+str(startPos)+" and end < "+str(endPos)+" and chromosome="+str(chr)+" and segment_type='gene' order by start")
		genes = []
		currTairID=""
		while(1):
			try:
				row = cursor.fetchone()
				if not row:
					break;
				gene = gwaResults.Gene()
				gene.startPos = int(row[1])
				gene.endPos = int(row[2])
				gene.tairID=row[0]
				gene.chromosome=chr
				genes.append(gene)
			except Warning:
				pass

		cursor.execute("select distinct t8.tair_id, t8_fd.short_description, t8_fd.description from t8_063008 t8, t8_func_desc t8_fd where t8.pub_locus=t8_fd.tair_id and t8.start > "+str(startPos)+" and t8.end < "+str(endPos)+" and t8.chromosome="+str(chr)+" order by t8.tair_id")

		functionDescriptions = []
		while(1):
			row = cursor.fetchone()
			if not row:
				break;
			functionDescriptions.append(row)   
		cursor.close ()
		conn.close ()

		for gene in genes:
			for fdesc in functionDescriptions:
				if gene.tairID==fdesc[0]:
					gene.shortDescriptions.append(fdesc[1])
					gene.functionDescriptions.append(fdesc[2])
			

		return genes


	def drawGenes(self,genes, gene_position_cycle=6):
		"""
		2008-09-27 More or less borrowed from Yu Huang..
		"""
		#print "\t Drawing gene model  ..."
		no_of_genes_drawn = 0
		for gene in genes:
			y_value = no_of_genes_drawn%gene_position_cycle+0.3 #cycling through the y position to avoid clogging
			self.plot_one_gene(gene, y_value=y_value)
			no_of_genes_drawn += 1
		#print "Done drawing genes..."



	def plot_one_gene(self,gene , y_value=1, buffer=0.4):  #ax: pylab.axis obj.
		"""
		2008-09-29: Code largely borrowed from Yu Huang..		  
		"""
		y_value = buffer+y_value/2.0
		pylab.text(gene.startPos, -y_value+0.08, gene.tairID, size=8)
		pylab.plot([gene.startPos,gene.endPos],[-y_value,-y_value],"k",linewidth=2)


	def drawGWPlot(phenotypeID):
		"""
		Draws all the GWA plots for 6 methods.
		"""
		import plotResults		
		results = self.results_map[phenotypeID]
		for m_i in range(0,len(results)): #For all methods 
			result = results[m_i]
			print "\nPlotting result",result.name,":"
			if m_i==0:
				plotResults.plotResult(result,ylab="KW: -log(p-value)")
			elif m_i==1:
				plotResults.plotResult(result,ylab="Emma: -log(p-value)")
			elif m_i==2:
				plotResults.plotResult(result,type="score",ylab="Margarita: ARG score")
			elif m_i==3:
				plotResults.plotResult(result,type="score",ylab="RF: Importance score")
			elif m_i==4:
				plotResults.plotResult(result,type="score",ylab="Simple composite rank score")





def fun1(rp):
	"""
	chr5 2.29 - 2.8 M (complex regions)
	chr5 5.17 - 6.8M (complex regions)
	chr5 3.5M (myb - emma peak)
	chr5 6.8 ( Emma peak -HHP1) 
	chr5 7.78M (Agamous-like)
	chr5 8.3 (Marg high - other methods not so)
	chr2 9.58 - 9.62M (SVP and Agamous like)
	chr1 19.79M (spa4 - nice peak)
	"""
	chr = 5
	startPos = 2250000
	endPos = 2810000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)
   
 

def fun3(rp):
	"""
	"""
	chr = 5
	startPos = 6480000
	endPos = 7020000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)



def fun4(rp):
	"""
	"""
	chr = 1
	startPos = 3800000
	endPos = 4010000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)



def fun5(rp):
	"""
	"""
	chr = 4
	startPos = 130000
	endPos = 700000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_noGenes"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",plotGenes=False)


def fri(rp):
	"""
	"""
	chr = 4
	startPos = 200000
	endPos = 350000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_ler"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=268809)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=268809)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_col"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=269962)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")

def fri2(rp):
	"""
	"""
	chr = 4
	startPos = 360000
	endPos = 510000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=409692)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=409692)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=429928)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig3"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=454542)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")

def fri_test():
	"""
	"""
	resDir = "/Network/Data/250k/tmp-bvilhjal/"
	chr = 4
	startPos = 150000
	endPos = 550000
	
	
	resultFiles=[("/Users/bjarni/tmp/FRI_FT_22.pvals","Emma","local"),("/Network/Data/250k/tmp-bvilhjal/emma_results/Emma_newDataset_7_FT_22C.pvals","Emma","global")]
	results = []		

	for (filename,method,name) in resultFiles:
		rt = gwaResults.ResultType(resultType=method,datasetName=name,mafCutoff=15,logTransform=True)
		result = gwaResults.Result(filename,name=method+"_"+name, resultType=rt, phenotypeID=7)
		result.negLogTransform()
		print "Log transformed the p-values"
		result.filterMAF(minMaf=15)
		results.append(result)

	results_map={7:results}
	rp = RegionPlotter(phenotypeIndices=[7],results_map=results_map)
	rp.plotRegion(chr,startPos,endPos,7)
	
	
def flc_test():
	"""
	"""
	resDir = "/Network/Data/250k/tmp-bvilhjal/"
	chr = 5
	startPos = 2250000
	endPos = 3350000
	
	
	resultFiles=[("/Users/bjarni/tmp/FLC_FT_22.pvals","Emma","local"),("/Network/Data/250k/tmp-bvilhjal/emma_results/Emma_newDataset_7_FT_22C.pvals","Emma","global")]
	results = []		

	for (filename,method,name) in resultFiles:
		rt = gwaResults.ResultType(resultType=method,datasetName=name,mafCutoff=15,logTransform=True)
		result = gwaResults.Result(filename,name=method+"_"+name, resultType=rt, phenotypeID=7)
		result.negLogTransform()
		print "Log transformed the p-values"
		result.filterMAF(minMaf=15)
		results.append(result)

	results_map={7:results}
	rp = RegionPlotter(phenotypeIndices=[7],results_map=results_map)
	rp.plotRegion(chr,startPos,endPos,7)
	


def svp(rp):
	"""
	"""
	chr = 2
	startPos = 9530000
	endPos = 9680000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9588500)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=9588500)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9611600)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")



def ie1(rp):
	"""
	Interesting example
	"""
	chr = 3
	startPos = 18875000
	endPos = 19025000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18923922)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=18923922)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")


def ie2(rp):
	"""
	Interesting example
	"""
	chr = 5
	startPos = 18560000
	endPos = 18660000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18625634)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=18625634)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=18582030)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")




def ie3(rp):
	"""
	Interesting example
	"""
	chr = 1
	startPos = 6340000
	endPos = 6410000
	for pi in ft:
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
		if pi != ft[0]:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=6369765)
		else:
			rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",snpPos=6369765)
		filename = resDir+"ft_chr"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf")


def lesioning(rp):
	"""
	Interesting example
	"""
	chr = 4
	startPos = 8220000
	endPos = 8360000
	pi=77
	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
	tairFileName = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig1"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=8274503,binary=True)
	filename = resDir+"Lesioning_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld_fig2"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=8297531,binary=True)

 


def avr(rp):
	"""
	Interesting example
	"""
	chr = 3
	startPos = 2150000
	endPos = 2300000
	pi=33
	phenName = "avrB"
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=2227823,binary=True)






def rp2(rp):
	"""
	Interesting example
	"""
	chr = 4
	startPos = 13190000
	endPos = 13270000
	pi=34
	phenName = "Rps2"
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=13226419,binary=True)
	



def rp2_b(rp):
	"""
	Interesting example
	"""
	chr = 4
	startPos = 9760000
	endPos = 9840000
	pi=34
	phenName = "Rps2"
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt",binary=True)
	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=9803559,binary=True)
	

def hyp(rp):
	"""
	Interesting example
	"""
	chr = 5
	startPos = 5015000
	endPos = 5120000
	pi=182
	phenName = "Hypocot_length"
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=5069348)
	




def ant(rp):
	"""
	Interesting example
	"""
	chr = 2
	startPos = 1890000
	endPos = 2000000
	pi=170
	phenName = "Antho_10"
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
	tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
	filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
	rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=1937020)
	



def na_23(rp):
	"""
	Interesting example
	"""
	chr = 4
	startPos = 6350000
	endPos = 6470000
	for pi in [16,85]:
		phenName = "Na_23_Soil_L"
		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",tairFile=filename+"_tair.txt")
		tairFileName = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)
		filename = resDir+phenName+"_chr_"+str(chr)+"_"+str(startPos)+"_"+str(endPos)+"_pi"+str(pi)+"_wld"
		rp.plotRegion(chr,startPos,endPos,pi,pdfFile=filename+".pdf",snpPos=6391200)
	


def _test_():
	rp = RegionPlotter([1,2,3,4,5,6,7])
	rp.plotRegion(5, 18550000, 18670000,1,pngFile='/Users/bjarni/tmp/DOG_region_1.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,2,pngFile='/Users/bjarni/tmp/DOG_region_2.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,3,pngFile='/Users/bjarni/tmp/DOG_region_3.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,4,pngFile='/Users/bjarni/tmp/DOG_region_4.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,5,pngFile='/Users/bjarni/tmp/DOG_region_5.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,6,pngFile='/Users/bjarni/tmp/DOG_region_6.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])
	rp.plotRegion(5, 18550000, 18670000,7,pngFile='/Users/bjarni/tmp/DOG_region_7.png',labels=["KW","KW SR","Emma","Emma SR"],plotColors=['b','g','r','c'],plotOrder=[2,4,1,3])


if __name__ == '__main__':
	_test_()
	pass


