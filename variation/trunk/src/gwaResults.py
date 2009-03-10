import pdb,gc

class ResultType(object):
	def __init__(self,resultType=None,fileType=None,datasetName=None,resultDir=None,mafCutoff=0,logTransform=None,name=None):
		self.resultType=resultType
		self.fileType=fileType
		self.datasetName=datasetName
		self.resultDir=resultDir
		self.mafCutoff = mafCutoff
		self.logTransform = False
		self.logTransform = logTransform
		if self.logTransform==None:
			if self.fileType==".pvals":
				self.logTransform = True
			else:
				self.logTransform = False
		self.name = name
		if not name and datasetName:
			self.name = resultType+"_"+datasetName
		if resultDir and datasetName:
			self.filenamePrefix = resultDir+resultType+"_"+datasetName+"_"
		
		
	def getFileName(self,phed,phenotypeID,secondRun=False):
		print phenotypeID
		if secondRun:
			filename = self.filenamePrefix+str(phed.getPhenotypeName(phenotypeID))+".sr"+self.fileType
		else:
			filename = self.filenamePrefix+str(phed.getPhenotypeName(phenotypeID))+self.fileType
		return filename

	def __str__(self):
		return self.resultType+"_"+self.datasetName
		

class Region(object):
	
	def __init__(self,chromosome,startPos,endPos,snps=None,snpsd_indices=None,snpsInfo=None,name=""):
		self.name=name
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.size = endPos-startPos
		self.snps = []
		if snps:
			self.snps = snps
		self.snpsd_indices = []  #Indicies in other data structures.
		if snpsd_indices:
			self.snpsd_indices = snpsd_indices
		self.snpsInfo = {}       #E.g. a dictionary of information about the snps.  
		if snpsInfo:
			self.snpsInfo = snpsInfo
			
	def __cmp__(self,other):
		return cmp((self.chromosome,self.startPos),(other.chromosome,other.startPos))

	def __str__(self):
		return "Chr.:"+str(self.chromosome)+", start pos.:"+str(self.startPos)+", end pos.:"+str(self.endPos)
		
	def get_chr_pos_str(self):
		return str(self.chromosome)+"_"+str(self.startPos)+"_"+str(self.endPos)


fri_region_small = Region(4,220000,320000,name="FRI_small")
flc_region_small = Region(5,3130000,3230000,name="FLC_small")

fri_region = Region(4,120000,550000,name="FRI_large")
flc_region = Region(5,2100000,3300000,name="FLC_large")

fri_flc_regions = [fri_region_small, fri_region, flc_region_small, flc_region]

def getRegions(regionSet,window=[25000,25000]):
	"""
	Converts a set of crh_pos into a list of regions objects.
	"""
	
	res_ls = list(regionSet)
	res_ls.sort()
	
	oldPos = 0
	countRegions = 0
	chr = -1
	regions = []
	curPosList = []
	for i in range(0,len(res_ls)):
		pos = res_ls[i][1]
		if chr!=res_ls[i][0]:
			if len(curPosList):
				regions.append(curPosList)
			curPosList=[]
			countRegions += 1
			chr = res_ls[i][0]
		elif pos-oldPos>sum(window):
			if len(curPosList):
				regions.append(curPosList)
			curPosList=[]
			countRegions += 1
		
		curPosList.append((chr,pos))
		oldPos = pos
	
	print countRegions,len(regions)
	
	regionList = []
	for region in regions:
		chr = region[0][0]
		positions = []
		for (chr,pos) in region:
			positions.append(pos)
		regionList.append(Region(chr,min(positions)-window[0],max(positions)+window[1]))
	regionList.sort()
	return regionList




class Result(object):
	"""
	Contains information on the result.
	"""
	def __init__(self,resultFile=None,snpsds=None,name=None,resultType=None,phenotypeID=None,interactionResult=False):

		self.phenotypeID =phenotypeID
		self.resultType=resultType
		self.name = name
		self.scores = [] #Scores or p-values
		self.positions = []
		self.chromosomes = []
		self.chromosomeEnds = []
		self.marfs = [] #Minor allele relative frequencies.
		self.mafs = [] #Minor allele frequencies.
		self.accessions = None

		self.secondaryScores = []
		self.secondaryPositions = []
		self.secondaryChromosomes = []
		self.orders = None
		self.ranks = None
		self.snps = None
		self.interactionResult=interactionResult
		if interactionResult:
			self.interactionPositions = []			
		if resultFile:
			self._loadResult_(resultFile)
		if snpsds:
			self._loadSnpsData_(snpsds)


	def _calcGlobalRanks_(self):
		"""
		This function should only be called right after loading the full results.
		"""
		self._rankScores_()
		self.globalRanks = self.ranks
		self.ranks = None


	def _loadResult_(self,resultFile):
		f = open(resultFile,"r")
		lines = f.readlines()
		line = lines[0].split(",")
		withMAF = len(line)>3
		start = 0
		currChrom = -1
		lastPos = 0
		if not (line[0].strip()).isdigit():
			start = 1
			#print "Header detected"
			#print line[0].strip()
		if not withMAF:
			for i in range(start,len(lines)):
				line = lines[i].split(",")
				newChrom = int(line[0].strip())
				if newChrom!=currChrom:
					currChrom = newChrom
					if lastPos:
						self.chromosomeEnds.append(lastPos)
				self.chromosomes.append(newChrom)
				if self.interactionResult:
					iPos = [int(line[-1].strip())]
					self.interactionPositions.append(iPos)
				lastPos = int(line[1].strip())
				self.positions.append(lastPos)
				self.scores.append(float(line[2].strip()))
		else:
			for i in range(start,len(lines)):
				line = lines[i].split(",")
				newChrom = int(line[0].strip())
				if newChrom!=currChrom:
					currChrom = newChrom
					if lastPos:
						self.chromosomeEnds.append(lastPos)
				self.chromosomes.append(newChrom)
				if self.interactionResult:
					iPos = [int(line[-1].strip())]
					self.interactionPositions.append(iPos)
				lastPos = int(line[1].strip())
				self.positions.append(lastPos)
				self.scores.append(float(line[2].strip()))
				self.marfs.append(float(line[3].strip()))
				self.mafs.append(float(line[4].strip()))
								
		self.chromosomeEnds.append(lastPos)
		self._calcGlobalRanks_()

	def _loadSnpsData_(self,snpsds):
		"""
		Loads the SNP data.
		"""		
		self.snps = []
		self.accessions=snpsds[0].accessions
		i = 0 #result index
		chr = -1 # chromosome (index)
		while i<len(self.scores):
			if chr != self.chromosomes[i]-1:
				chr = self.chromosomes[i]-1
				j = 0 #snpsdata index
			pos = self.positions[i]
			#print i,chr, j,len(snpsds[chr].positions)
			while j<len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j<len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])
			i += 1
		
		if pos > snpsds[chr].positions[j]:
			while j<len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j<len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])

		if i!= len( self.snps):
			print "Problems with loading SNPs",i,len( self.snps)
		else:
			print "Loading SNPs appears to have worked?"
			
		print "Loaded",len(self.snps)," SNPs and",len(self.accessions),"accessions."


	def _rankScores_(self):
		
		rank_ls = zip(self.scores, range(0,len(self.scores)))
		rank_ls.sort()
		rank_ls.reverse()
		self.orders = []
		for j in range(0,len(rank_ls)):
			(s,i) = rank_ls[j]
			self.orders.append(i)
			rank_ls[j] = (i,j)

		rank_ls.sort()
		self.ranks = []
		for (i,j) in rank_ls:
			self.ranks.append(j+1)


	def getChromosomeSplit(self):
		"""
		"""
		oldChrom = 0
		chromosomeSplits = []
		for i in range(0,len(self.scores)):
			newChrom = self.chromosomes[i]
			if oldChrom != newChrom:
				while oldChrom < newChrom:
					oldChrom += 1
					chromosomeSplits.append((i,oldChrom))
		chromosomeSplits.append((i,-1))
		return chromosomeSplits

			
	def getSecondaryChromosomeSplit(self):
		"""
		"""
		oldChrom = 0
		chromosomeSplits = []
		for i in range(0,len(self.secondaryScores)):
			newChrom = self.secondaryChromosomes[i]
			if oldChrom != newChrom:
				while oldChrom < newChrom:
					oldChrom += 1
					chromosomeSplits.append((i,oldChrom))
		chromosomeSplits.append((i,-1))
		return chromosomeSplits
 

	def negLogTransform(self):
		"""
		apply -log(x) to the pvalues (scores)
		"""
		import math
		newScores = []
		for score in self.scores:
			if score != 0.0:
				newScore = -math.log(score,10)
				newScores.append(newScore)
			else:
				newScores.append(50)  				
		self.scores = newScores
		

	def filterNicePeaks(self,scoreThreshold,singletonScoreThreshold,window=[20000,20000], method=1):
		currScoreWeight=0.2
		currChrom = -1
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		singletonCount = 0
		lastSecondaryEnd = 0
		if method == 1:
			for i in range(0,len(self.positions)):
				if currChrom != self.chromosomes[i]: #Restart
					currChrom = self.chromosomes[i]
					startIndex = i #windowIndices 
					stopIndex = i #windowIndices 
					curScoreSum = currScoreWeight*self.scores[i]				
					oldScore=0
					numSNPs = 1
				currPos = self.positions[i]

				while currPos - self.positions[startIndex]>window[0]:
					curScoreSum -= self.scores[startIndex]
					startIndex += 1
					numSNPs -= 1
				
				while stopIndex+1 < len(self.positions) and self.positions[stopIndex+1]-currPos<window[1]:
					stopIndex += 1
					curScoreSum += self.scores[stopIndex]
					numSNPs += 1

				curScoreSum -= oldScore				
				oldScore = currScoreWeight*numSNPs*self.scores[i]
				curScoreSum += oldScore

				if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold) or  (numSNPs > 5 and curScoreSum/float(numSNPs) > scoreThreshold):
					if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold):
						singletonCount +=1 
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.mafs[i])

		elif method==2:
			for i in range(0,len(self.scores)):
				if self.scores[i]>=singletonScoreThreshold:
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.mafs[i])

		# The following code locates the regions before and after the "nice" SNPs.
		j = 0
		for i in range(0,len(self.positions)):
			pos = self.positions[i]
			chr = self.chromosomes[i]
			if j < len(newPositions) and pos == newPositions[j] and chr==newChromosomes[j]:
				k = 0				
				while  i+k > lastSecondaryEnd and pos-self.positions[i+k-1]<window[0] and self.chromosomes[i+k-1]==chr:
					k -= 1
				while i+k < len(self.positions)-1 and self.positions[i+k]-pos < window[1] and self.chromosomes[i+k]==chr:
					if i+k > lastSecondaryEnd:
						self.secondaryScores.append(self.scores[i+k])
						self.secondaryPositions.append(self.positions[i+k])
						self.secondaryChromosomes.append(self.chromosomes[i+k])						
					k+= 1 
				lastSecondaryEnd = i+k-1
				j += 1 
				
		
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

		print singletonCount,"singletons were added"


	def filterPercentile(self,percentile):
		newScores = []
		for score in self.scores:
			newScores.append(score)
		newScores.sort() 
		scoreCutoff = newScores[int(len(newScores)*percentile)]
		self.filterScoreCutoff(scoreCutoff)


	def filter(self, quantile=0.98, window=[25000,25000], singletonScoreCutoff=None, nicePeaks = False, method=1):
		"""
		Filter out scores/pvalues.
				
		"""
		originalSize = len(self.scores)
		newScores = []
		for score in self.scores:
			newScores.append(score)
		newScores.sort() 
		
		if nicePeaks: 
			basic_quantile = 0.90
			top_quantile = 0.998
			singleton_quantile=0.9998
			
			scoreCutoff = newScores[int(len(newScores)*basic_quantile)]
			self._filterScoreCutoff_(scoreCutoff)
			scoreCutoff = newScores[int(len(newScores)*top_quantile)]				
			if not singletonScoreCutoff:
				singletonScoreCutoff = newScores[int(len(newScores)*singleton_quantile)]
			self.filterNicePeaks(scoreCutoff,singletonScoreCutoff,window,method)

			
		scoreCutoff = newScores[int(len(newScores)*quantile)]
		self._filterScoreCutoff_(scoreCutoff)

		finalSize = len(self.scores)

		print "Original results size =",originalSize
		print "results size after filtration =",finalSize


	def filterMAF(self, minMaf=15):
		"""
		Filter out scores/pvalues which have maf<minMaf.		
		"""
		if self.interactionResult:
			newIPos = []
			for i in range(0,len(self.scores)):
				if self.mafs[i]>= minMaf:
					newIPos.append(self.interactionPositions[i])
			del self.interactionPositions
			self.interactionPositions = newIPos

		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		for i in range(0,len(self.scores)):
			if self.mafs[i]>= minMaf:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])

		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs
		
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def filterNonSegregatingSnps(self, ecotype1, ecotype2, snpsd):
		"""
		Filter out all SNPs which are not segregating in the two accessions.		
		"""
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []

		ecotypes = snpsd.accessions
		
		e_i1 = ecotypes.index(ecotype1)
		e_i2 = ecotypes.index(ecotype2)

		res_chr_pos_list = zip(self.chromosomes,self.positions)
		snpsd_chr_pos_snp_list = snpsd.getChrPosSNPList()

		i = 0 #index in snpsd
		j = 0 #index in result
		
		while i < len(snpsd_chr_pos_snp_list):
			(chr,pos,snp) = snpsd_chr_pos_snp_list[i]
			if j<len(res_chr_pos_list) and (chr,pos)==res_chr_pos_list[j]:
				if snp[e_i1]!=snp[e_i2]:
					newScores.append(self.scores[j])
					newPositions.append(self.positions[j])
					newChromosomes.append(self.chromosomes[j])
					newMafs.append(self.mafs[j])
					newMarfs.append(self.marfs[j])					
				j += 1				
			elif j<len(res_chr_pos_list) and (chr,pos)>res_chr_pos_list[j]:
				import pdb;pdb.set_trace()
				
				print "ERROR!!!"
			i += 1
			
		n = len(self.scores)
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs
		
		print n-len(newScores),"SNPs were removed, with",len(newScores),"remaining."

			
	def filterScoreCutoff(self, scoreCutoff):

		if self.interactionResult:
			newIPos = []
			for i in range(0,len(self.scores)):
				if self.scores[i]>scoreCutoff:
					newIPos.append(self.interactionPositions[i])
			del self.interactionPositions
			self.interactionPositions = newIPos

		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		for i in range(0,len(self.scores)):
			if self.scores[i]>scoreCutoff:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs
		
	def getRegions(self,window=[25000,25000]):
		self._rankScores_()
		snpIndex = 0
		startPos = startPos = max(0,self.positions[0]-window[0])
		endPos = self.positions[0]+window[1]
		regions = []
		snpRank = self.ranks[0]
		chr = self.chromosomes[0]
		maxScore = self.scores[0]
		maxPos = self.positions[0]
		oldPos = self.positions[0]
		for i in range(1,len(self.positions)):
			pos = self.positions[i]
			if chr!=self.chromosomes[i]:
				size = endPos-startPos
				regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
				snpIndex = -1
				maxScore = 0
				startPos = max(0,pos-window[0])
				chr = self.chromosomes[i]
			elif pos-oldPos>sum(window):
				size = endPos-startPos
				regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
				maxScore = 0
				startPos = pos-window[0]
			if self.scores[i]>maxScore:
				snpIndex = i 
				snpRank = self.ranks[snpIndex]
				maxScore = self.scores[i]
				maxPos = self.positions[i]
			
			endPos = pos+window[1]
			oldPos = pos
		
		regions.append((snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos))
		
		regions.sort()
		
		self.regions = []
		for i in range(0,len(regions)):
			(snpRank,snpIndex,startPos,endPos,chr,size,maxScore,maxPos) = regions[i]
			#self.regions.append((snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,i+1))
			#def __init__(self,chromosome,startPos,endPos,snps=None,snpsd_indices=None,maxScores=None,maxPositions=None,maxSnpIndices=None,ranks=None):
			self.regions.append(Region(chr,startPos,endPos))#,maxScores=[maxScore],maxPositions=[maxPos],maxSnpIndices=[snpIndex],ranks=[snpRank]))
		
		#pdb.set_trace()
		self.regions.sort()
		#print self.regions		
		return self.regions
		

	def getTopSnps(self,n):
		"""
		returns top n SNPs
		"""
		result = self.clone()
		result.filterTopSNPs(n) #FIXME: implement
		return result
			

	def _countRegions_(self,res_ls,window=[25000,25000]):
		oldPos = 0
		countRegions = 0
		chr = -1
		for i in range(0,len(res_ls)):
			pos = res_ls[i][1]
			if chr!=res_ls[i][0]:
				countRegions += 1
				chr = res_ls[i][0]
			elif pos-oldPos>sum(window):
				countRegions += 1
			oldPos = pos
		
		return countRegions


	def getChrScorePos(self,chromosome):
		"""
		returns a list of (score,pos) tuples.
		"""
		i = 0
		while i<len(self.chromosomes) and self.chromosomes[i]<chromosome:
			i += 1
		
		posList = []
		scoreList = []
		while i<len(self.chromosomes) and self.chromosomes[i]==chromosome:
			posList.append(self.positions[i])
			scoreList.append(self.scores[i])
			i += 1

		return zip(scoreList,posList)

	def getChrPos(self):
		"""
		returns a list of (chr,pos) tuples
		"""
		posList = []
		chrList = []
		for i in range(0,len(self.positions)):
			posList.append(self.positions[i])
			chrList.append(self.chromosomes[i])
		return zip(chrList,posList)



	def getTopRegions(self,n,window=[25000,25000],minScore=None):
		"""
		returns a regionSet top n SNPs
		"""
		self._rankScores_() 
		
		regionCount = 0 
		i = 0 
		
		res_ls = []
		while regionCount < n:
			res_ls.append((self.chromosomes[self.orders[i]],self.positions[self.orders[i]]))
			res_ls.sort()
			if minScore and self.scores[self.orders[i]]<minScore:
				break
			regionCount = self._countRegions_(res_ls,window=window)
			i += 1 
		
		#print res_ls
		return set(res_ls)


	def filterTopSNPs(self,n):
		self._rankScores_() #Making sure the ranks are updated
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		l = self.orders[0:n]
		l.sort()
	
		for i in l:
			newScores.append(self.scores[i])
			newPositions.append(self.positions[i])
			newChromosomes.append(self.chromosomes[i])
			newMafs.append(self.mafs[i])
			newMarfs.append(self.marfs[i])

		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def _sortByChrPos_(self):
		res_ls = zip(self.chromosomes,self.positions,self.scores,self.mafs,self.marfs)
		res_ls.sort()
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		for (chr,pos,score,maf,marf) in res_ls:
			newScores.append(score)
			newPositions.append(pos)
			newChromosomes.append(chr)
			newMafs.append(maf)
			newMarfs.append(marf)

		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def filterTopRegions(self,n,window=[25000,25000],minScore=None):
		self._rankScores_() 
		oldScores = self.scores
		oldPositions = self.positions
		oldChromosomes = self.chromosomes
		oldMafs = self.mafs
		oldMarfs = self.marfs
		oldGlobalRanks = self.globalRanks

		self.scores = []
		self.positions = []
		self.chromosomes = []
		self.mafs = []
		self.marfs = []
		elf.globalRanks = []
		regionCount = 0 
		i = 0 
		
		res_ls = []
		while regionCount < n:
			res_ls.append((oldChromosomes[self.orders[i]],oldPositions[self.orders[i]]))
			res_ls.sort()
			if minScore and oldScores[self.orders[i]]<minScore:
				break
			self.scores.append(oldScores[self.orders[i]])
			self.positions.append(oldPositions[self.orders[i]])
			self.chromosomes.append(oldChromosomes[self.orders[i]])
			self.mafs.append(oldMafs[self.orders[i]])
			self.marfs.append(oldMarfs[self.orders[i]])
			self.globalRanks.append(oldGlobalRanks[self.orders[i]])
			regionCount = self._countRegions_(res_ls,window=window)
			i += 1 
		
		del oldScores
		del oldPositions
		del oldChromosomes
		del oldMafs
		del oldMarfs
		del oldGlobalRanks

		self._sortByChrPos_()


		
			

	def mergeWith(self,snpResult):
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		
		i1 = 0
		
		pos1 = self.positions[i1] 
		chr1 = self.chromosomes[i1]

		for i2 in range(0,len(snpResult.positions)):
			pos2 = snpResult.positions[i2]
			chr2 = snpResult.chromosomes[i2]
			while i1 <len(self.positions)-1 and (chr1,pos1)<(chr2,pos2):
				newPositions.append(self.positions[i1])
				newChromosomes.append(self.chromosomes[i1])
				newScores.append(self.scores[i1])
				newMafs.append(self.mafs[i1])
				newMarfs.append(self.marfs[i1])
				i1 += 1
				pos1 = self.positions[i1]
				chr1 = self.chromosomes[i1]

			if i1 <len(self.positions)-1 and (chr1,pos1)==(chr2,pos2):
				i1 += 1
				pos1 = self.positions[i1]
				chr1 = self.chromosomes[i1]

			newPositions.append(snpResult.positions[i2])
			newChromosomes.append(snpResult.chromosomes[i2])
			newScores.append(snpResult.scores[i2])
			newMafs.append(snpResult.mafs[i2])
			newMarfs.append(snpResult.marfs[i2])

		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

	def clone(self):
		if self.snps:
			result = SNPResult(name=self.name)
			result.snps = self.snps + []
		else:
			result = Result(name=self.name)
		result.resultType = self.resultType
		result.name = self.name
		result.scores = self.scores + []
		result.positions = self.positions + []
		result.chromosomes = self.chromosomes + []
		result.chromosomeEnds = self.chromosomeEnds + []
		result.marfs = self.marfs + [] #Minor allele relative frequencies.
		result.mafs = self.mafs + [] #Minor allele frequencies.
		if self.accessions:
			result.accessions = self.accessions + []
		if self.orders:
		   result.orders = self.orders + []
		if self.ranks:
		   result.ranks = self.ranks + []
		return result

	
	def naMAF(self, minMaf=10):
		"""
		NA scores/pvalues which have maf<minMaf.		
		"""
		for i in range(0,len(self.scores)):
			if self.mafs[i]< minMaf:
				self.scores[i]= "NA"


	def alexFiltering(self,emmaScores,cutoff=6,window=[50000,50000]): 
		"""
		
		"""
		scoreList = zip(emmaScores,self.scores)
		
		scoreList.sort()
		scoreList.reverse()
		
		cScores = []
		i=0
		emmaScore = scoreList[0][0]
		while emmaScore>cutoff: 
			cScores.append(scoreList[i][1])
			i += 1
			emmaScore = scoreList[i][0]


		#Always top 5 regions, and at most 20 regions,  continue until cutoff is reached.


		if len(cScores):
			first5Regions = self.clone()
			first5Regions.filterTopRegions(5,window=window)
			cScoreCutoff = min(cScores)
			self.filterTopRegions(20,window=window)
			self.filterScoreCutoff(cScoreCutoff)
			self.mergeWith(first5Regions)
		else:
			self.filterTopRegions(5,window=window)
			
	def updateRegions(self,regionList):
		self._rankScores_()
		i=0
		rl_i=0 #region list index
		while rl_i<len(regionList):
			region = regionList[rl_i]
			cp1=(self.chromosomes[i],self.positions[i])
			cp_start=(region.chromosome,region.startPos)
			while cp1 < cp_start:
				i += 1
				cp1=(self.chromosomes[i],self.positions[i])
			cp_end=(region.chromosome,region.endPos)

			maxScore = 0
			maxRank = 0
			maxPos = 0 
			while cp1<=cp_end:
				"""Update current region!"""
				if maxScore < self.scores[i]:
					maxScore = self.scores[i]
					maxRank = self.ranks[i]
					maxPos = self.positions[i]
				i += 1
				cp1=(self.chromosomes[i],self.positions[i])
						
			region.snpsInfo[self.resultType.name] = {"maxScore":maxScore,"maxRank":maxRank,"maxPos":maxPos}
			rl_i += 1
			
	
	def writeToFile(self,filename,format="simple"):
		f = open(filename,"w")
		f.write("Chromosome,Position,Score,MARF,MAF \n")
		for (ch,pos,score,marf,maf) in zip(self.chromosomes,self.positions,self.scores,self.marfs,self.mafs):
			l = map(str,[ch,pos,score,marf,maf])
			f.write(",".join(l)+"\n")
		f.close()
		
		
		
			
		
			

class SNPResult(Result):
	"""
	Contains information on the result.
	"""

	def _loadSnpsData_(self,snpsds):
		"""
		Loads the SNP data.
		"""		
		self.snps = []
		self.accessions=snpsds[0].accessions
		i = 0 #result index
		chr = -1 # chromosome (index)
		while i<len(self.scores):
			if chr != self.chromosomes[i]-1:
				chr = self.chromosomes[i]-1
				j = 0 #snpsdata index
			pos = self.positions[i]
			#print i,chr, j,len(snpsds[chr].positions)
			while j<len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j<len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])
			i += 1
		
		if pos > snpsds[chr].positions[j]:
			while j<len(snpsds[chr].positions) and pos > snpsds[chr].positions[j]:
				j += 1
			if j<len(snpsds[chr].positions) and pos == snpsds[chr].positions[j]:
				self.snps.append(snpsds[chr].snps[j])

		if i!= len( self.snps):
			print "Problems with loading SNPs",i,len( self.snps)
			
		print "Loaded",len(self.snps)," SNPs and",len(self.accessions),"accessions."



	def filterNicePeaks(self,scoreThreshold,singletonScoreThreshold,window=[20000,20000], method=1):
		currScoreWeight=0.2
		currChrom = -1
		newScores = []
		newSnps = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		singletonCount = 0
		lastSecondaryEnd = 0
		if method == 1:
			for i in range(0,len(self.positions)):
				if currChrom != self.chromosomes[i]: #Restart
					currChrom = self.chromosomes[i]
					startIndex = i #windowIndices 
					stopIndex = i #windowIndices 
					curScoreSum = currScoreWeight*self.scores[i]				
					oldScore=0
					numSNPs = 1
				currPos = self.positions[i]

				while currPos - self.positions[startIndex]>window[0]:
					curScoreSum -= self.scores[startIndex]
					startIndex += 1
					numSNPs -= 1
				
				while stopIndex+1 < len(self.positions) and self.positions[stopIndex+1]-currPos<window[1]:
					stopIndex += 1
					curScoreSum += self.scores[stopIndex]
					numSNPs += 1

				curScoreSum -= oldScore				
				oldScore = currScoreWeight*numSNPs*self.scores[i]
				curScoreSum += oldScore

				if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold) or  (numSNPs > 5 and curScoreSum/float(numSNPs) > scoreThreshold):
					if (numSNPs < 5 and self.scores[i] > singletonScoreThreshold):
						singletonCount +=1 
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.marfs[i])
					newSnps.append(self.snps[i])

		elif method==2:
			for i in range(0,len(self.scores)):
				if self.scores[i]>=singletonScoreThreshold:
					newScores.append(self.scores[i])
					newPositions.append(self.positions[i])
					newChromosomes.append(self.chromosomes[i])
					newMafs.append(self.mafs[i])
					newMarfs.append(self.marfs[i])
					newSnps.append(self.snps[i])

		# The following code locates the regions before and after the "nice" SNPs.
		j = 0
		for i in range(0,len(self.positions)):
			pos = self.positions[i]
			chr = self.chromosomes[i]
			if j < len(newPositions) and pos == newPositions[j] and chr==newChromosomes[j]:
				k = 0				
				while  i+k > lastSecondaryEnd and pos-self.positions[i+k-1]<window[0] and self.chromosomes[i+k-1]==chr:
					k -= 1
				while i+k < len(self.positions)-1 and self.positions[i+k]-pos < window[1] and self.chromosomes[i+k]==chr:
					if i+k > lastSecondaryEnd:
						self.secondaryScores.append(self.scores[i+k])
						self.secondaryPositions.append(self.positions[i+k])
						self.secondaryChromosomes.append(self.chromosomes[i+k])						
					k+= 1 
				lastSecondaryEnd = i+k-1
				j += 1 
				
		
		self.snps = newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

		print singletonCount,"singletons were added"



	def filterMAF(self, minMaf=15):
		"""
		Filter out scores/pvalues which have maf<minMaf.		
		"""
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps=[]
		for i in range(0,len(self.scores)):
			if self.mafs[i]>= minMaf:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps=newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

	
		
	def filterScoreCutoff(self, scoreCutoff):
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps=[]
		for i in range(0,len(self.scores)):
			if self.scores[i]>scoreCutoff:
				newScores.append(self.scores[i])
				newPositions.append(self.positions[i])
				newChromosomes.append(self.chromosomes[i])
				newMafs.append(self.mafs[i])
				newMarfs.append(self.marfs[i])
				newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps=newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def getRegionSNPs(self,window=[25000,25000],snpsDataFile=None):
		if self.snps:
			regionSNPs = []
			self.getRegions(window=window)
			for (snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,regionRank) in self.regions:
				snp = SNP(self.positions[snpIndex],self.chromosomes[snpIndex],alleles=self.snps[snpIndex],accessions=self.accessions,score=maxScore)
				regionSNPs.append(snp)
			return regionSNPs
			
	
	def filterTopSNPs(self,n):
		self._rankScores_() #Making sure the ranks are updated
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps=[]
		l = self.orders[0:n]
		l.sort()
	
		for i in l:
			newScores.append(self.scores[i])
			newPositions.append(self.positions[i])
			newChromosomes.append(self.chromosomes[i])
			newMafs.append(self.mafs[i])
			newMarfs.append(self.marfs[i])
			newSnps.append(self.snps[i])
		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps=newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def _sortByChrPos_(self):
		res_ls = zip(self.chromosomes,self.positions,self.scores,self.mafs,self.marfs,self.snps)
		res_ls.sort()
		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps=[]
		for (chr,pos,score,maf,marf,snp) in res_ls:
			newScores.append(score)
			newPositions.append(pos)
			newChromosomes.append(chr)
			newMafs.append(maf)
			newMarfs.append(marf)
			newSnps.append(snp)

		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs

		self.snps=newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs


	def filterTopRegions(self,n,window=[25000,25000],minScore=None):
		self._rankScores_() 
		oldScores = self.scores
		oldPositions = self.positions
		oldChromosomes = self.chromosomes
		oldMafs = self.mafs
		oldMarfs = self.marfs
		oldSnps = self.snps
		oldGlobalRanks = self.globalRanks
		self.snps=[]
		self.scores = []
		self.positions = []
		self.chromosomes = []
		self.mafs = []
		self.marfs = []
		self.globalRanks = []
		
		regionCount = 0 
		i = 0 
		
		res_ls = []
		if minScore:
			while regionCount < n:
				res_ls.append((oldChromosomes[self.orders[i]],oldPositions[self.orders[i]]))
				res_ls.sort()
				if oldScores[self.orders[i]]<minScore:
					break
				self.scores.append(oldScores[self.orders[i]])
				self.positions.append(oldPositions[self.orders[i]])
				self.chromosomes.append(oldChromosomes[self.orders[i]])
				self.mafs.append(oldMafs[self.orders[i]])
				self.marfs.append(oldMarfs[self.orders[i]])
				self.snps.append(oldSnps[self.orders[i]])		   
				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
				regionCount = self._countRegions_(res_ls,window=window)
				i += 1 
		else:
			while regionCount < n:
				res_ls.append((oldChromosomes[self.orders[i]],oldPositions[self.orders[i]]))
				res_ls.sort()
				self.scores.append(oldScores[self.orders[i]])
				self.positions.append(oldPositions[self.orders[i]])
				self.chromosomes.append(oldChromosomes[self.orders[i]])
				self.mafs.append(oldMafs[self.orders[i]])
				self.marfs.append(oldMarfs[self.orders[i]])
				self.snps.append(oldSnps[self.orders[i]])		   
				self.globalRanks.append(oldGlobalRanks[self.orders[i]])
				regionCount = self._countRegions_(res_ls,window=window)
				i += 1 
			
		del oldScores
		del oldPositions
		del oldChromosomes
		del oldMafs
		del oldMarfs
		del oldSnps
		del oldGlobalRanks

		self._sortByChrPos_()


	def mergeWith(self,snpResult):
		#pdb.set_trace()

		newScores = []
		newPositions = []
		newChromosomes = []
		newMafs = []
		newMarfs = []
		newSnps=[]	   
		
		if len(snpResult.scores)==0:
			return 
		elif len(self.scores)==0:
			newPositions = snpResult.positions
			newChromosomes = snpResult.chromosomes
			newScores = snpResult.scores
			newMafs = snpResult.mafs
			newMarfs = snpResult.marfs
			newSnps = snpResult.snps
		else:

			i1 = 0
		
			pos1 = self.positions[i1] 
			chr1 = self.chromosomes[i1]

			for i2 in range(0,len(snpResult.positions)):
				pos2 = snpResult.positions[i2]
				chr2 = snpResult.chromosomes[i2]
				while i1 <len(self.positions)-1 and (chr1,pos1)<(chr2,pos2):
					newPositions.append(self.positions[i1])
					newChromosomes.append(self.chromosomes[i1])
					newScores.append(self.scores[i1])
					newMafs.append(self.mafs[i1])
					newMarfs.append(self.marfs[i1])
					newSnps.append(self.snps[i1])	
					i1 += 1
					pos1 = self.positions[i1]
					chr1 = self.chromosomes[i1]

				if i1 <len(self.positions)-1 and (chr1,pos1)==(chr2,pos2):
					i1 += 1
					pos1 = self.positions[i1]
					chr1 = self.chromosomes[i1]

				newPositions.append(snpResult.positions[i2])
				newChromosomes.append(snpResult.chromosomes[i2])
				newScores.append(snpResult.scores[i2])
				newMafs.append(snpResult.mafs[i2])
				newMarfs.append(snpResult.marfs[i2])
				newSnps.append(snpResult.snps[i2])
				#pdb.set_trace()

		del self.snps
		del self.scores
		del self.positions
		del self.chromosomes
		del self.mafs
		del self.marfs
		del snpResult

		self.snps=newSnps
		self.scores = newScores
		self.positions = newPositions
		self.chromosomes = newChromosomes
		self.mafs = newMafs
		self.marfs = newMarfs

			
		

class RegionsTable(object):
	"""
	A table or regions X methods/phenotypes.
	"""
	def __init__(self,result_ls,window=[25000,25000]):
		merged_result = result_ls[0].clone()
		for i in range(1,len(result_ls)):
			result = result_ls[i]
			merged_result.mergeWith(result)
		merged_result.getRegions(window=window)
		totalRegSize = 0
		if len(merged_result.positions):
			for reg in merged_result.regions:
				totalRegSize += reg[4]
		#print "The union of the results: Number of sign. SNPs: "+str(len(merged_result.scores))+", number of sign. regions: "+str(len(merged_result.regions))+", ave. region size: "+str(totalRegSize/float(len(merged_result.regions)))+".\n"

		self.regions = []
		for (snpIndex,startPos,endPos,chr,size,maxScore,maxPos,snpRank,regionRank) in merged_result.regions:
			self.regions.append((chr,startPos,endPos,size))
		del merged_result

		self.region_by_methods_table = [] # list[region_index][method/phenotype_index][snp_index]
		for (chr,startPos,endPos,size) in self.regions:  #For all regions
			methods_snps_ls = []
			for m_i in range(0,len(result_ls)):      #for all results
				result = result_ls[m_i]
				snps_ls = []				
				
				if result.snps:
					##Finding all SNPs in region of interest
					regionRank = None
					snpRank = None
					result.getRegions(window=window)
					for (si,startPos2,endPos2,chr2,size2,mScore,mPos,sRank,rRank) in result.regions:
						if chr2==chr and startPos<=startPos2 and endPos2<=endPos:
							#print startPos,startPos2,endPos2,endPos
							regionRank = rRank
							snpRank = sRank
					#if not regionRank:
						#pdb.set_trace()

				##Finding all SNPs in region of interest
				for i in range(0,len(result.positions)):
					if result.chromosomes[i]==chr and startPos<result.positions[i]<endPos:
						if result.snps:
							#print result.snps[i]
							snp = SNP(result.positions[i],chr,alleles=result.snps[i],score=result.scores[i],rank=snpRank,regionRank=regionRank)
						else:
							snp = (chr,result.positions[i],result.scores[i],i)  #(chr,pos,score,index)
						snps_ls.append(snp)  
				methods_snps_ls.append(snps_ls)
			self.region_by_methods_table.append(methods_snps_ls)

		self.resultNames = []
		for result in result_ls:
			self.resultNames.append(result.name)




class SNP(object):
	"""
	A class to represent a SNP. 

	It's used only when analysing a SNP.
	"""

	def __init__(self,position,chromosome,accessions=None,alleles=None,snpsds=None,score=None,rank=None,regionRank=None):
		self.position = position
		self.chromosome = chromosome

		self.alleles = None
		self.accessions=None
		self.score=None
		self.rank=None
		self.regionRank=None

		if not alleles and snpsds:
			self._getAllele_(snpsds)
		else:
			self.alleles = alleles
		if accessions:
			self.accessions=accessions
		if score:
			self.score=score
		if rank:
			self.rank=rank
		if regionRank:
			self.regionRank=regionRank

	def _getAllele_(self,snpsds):
		chr = self.chromosome - 1		
		snpsd = snpsds[chr]
		self.snpsdIndex =-1
		self.alleles = None
		for i in range(0,len(snpsd.positions)):
			if snpsd.position==snpsd.positions[i]:
				self.alleles = snpsd.snps[i]
				self.snpsdIndex = i
				break
		if not self.alleles:
			print "The corresponding allele was not found in the data."

	
		
 
class Gene(object):
	"""
	A class which encompasses basic information about a gene.
	"""
	def __init__(self,chromosome=None,startPos=None,endPos=None,name="",description=None, dbRef=""):
		self.chromosome = chromosome
		self.startPos = startPos
		self.endPos = endPos
		self.exons = []
		self.introns = []
		self.tairID = ""
		self.dbRef = dbRef
		self.name = name
		self.description = description
		self.functionDescriptions = []
		self.shortDescriptions = []
		self.direction = None


	def __str__(self):
		if not self.description:
			return "Chromosome="+str(self.chromosome)+", position=("+str(self.startPos)+","+str(self.endPos)+"), tair ID="+self.tairID+", short descriptions="+str(self.shortDescriptions)+", function descriptions="+str(self.functionDescriptions)+"."
		else:
			return "Chromosome="+str(self.chromosome)+", position=("+str(self.startPos)+","+str(self.endPos)+"), tdbRef="+self.dbRef+", name="+str(self.name)+", description="+str(self.description)+"."



def getCandidateGeneList(cgl_id,host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123",db="stock_250k"):
	import MySQLdb
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

	numRows = int(cursor.execute("select distinct gm.chromosome, gm.start, gm.stop, g.locustag, g.gene_symbol, g.description, g.dbxrefs from genome.entrezgene_mapping gm, genome.gene g, stock_250k.candidate_gene_list cgl where cgl.list_type_id="+str(cgl_id)+" and gm.gene_id = g.gene_id and cgl.gene_id=g.gene_id order by gm.chromosome, gm.start, gm.stop"))
	candGenes = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		gene = Gene(int(row[0]),int(row[1]),int(row[2]),name=row[4],description=row[5],dbRef=row[6])
		candGenes.append(gene)
	cursor.close ()
	conn.close ()
	print "Candiate genelists fetched"
	return candGenes


def getResultsFilename(host,user,passwd,callMethodID,phenotypeMethodID,analysisMethodID):
	"""
	Retrieve the filename with the results.
	"""
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "stock_250k")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	#Retrieve the filenames
	print "Fetching data"
	numRows = int(cursor.execute("select rm.filename from stock_250k.results_method rm where rm.call_method_id="+str(callMethodID)+" and rm.phenotype_method_id="+str(phenotypeMethodID)+" and analysis_method_id="+str(analysisMethodID)+" "))
	filenames = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		filenames.append(row[0])
	cursor.close ()
	conn.close ()
	return filenames


	

def _getStandardResultTypes_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","newDataset",res_path+"kw_results/"))
	resultTypes.append(ResultType("Emma",".pvals","newDataset",res_path+"emma_results/",mafCutoff=15))
	resultTypes.append(ResultType("Marg",".score","newDataset",res_path+"marg_results/"))
	resultTypes.append(ResultType("RF",".imp","newDataset",res_path+"rf_results/",mafCutoff=15))
	return resultTypes

def _getStandardResultTypes2_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","raw",res_path+"kw_results/"))
	resultTypes.append(ResultType("Emma",".pvals","newDataset",res_path+"emma_results/",mafCutoff=15))
	resultTypes.append(ResultType("Marg",".score","newDataset",res_path+"marg_results/"))
	resultTypes.append(ResultType("RF",".imp","newDataset",res_path+"rf_results/",mafCutoff=15))
	return resultTypes

def _getStandardResultTypes3_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","raw",res_path+"kw_results/"))
	resultTypes.append(ResultType("Emma",".pvals","newDataset",res_path+"emma_results/",mafCutoff=15))
	return resultTypes

def _getStandardResultTypes4_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","raw",res_path+"kw_results/"))
	resultTypes.append(ResultType("Emma",".pvals","new_trans_",res_path+"emma_results/",mafCutoff=10))
	return resultTypes


def _getStandardBinaryResultTypes_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","newDataset",res_path+"kw_results/"))
	resultTypes.append(ResultType("Marg",".score","newDataset",res_path+"marg_results/"))
	resultTypes.append(ResultType("RF",".imp","newDataset",res_path+"rf_results/",mafCutoff=15))
	return resultTypes

def _getStandardSecondRunResultTypes_():
	res_path="/Network/Data/250k/tmp-bvilhjal/"	
	resultTypes = []
	resultTypes.append(ResultType("KW",".pvals","raw",res_path+"kw_results/"))
	resultTypes.append(ResultType("KW",".pvals","raw",res_path+"kw_results/"))
#	resultTypes.append(ResultType("Emma",".pvals","logTransform",res_path+"emma_results/",mafCutoff=20))
#	resultTypes.append(ResultType("Emma",".pvals","logTransform",res_path+"emma_results/",mafCutoff=20))
#	resultTypes.append(ResultType("Emma",".pvals","raw",res_path+"emma_results/",mafCutoff=15))
#	resultTypes.append(ResultType("Emma",".pvals","raw",res_path+"emma_results/",mafCutoff=15))
	return resultTypes


def loadResults(phenotypeIndices,resultTypes=None,phed=None,snpsds=None,filterPercentile=None,filterCutoffs=None,phenotypeFile="/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_111008.tsv",secondRun=False):
	
	if not phed:
		phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter='\t')
		
	if not resultTypes:
		if secondRun:
			resultTypes = _getStandardSecondRunResultTypes_()
		else:
			resultTypes = _getStandardResultTypes4_() 
		
	results_map = {}
	for i in phenotypeIndices:
		
		results = []
		for j in range(0,len(resultTypes)):
			resultType = resultTypes[j]
			phenName = phed.getPhenotypeName(i)
			if phenName:
				resultFile=resultType.getFileName(phed,i,secondRun=(secondRun and j%2==1))  #Modify back to get old results 120708
				try:
					print "Loading result file",resultFile
					if snpsds:
						result = SNPResult(resultFile,snpsds=snpsds,name=str(resultType)+"_"+phenName, resultType=resultType, phenotypeID=i)
					else:
						result = Result(resultFile,name=str(resultType)+"_"+phenName, resultType=resultType, phenotypeID=i)					
					if resultType.logTransform:
						print "Log transformed the p-values"
						result.negLogTransform()
	
					result.filterMAF(minMaf=resultType.mafCutoff)
					if filterPercentile:
						result.filterPercentile(filterPercentile)
					elif filterCutoffs:
						result.filterScoreCutoff(filterCutoffs[j])
						
					results.append(result)
				except Exception, e:
					print e.message
					print "Couldn't load",resultFile
				
		results_map[i] = results
		gc.collect()  #Calling garbage collector, in an attempt to clean up memory..
	return results_map
