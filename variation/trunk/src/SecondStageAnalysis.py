import sys
	
def _snpAnd_(snp1,snp2):
	newSnp = []
	for a1,a2 in zip(snp1,snp2):
		newSnp.append(a1 and a2)
	return newSnp

def _snpOr_(snp1,snp2):
	newSnp = []
	for a1,a2 in zip(snp1,snp2):
		newSnp.append(a1 or a2)
	return newSnp
	

def retrieveSecondRunSNPs(result,srTopQuantile,srWindowSize):
	result.filterPercentile(srTopQuantile)
	print "Generating pseudo SNPs"
	sys.stdout.flush()
	newSNPs = []
	newSNPsTypes = []
	newSNPsPositions = []
	newChromosomes = []
	#Generate new SNPs for all SNPs within WindowSize
	for i in range(0,len(result.positions)-1):
		chr = result.chromosomes[i]
		pos = result.positions[i]
		j = 1
		nextChr = result.chromosomes[i+j]
		nextPos = result.positions[i+j]
		while (i+j)<len(result.positions) and chr == nextChr and (pos + srWindowSize) > nextPos:
			if result.snps[i]!=result.snps[i+j]:
				newSNP = _snpAnd_(result.snps[i], result.snps[i+j])
				if newSNP!=result.snps[i] and newSNP!=result.snps[i+j] and sum(newSNP) != 0 and sum(newSNP) != len(newSNP):
					newSNPs.append(newSNP)
					newSNPsTypes.append("and1")
					newSNPsPositions.append((result.positions[i],result.positions[i+j]))
					newChromosomes.append(nextChr)
				
				newSNP = _snpOr_(result.snps[i], result.snps[i+j])
				if newSNP!=result.snps[i] and newSNP!=result.snps[i+j] and sum(newSNP) != 0 and sum(newSNP) != len(newSNP):
					newSNPs.append(newSNP)
					newSNPsTypes.append("or1")
					newSNPsPositions.append((result.positions[i],result.positions[i+j]))
					newChromosomes.append(nextChr)
				snp = []
				for allele in result.snps[i]:
					snp.append(abs(allele-1))
				
				newSNP = _snpAnd_(snp, result.snps[i+j])
				andSNP = newSNP
				if newSNP!=result.snps[i] and newSNP!=result.snps[i+j] and sum(newSNP) != 0 and sum(newSNP) != len(newSNP):
					newSNPs.append(newSNP)
					newSNPsTypes.append("and2")
					newSNPsPositions.append((result.positions[i],result.positions[i+j]))
					newChromosomes.append(nextChr)

				newSNP = _snpOr_(snp, result.snps[i+j])
				if newSNP!=result.snps[i] and newSNP!=result.snps[i+j] and newSNP!=andSNP and sum(newSNP) != 0 and sum(newSNP) != len(newSNP):
					newSNPs.append(newSNP)
					newSNPsTypes.append("or2")
					newSNPsPositions.append((result.positions[i],result.positions[i+j]))
					newChromosomes.append(nextChr)
			j += 1
			if (i+j)<len(result.positions):
				nextChr = result.chromosomes[i+j]
				nextPos = result.positions[i+j]

	print "Total # of second run SNPs pairs:",len(newSNPs)
	#print newSNPs[0:10]

	#Calculate MAFs and and Allele Counts for new SNPs
	newMafs = []
	newMarfs = []
	for snp in newSNPs:
		c = snp.count(0)
		l = len(snp)
		maf = min(c,l-c)
		marf = maf/float(l)
		newMafs.append(maf)
		newMarfs.append(marf)

	snps = []
	chromosomes = []
	positions = []
	snpTypes = []
	mafs = []
	marfs = []
	for i in range(0,len(newSNPs)):
		if not newMarfs[i]==0.0:
			chromosomes.append(newChromosomes[i])
			snps.append(newSNPs[i])
			positions.append(newSNPsPositions[i])
			snpTypes.append(newSNPsTypes[i])
			mafs.append(newMafs[i])
			marfs.append(newMarfs[i])
	print "Total # of second run SNPs pairs after monomorphic filtering:",len(snps)
	return {"snps":snps, "snpTypes":snpTypes, "marfs":marfs, "mafs":mafs, "positions":positions, "chromosomes":chromosomes}
