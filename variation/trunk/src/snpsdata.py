"""
This python library aims to do two things.
1. Offer general wrapper classes around SNPs datasets.
2. Offer basic functions which can aid analysis of the SNPs.

Bjarni Vilhjalmsson, bvilhjal@usc.edu
"""

class _SnpsData_:
    """
    An abstract superclass.
    """
    accessions = [] #list[accession_index]
    arrayIds = None #list[accession_index]
    snps = []  #list[position_index][accession_index]
    positions = [] #list[position_index]
    def __init__(self,snps,positions,baseScale=None,accessions=None,arrayIds=None):
        self.snps = snps
        self.positions = positions
        if accessions: 
            self.accessions=accessions
        if arrayIds: 
            self.arrayIds=arrayIds

    def scalePositions(self,baseScale):
        for i in range(0,len(self.positions)):
            self.positions[i] = int(self.positions[i]*baseScale)
        self.baseScale = baseScale

    def addSnp(self,snp):
        self.snps.append(snp)
	
    def addPos(self,position):
        self.positions.append(position)


class RawDecoder(dict):
    def __missing__(self, key):
        return 'NA'
    def __init__(self):
        for letter in ['A','C','G','T']:
            self[letter]=letter
    

class RawSnpsData(_SnpsData_):
    """
    Similar to the SnpsData class, except it deals with bases (A,C,G,T), instead of allele numbers (0s and 1s).
    
    Alphabet: A,C,G,T, and NA

    """
    
    def writeToFile(self,filename,chromosome, withArrayId=False):
        """
        Writes data to a file.  1 file for each chromosome.

        WARNING OLD, outdated!
        """
        outStr = ""
        fieldStrings = ["Chromosome", "Positions"]
        for acc in self.accessions:
            fieldStrings.append(str(acc))
        outStr += ", ".join(fieldStrings)+"\n"
        for i in range(0,len(self.positions)):
            outStr +=str(chromosome)+", "+str(self.positions[i])+", "+", ".join(self.snps[i])+"\n"
        f = open(filename,"w")
        f.write(outStr)
        f.close()


    def mergeData(self,snpsd, multIdenticalAcc=False):
        """
        Merges two RawSnpsData objects.

        If snps disagree, then the snps from the object called from is used.        
        """
        print "Merging datas"
        print "Number of snps:",len(self.snps),"and",len(snpsd.snps)
	# Find common accession indices
        accessionsIndices = []
        commonAccessions = []
        allAccessions = []
        if not multIdenticalAcc:            
            commonAccessions = list(set(self.accessions).intersection(set(snpsd.accessions)))
            allAccessions = list(set(self.accessions).union(set(snpsd.accessions)))
            for i in range(0,len(allAccessions)):
                acc = allAccessions[i]
                index1 = -1
                index2 = -1
                if self.accessions.count(acc):
                    index1 = self.accessions.index(acc)
                if snpsd.accessions.count(acc):
                    index2 = snpsd.accessions.index(acc)
                accessionsIndices.append([i,index1,index2])
        else:
            i = 0
            for j in range(0,len(self.accessions)):
                acc1 = self.accessions[j]
                count = 0
                for k in range(0,len(snpsd.accessions)):
                    acc2 = snpsd.accessions[k]
                    if acc1==acc2:
                        count += 1
                        accessionsIndices.append([i,j,k])
                        acc = acc1
                        if count >1:
                            acc = acc1+"_"+str(count)
                        commonAccessions.append(acc)
                        allAccessions.append(acc)
                        i += 1

            for j in range(0,len(self.accessions)):
                acc = self.accessions[j]
                if not acc in snpsd.accessions:
                    allAccessions.append(acc)
                    accessionsIndices.append([i,j,-1])
                    i += 1

            for j in range(0,len(snpsd.accessions)):
                acc = snpsd.accessions[j]
                if not acc in self.accessions:
                    allAccessions.append(acc)
                    accessionsIndices.append([i,-1,j])                        
                    i += 1

        print "Common accessions:", len(commonAccessions)
        print "Only in 1st data set", list(set(self.accessions).difference(set(commonAccessions)))
        print "Only in 2st data set", list(set(snpsd.accessions).difference(set(commonAccessions)))
        print "All accessions:",len(allAccessions)
        print snpsd.accessions
        print len(snpsd.accessions), len(list(set(snpsd.accessions)))
            
        commonSnpsPos = []
        snpErrorRate = []
        goodSnpsCounts = []
        newSnps = []
        newPositions = []
        i = 0
        j = 0 
        action = ""
        while i <= len(self.positions) and j <= len(snpsd.positions):
            if i < len(self.positions):
                pos1 = self.positions[i]
            if j < len(snpsd.positions):
                pos2 = snpsd.positions[j] 
            if i < len(self.positions) and pos1 < pos2:
                newPositions.append(pos1)
                newSnp = ['NA']*len(allAccessions)
                oldSnp = self.snps[i]
                for index in accessionsIndices:
                    newSnp[index[0]]=oldSnp[index[1]]                   
                newSnps.append(newSnp)
                i = i+1
                action ="pos1<pos2"
            elif j < len(snpsd.positions) and pos2 < pos1:
                newPositions.append(pos2)
                newSnp = ['NA']*len(allAccessions)
                oldSnp = snpsd.snps[j]
                for index in accessionsIndices:
                    newSnp[index[0]]=oldSnp[index[2]]                   
                newSnps.append(newSnp)
                j = j+1
                action ="pos2<pos1"
            elif i < len(self.positions) and j < len(snpsd.positions) and pos1==pos2:
                newPositions.append(pos1)
                newSnp = ['NA']*len(allAccessions)
                for index in accessionsIndices:
                    if index[1] !=-1 and self.snps[index[1]] != 'NA':
                        newSnp[index[0]]=self.snps[i][index[1]]            
                    elif index[2] !=-1:
                        newSnp[index[0]]=snpsd.snps[j][index[2]]            
                newSnps.append(newSnp)
                commonSnpsPos.append(pos1)
                fails = 0
                counts = 0
                for k in range(0,len(accessionsIndices)):
                    accIndices = accessionsIndices[k]
                    snp1 = self.snps[i][accIndices[1]]
                    snp2 = snpsd.snps[j][accIndices[2]]
                    if snp1!='NA' and snp2!='NA':
                        counts += 1
                        if snp1!=snp2:
                            fails = fails+1
                goodSnpsCounts.append(counts)
                error = 0
                if counts>0:
                    error = float(fails)/float(counts)
                snpErrorRate.append(error)                                        
                i = i+1
                j = j+1
                action ="pos2=pos1"
            else: 
                # One pointer has reached the end and the end and the other surpassed it, i.e. we only need to copy the remaining one..
                #print len(self.positions) -i-1, len(snpsd.positions) -j-1
                #print "last action:",action
                while i<len(self.positions):
                    newPositions.append(self.positions[i])
                    newSnp = ['NA']*len(allAccessions)
                    oldSnp = self.snps[i]
                    for index in accessionsIndices:
                        newSnp[index[0]]=oldSnp[index[1]]                   
                    newSnps.append(newSnp)
                    i = i+1
                while j<len(snpsd.positions):
                    newPositions.append(snpsd.positions[j])
                    newSnp = ['NA']*len(allAccessions)
                    oldSnp = snpsd.snps[j]
                    for index in accessionsIndices:
                        newSnp[index[0]]=oldSnp[index[2]]                   
                    newSnps.append(newSnp)
                    j = j+1
                 
                break
        
        #print newPositions[-10:]
        #print len(newPositions)
        #print len(set(newPositions))
        newSnpsd = RawSnpsData(newSnps,newPositions,accessions=allAccessions)


        print "In all",len(commonSnpsPos),"common snps found"
        print "In all",len(commonAccessions),"common accessions found"
        print "Mean Snp Error:",sum(snpErrorRate)/float(len(snpErrorRate))
        print "Number of SNPs in merged data:",len(newPositions)
        print "Number of SNPs in merged data:",len(newSnpsd.snps)
        print "Number of accessions in merged data:",len(newSnpsd.accessions)
        return newSnpsd



    def compareWith(self,snpsd, withArrayIds=0):
        """
        This function performs QC on two datasets.

        Requires accessions to be defined.

        withArrayIds = 0 (no array IDs), =1 the object called from has array IDs, =2 both objects have array IDs 
        """
        print "Comparing datas"
        print "Number of snps:",len(self.snps),"and",len(snpsd.snps)
	# Find common accession indices
        accessionsIndices = []
        accessionErrorRate = []
        accessionCallRates = [[],[]]
        accessionCounts = []
        commonAccessions = []
        arrayIds = []
        for i in range(0,len(self.accessions)):
            acc = self.accessions[i]
            if snpsd.accessions.count(acc):
                j = snpsd.accessions.index(acc)
                accessionsIndices.append([i,j])
                accessionErrorRate.append(0)
                accessionCounts.append(0)
                accessionCallRates[0].append(0)
                accessionCallRates[1].append(0)
                commonAccessions.append(acc)
                if withArrayIds>0:
                    if withArrayIds==1:
                        aId = self.arrayIds[i]
                    elif withArrayIds==2:
                        aId = (self.arrayIds[i], snpsd.arrayIds[j])
                    arrayIds.append(aId)                


        commonSnpsPos = []
        snpErrorRate = []
        snpCallRate = [[],[]]
        goodSnpsCounts = []
        i = 0
        j = 0
        while i <= len(self.positions) and j <= len(snpsd.positions):
            if i < len(self.positions):
                pos1 = self.positions[i]
            if j < len(snpsd.positions):
                pos2 = snpsd.positions[j] 
            if i < len(self.positions) and pos1 < pos2:
                i = i+1
            elif j < len(snpsd.positions) and pos2 < pos1:
                j = j+1
            elif i < len(self.positions) and j < len(snpsd.positions) and pos1==pos2:
                commonSnpsPos.append(pos1)
                fails = 0
                counts = 0
                missing1 = 0
                missing2 = 0
                for k in range(0,len(accessionsIndices)):
                    accIndices = accessionsIndices[k]
                    snp1 = self.snps[i][accIndices[0]]
                    snp2 = snpsd.snps[j][accIndices[1]]
                    if snp1!='NA' and snp2!='NA':
                        accessionCounts[k] += 1
                        counts += 1
                        if snp1!=snp2:
                            fails = fails+1
                            accessionErrorRate[k] += 1
                    elif snp1=='NA':
                        accessionCallRates[0][k]+=1
                        missing1 += 1
                    elif snp2=='NA':
                        accessionCallRates[1][k]+=1
                        missing2 += 1
                goodSnpsCounts.append(counts)
                error = 0
                if counts>0:
                    error = float(fails)/float(counts)
                snpErrorRate.append(error) 
                snpCallRate[0].append(missing1/float(len(accessionsIndices)))                                       
                snpCallRate[1].append(missing2/float(len(accessionsIndices)))                                       
                i = i+1
                j = j+1
            else: 
                # One pointer has reached and the end and the other surpassed it.
                break
        

        for i in range(0,len(accessionErrorRate)):
            if accessionCounts[i]>0:
                accessionErrorRate[i] = accessionErrorRate[i]/float(accessionCounts[i])
            accessionCallRates[0][i] = accessionCallRates[0][i]/float(len(commonSnpsPos))
            accessionCallRates[1][i] = accessionCallRates[1][i]/float(len(commonSnpsPos))
        
        """
        print "In all",len(snpErrorRate),"common snps found"
        print "In all",len(commonAccessions),"common accessions found"
        print "Common accessions IDs :",commonAccessions
        print "Common SNPs positions :", commonSnpsPos
        print "Accessions error rates",accessionErrorRate 
        print "Average Accession SNP Error:",sum(accessionErrorRate)/float(len(accessionErrorRate))
        print "SNP error rates",snpErrorRate
        print "Average Snp Error:",sum(snpErrorRate)/float(len(snpErrorRate))
        """
        
        return [commonSnpsPos, snpErrorRate, commonAccessions, accessionErrorRate, accessionCallRates, arrayIds, accessionCounts, snpCallRate]


    def getSnpsData(self):
        """
        Returns a SnpsData object correspoding to this RawSnpsData object.

        Note that some of the SnpsData attributes are a shallow copy of the RawSnpsData obj.
        """
        decoder = {'NA':-1}
        snps = []
        for i in range(0,len(self.snps)):
            k = 0
            for nt in ['A','C','G','T']:
                if nt in self.snps[i]:
                    decoder[nt]=k
                    k = k+1
            snp = []
            for nt in self.snps[i]:
                snp.append(decoder[nt])
            snps.append(snp)

        accessions = self.accessions
        positions = self.positions
        baseScale = self.baseScale
        return SnpsData(snps,positions,baseScale=baseScale,accessions=accessions)


    def filterBadSnps(self,snpsd,maxNumError=0):
        """
        Filters snps with high rate mismatches, when compared against another snpsd.
        """

        newSnps = []
        newPositions = []
        print "Comparing datas"
        print "Number of snps:",len(self.snps),"and",len(snpsd.snps)
	# Find common accession indices
        accessionsIndices = []
        commonAccessions = []
        for i in range(0,len(self.accessions)):
            acc = self.accessions[i]
            if snpsd.accessions.count(acc):
                j = snpsd.accessions.index(acc)
                accessionsIndices.append([i,j])
                commonAccessions.append(acc)         


        commonSnpsPos = []
        i = 0
        j = 0
        while i <= len(self.positions) and j <= len(snpsd.positions):
            if i < len(self.positions):
                pos1 = self.positions[i]
            if j < len(snpsd.positions):
                pos2 = snpsd.positions[j] 
            if i < len(self.positions) and pos1 < pos2:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
                i = i+1
            elif j < len(snpsd.positions) and pos2 < pos1:
                j = j+1
            elif i < len(self.positions) and j < len(snpsd.positions) and pos1==pos2:
                commonSnpsPos.append(pos1)
                fails = 0
                counts = 0
                for k in range(0,len(accessionsIndices)):
                    accIndices = accessionsIndices[k]
                    snp1 = self.snps[i][accIndices[0]]
                    snp2 = snpsd.snps[j][accIndices[1]]
                    if snp1!='NA' and snp2!='NA':
                        counts += 1
                        if snp1!=snp2:
                            fails = fails+1
                error = 0
                if counts>0:
                    error = float(fails)/float(counts)
                if error<=maxNumError:
                    newSnps.append(self.snps[i])
                    newPositions.append(self.positions[i])                
                i = i+1
                j = j+1
            else: 
                # One pointer has reached and the end and the other surpassed it.
                break
        

        print len(self.snps)-len(newSnps),"were filtered"
        self.snps = newSnps
        self.positions = newPositions
        




    def filterMissingSnps(self,maxNumMissing=0):
        """
        Removes SNPs from the data which have more than maxNumMissing missing values.
        """
        newPositions = []
        newSnps = []
        for i in range(0,len(self.positions)):
            missingCount = 0
            for nt in self.snps[i]:
                if nt =='NA':
                    missingCount += 1
            if missingCount<=maxNumMissing:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
        numRemoved = len(self.positions)-len(newPositions)
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved

    def filterMonoMorphicSnps(self):
        """
        Removes SNPs from the data which are monomorphic.
        """
        newPositions = []
        newSnps = []
        for i in range(0,len(self.positions)):
            count = 0
            for nt in ['A','C','G','T']:
                if nt in self.snps[i]:
                    count += 1
            if count>1:
                newSnps.append(self.snps[i])
                newPositions.append(self.positions[i])
        numRemoved = len(self.positions)-len(newPositions)
        self.snps = newSnps
        self.positions = newPositions
        return numRemoved


    def removeAccession(self,accession):
        """
        Removes an accession from the data.
        """
        pass


    def mergeIdenticalAccessions(self,accessionIndexList, priority):
        """
        The priority argument gives particular accessions in the list priority
        if priority is set to 0, then majority (if any) rules.
        """
        pass


    def filterMissingAccessions(self,maxMissingFraction=0.8):
        """
        Removes Accessions from the data which have are monomorphic.
        """
        pass


    def getStatString(self):
        st = "Number of accessions: "+str(len(self.accessions))+"\n"
        st += "Number of SNPs: "+str(len(self.snps))+"\n"
        uniqueAccessions = list(set(self.accessions))
        if len(uniqueAccessions)<len(self.accessions):
            st += "Not all accessions are unique. \n"+"Number of unique accessions: "+str(len(uniqueAccessions))+"\n"
            for acc in uniqueAccessions:
                count = self.accessions.count(acc)
                if count>1:
                    st += acc+" has "+str(count)+" occurrences.\n\n"
        return st

 
class SnpsData(_SnpsData_):
    """
    A class for SNPs data.  It uses 0, 1 (and 2, 3 if there are more than 2 alleles) to represent SNPs.  
    -1 is used if the allele data is missing.

    Contains various functions that aid the analysis of the data.
    """
    freqs = [] #list[position_index1][position_index2-position_index1+1] Linkage frequencies. 
    baseScale = 1 #Scaling for positions
    def __init__(self,snps,positions,baseScale=None,accessions=None):
        self.snps = snps
        self.positions = positions
        if baseScale:
            self.scalePositions(baseScale)
        if accessions: 
            self.accessions=accessions


    def calcFreqs(self,windowSize, innerWindowSize = 0): # Returns a list of two loci comparison frequencies with in a window.
        freqs =	[]
        delta =0
        if len(self.snps)>0:
            delta = 1.0/float(len(self.snps[0]))
        for i in xrange(0,len(self.snps)-1):
            l = []
            j = i+1
            while j<len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                j = j+1
            while j<len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                jfreqs = [0.0]*4
                snp1 = self.snps[i]
                snp2 = self.snps[j]
                count =	0
                for k in xrange(0,len(snp1)):
                    val = snp1[k]*2+snp2[k]
                    jfreqs[val] = jfreqs[val]+delta
                l.append(jfreqs)
                j = j+1
            freqs.append(l)
        self.freqs = freqs
        return self.freqs

    def calcFreqsUnbiased(self,windowSize, innerWindowSize = 0): # Returns a list of two loci comparison frequencies with in a window.
        """ 
        Uses a distribution of frequencies that is not dependent on the lenght of the sequence.  (Alot of data is disregarded.)
        """
        numPairs = 0  #Counts the number of pairs
        freqs =	[]
        delta =0
        if len(self.snps)>0:
            delta = 1.0/float(len(self.snps[0]))
        for i in xrange(0,len(self.snps)-1):
            if (self.positions[len(self.snps)-1] - self.positions[i]) >= windowSize:
                j = i+1
                l = []
                while j<len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                    j = j+1
                while  j<len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                    jfreqs = [0.0]*4
                    snp1 = self.snps[i]
                    snp2 = self.snps[j]
                    count =	0
                    for k in xrange(0,len(snp1)):
                        val = snp1[k]*2+snp2[k]
                        jfreqs[val] = jfreqs[val]+delta
                    l.append(jfreqs)
                    j = j+1
                    numPairs = numPairs+1
                freqs.append(l)
            else:
                break
        self.freqs = freqs
        #return self.freqs
        return numPairs
		



    def calcFreqsSimple(self):
        freqs =	[]
        delta =0
        if len(self.snps)>0:
            delta = 1.0/float(len(self.snps[0]))
        for i in xrange(0,len(self.snps)-1):
            l = []
            for j in xrange(i+1,len(self.snps)):
                jfreqs = [0.0]*4
                snp1 = self.snps[i]
                snp2 = self.snps[j]
                count =	0
                for k in xrange(0,len(snp1)):
                    val = snp1[k]*2+snp2[k]
                    jfreqs[val] = jfreqs[val]+delta
                l.append(jfreqs)
            freqs.append(l)
        self.freqs = freqs
        return self.freqs
	

    def snpsFilter(self):
        """
        Splits the SNPs up after their allele frequency ..
        """
        snpsDatas = [SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[]), SnpsData([],[])]

        l = len(self.snps[0])
        for j in xrange(0,len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0)/float(l)
            """
            if c>0.5:
            c = 1-c
            if c ==	0.5:
            c =0.499
            """
            if c ==	1:
                c =0.99999
            i = int(c*10)
            snpsDatas[i].addSnp(snp)
            snpsDatas[i].addPos(self.positions[j])
        return snpsDatas
		
    def snpsFilterMAF(self,mafs):
        """
        Filters all snps with MAF not in the interval out of dataset.
        """
        newsnps = []
        newpos = []
        #print self.snps
        l = len(self.snps[0])
        if l == 0:
            print self.snps
        for j in xrange(0,len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0)/float(l)
            if c>0.5:
                c = 1-c
            if c >mafs[0] and c<=mafs[1]:
                newsnps.append(snp)
                newpos.append(self.positions[j])
        if len(newsnps)==0:  
            print "Filtered out all snps from",len(self.snps)," what to do what to do?"
        del self.snps
        self.snps = newsnps
        del self.positions
        self.positions = newpos


    def snpsFilterRare(self,threshold=0.1):
        """Filters all snps with MAF of less than threshold out of dataset."""
        newsnps = []
        newpos = []
        #print self.snps
        l = len(self.snps[0])
        if l == 0:
            print self.snps
        for j in xrange(0,len(self.snps)):
            snp = self.snps[j]
            c = snp.count(0)/float(l)
            if c>0.5:
                c = 1-c
            if c >threshold:
                newsnps.append(snp)
                newpos.append(self.positions[j])
        #if len(newsnps)==0:  
            #print "Filtered out all snps from",len(self.snps)," what to do what to do?"
        del self.snps
        self.snps = newsnps
        del self.positions
        self.positions = newpos


    def _genRecombFile(self,filename,windowSize,maxNumPairs):
        n = len(self.snps[0]) #number of individuals/accessions
        numPairs = self.calcFreqsUnbiased(windowSize)
        if numPairs!= 0:
            filterProb = float(maxNumPairs)/numPairs
        else:
            return numPairs
        f = open(filename, 'w')
        numPairs = 0
        for i in range(0,len(self.freqs)):
            for j in range(0,len(self.freqs[i])):
                if random.random()<=filterProb:
                    numPairs = numPairs +1
                    st = str(i+1)+" "+str(j+i+2)+" "+str(self.positions[j+i+1]-self.positions[i])+" u "  # u denotes unknown as opposed to ad ancestral derived
                    st = st+str(int(self.freqs[i][j][0]*n+0.5))+" "+str(int(self.freqs[i][j][1]*n+0.5))+" "
                    st = st+str(int(self.freqs[i][j][2]*n+0.5))+" "+str(int(self.freqs[i][j][3]*n+0.5))+" 0 0 0 0 "+str(100-n)+"\n"
                    f.write(st)
        f.close()
        return numPairs

    def estimateRecomb(self, windowSize, maxNumPairs = 10000,  tempfile1 = "tmp1", tempfile2="tmp2", meanTract=200, numPoints=50):
        num = self._genRecombFile(tempfile1,windowSize,maxNumPairs)
        if num < 1:
            return [0,0,0]
        os.system(homedir+"Projects/programs/Maxhap/maxhap 1 "+homedir+"Projects/programs/Maxhap/h100rho  .0008 10 .1 0.01 500 "+str(numPoints)+" "+str(meanTract)+" < "+tempfile1+" > "+tempfile2)
        f = open(tempfile2,'r')
        lines = f.readlines()
        npairs = float(lines[1].split()[2])
        i = 2
        while(lines[i].split()[0]=="Warning:"):
            i = i+1
        rho = float(lines[i].split()[1])
        ratio = float(lines[i].split()[2])
        f.close()
        return [rho,ratio,npairs]


    def meanAF(self):
        """ Mean allele frequency. """
        if len(self.snps):
            l = float(len(self.snps[0]))
            c = 0
            for j in xrange(0,len(self.snps)):
                snp = self.snps[j]
                snpsc = snp.count(0)
                if snpsc<(l/2.0):
                    c = c+snpsc/l
                else:
                    c = c+abs((l/2.0)-snpsc)/l
            return c/len(self.snps)
        else: 
            return 0
		
    def EHH(self,snp1,snp2):
        """ Calculates the EHH between two SNPs"""
        data = self.snps[snp1:snp2]
        haplotypes = []
        haplotypecount = []
        for i in range(0,len(self.snps[0])):
            haplotype = []
            for j in range(snp1,snp2+1):
                haplotype.append(self.snps[j][i])
            if not haplotype in haplotypes:
                haplotypes.append(haplotype)
                haplotypecount.append(1.0)
            else:
                k = haplotypes.index(haplotype)
                haplotypecount[k] = haplotypecount[k] + 1.0
        s = 0.0
        for i in range(0,len(haplotypes)):
            if haplotypecount[i]>1:
                s = s+haplotypecount[i]*(haplotypecount[i]-1)
        s = s/(len(self.snps[0])*(len(self.snps[0])-1))
        return s

    def totalEHH(self,windowSize,innerWindowSize):
        """ 
        Lenght indep mean EHH statistics.. (Note: no data filtering!)
        """
        ehhcount = 0
        ehh = 0.0
        for i in xrange(0,len(self.snps)-1):
            if (self.positions[len(self.snps)-1] - self.positions[i]) >= windowSize:
                j = i+1
                l = []
                while j<len(self.snps) and (self.positions[j] - self.positions[i]) < innerWindowSize :
                    j = j+1
                while  j<len(self.snps) and (self.positions[j] - self.positions[i]) <= windowSize:
                    ehh = ehh+self.EHH(i,j)
                    ehhcount = ehhcount + 1
                    j = j+1
            else:
                break
        return [ehh,ehhcount]


def writeRawSnpsDatasToFile(filename,snpsds,chromosomes=[1,2,3,4,5], deliminator=", ", missingVal = "NA", accDecoder=None, withArrayIds = False):
    """
    Writes data to a file. 
    """
    
    print "Writing data to a file."
    numSnps = 0
    for i in range(0,len(chromosomes)):
        numSnps += len(snpsds[i].positions)
        
    accessions = snpsds[0].accessions
    for i in range(1,len(chromosomes)):
        if accessions != snpsds[i].accessions:
            raise Exception("Accessions are different between SNPs datas")


    decoder = RawDecoder()
    decoder['NA']=missingVal

    #outStr = "NumSnps: "+str(numSnps)+", NumAcc: "+str(len(accessions))+"\n"
    if withArrayIds:
        outStr = "-, -, "+", ".join(snpsds[0].arrayIds)+"\n"
    else:
        outStr = ""
    fieldStrings = ["Chromosome", "Positions"]
    if accDecoder:
        for acc in snpsds[i].accessions:
            fieldStrings.append(str(accDecoder[acc]))
    else:
        for acc in snpsds[i].accessions:
            fieldStrings.append(str(acc))
    outStr += deliminator.join(fieldStrings)+"\n"
    for i in range(0,len(chromosomes)):
        for j in range(0,len(snpsds[i].positions)):
            outStr += str(chromosomes[i])+deliminator+str(snpsds[i].positions[j])
            for k in range(0, len(snpsds[0].accessions)):
                outStr += deliminator+decoder[snpsds[i].snps[j][k]]
            outStr +="\n"
    f = open(filename,"w")
    f.write(outStr)
    f.close()


def estimateRecomb(snpsdList,baseNum,filterProb,id):
    rho = 0
    npairs = 0
    for i in range(0,len(snpsdList)):
        snpsd = snpsdList[i]
        tmp1 = "tmp"+id+"1"
        tmp2 = "tmp"+id+"2"
        (rho2,npairs2) = snpsd.estimateRecomb(baseNum,filterProb,tmp1,tmp2)
        rho = rho + rho2*npairs2
        npairs = npairs + npairs2
    rho = rho/float(npairs)
    print "rho: "+str(rho)+", npairs: "+str(npairs)
    return rho
		
def D(freqs):
    """ Returns the D' LD measure """
    p1 = freqs[1]+freqs[3]  #Always positive (if no trivial SNPs are allowed).
    p2 = freqs[2]+freqs[3]
    Dmax = min(p1*(1-p2),p2*(1-p1))
    Dmin = -min(p1*p2,(1-p2)*(1-p1))
    D = freqs[3]-p1*p2
    if D >=0.0:
        return D/Dmax
    else:
        return D/Dmin
		
def r2(freqs):
    f1 = freqs[1]+freqs[3]
    f2 = freqs[2]+freqs[3]
    D = freqs[3]-f1*f2
    divisor = f1*f2*(1-f1)*(1-f2)
    if divisor != 0:
        return D*D/divisor
    else:
        return -1


	
def r2listAll(data,windowSize,nbins=10):
    sums = {1:[0.0]*nbins,2:[0.0]*nbins,3:[0.0]*nbins,4:[0.0]*nbins,5:[0.0]*nbins,6:[0.0]*nbins,7:[0.0]*nbins,8:[0.0]*nbins,"mean":[0.0]*nbins}	    
    counts = {1:[0]*nbins,2:[0]*nbins,3:[0]*nbins,4:[0]*nbins,5:[0]*nbins,6:[0]*nbins,7:[0]*nbins,8:[0]*nbins,"tot":[0]*nbins}
    for snpsd in data:
        fsnpd = snpsd.snpsFilter()
        for k in [1,2,3,4,5,6,7,8]:
            freqs =	fsnpd[k].calcFreqs(windowSize)
            for i in xrange(0,len(freqs)):
                for j in xrange(0,len(freqs[i])):
                    r = r2(freqs[i][j])
                    if r !=	-1:
                        bin = int((fsnpd[k].positions[j+i+1]-fsnpd[k].positions[i])*nbins/(windowSize+0.01))
                        #print fsnpd[k].positions[j+i+1], fsnpd[k].positions[i]
			#print bin						
			counts[k][bin] = counts[k][bin]+1
                        sums[k][bin] = sums[k][bin] +	r
						
    for k in [1,2,3,4,5,6,7,8]:
        for i in xrange(0,nbins):
            if counts[k][i] != 0:
                sums["mean"][i]	= sums["mean"][i]+sums[k][i] 
                sums[k][i] = float(sums[k][i])/float(counts[k][i])
                counts["tot"][i] = counts["tot"][i]+counts[k][i]
    for i in xrange(0,nbins):
        if counts["tot"][i] != 0:
            sums["mean"][i] = float(sums["mean"][i])/float(counts["tot"][i])
			
    return (sums)



def DlistAll(snpsdlist,windowSize,nbins=10):
    sums = {1:[0.0]*nbins,2:[0.0]*nbins,3:[0.0]*nbins,4:[0.0]*nbins,5:[0.0]*nbins,6:[0.0]*nbins,7:[0.0]*nbins,8:[0.0]*nbins,"mean":[0.0]*nbins}	    
    counts = {1:[0]*nbins,2:[0]*nbins,3:[0]*nbins,4:[0]*nbins,5:[0]*nbins,6:[0]*nbins,7:[0]*nbins,8:[0]*nbins,"tot":[0]*nbins}
    for snpsd in snpsdlist:
        fsnpd = snpsd.snpsFilter()
        for k in [1,2,3,4,5,6,7,8]:
            freqs =	fsnpd[k].calcFreqs(windowSize)
            for i in xrange(0,len(freqs)):
                for j in xrange(0,len(freqs[i])):
                    r = D(freqs[i][j])
                    if r !=	-1:
                        bin = int((fsnpd[k].positions[j+i+1]-fsnpd[k].positions[i])*nbins/(windowSize+0.01))
                        counts[k][bin] = counts[k][bin]+1
                        sums[k][bin] = sums[k][bin] + r
						
    for k in [1,2,3,4,5,6,7,8]:
        for i in xrange(0,nbins):
            if counts[k][i]	!= 0:
                sums["mean"][i]	= sums["mean"][i]+sums[k][i] 
                sums[k][i] = float(sums[k][i])/float(counts[k][i])
                counts["tot"][i] = counts["tot"][i]+counts[k][i]
    for i in xrange(0,nbins):
        if counts["tot"][i] != 0:
            sums["mean"][i]	= float(sums["mean"][i])/float(counts["tot"][i])
	
    return (sums)


def DlistAll2(snpsdlist,windowSize,nbins=10):
    sums = {0:[0.0]*nbins,1:[0.0]*nbins,2:[0.0]*nbins,3:[0.0]*nbins,4:[0.0]*nbins,"mean":[0.0]*nbins}	    
    counts = {0:[0]*nbins,1:[0]*nbins,2:[0]*nbins,3:[0]*nbins,4:[0]*nbins,"tot":[0]*nbins}
    for snpsd in snpsdlist:
        fsnpdata = snpsd.snpsFilter()
        for k in [0,1,2,3,4]:
            fsnpsd1 = fsnpdata[k]
            fsnpsd2 = fsnpdata[9-k]
            for fsnpd in [fsnpsd1,fsnpsd2]:
                freqs =	fsnpd.calcFreqs(windowSize)
                for i in xrange(0,len(freqs)):
                    for j in xrange(0,len(freqs[i])):
                        r = D(freqs[i][j])
                        if r !=	-1:
                            bin = int((fsnpd.positions[j+i+1]-fsnpd.positions[i])*nbins/(windowSize+0.01))
                            counts[k][bin] = counts[k][bin]+1
                            sums[k][bin] = sums[k][bin] + r
			
						
    for k in [0,1,2,3,4]:
        for i in xrange(0,nbins):
            if counts[k][i] != 0:
                sums["mean"][i]	= sums["mean"][i]+sums[k][i] 
                sums[k][i] = float(sums[k][i])/float(counts[k][i])
                counts["tot"][i] = counts["tot"][i]+counts[k][i]
    for i in xrange(0,nbins):
        if counts["tot"][i] != 0:
            sums["mean"][i]	= float(sums["mean"][i])/float(counts["tot"][i])
	
    return (sums)

def r2listAll2(snpsdlist,windowSize,nbins=10):
    sums = {0:[0.0]*nbins,1:[0.0]*nbins,2:[0.0]*nbins,3:[0.0]*nbins,4:[0.0]*nbins,"mean":[0.0]*nbins}	    
    counts = {0:[0]*nbins,1:[0]*nbins,2:[0]*nbins,3:[0]*nbins,4:[0]*nbins,"tot":[0]*nbins}
    for snpsd in snpsdlist:
        fsnpdata = snpsd.snpsFilter()
        for k in [0,1,2,3,4]:
            fsnpsd1 = fsnpdata[k]
            fsnpsd2 = fsnpdata[9-k]
            for fsnpd in [fsnpsd1,fsnpsd2]:
                freqs =	fsnpd.calcFreqs(windowSize)
                for i in xrange(0,len(freqs)):
                    for j in xrange(0,len(freqs[i])):
                        r = r2(freqs[i][j])
                        if r !=	-1:
                            bin = int((fsnpd.positions[j+i+1]-fsnpd.positions[i])*nbins/(windowSize+0.01))
                            counts[k][bin] = counts[k][bin]+1
                            sums[k][bin] = sums[k][bin] + r
			
						
    for k in [0,1,2,3,4]:
        for i in xrange(0,nbins):
            if counts[k][i] != 0:
                sums["mean"][i]	= sums["mean"][i]+sums[k][i] 
                sums[k][i] = float(sums[k][i])/float(counts[k][i])
                counts["tot"][i] = counts["tot"][i]+counts[k][i]
    for i in xrange(0,nbins):
        if counts["tot"][i] != 0:
            sums["mean"][i]	= float(sums["mean"][i])/float(counts["tot"][i])
	
    return (sums)


def writeRFormat(plotData, list = [0,1,2,3,4,"mean"],windowSize=20000,ylab="",lab=""):
    delta = float(windowSize)/float(len(plotData[list[0]]))
    st = "xv <- seq("+str(delta/2.0)+","+str(windowSize)+","+str(delta)+")\n"
    maxlim = 0.0
    minlim = 1.0
    for k in list:
        data = plotData[k]
        st = st+"d"+str(k)+" <- c("
        for i in range(0,len(data)-1):
            st = st + str(data[i])+", "
            maxlim = max(maxlim,data[i])
            minlim = min(minlim,data[i])
        st = st+str(data[len(data)-1])+")\n"
        maxlim = max(maxlim,data[len(data)-1])
        minlim = min(minlim,data[len(data)-1])
    i = 1
    for k in list[0:len(list)-1]:
        st = st+"plot(y=d"+str(k)+", x=xv, ylim = c("+str(minlim)+","+str(maxlim)+"), col = "+str(i)+', type = "b",ylab="",xlab="")\npar(new=T)\n'
        i = i + 1	
    st = st+"plot(y=d"+str(list[len(list)-1])+", x=xv, ylim = c("+str(minlim)+","+str(maxlim)+'), type = "b",ylab="'+ylab+'",xlab="Bases", col = '+str(i)+', main="'+lab+'")\n'
    return st


if __name__ == "__main__":
    pass 
