class PhenotypeData:
    """
    A class that knows how to read simple tsv or csv phenotypefiles and facilitates their interactions with SnpsData objects.    
    """
    
    def __init__(self, accessions, phenotypeNames, phenotypeValues, accessionNames=None):
        self.accessions = accessions
        self.phenotypeNames = phenotypeNames
        self.phenotypeValues=phenotypeValues # list[accession_index][phenotype_index]
        self.accessionNames = accessionNames
        self.phenIds = []
        for phenName in self.phenotypeNames:
            id = phenName.split("_")
            self.phenIds.append(int(id[0]))
        

    def _getIndexMapping_(self):
        indexMapping = dict()
        for i in range(0,len(self.phenIds)):
            indexMapping[self.phenIds[i]]=i
        return indexMapping


    def getPhenIndex(self,phenId):
        index=None
        i = 0
        while self.phenIds[i]!=phenId:
            i += 1
        if self.phenIds[i]==phenId:
            index = i
        return index
                
    
    def onlyBiologyCategory(self,phenotypeCategory,host=None,user=None,passwd=None):
        import MySQLdb

	#print "Connecting to db, host="+host
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
	#print "Fetching biological info on phenotypes."  

	numRows = int(cursor.execute("select distinct id, short_name, biology_category_id from stock_250k.phenotype_method order by id"))
    
	bioinfo = []
	currTairID=""
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		bioinfo.append([int(row[0]),row[1],int(row[2])])

	cursor.close ()
	conn.close ()
	#print "Biol. info fetched"

        phenIds = []
        for phenName in self.phenotypeNames:
            id = phenName.split("_")
            phenIds.append(int(id[0]))

        indicesToKeep = []
        j = 0
        for i in range(0, len(phenIds)):
            while phenIds[i]>bioinfo[j][0]:
                j += 1
            if phenIds[i]==bioinfo[j][0]:
                if bioinfo[j][2]==phenotypeCategory:
                    indicesToKeep.append(i)
                j += 1

        self.removePhenotypes(indicesToKeep)
        
    def removePhenotypes(self, indicesToKeep):
        """
        Removes phenotypes from the data.
        """
        numRemoved = len(self.phenotypeNames)-len(indicesToKeep)
        #print "Removing",numRemoved,"phenotypes in phenotype data, out of",len(self.phenotypeNames), "phenotypes."
        newPhenotypeNames = []
        newPhenotVals = []
        newPhenIds = []
        for j in range(0,len(self.accessions)):
            newPhenotVals.append([])
        for i in indicesToKeep:
            newPhenotypeNames.append(self.phenotypeNames[i])
            newPhenIds.append(self.phenIds[i])
            for j in range(0,len(self.accessions)):
                newPhenotVals[j].append(self.phenotypeValues[j][i])
        self.phenIds = newPhenIds
        self.phenotypeValues = newPhenotVals
        self.phenotypeNames = newPhenotypeNames

        #print "len(self.phenotypeNames):",len(self.phenotypeNames)
        #print "len(self.phenotypeValues[0]):",len(self.phenotypeValues[0])



    def getPhenotypeName(self, phenotypeIndex):
        indexMap = self._getIndexMapping_()
        return self.phenotypeNames[indexMap[phenotypeIndex]]


    def logTransform(self, pIndex):
        indexMap = self._getIndexMapping_() #Modified 9/26/08
        phenotypeIndex = indexMap[pIndex]

        if not self.isBinary(pIndex) and not self._lessOrEqualZero_(phenotypeIndex):
            import math
            for i in range(0,len(self.accessions)):
                self.phenotypeValues[i][phenotypeIndex] = str(math.log(float(self.phenotypeValues[i][phenotypeIndex])))
        else:
            print "Can't log-transform, since phenotype is binary OR values are out of logarithm range!"

    def _lessOrEqualZero_(self,phenotypeIndex):
        lessOrEqualZero = False
        for i in range(0,len(self.accessions)):
            if self.phenotypeValues[i][phenotypeIndex] <= 0:
                lessOrEqualZero = True
                break
        return lessOrEqualZero
            

    def isBinary(self, phenotypeIndex):
        indexMap = self._getIndexMapping_() #Modified 9/26/08
        phenotypeIndex = indexMap[phenotypeIndex]

        l = []
        for i in range(0,len(self.accessions)):
            val = self.phenotypeValues[i][phenotypeIndex]
            if val != 'NA':
                if not val in l:
                    l.append(val)
                
        if 1<len(l)<3:
            return True
        elif 1==len(l):
            raise Exception("Only one phenotype value")
        return False

    def orderAccessions(self, accessionMapping=None):
        """
        Orders the accession alphabetically if no mapping is given.
        """
        print "Ordering phenotype data accessions."
        newAccessions = []
        newPhenotVals = []
        for acc in self.accessions:
            newAccessions.append("")
            newPhenotVals.append([])

        if not accessionMapping:
            accessionMapping = []
            l = range(0,len(self.accessions))
            l1 = zip(self.accessions,l)
            l1.sort()
            j = 0
            for (acc,i) in l1:
                accessionMapping.append((i,j))
                j += 1

        for (i,j) in accessionMapping:
            #print j, len(newAccessions)
            newAccessions[j] = self.accessions[i]
            newPhenotVals[j] = self.phenotypeValues[i]
        self.accessions = newAccessions
        self.phenotypeValues = newPhenotVals
        
        if self.accessionNames:
            newAccessionNames = []
            for acc in self.accessions:
                newAccessionNames.append("")
            for (i,j) in accessionMapping:
                newAccessionNames[j] = self.accessionNames[i]
            self.accessionNames = newAccessionNames
                        
        
    def removeAccessions(self, indicesToKeep):
        """
        Removes accessions from the data.
        """
        numAccessionsRemoved = len(self.accessions)-len(indicesToKeep)
        print "Removing",numAccessionsRemoved,"accessions in phenotype data, out of",len(self.accessions), "accessions."
        newAccessions = []
        newPhenotVals = []
        for i in indicesToKeep:
            newAccessions.append(self.accessions[i])
            newPhenotVals.append(self.phenotypeValues[i])
        self.accessions = newAccessions
        self.phenotypeValues = newPhenotVals
        print "len(self.accessions):",len(self.accessions)
        print "len(self.phenotypeValues):",len(self.phenotypeValues)
        if self.accessionNames:
            newAccessionNames = []
            for i in indicesToKeep:
                newAccessionNames.append(self.accessionNames[i])
            self.accessionNames = newAccessionNames
            print "len(self.accessionNames):",len(self.accessionNames)

    def writeToFile(self, outputFile, phenotypes=None, delimiter=','):
        print "Writing out phenotype file:",outputFile
        outStr = "ecotype_id"
        if self.accessionNames:
            outStr += delimiter+"accession_name"
        if phenotypes:
            for i in phenotypes:
                name = self.phenotypeNames[i]
                outStr += delimiter+name
            outStr += '\n'
            for i in range(0,len(self.accessions)):
                outStr += str(self.accessions[i])
                for j in phenotypes:
                    outStr += delimiter+str(self.phenotypeValues[i][j])
                outStr +="\n"
        else:
            for name in self.phenotypeNames:
                outStr += delimiter+name
            outStr += '\n'
            for i in range(0,len(self.accessions)):
                outStr += str(self.accessions[i])
                if self.accessionNames:
                    outStr += delimiter+str(self.accessionNames[i])
                for j in range(0, len(self.phenotypeNames)):
                    outStr += delimiter+str(self.phenotypeValues[i][j])
                outStr +="\n"
                #print outStr

        f = open(outputFile,"w")
        f.write(outStr)
        f.close()

def readPhenotypeFile(filename, delimiter=',', missingVal='NA', accessionDecoder=None, type=1):
    """
    Reads a phenotype file and returns a phenotype object.
    """
    f = open(filename,"r")
    lines = f.readlines()
    f.close()
    shift = 2
    if type==2:
        shift = 1
    line = (lines[0].rstrip()).split(delimiter)
    phenotypeNames = line[shift:]
    accessions = []
    phenotypeValues = []
    for i in range(1, len(lines)):
        line = (lines[i].rstrip()).split(delimiter)
        if accessionDecoder:
            accessions.append(accessionDecoder[line[0]])
        else:
            accessions.append(line[0])
        phenotypeValues.append(line[shift:])

    return PhenotypeData(accessions,phenotypeNames,phenotypeValues)
    




def _runTest_():
    import dataParsers
    filename = "/Network/Data/250k/finalData_051808/phenotypes_052208.tsv"
    phed = readPhenotypeFile(filename,delimiter='\t')    
    #print phed.accessions
    #print phed.phenotypeNames
    #print phed.phenotypeValues
    for i in range(0,len(phed.phenotypeNames)):
        binary = phed.isBinary(i)
        print i, phed.phenotypeNames[i], binary
        
if __name__ == '__main__':
	_runTest_()

