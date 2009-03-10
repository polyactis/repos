import pdb

# Phenotype categories: (category,order)
phenotypeCategories = {
				1: (1,1), 2: (1,2), 3: (1,3), 4: (1,4), 39: (1,5), 40: (1,6), 41: (1,7), 42: (1,8), 
				5: (1,9), 6: (1,10), 7: (1,11), 80: (1,12), 81: (1,13), 82: (1,14),	45: (1,15), 
				46: (1,16), 47: (1,17), 48: (1,18), 57: (1,19), 59: (1,20), 58: (1,21), 43: (1,22), 44: (1,23),
					
				60: (2,1), 61: (2,2), 62: (2,3), 63: (2,4), 64: (2,5), 8: (2,6), 161: (2,7), 162: (2,8), 163: (2,9), 182: (2,10),
				164: (2,11), 165: (2,12), 166: (2,13), 173: (2,14), 174: (2,15), 175: (2,16), 176: (2,17), 177: (2,18), 178: (2,19), 179: (2,19),
				158: (2,19), 159: (2,19), 75: (2,20), 76: (2,21), 77: (2,22), 78: (2,23), 79: (2,24), 
				167: (2,25), 180: (2,28), 181: (2,29), 170: (2,30), 171: (2,31), 172: (2,32), 183: (2,33), 184: (2,34),
					
				32: (3,1), 36: (3,2), 33: (3,3), 35: (3,4), 37: (3,5), 34: (3,6), 38: (3,7),
				9: (3,8), 10: (3,9), 11: (3,10), 12: (3,11), 13: (3,12),
				65: (3,13), 67: (3,14), 69: (3,15), 71: (3,16), 73: (3,17), 
				66: (3,18), 68: (3,19), 70: (3,20), 72: (3,21), 74: (3,22),
				186: (3,23), 185: (3,24),
					
				14: (4,1), 15: (4,2), 16: (4,3), 17: (4,4), 18: (4,5), 19: (4,6), 20: (4,7), 21: (4,8), 22: (4,9), 23: (4,10), 
				24: (4,11), 25: (4,12), 26: (4,13), 27: (4,14), 28: (4,15), 29: (4,16), 30: (4,17), 31: (4,18)
			}

categories_2_phenotypes = {1:[1,2,3,4,39,40,41,42,5,6,7,80,81,82,45,46,47,48,57,59,58,43,44],
				2:[60,61,62,63,64,8,161,162,163,182,164,165,166,173,174,175,176,177,178,179,158,159,75,76,77,78,79,167,180,181,170,171,172],
				3:[32,36,33,35,37,34,38,9,10,11,12,13,65,67,69,71,73,66,68,70,72,74,186,185,183,184],
				4:[14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]}


#Log 02/25/09 - Chlor_16,Chlor_22 were removed from cat. 2: 168,169, { 168: (2,26), 169: (2,27),} 
 
 
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

	def getVariance(self, phenId):
		import util
		vals = self.getPhenVals(phenId)
		variance = util.calcVar(vals)
		return variance
		
	def getMinValue(self, phenId):
		import util
		vals = self.getPhenVals(phenId)
		return min(vals)
		

	def getPhenVals(self,phenId,asString=False, noNAs = True):
		p_i = self.getPhenIndex(phenId)
		vals = [] 
		for phenVals in self.phenotypeValues:
			if not (noNAs and phenVals[p_i]=='NA'):
				if asString or phenVals[p_i]=='NA':
					vals.append(phenVals[p_i])
				else:
					vals.append(float(phenVals[p_i]))
		return vals

	def getAccessionsWithValues(self,phenId):
		p_i = self.getPhenIndex(phenId)
		accessions = [] 
		for i in range(0,len(self.phenotypeValues)):
			phenVals = self.phenotypeValues[i]
			acc = self.accessions[i]
			if phenVals[p_i]!='NA':
				accessions.append(acc)
		return accessions


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


	def removePhenotypeIDs(self, idsToKeep):
		"""
		Removes phenotypes from the data.
		"""
		idMap = self._getIndexMapping_()
		indicesToKeep = []
		for p_i in idsToKeep:
			indicesToKeep.append(idMap[p_i])
		self.removePhenotypes(indicesToKeep)

	def naOutliers(self, phenotypeID,iqrThreshold=10.0):
		"""
		Removes outliers from the data
		"""
		numNA = 0
		idMap = self._getIndexMapping_()
		#pdb.set_trace()
		vals = zip(*self.phenotypeValues[:])[idMap[phenotypeID]]
		vals = list(vals)
		indices = []
		values = []
		for i in range(0,len(vals)):
			if vals[i]!="NA":
				indices.append(i)
				values.append(float(vals[i]))
		import util
		quantiles = util.calcQuantiles(values[:])
		print "Quantiles:",quantiles
		median = quantiles[1]
		iqr = abs(quantiles[2]-quantiles[0])
		print iqr
		for i in range(0,len(values)):
			if abs(values[i]-median)>iqrThreshold*iqr:
				print "removed",values[i]
				vals[indices[i]]= "NA"
				numNA += 1
		for i in range(0,len(self.accessions)):
			#print vals[i],idMap,phenotypeID,i
			#pdb.set_trace()
			#print self.phenotypeValues
			self.phenotypeValues[i][idMap[phenotypeID]]=vals[i]
		print "NAed",numNA,"values."		

	def getPhenotypeName(self, phenotypeIndex,rawName=False):
		indexMap = self._getIndexMapping_()
		if indexMap.has_key(phenotypeIndex):
			phenName = self.phenotypeNames[indexMap[phenotypeIndex]]
			if not rawName:
				phenName = phenName.replace("/","_div_")
				phenName = phenName.replace("*","_star_")
			return phenName
		else:
			print "Phenotype with ID",phenotypeIndex,"not found"
			return False

	def logTransform(self, pIndex):
		indexMap = self._getIndexMapping_() #Modified 9/26/08
		phenotypeIndex = indexMap[pIndex]

		if not self.isBinary(pIndex) and not self._lessOrEqualZero_(phenotypeIndex):
			import math
			for i in range(0,len(self.accessions)):
				if self.phenotypeValues[i][phenotypeIndex] !='NA':
					self.phenotypeValues[i][phenotypeIndex] = str(math.log(float(self.phenotypeValues[i][phenotypeIndex])))
		else:
			print "Can't log-transform, since phenotype is binary OR values are out of logarithm range!"
			return False
		return True


	def transformToRanks(self, pIndex):
		"""
		Transformes the phenotypic values to ranks.
		"""
		indexMap = self._getIndexMapping_() 
		phenotypeIndex = indexMap[pIndex]

		values = []
		for i in range(0,len(self.accessions)):
			if self.phenotypeValues[i][phenotypeIndex] !='NA':
				values.append((self.phenotypeValues[i][phenotypeIndex],i))

		values.sort()
		ranks = []
		for i in range(0,len(values)):
			vals = values[i]
			ranks.append(vals[1])
		

		j = 1
		for r_i in ranks:
			self.phenotypeValues[r_i][phenotypeIndex] = float(j)
			j += 1
			
		


	def _lessOrEqualZero_(self,phenotypeIndex):
		lessOrEqualZero = False
		minVal = None
		for i in range(0,len(self.accessions)):
			if self.phenotypeValues[i][phenotypeIndex] != 'NA' and float(self.phenotypeValues[i][phenotypeIndex]) <= 0.0:
				lessOrEqualZero = True
				minVal = float(self.phenotypeValues[i][phenotypeIndex])
				break
			if self.phenotypeValues[i][phenotypeIndex] != 'NA' and (minVal==None or float(self.phenotypeValues[i][phenotypeIndex])<minVal):
				minVal = float(self.phenotypeValues[i][phenotypeIndex])
		print "Minimum value =",minVal
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

	def countValues(self, phenotypeIndex):
		indexMap = self._getIndexMapping_() #Modified 9/26/08
		phenotypeIndex = indexMap[phenotypeIndex]

		valCount = 0
		for i in range(0,len(self.accessions)):
			val = self.phenotypeValues[i][phenotypeIndex]
			if val != 'NA':
				valCount += 1
		return valCount
		
	def countPhenotypes(self):
		return len(self.phenotypeNames)
		


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


	def addConstant(self, phenotypeID, constant):
		i = self.getPhenIndex(phenotypeID)
		for j in range(0,len(self.accessions)):
			if self.phenotypeValues[j][i] != "NA":
				self.phenotypeValues[j][i] = str(float(self.phenotypeValues[j][i])+constant) 
		
		
	def addSDscaledConstant(self,p_i,scale = 0.1):
		import math
		addConstant = math.sqrt(self.getVariance(p_i))*scale
		addConstant = addConstant - self.getMinValue(p_i)			
		print "Adding a constant to phenotype",p_i,":",addConstant
		self.addConstant(p_i,addConstant)

	def negateValues(self, phenotypeID):
		i = self.getPhenIndex(phenotypeID)
		for j in range(0,len(self.accessions)):
			self.phenotypeValues[j][i] = str(-float(self.phenotypeValues[j][i])) 
		

	def removeAccessionsNotInSNPsData(self,snpsd):
		indicesToKeep = []
		for i in range(0,len(self.accessions)):
			acc = self.accessions[i]
			if acc in snpsd.accessions:
				indicesToKeep.append(i)
		self.removeAccessions(indicesToKeep)

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


	def filterAccessions(self, accessionsToKeep):
		"""
		Filter accessions in the data.
		"""
		disregardedEcotypes = set(self.accessions)
 		disregardedEcotypes = disregardedEcotypes.difference(set(accessionsToKeep))

		indicesToKeep = []
		for i in range(0,len(self.accessions)):
		 	acc = self.accessions[i]
		 	if acc in accessionsToKeep:
		 		indicesToKeep.append(i)
		self.removeAccessions(indicesToKeep)

		return disregardedEcotypes
		

	def getNonNAEcotypes(self,phenotypeID):
		i = self.getPhenIndex(phenotypeID)
		ecotypes = []
		for j in range(0,len(self.accessions)):
			if self.phenotypeValues[j][i]!='NA':
				ecotypes.append(self.accessions[j])
		return ecotypes
			


		
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


def readPhenotypeFile(filename, delimiter='\t', missingVal='NA', accessionDecoder=None, type=1):
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


# A list of publishable phenotypes and their group.
_publishablePhenotypes_ = []  



def getPhenotypes(host="papaya.usc.edu", user=None, passwd=None, onlyBinary=False, onlyQuantitative=False, onlyCategorical=False, onlyReplicates=False, includeSD=False, rawPhenotypes=False, onlyPublishable=False):
	print "onlyPublishable:",onlyPublishable
	import dataParsers
	e2a = dataParsers.getEcotypeToAccessionDictionary(host,user=user,passwd=passwd,defaultValue='100')

	import MySQLdb
	print "Connecting to db, host="+host
	try:
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

		#"at.accession2tg_ecotypeid eva "
		
		#Retrieve the ecotypes
	print "Fetching ecotypes"
	if onlyPublishable:
		numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei where ei.ecotypeid=pa.ecotype_id and pa.ready_for_publication=1 order by ei.tg_ecotypeid"))
	else:
		numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei where ei.ecotypeid=pa.ecotype_id order by ei.tg_ecotypeid"))
	
	
	ecotypes = []
	accessions = []
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(str(int(row[0])))
		accessions.append(row[1])
	print len(ecotypes), "ecotypes (accessions) were found."

	#Get the phenotypic values
	phenotypeValues = []
	phenotypeNames = []
	valueColumn = "pa.transformed_value"
	if rawPhenotypes:
		print "Fetching raw phenotypic values"
		valueColumn = "pa.value"
	else:
		print "Fetching phenotypic values"
		
	if onlyPublishable:
		numRows = int(cursor.execute("select distinct pa.method_id, ei.tg_ecotypeid, "+valueColumn+", pm.short_name from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id and pa.ready_for_publication=1 order by pa.method_id, ei.tg_ecotypeid"))
	else:
		numRows = int(cursor.execute("select distinct pa.method_id, ei.tg_ecotypeid, "+valueColumn+", pm.short_name from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id order by pa.method_id, ei.tg_ecotypeid"))
		
	pvalues = [[] for j in range(0,len(ecotypes))]
	row = cursor.fetchone()
	currentMethod = int(row[0])
	phenName = str(int(row[0]))+"_"+row[3]
	#print phenName
	phenotypeNames.append(phenName)
	i=ecotypes.index(str(int(row[1])))
	pvalues[i]+=[float(row[2])]
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		nextMethod = int(row[0])
		if currentMethod != nextMethod:
			phenName = str(int(row[0]))+"_"+row[3]
			#print phenName
			phenotypeNames.append(phenName)
			valCount = 0
			for j in range(0,len(pvalues)):
				if pvalues[j] == []:
					pvalues[j] = "NA"
				else:
					valCount += len(pvalues[j])
					pvalues[j] = sum(pvalues[j])/float(len(pvalues[j]))
			phenotypeValues.append(pvalues)
			pvalues = [[] for j in range(0,len(ecotypes))]
			currentMethod = nextMethod
		
		i=ecotypes.index(str(int(row[1])))
		if row[2]!=None:
			pvalues[i]+=[float(row[2])]

	for j in range(0,len(pvalues)):
		if pvalues[j] == []:
			pvalues[j] = "NA"
		else:
			pvalues[j] = sum(pvalues[j])/float(len(pvalues[j]))
	phenotypeValues.append(pvalues)
	print len(phenotypeValues), "phenotypes were found."
	print len(phenotypeNames)
	
	cursor.close ()
	conn.close ()
	
	phenotypeValues = zip(*phenotypeValues)  
	phenotypeValues = map(list,phenotypeValues)

	phenDat = PhenotypeData(ecotypes, phenotypeNames, phenotypeValues, accessionNames=accessions)
	return phenDat


def _getFirst96Ecotypes_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123"):
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	print "Fetching data"
	ecotypes = []
	numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, c2010.accession_id, e2te.nativename, e2te.stockparent from at.complete_2010_strains_in_stock c2010, stock.ecotypeid2tg_ecotypeid e2te where e2te.ecotypeid=c2010.ecotypeid and e2te.stockparent = c2010.stockparent and c2010.accession_id<97 order by c2010.accession_id"))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(int(row[0]))
	cursor.close ()
	conn.close ()
	return ecotypes

def _getFirst192Ecotypes_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123"):
	return map(int,['100000', '5837', '6008', '6009', '6016', '6024', '6039', '6040', '6042', '6043', '6046', '6064', '6074', '6088', '6243', '6709', '6830', '6897', '6898', '6899', '6900', '6901', '6903', '6904', '6905', '6906', '6907', '6908', '6909', '6910', '6911', '6913', '6914', '6915', '6916', '6917', '6918', '6919', '6920', '6921', '6922', '6923', '6924', '6926', '6927', '6928', '6929', '6930', '6931', '6932', '6933', '6936', '6937', '6938', '6939', '6940', '6942', '6943', '6944', '6945', '6946', '6951', '6956', '6957', '6958', '6959', '6960', '6961', '6962', '6963', '6964', '6965', '6966', '6967', '6968', '6969', '6970', '6971', '6972', '6973', '6974', '6975', '6976', '6977', '6978', '6979', '6980', '6981', '6982', '6983', '6984', '6985', '6988', '7000', '7014', '7033', '7062', '7064', '7081', '7094', '7123', '7147', '7163', '7231', '7255', '7275', '7282', '7296', '7306', '7323', '7327', '7329', '7346', '7418', '7424', '7438', '7460', '7461', '7477', '7514', '7515', '7516', '7517', '7518', '7519', '7520', '7521', '7522', '7523', '7524', '7525', '7526', '8213', '8214', '8215', '8222', '8230', '8231', '8233', '8234', '8235', '8236', '8237', '8239', '8240', '8241', '8242', '8243', '8244', '8245', '8246', '8247', '8249', '8254', '8256', '8257', '8258', '8259', '8264', '8265', '8266', '8270', '8271', '8274', '8275', '8283', '8284', '8285', '8290', '8296', '8297', '8300', '8304', '8306', '8307', '8310', '8311', '8312', '8313', '8314', '8323', '8325', '8326', '8329', '8334', '8335', '8337', '8343', '8351', '8353', '8354', '8356', '8357', '8365', '8366', '8369', '8374', '8376', '8378', '8387', '8388', '8389', '8395', '8411', '8412', '8419', '8420', '8422', '8423', '8424', '8426', '8428', '8430', '9057', '9058'])
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	print "Fetching data"
	ecotypes = []
	#numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, c2010.accession_id, e2te.nativename, e2te.stockparent from at.complete_2010_strains_in_stock c2010, stock.ecotypeid2tg_ecotypeid e2te, at.ecotype_192_vs_accession_192 ea192 where e2te.ecotypeid=c2010.ecotypeid and e2te.stockparent = c2010.stockparent and ea192.accession_id =c2010.accession_id order by c2010.accession_id"))
	numRows = int(cursor.execute("select distinct e2te.tg_ecotypeid, ai2ei.accession_id, e2te.nativename, e2te.stockparent from at.accession2tg_ecotypeid ai2ei, stock.ecotypeid2tg_ecotypeid e2te where e2te.tg_ecotypeid = ai2ei.ecotype_id order by ai2ei.accession_id"))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypes.append(int(row[0]))
	cursor.close ()
	conn.close ()
	return ecotypes
"""

def _getEcotypeIdToStockParentDict_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123", defaultValue=None):
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	ecotypeDict = {}
	numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid "))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		ecotypeDict[int(row[0])] = (row[1],row[2])
		sp = "NA"
		if row[2]:
			sp = row[2]
		#print str(int(row[0]))+","+str(row[1])+","+sp
	cursor.close ()
	conn.close ()
	return ecotypeDict
	
	

def _getEcotypeIdInfoDict_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123", defaultValue=None):
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	ecotypeDict = {}
	numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename, ei.stockparent, e.latitude from stock.ecotype e, stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and e.id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid"))
	i = 0
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		latitude = None
		if row[3]:
			latitude = float(row[3])
		ecotypeDict[int(row[0])] = (row[1],row[2],latitude)
		sp = "NA"
		if row[2]:
			sp = row[2]
		#print str(int(row[0]))+","+str(row[1])+","+sp
	cursor.close()
	conn.close()
	return ecotypeDict
	


def _getEcotype2TgEcotypeDict_(host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123", defaultValue=None):
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	eDict = {}
	numRows = int(cursor.execute("select distinct ecotypeid, tg_ecotypeid from stock.ecotypeid2tg_ecotypeid "))
	while(1):
		row = cursor.fetchone()
		if row:
			eDict[int(row[0])] = int(row[1])
		else:
			break
	cursor.close ()
	conn.close ()
	return eDict



def _getAccessionToEcotypeIdDict_(accessions,stockParents=None,host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123", defaultValue=None):
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Fetching data"
	accDict = {}
	for i in range(0,len(accessions)):
		acc = accessions[i]
		if stockParents!=None and stockParents[i]!="NA":
			sp = stockParents[i]
			sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and ei.nativename like '"+str(acc)+"' and ei.stockparent like '"+str(sp)+"'"
		else:
			sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and ei.nativename like '"+str(acc)+"'"
		#sql_statement = "select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ei.nativename like '"+str(acc)+"'"
		#print sql_statement
		numRows = int(cursor.execute(sql_statement))
		i = 0
		while(1):
			row = cursor.fetchone()
			#print row
			if not row:
				if i==0:
					print "Accession",acc,"wasn't found."
				if i>1:
					print i, acc
				break;
			else: 
				i += 1
			if not accDict.has_key(acc.lower()):
				accDict[acc.lower()] = int(row[0])
			elif  accDict[acc.lower()]>int(row[0]):
				accDict[acc.lower()] = int(row[0])
			#print acc.lower(),":", row[1], int(row[0])
	cursor.close ()
	conn.close ()
	return accDict



#def _insertPhenotypesIntoDb_(host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123"):
#	name = "male"
#	raw_phen_dir = "/Network/Data/250k/dataFreeze_080608/raw_phenotypes/"
#	filename = raw_phen_dir+"FT_ALL_AC_16.csv"
#	f = open(filename,"r")
#	lines = f.readlines()
#	print lines
#	pNames = lines[0].split(",")[1:]
#	pNames[-1] = pNames[-1].strip()
#	phenValues = []
#	accessions = []
#	for i in range(1,len(lines)):
#		line = lines[i]
#		values = lines[i].split(",")
#		accessions.append(values[0])
#		values = values[1:]
#		for i in range(0,len(values)):
#			values[i] = float(values[i])
#		phenValues.append(values)
#
#	print accessions
#	#Retrieve acc_2_ei mapping 
#	accDict = _getAccessionToEcotypeIdDict_(accessions)
#	print accDict
#
#
#	#Insert data into db.
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
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Inserting data"
#	p_i = 6
#	m_i = 209
#	for i in range(0,len(accessions)):
#		val = phenValues[i][p_i]
#		e_i = accDict[accessions[i]]
#		numRows = int(cursor.execute("INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( "+str(e_i)+", "+str(val)+", 0, "+str(m_i)+", "+str(val)+" )"))
#		row = cursor.fetchone()
#		print row
#	conn.commit()
#	cursor.close ()
#	conn.close ()


#def _insertPhenotypesIntoDb_(host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123"):
#	name = "male"
#	raw_phen_dir = "/Network/Data/250k/dataFreeze_080608/raw_phenotypes/"
#	filename = raw_phen_dir+"FT_ALL_AC_16.csv"
#	f = open(filename,"r")
#	lines = f.readlines()
#	print lines
#	phenValues = []
#	accessions = []
#	ecotypes = []
#	for i in range(1,len(lines)):
#		line = lines[i].strip()
#		values = line.split(",")
#		print values
#		if values[2]!='na' and values[2]!='NA' :
#			accessions.append(values[0])
#			ecotypes.append(int(values[1]))
#			phenValues.append(float(values[2]))
#
#
#	print accessions
#	#Retrieve acc_2_ei mapping 
#	accDict = _getAccessionToEcotypeIdDict_(accessions)
#	eDict = _getEcotype2TgEcotypeDict_()
#	
#	#Verify accessions:
#	discrepancy=False
#	for i in range(0,len(accessions)):
#		if accDict.has_key(accessions[i].lower()):
#			if accDict[accessions[i].lower()]!=eDict[ecotypes[i]]:
#				print "Discrepancies!!!!",accessions[i],accDict[accessions[i].lower()],ecotypes[i],eDict[ecotypes[i]],phenValues[i]
#				accDict[accessions[i].lower()] = ecotypes[i]  #Nasty fix.
#				discrepancy = True
#		else:
#			print "Key not found",accessions[i],ecotypes[i],eDict[ecotypes[i]]
#	#if discrepancy:
#	#	import pdb;pdb.set_trace()
#
#	#Insert data into db.
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
#		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
#	except MySQLdb.Error, e:
#		print "Error %d: %s" % (e.args[0], e.args[1])
#		sys.exit (1)
#	cursor = conn.cursor ()
#	#Retrieve the filenames
#	print "Inserting data"
#	
#	m_i = 6  #Phenotype ID
#	for i in range(0,len(accessions)):
#		val = phenValues[i]
#		if accDict.has_key(accessions[i].lower()):
#			e_i = accDict[accessions[i].lower()]
#		else:
#			e_i = ecotypes[i]
#		sql_statement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( "+str(e_i)+", "+str(val)+", 0, "+str(m_i)+", "+str(val)+" )"
#		print sql_statement
#		numRows = int(cursor.execute(sql_statement))
#		row = cursor.fetchone()
#		print row
#	conn.commit()
#	cursor.close ()
#	conn.close ()
	
def _simpleInsertPhenotypesIntoDb_(host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123"):
	ion_number=19
	method_description = "Ratio: Delta_wet/Delta_dry. (Thomas Juenger)"
	
	p_i = 3  #FIXME
	m_i = 262  #FIXME
	data_description = ""
	raw_phen_dir = "/Network/Data/250k/tmp-bvilhjal/thomas_juenger/"
	filename = raw_phen_dir+"Delta_Magnus_Set_DEC08.csv"
	f = open(filename,"r")
	lines = f.readlines()
	#lines = lines[0].split("\r")
	print lines
	phenValues = []
	accessions = []
	stockParents = []
	line = lines[0].strip()
	values = line.split(",")
	phenName ="Delta_plasticity"
	for i in range(1,len(lines)):
		line = lines[i].strip()
		values = line.split(",")
		#print values
		if values[p_i]!='' :
			accName = values[0].strip()
			accessions.append(accName)
			phenValues.append(float(values[p_i].strip()))
		

	print accessions
	print phenValues
	

	#Insert data into db.
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Inserting data"
	
	print phenName
	
	sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, biology_category_id, method_description, data_description, data_type) VALUES ("+str(m_i)+", '"+phenName+"', true, 3, '"+method_description+"', '"+data_description+"', 'quantitative')"
	numRows = int(cursor.execute(sql_statement))
	row = cursor.fetchone()
	if row:
		print row
	
	for i in range(0,len(accessions)):
		val = phenValues[i]
		e_i = accessions[i]
		sql_statement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( "+str(e_i)+", "+str(val)+", 0, "+str(m_i)+", "+str(val)+" )"
		#print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			print row
	conn.commit()
	cursor.close ()
	conn.close ()
		


def _insertPhenotypesIntoDb_(host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123"):
	ion_number=19
	method_description = "Cadmium concentrations in leaves, grown in soil, version 3, obtained at Purdue"
	
	p_i = 3+ion_number  #FIXME
	m_i = 223+ion_number  #FIXME
	data_description = ""
	raw_phen_dir = "/Network/Data/250k/dataFreeze_080608/raw_phenotypes/"
	filename = raw_phen_dir+"Full_pop_matchedLines.csv"
	f = open(filename,"r")
	lines = f.readlines()
	lines = lines[0].split("\r")
	#print lines
	phenValues = []
	accessions = []
	stockParents = []
	line = lines[0].strip()
	values = line.split(",")
	phenName = values[p_i]+"_Soil_3"
	for i in range(1,len(lines)):
		line = lines[i].strip()
		values = line.split(",")
		#print values
		if values[p_i]!='na' and values[p_i]!='NA' :
			accName = values[0].strip()
			accessions.append(accName)
			stockParents.append(values[1].strip())
			phenValues.append(float(values[p_i].strip()))


	print accessions
	#Retrieve acc_2_ei mapping 
	accDict = _getAccessionToEcotypeIdDict_(accessions,stockParents)

	#Insert data into db.
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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()
	#Retrieve the filenames
	print "Inserting data"
	
	print phenName
	
	sql_statement = "INSERT INTO stock_250k.phenotype_method  (id, short_name, only_first_96, biology_category_id, method_description, data_description, data_type) VALUES ("+str(m_i)+", '"+phenName+"', true, 3, '"+method_description+"', '"+data_description+"', 'quantitative')"
	numRows = int(cursor.execute(sql_statement))
	row = cursor.fetchone()
	if row:
		print row
	
	for i in range(0,len(accessions)):
		val = phenValues[i]
		if accDict.has_key(accessions[i].lower()):
			e_i = accDict[accessions[i].lower()]
		else:
			e_i = ecotypes[i]
		sql_statement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( "+str(e_i)+", "+str(val)+", 0, "+str(m_i)+", "+str(val)+" )"
		#print sql_statement
		numRows = int(cursor.execute(sql_statement))
		row = cursor.fetchone()
		if row:
			print row
	conn.commit()
	cursor.close ()
	conn.close ()
		

def _createRatioPhenotype_(pi1,pi2,methodID,short_name,method_description,data_description="",data_type="quantitative",created_by="bvilhjal",comment="",biology_category_id=1,only_first_96=0,readyForPublication = 0,user = "bvilhjal",host="papaya.usc.edu",passwd="bamboo123"):
	phed = getPhenotypes(user=user, passwd=passwd,rawPhenotypes=True) 
	pi1 = phed.getPhenIndex(pi1)
	pi2 = phed.getPhenIndex(pi2)
	print len(phed.accessions)
	
	newVals = []
	for i in range(0,len(phed.accessions)):
		newVal = 'NA'
		if phed.phenotypeValues[i][pi1] != 'NA' and phed.phenotypeValues[i][pi2]!='NA':
			v1 = float(phed.phenotypeValues[i][pi1])
			v2 = float(phed.phenotypeValues[i][pi2])
			#print v1, v2
			if v2 != 0.0:
				newVal = v1/v2
		newVals.append(newVal)

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
		conn = MySQLdb.connect (host = host, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor ()

	sqlStat = """INSERT INTO stock_250k.phenotype_method 
    (id,short_name,only_first_96,biology_category_id,method_description, data_description, comment, created_by, data_type) 
    VALUES """
	sqlStat += "("+str(methodID)+",'"+short_name+"',"+str(only_first_96)+","+str(biology_category_id)+",'"+method_description+"','"+data_description+"','"+comment+"','"+created_by+"','"+data_type+"')"	
	#print sqlStat	
	#numRows = int(cursor.execute(sqlStat))
	#print "Inserted data into stock_250k.phenotype_method:",numRows
	#row = cursor.fetchone()
	#if row:
	#	print row
	
	for i in range(0,len(phed.accessions)):
		
		val = newVals[i]
		if val !='NA':
			e_i = phed.accessions[i]
			sqlStatement = "INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id) VALUES ( "+str(e_i)+", "+str(val)+", "+str(readyForPublication)+", "+str(methodID)+")"
			#print sqlStatement
			numRows = int(cursor.execute(sqlStatement))
			row = cursor.fetchone()
			if row:
				print row
	print "Done inserting data into table!"
	conn.commit()
	cursor.close ()
	conn.close ()


	

#		
#def _runTest_():
#	filename = "/Network/Data/250k/dataFreeze_080608/phenotypes_all_raw_120308.tsv"
#	phed = readPhenotypeFile(filename)
#	print	phed._getIndexMapping_() #Modified 9/26/08
#	p_i = 39
#	print "phentoype:",p_i
#	print "Number of non-zero values:",phed.countValues(p_i)
#	phenIndex = phed.getPhenIndex(p_i)
#	vals = [] 
#	accessions = []
#	for i in range(0,len(phed.phenotypeValues)):
#		phenVals = phed.phenotypeValues[i]
#		if phenVals[phenIndex]!='NA':
#			vals.append(phenVals[phenIndex])
#			accessions.append(phed.accessions[i])
#	
#	print "len(phenVals):",len(vals)
#	print "len(accessions):",len(accessions)
#	print "len(set(accessions)):",len(set(accessions))
#	print vals
#	print accessions
#	
	
def _runTest_():
	filename = "/Network/Data/250k/dataFreeze_011209/250K_f13_012509.csv"
	import dataParsers,snpsdata
	snpsds = dataParsers.parseCSVData(filename, format=1, deliminator=",")#,debug=True)
	snpsd = snpsdata.SNPsDataSet(snpsds,[1,2,3,4,5])
	eDict = _getEcotypeIdToStockParentDict_()
	accessions = map(int,snpsd.accessions)
	accessions.sort()
	print "ecotype_id, native_name, stock_parent"
	i = 0
	for et in accessions:
		et = int(et)
		print str(et)+", "+str(eDict[et][0])+", "+str(eDict[et][1])
	
		
		
	

if __name__ == '__main__':
	#_createRatioPhenotype_(184,183,222,"Trich_avg_JA_div_Trich_avg_C","Ratio: Trich_avg_JA/Trich_avg_C",biology_category_id=7)
	#pass
	#_getEcotypeIdToStockParentDict_()
	#eDict = _getEcotype2TgEcotypeDict_()
	#_insertPhenotypesIntoDb_()
	#_simpleInsertPhenotypesIntoDb_()
	_runTest_()
