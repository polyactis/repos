import pdb

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
		phenName = self.phenotypeNames[indexMap[phenotypeIndex]]
		if not rawName:
			phenName = phenName.replace("/","_div_")
			phenName = phenName.replace("*","_star_")
		return phenName


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
			self.phenotypeValues[j][i] = str(float(self.phenotypeValues[j][i])+constant) 
		

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
	
	cursor.close ()
	conn.close ()
	
	phenotypeValues = zip(*phenotypeValues)  
	phenotypeValues = map(list,phenotypeValues)

	phenDat = PhenotypeData(ecotypes, phenotypeNames, phenotypeValues, accessionNames=accessions)
	return phenDat



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
		print str(int(row[0]))+","+str(row[1])+","+sp
	cursor.close ()
	conn.close ()
	return ecotypeDict
	



def _getAccessionToEcotypeIdDict_(accessions,host="papaya.usc.edu", user="bvilhjal", passwd="bamboo123", defaultValue=None):
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
	for acc in accessions:
		numRows = int(cursor.execute("select distinct ei.tg_ecotypeid, ei.nativename from stock.ecotypeid2tg_ecotypeid ei, stock_250k.array_info ai where ai.paternal_ecotype_id=ei.tg_ecotypeid and ai.maternal_ecotype_id=ei.tg_ecotypeid and ei.nativename like '"+str(acc)+"'"))
		i = 0
		while(1):
			row = cursor.fetchone()
			if not row:
				if i==0:
					print "Accession",acc,"wasn't found."
				if i>1:
					print i
				break;
			else: 
				i += 1
			accDict[acc] = int(row[0])
			print acc,":", row[1], int(row[0])
	cursor.close ()
	conn.close ()
	return accDict



def _insertPhenotypesIntoDb_(host="papaya.usc.edu",user="bvilhjal",passwd="bamboo123"):
	name = "male"
	filename = "/Users/bjarni/Projects/gwa2008/ThomasJunger/FEMALE_DATA.csv"
	f = open(filename,"r")
	lines = f.readlines()[0].rstrip().split("\r")
	print lines
	pNames = lines[0].split(",")[1:]
	pNames[-1] = pNames[-1].strip()
	phenValues = []
	accessions = []
	for i in range(1,len(lines)):
		line = lines[i]
		values = lines[i].split(",")
		accessions.append(values[0])
		values = values[1:]
		for i in range(0,len(values)):
			values[i] = float(values[i])
		phenValues.append(values)

	print accessions
	#Retrieve acc_2_ei mapping 
	accDict = _getAccessionToEcotypeIdDict_(accessions)
	print accDict


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
	p_i = 6
	m_i = 209
	for i in range(0,len(accessions)):
		val = phenValues[i][p_i]
		e_i = accDict[accessions[i]]
		numRows = int(cursor.execute("INSERT INTO stock_250k.phenotype_avg (ecotype_id, value, ready_for_publication, method_id, transformed_value) VALUES ( "+str(e_i)+", "+str(val)+", 0, "+str(m_i)+", "+str(val)+" )"))
		row = cursor.fetchone()
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
	numRows = int(cursor.execute(sqlStat))
	print "Inserted data into stock_250k.phenotype_method:",numRows
	row = cursor.fetchone()
	if row:
		print row
	
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
	_createRatioPhenotype_(184,183,222,"Trich_avg_JA_div_Trich_avg_C","Ratio: Trich_avg_JA/Trich_avg_C",biology_category_id=7)
	pass
	#_getEcotypeIdToStockParentDict_()
	#_insertPhenotypesIntoDb_()
	#_runTest_()
