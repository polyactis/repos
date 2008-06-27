#!/usr/bin/env python2.5
"""
Usage: outputPhenotypeData.py [OPTIONS] -o OUTPUT_FILE 

Option:

        -o ...,	output file
	-d ..., --delim=...         default is \", \"      
        -m ..., --missingval=...    default is \"NA\"
	-p ...                      Password
	-h ...                      Host (papaya.usc.edu is default)
	-u ...                      Username
	--rawPhenotypes             Outputs the raw phenotype values (otherwise it outputs the transformed ones).
	--onlyBinary                Outputs only binary phenotypes
	--onlyCategorical           Outputs only categorical phenotypes
	--onlyQuantitative          Outputs only quantitative phenotypes
	--onlyReplicates            Outputs phenotypes which have replicates, with the replicates
	--includeSD                 Include standard deviations (if available) as phenotypes.
	--orderByGenotypeFile=...   Order the accessions so that they match the accessions in the given genotype file.
	--onlyPublishable           Outputs only the publishable phentotypic values.
	-h, --help	            show this help

Examples:
	
Description:

"""

import sys, getopt, traceback
import phenotypeData, dataParsers

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["rawPhenotypes","onlyBinary", "onlyCategorical", "onlyQuantitative", "onlyReplicates", "delim=", "missingval=", "help","includeSD","orderByGenotypeFile=", "onlyPublishable"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "p:u:h:o:d:m:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	output_fname = None
	delim = ","
	missingVal = "NA"
	rawPhenotypes = False
	onlyBinary = False
	onlyCategorical = False
	onlyQuantitative = False
	onlyReplicates = False
	includeSD = False
	onlyPublishable = False
	genotypeFile = None
	help = 0
	passwd = None
	user=None
	host="papaya.usc.edu"
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-u",):
			user = arg
		elif opt in ("-p",):
			passwd = arg
		elif opt in ("-h",):
			host = arg
		elif opt in ("-m","--missingval"):
			missingVal = arg
		elif opt in ("-d","--delim"):
			delim = arg
		elif opt in ("--rawPhenotypes"):
			rawPhenotypes = True
		elif opt in ("--onlyBinary"):
			onlyBinary = True
		elif opt in ("--onlyCategorical"):
			onlyCategorical = True
		elif opt in ("--onlyQuantitative"):
			onlyQuantitative = True
		elif opt in ("--onlyReplicates"):
			onlyReplicates = True
		elif opt in ("--includeSD"):
			includeSD = True
		elif opt in ("--onlyPublishable"):
			onlyPublishable = True
		elif opt in ("--orderByGenotypeFile"):
			genotypeFile = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not output_fname:
		output_fname
		if help==0:
			print "Output file missing!!\n"
			print __doc__
		sys.exit(2)


        #Retrieve phenotype data.
	phenData = getPhenotypes(host=host,user=user,passwd=passwd,onlyBinary=onlyBinary, onlyQuantitative=onlyQuantitative, onlyCategorical=onlyCategorical, onlyReplicates=onlyReplicates, includeSD=includeSD, rawPhenotypes=rawPhenotypes, onlyPublishable=onlyPublishable)

        #Sort in correct order.
	if genotypeFile:
		snpsds = dataParsers.parseCSVData(genotypeFile, format=1, deliminator=delim, missingVal=missingVal)
		print "Removing accessions which are not in the genotype file."
		indicesToKeep = []
		for i in range(0,len(phenData.accessions)):
			if phenData.accessions[i] in snpsds[0].accessions:
				indicesToKeep.append(i)
		phenData.removeAccessions(indicesToKeep)

		print "Ordering accessions to match the genotype file order"
		associationMapping = []
		j = 0
		for acc in snpsds[0].accessions:
			if acc in phenData.accessions:
				associationMapping.append((phenData.accessions.index(acc),j))
				j += 1
	
		phenData.orderAccessions(associationMapping)
			
        #Output phenotypes to file.
	phenData.writeToFile(output_fname, delimiter='\t')

def getPhenotypes(host="papaya.usc.edu", user=None, passwd=None, onlyBinary=False, onlyQuantitative=False, onlyCategorical=False, onlyReplicates=False, includeSD=False, rawPhenotypes=False, onlyPublishable=False):
	import dataParsers
	e2a = dataParsers.getEcotypeToAccessionDictionary(host,user=user,passwd=passwd,defaultValue=100)

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

        #"at.accession2tg_ecotypeid eva "
        
        #Retrieve the ecotypes
	print "Fetching ecotypes"
	#numRows = int(cursor.execute("select distinct pa.ecotype_id, ei.nativename from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei where ei.ecotypeid=pa.ecotype_id order by pa.ecotype_id"))
	
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
	#numRows = int(cursor.execute("select distinct pa.method_id, pa.ecotype_id, pa.value, pm.short_name from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id order by pa.method_id, pa.ecotype_id"))
	
	numRows = int(cursor.execute("select distinct pa.method_id, ei.tg_ecotypeid, "+valueColumn+", pm.short_name, pm.only_first_96 from stock_250k.phenotype_avg pa, stock.ecotypeid2tg_ecotypeid ei, stock_250k.phenotype_method pm where ei.ecotypeid=pa.ecotype_id and pa.method_id=pm.id order by pa.method_id, ei.tg_ecotypeid"))
		
	onlyFirst96 = []
	pvalues = ['NA' for j in range(0,len(ecotypes))]
	row = cursor.fetchone()
	currentMethod = int(row[0])
	phenotypeNames.append(str(int(row[0]))+"_"+row[3])
	i=ecotypes.index(str(int(row[1])))
	pvalues[i]=float(row[2])
	onlyFirst96.append(bool(int(row[4])))
	while(1):
		row = cursor.fetchone()
		if not row:
			break;
		nextMethod = int(row[0])
		if currentMethod != nextMethod:
			onlyFirst96.append(bool(int(row[4])))
			phenotypeNames.append(str(int(row[0]))+"_"+row[3])
			phenotypeValues.append(pvalues)
			pvalues = ['NA' for j in range(0,len(ecotypes))]
			currentMethod = nextMethod
		
		i=ecotypes.index(str(int(row[1])))
		if row[2]!=None:
			pvalues[i]=float(row[2])

	phenotypeValues.append(pvalues)
	print len(phenotypeValues), "phenotypes were found."
    
	cursor.close ()
	conn.close ()
	import util
	phenotypeValues = util.transposeDoubleLists(phenotypeValues)

	if onlyPublishable:  #Censor the data.
		for i in range(0,len(phenotypeValues)):
			for j in range(0,len(ecotypes)):
				if onlyFirst96[i]:
					if e2a[ecotypes[j]][0]>97:
						phenotypeValues[i][j]='NA'

	phenDat = phenotypeData.PhenotypeData(ecotypes, phenotypeNames, phenotypeValues, accessionNames=accessions)
	return phenDat
	


def _transformPhenotypes_():
	"""
	Transforms some of the phenotypes.
	
	This function should ideally only be executed once.
	"""
	
	transformationDict = {0:"None", 1:"Log(x)", 2:"Log(0.5+x)", 3:"Log(5+x)", 4:"(x-3)"}
	#                                                                                        74                     116                       134        152
	transformationList = 7*[1]+[3]+5*[0]+18*[1]+7*[0]+24*[1]+3*[0]+[1]+[2]+3*[0]+[2]+[0]+[2]+[0]+2*[1]+3*[0]+35*[1]+[2]+7*[1]+2*[0]+7*[1]+[0]+[2]+17*[1]+[0]+5*[1]+[0]+2*[1]+[0]+[4]+[0]+2*[1]+18*[0]+2*[1]+[0]

	#for i in range(0,len(transformationList)):
	#	print (i+1), transformationList[i]
	
	datatypeDict = {1:"quantitative",2:"ordered_categorical",3:"binary"}
	datatypeList = 8*[1]+5*[3]+18*[1]+7*[3]+38*[1]+3*[3]+52*[1]+[2]+24*[1]+2*[2]+[1]+[2]+[3]+[2]+3*[1]+6*[3]+3*[2]+3*[3]+3*[2]+5*[1]

	
	phenData = getPhenotypes(user="bvilhjal", passwd="bamboo123")
	print phenData.phenotypeNames
        hostname="papaya.usc.edu"
	user="bvilhjal"
	passwd="bamboo123"
        #Connect to DB
	import MySQLdb
	print "Connecting to db, host="+hostname
	if not user:
		import sys
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()
	try:
		conn = MySQLdb.connect (host = hostname, user = user, passwd = passwd, db = "at")
	except MySQLdb.Error, e:
		print "Error %d: %s" % (e.args[0], e.args[1])
		sys.exit (1)
	cursor = conn.cursor()
	print "Inserting results into DB."
	for i in range(0,len(phenData.phenotypeNames)):
		for j in range(0,len(phenData.accessions)):
			transFormedValue = _transformValue_(phenData.phenotypeValues[j][i],transformationList[i])
			methodID = phenData.phenotypeNames[i].split("_")[0]
			acc=phenData.accessions[j]
			#Insert info
			print phenData.phenotypeValues[j][i], transFormedValue, methodID, acc
			sqlStatm = "UPDATE stock_250k.phenotype_avg SET transformed_value="+str(transFormedValue)+" WHERE ecotype_id="+acc+" and method_id="+methodID+""
			numRows = int(cursor.execute(sqlStatm))	
			conn.commit()
	for i in range(0,len(phenData.phenotypeNames)):
		methodID = phenData.phenotypeNames[i].split("_")[0]
		sqlStatm = "UPDATE stock_250k.phenotype_method SET data_type='"+datatypeDict[datatypeList[i]]+"' WHERE id="+methodID+""
		print sqlStatm
		numRows = int(cursor.execute(sqlStatm))	
		sqlStatm = "UPDATE stock_250k.phenotype_method SET transformation_description='"+transformationDict[transformationList[i]]+"' WHERE id="+methodID+""
		print sqlStatm
		numRows = int(cursor.execute(sqlStatm))	
		conn.commit()

	print "Committing transaction (making changes permament)."
	print "Closing connection."
	#Close connection
	cursor.close()
	conn.close()


def _transformValue_(value,transformType):
	if value!="NA":
		import math
		if transformType == 1:
			return math.log(value)
		elif transformType == 2:
			return math.log(0.5+value)
		elif transformType == 3:
			return math.log(5+value)
		elif transformType == 4:
			return value-3
	else:
		value='NULL'
	return value

def _test1_():
	pass


if __name__ == '__main__':
	#_test1_()
	#_transformPhenotypes_()
	_run_()


