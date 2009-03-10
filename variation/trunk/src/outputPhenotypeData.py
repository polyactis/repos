#!/usr/bin/env python2.5
"""
Usage: outputPhenotypeData.py [OPTIONS] -o OUTPUT_FILE 

Option:

	-o ...,	output file
	-d ..., --delim=...		  default is \", \"	  
	-m ..., --missingval=...  default is \"NA\"
	-p ...					  Password
	-h ...					  Host (papaya.usc.edu is default)
	-u ...					  Username
	--rawPhenotypes			  Outputs the raw phenotype values (otherwise it outputs the transformed ones).
	--phenotypeIDs=n1,n2,...  Outputs only the phenotypes with the given IDs.
	--onlyBinary			  Outputs only binary phenotypes
	--onlyCategorical		  Outputs only categorical phenotypes
	--onlyQuantitative		  Outputs only quantitative phenotypes
	--onlyReplicates		  Outputs phenotypes which have replicates, with the replicates
	--includeSD				  Include standard deviations (if available) as phenotypes.
	--orderByGenotypeFile=... Order the accessions so that they match the accessions in the given genotype file. 
							  Phenotypes which are not in the genotype file are not included in the phenotype file.
	--onlyPublishable		  Outputs only the publishable phentotypic values.
	-h, --help				  show this help

Examples:
	
Description:

"""
from doctest import SKIP

import sys, getopt, traceback
import phenotypeData, dataParsers

from phenotypeData import *

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["phenotypeIDs=","rawPhenotypes","onlyBinary", "onlyCategorical", "onlyQuantitative", "onlyReplicates", "delim=", "missingval=", "help","includeSD","orderByGenotypeFile=", "onlyPublishable"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "p:u:h:o:d:m:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)

	output_fname = None
	phenotypeIDs = None
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
		elif opt in ("--phenotypeIDs"):
			phenotypeIDs = []
			for num in arg.split(","):
				phenotypeIDs.append(int(num))
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


	if not user:
		sys.stdout.write("Username: ")
		user = sys.stdin.readline().rstrip()
	if not passwd:
		import getpass
		passwd = getpass.getpass()

		#Retrieve phenotype data.
	phenData = getPhenotypes(host=host,user=user,passwd=passwd,onlyBinary=onlyBinary, onlyQuantitative=onlyQuantitative, onlyCategorical=onlyCategorical, onlyReplicates=onlyReplicates, includeSD=includeSD, rawPhenotypes=rawPhenotypes, onlyPublishable=onlyPublishable)

	if phenotypeIDs:
		print "Remoing phenotypes."
		phenData.removePhenotypeIDs(phenotypeIDs)
	
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


#transformationMap = {210: (1,1), 211: (1,1), 212: (1,1), 213: (5,1), 214: (5,1), 215: (5,1)}
#transformationMap = {253: (0,1), 254: (0,1), 255: (0,1), 256: (8,1), 257: (0,1), 258: (1,1), 259: (0,1)}
transformationMap = {5: (9,1), 6: (9,1), 7: (9,1), 65: (6,1), 67: (6,1), 71: (6,1), 73: (6,1)}

#transformationMap = {1: (1,1), 2: (1,1), 3: (1,1), 4: (1,1), 5: (1,1), 6: (1,1), 7: (1,1), 
#					8: (0,1),
#					9: (0,3), 10: (0,3), 11: (0,3), 12: (0,3), 13: (0,3),
#					14: (0,1), 15: (0,1), 16: (1,1), 17: (0,1), 18: (0,1), 19: (0,1), 20: (0,1), 
#					21: (0,1), 22: (1,1), 23: (0,1), 24: (1,1), 25: (0,1), 26: (0,1), 27: (0,1), 28: (1,1), 29: (0,1), 30: (1,1), 31: (1,1),
#					32: (0,3), 33: (0,3), 34: (0,3), 35: (0,3), 36: (0,3), 37: (0,3), 38: (0,3),
#					39: (1,1), 40: (1,1), 41: (1,1), 42: (1,1), 43: (1,1), 44: (1,1), 45: (1,1), 46: (1,1), 47: (1,1), 48: (7,1),
#					57: (0,1), 58: (0,1), 59: (1,1), 60: (1,1), 61: (1,1), 62: (1,1), 
#					63: (0,1), 64: (0,1), 65: (6,1), 66: (0,1), 67: (6,1), 68: (0,1), 69: (0,1), 70: (0,1), 71: (6,1), 72: (0,1), 73: (6,1), 74: (0,1), 75: (1,1), 76: (1,1), 
#					77: (0,3), 78: (0,3), 79: (0,3), 
#					80: (1,1), 81: (1,1), 82: (1,1), 
#					158: (0,2), 159: (0,2), 160: (1,2), 161: (0,2), 162: (4,3), 163: (1,2), 164: (1,1), 165: (0,1), 166: (0,1), 167: (0,3), 168: (0,3), 169: (0,3), 170: (0,3), 
#					171: (0,3), 172: (0,3), 173: (0,2), 174: (0,2), 175: (0,2), 176: (0,3), 177: (0,3), 178: (0,3), 179: (0,2), 180: (0,2), 181: (0,2), 182: (1,1), 
#					183: (1,1), 184: (1,1), 185: (0,1), 186: (0,1)}


transformation2phenotype = {0: [8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 23, 25, 26, 27, 29, 32, 33, 34, 35, 36, 37, 
							38, 57,	58, 63, 64, 66, 68, 69, 70, 72, 74, 77, 78, 79, 158, 159, 161, 165, 166, 167, 168, 169, 
							170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 185, 186], 
							1: [1, 2, 3, 4, 5, 6, 7, 16, 22, 24, 28, 30, 31, 39, 40, 41, 42, 43, 44, 45, 46, 47, 59, 60, 61, 62,
							75, 76, 80, 81, 82, 160, 163, 164, 182, 183, 184], 
							4: [162], 
							6: [ 65, 67, 71, 73], 
							7: [48]}


datatypeDict = {1:"quantitative",2:"ordered_categorical",3:"binary"}
transformationTypes = {0:"None", 1:"Log(x)", 2:"Log(0.5+x)", 3:"Log(5+x)", 4:"(x-3)", 5:"Log(x) w. outlier cutoff = IQR*10", 6: "Log(SD/10+x)", 7:"Ranks", 8:"Log(-x-minVal+SD/10)", 9:"Log(x) w. 192 only"}


def _fakeTransformPhenotypes_():
	import os
	res_dir = "/Network/Data/250k/tmp-bvilhjal/emma_results/"
	phenotypeFile = "/Network/Data/250k/dataFreeze_011209/phenotypes_all_raw_012509.tsv"
	phed = phenotypeData.readPhenotypeFile(phenotypeFile, delimiter = '\t')
	for p_i in transformationMap: 
		(trans_type,data_type) = transformationMap[p_i]
		print trans_type
		if trans_type==1 or trans_type==6:
			t_type = "new_logTransform" #m"
		elif trans_type==9 :
			t_type = "new_logTransform_f192" #m"
		elif trans_type==5:
			t_type = "newDataset_logTransform_noOutliers"
		elif trans_type==7:
			t_type = "new_ranks"
		elif trans_type==8:
			t_type = "neg_const_logTransform"
		else: #if trans_type==0:
			t_type = "new_raw"
		phenName = phed.getPhenotypeName(p_i)
#		oldName = res_dir+"Emma_"+t_type+"_"+phenName+".pvals"
#		newName = res_dir+"Emma_new_trans_"+phenName+".pvals"
#		cp_command = "sudo cp "+oldName+" "+newName
#		print cp_command
#		os.system(cp_command)
		oldName = res_dir+"Emma_"+t_type+"_"+phenName+".sr.pvals"
		newName = res_dir+"Emma_new_trans_"+phenName+".sr.pvals"
		cp_command = "sudo cp "+oldName+" "+newName
		print cp_command
		os.system(cp_command)


def _transformPhenotypes_(phenotypeTransformations = {224:(0,1),225:(0,1),226:(1,1),227:(0,1),228:(1,1),229:(0,1),230:(1,1),231:(0,1),232:(1,1),233:(1,1),234:(1,1),235:(1,1),236:(0,1),237:(1,1),238:(1,1),239:(0,1),240:(1,1),241:(0,1),242:(1,1)}):
	"""
	Transforms phenotypes.  {phenotype_id: (type_transform,data_type)}
	"""
	
	#																			 32			   58						66							74								  90				   110 
	
	#transformationList = 4*[1]+[0]+2*[1]+[3]+6*[0]+6*[1]+2*[0]+[1]+5*[0]+3*[1]+7*[0]+13*[1]+[0]+4*[1]+[0]+5*[1]+[2]+[0]+[2]+[0]+[2]+2*[0]+[1]+[2]+[1]+[2]+[0]+2*[1]+3*[0]+3*[1]+2*[0]+4*[1]+2*[0]+7*[1]+9*[0]+2*[1]+[0]

	#							  116	   123		 126		   130							 146				   161	 163					  
	#transformationList +=4*[1]+[0]+[2]+6*[1]+[0]+[1]+[0]+[1]+2*[0]+[1]+[0]+2*[1]+2*[2]+2*[1]+[0]+8*[1]+[0]+5*[1]+2*[0]+7*[1]+[0]+[4]+[0]+2*[1]+17*[0]+2*[1]+2*[0]

	#					 187
	#transformationList +=4*[1]+4*[0]


	for p_i in phenotypeTransformations:
		print p_i, transformationTypes[phenotypeTransformations[p_i][0]]
	
	#datatypeList = 8*[1]+5*[3]+18*[1]+7*[3]+38*[1]+3*[3]+52*[1]+[2]+24*[1]+2*[2]+[1]+[2]+[3]+[2]+3*[1]+6*[3]+3*[2]+3*[3]+3*[2]+5*[1]+8*[1]

	
	phenData = getPhenotypes(user="bvilhjal", passwd="bamboo123", rawPhenotypes=True)
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
	for p_i in phenotypeTransformations: #range(0,len(phenData.phenotypeNames)):
		i = phenData.getPhenIndex(p_i)
		
		#Remove outliers...
		if phenotypeTransformations[p_i][0]==5:
			phenData.naOutliers(p_i)
			
		for j in range(0,len(phenData.accessions)):
			transFormedValue = _transformValue_(phenData.phenotypeValues[j][i],phenotypeTransformations[p_i][0])
			acc=phenData.accessions[j]
			#Insert info
			print phenData.phenotypeValues[j][i], transFormedValue, p_i, acc
			sqlStatm = "UPDATE stock_250k.phenotype_avg SET transformed_value="+str(transFormedValue)+" WHERE ecotype_id="+acc+" and method_id="+str(p_i)+""
			numRows = int(cursor.execute(sqlStatm))	
			conn.commit()
			
	for p_i in phenotypeTransformations:		
		sqlStatm = "UPDATE stock_250k.phenotype_method SET data_type='"+datatypeDict[phenotypeTransformations[p_i][1]]+"' WHERE id="+str(p_i)+""
		print sqlStatm
		numRows = int(cursor.execute(sqlStatm))	
		sqlStatm = "UPDATE stock_250k.phenotype_method SET transformation_description='"+transformationTypes[phenotypeTransformations[p_i][0]]+"' WHERE id="+str(p_i)+""
		print sqlStatm
		numRows = int(cursor.execute(sqlStatm))	
		conn.commit()

	print "Committing transaction (making changes permament)."
	print "Closing connection."
	#Close connection
	cursor.close()
	conn.close()


	#UPDATE stock_250k.phenotype_avg p SET p.value=250 WHERE p.method_id=7 and p.value=

def _transformValue_(value,transformType):
	if value!="NA":
		import math
		if transformType == 1 or transformType == 5:
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
	_fakeTransformPhenotypes_()
	#_run_()

