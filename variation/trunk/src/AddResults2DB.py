#!/usr/bin/env python2.5
"""
Usage: AddResults2DB.py [OPTIONS] RESULTS_FILE

Option:

	-z ..., --hostname=...	    the hostname, (papaya.usc.edu is default).
	-u ..., --user=...          the username, (otherwise it will ask for it).
	-p ..., --passwd=...	    the password, (otherwise it will ask for it).
	--callMethodID=...          * Use the call_method_id in the DB. 
        --resultsMethodID=...       The results_method_type in the DB. (Default is 1, i.e. association).
        --phenotypeMethodID=...     * The phenotype_method_id from the DB.
        --analysisMethodID=...      * The analysis_method_id from the DB.
	--createdBy=...             * Who created the data (created_by in the DB.)
	--shortName=...             * Short name for the result.
	--methodDescription=...     Description of the method (method_description in the DB).
	--dataDescription=...       Description of the data (data_description in the DB).
	--comment=...               Comments (comment in the DB.) 
	-h, --help	            Show this help

	* these options are required!
Examples:

	
Description:
        Add a result file to the DB.


Important method IDs: 
callMethodID:
        6 : The full dataset.  (This dataset should always be used.)
        7 : First 96 accessions.
	10 : phenotyped (data for which we have phenotypes). (Dataset is possibly outdated)
	12 : phenotyped and merged with 2010 data and MAF <5% removed. (Dataset is possibly outdated)
	13 : phenotyped and merged with 2010 data.  (Dataset is possibly outdated)

phenotypeMethodID:
        The number that the phenotype had e.g. LD is number 1, since it's called 1_LD in the phenotype file.

analysisMethodID:
        1: Kruskal-Wallis.
	2: Fisher's test (binary phenotypes)
	3: Chi-square test.
	4: Emma
	5: Margarita 
	6: Random Forest
"""

import sys, getopt, traceback, env

def _run_():
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["callMethodID=", "resultsMethodID=", "phenotypeMethodID=","analysisMethodID=", "createdBy=", "shortName=", "methodDescription=", "dataDescription=", "comment=", "help","hostname=", "user=", "passwd="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:u:p:h", long_options_list)

	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = 'papaya.usc.edu'
	user = 'bvilhjal'
	passwd = 'bamboo123'
	resultsFile = args[0]
	help = 0
	callMethodID = None
	resultsMethodID = 1
	phenotypeMethodID = None
	analysisMethodID = None
	createdBy  = None
	shortName = None
	methodDescription = ""
	dataDescription = ""
	comment = ""
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("--callMethodID"):
			callMethodID = int(arg)
		elif opt in ("--resultsMethodID"):
			resultsMethodID = int(arg)
		elif opt in ("--phenotypeMethodID"):
			phenotypeMethodID = int(arg)
		elif opt in ("--analysisMethodID"):
			analysisMethodID = int(arg)
		elif opt in ("--createdBy"):
			createdBy = arg
		elif opt in ("--shortName"):
			shortName = arg
		elif opt in ("--methodDescription"):
			methodDescription = arg
		elif opt in ("--dataDescription"):
			dataDescription = arg
		elif opt in ("--comment"):
			comment = arg
		else:
			if help==0:
				print "Unkown option!!\n"
				print __doc__
			sys.exit(2)

	if not resultsFile:
		if help==0:
			print "Result file missing!!\n"
			print __doc__
		sys.exit(2)

	addResultsToDB(resultsFile,hostname,user,passwd,callMethodID,phenotypeMethodID,analysisMethodID,createdBy,shortName,resultsMethodID,methodDescription,dataDescription,comment,commit=commit)

def addResultsToDB(resultsFile,hostname,user,passwd,callMethodID,phenotypeMethodID,analysisMethodID,createdBy,shortName,resultsMethodID=1,methodDescription="", dataDescription="" ,comment="", commit=False):
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
	print "Fetching result ID from the DB"
	numRows = int(cursor.execute("select max(id) from stock_250k.results_method"))	
	row = cursor.fetchone()
	resultsID = int(row[0])

	print "Reading results file:",resultsFile
#Convert resultsfile to a tsv file.
	f = open(resultsFile,"r")
	lines = f.readlines()
	f.close()
	
        #Write to a designated place.
	resultsDir = "/Network/Data/250k/db/results/type_1/"		
	tempNewFile = resultsDir+str(resultsID)+"_results.tsv"
	if env.user=="bvilhjal": #Working on cluster.
		tempResultsDir = "/home/cmb-01/bvilhjal/tmp/"#FIXME
		tempNewFile = tempResultsDir+str(resultsID)+"_results.tsv"
	newFile=resultsDir+str(resultsID)+"_results.tsv"
	
	start=0
	try:
		int(lines[0].split(",")[0])
	except Exception:
		start=1 #The file contains a header.

	print "Writing tsv file:",newFile
	f = open(tempNewFile,"w")
	for i in range(start,len(lines)):
		line = lines[i]
		line = line.replace(",","\t")
		f.write(line)
	f.close()
	
	print "Inserting results into DB."
	#Insert info
	sqlStatm = "insert into stock_250k.results_method (short_name, filename, method_description, data_description, phenotype_method_id, call_method_id, results_method_type_id, comment, created_by, analysis_method_id) values"
	sqlStatm += "('"+shortName+"','"+newFile+"','"+methodDescription+"','"+dataDescription+"',"+str(phenotypeMethodID)+","+str(callMethodID)+","+str(resultsMethodID)+",'"+comment+"','"+createdBy+"',"+str(analysisMethodID)+");"
	numRows = int(cursor.execute(sqlStatm))	
	print "Committing transaction (making changes permament)."
	conn.commit()
	
	print "Closing connection."
        #Close connection
        cursor.close()
	conn.close()
	if env.user=="bvilhjal":
		print "Remember to transfer files !!!!"

if __name__ == '__main__':
	#_test1_()
	_run_()


