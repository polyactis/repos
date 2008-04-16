#!/usr/bin/env python
"""
Usage:	Results2DB_250k.py [OPTIONS] -i INPUT_FNAME -e PHENOTYPE_METHOD_ID

Argument list:
	-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
	-d ..., --dbname=...	the database name, stock_250k(default)
	-u ..., --user=...	the db username, (otherwise it will ask for it).
	-p ..., --passwd=...	the db password, (otherwise it will ask for it).
	-i ...,	input_fname*	File containing association results
	-e ...,	phenotype_method_id*	which phenotype you used, check table phenotype_method
	-c,	commit db transaction
	-b,	toggle debug
	-r, toggle report
Examples:
	#test run without commiting database (no records in database)
	Results2DB_250k.py -i pvalue.log -e 1
	
	
	Results2DB_250k.py -i pvalue.log -e 1 -c
	
Description:
	This program would submit simple association results into database.

	The input file format is 3-column, tab-delimited:
		chromosome	position	score/pvalue

	The program will ask for following things.
		short_name:	give a short name of what you did. try to incorporate method, genotype data and phenotype data.
			It has to be <=30 characters and UNIQUE in table results_method. Check the table beforehand.
			Like 'kruskal_wallis_96_LD'.
		method_description:	a longer description of your method
		data_description:	which data you used
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv, stat, getopt
import traceback, gc, subprocess

"""
2008-04-16 temporarily put here
	-s ...,	short_name*	give a short name of what you did. try to incorporate method, genotype data and phenotype data
	-m ...,	method_description*	a longer description of your method
	-a ...,	data_description*	which data you used
"""

class Results2DB_250k(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	def __init__(self, **keywords):
		"""
		2008-04-16
		"""
		from pymodule import process_function_arguments
		argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('input_fname',1, ):None,\
								('short_name',1,):None,\
								('method_description',1,):None,\
								('data_description',1, ):None,\
								('phenotype_method_id',1,int):None,\
								('results_table',1, ):'results',\
								('results_method_table',1, ):'results_method',\
								('phenotype_method_table',1, ):'phenotype_method',\
								('commit',0, int):0,\
								('debug',0, int):0,\
								('report',0, int):0}
		"""
		2008-02-28
			argument_default_dict is a dictionary of default arguments, the key is a tuple, ('argument_name', is_argument_required, argument_type)
			argument_type is optional
		"""
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def check_if_phenotype_method_id_in_db(self, curs, phenotype_method_table, phenotype_method_id):
		"""
		"""
		curs.execute("select id from %s where id=%s"%(phenotype_method_table, phenotype_method_id))
		rows = curs.fetchall()
		if len(rows)>0:
			return 1
		else:
			return 0
	
	def submit_results_method(self, curs, results_method_table, short_name, method_description, data_description):
		"""
		2008-04-16
			submit the method part into db and return the id
		"""
		sys.stderr.write("Submitting results method ...")
		curs.execute("insert into %s(short_name, method_description, data_description) values ('%s', '%s', '%s')"%\
					(results_method_table, short_name, method_description, data_description))
		curs.execute("select id from %s where short_name='%s'"%\
					(results_method_table, short_name))
		rows = curs.fetchall()
		sys.stderr.write("Done.\n")
		return rows[0][0]
	
	def submit_results(self, curs, results_table, input_fname, results_method_id, phenotype_method_id):
		"""
		"""
		sys.stderr.write("Submitting results from %s ..."%(os.path.basename(input_fname)))
		reader = csv.reader(open(input_fname), delimiter='\t')
		for row in reader:
			if len(row)==3:
				chr, start_pos, score = row
				stop_pos = 'NULL'
			elif len(row)==4:
				chr, start_pos, stop_pos, score = row
			else:
				sys.stderr.write("ERROR: Found %s columns.\n"%(len(row)))
				sys.exit(3)
			curs.execute("insert into %s(chr, start_pos, stop_pos, score, method_id, phenotype_method_id) values (%s, %s, %s, %s, %s, %s)"%\
					(results_table, chr, start_pos, stop_pos, score, results_method_id, phenotype_method_id))
		del reader
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		2008-04-16
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		if not self.check_if_phenotype_method_id_in_db(curs, self.phenotype_method_table, self.phenotype_method_id):
			sys.stderr.write("Error: phenotype_method_id %s doesn't exist in table %s.\n"%(self.phenotype_method_id,self.phenotype_method_table))
			sys.stderr.exit(2)
		if not self.short_name:
			sys.stderr.write("Error: short_name unspecified.\n"%(self.short_name))
			sys.stderr.exit(2)
		results_method_id = self.submit_results_method(curs, self.results_method_table, self.short_name, self.method_description, self.data_description)
		self.submit_results(curs, self.results_table, self.input_fname, results_method_id, self.phenotype_method_id)
		if self.commit:
			curs.execute("commit")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "help", "type=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:u:p:i:e:cbr", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print __doc__
		sys.exit(2)
	
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	input_fname =None
	phenotype_method_id = None
	commit = None
	debug = None
	report = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			help = 1
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-u", "--user"):
			user = arg
		elif opt in ("-p", "--passwd"):
			passwd = arg
		elif opt in ("-e", ):
			phenotype_method_id = int(arg)
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = Results2DB_250k(hostname=hostname, dbname=dbname, user=user, passwd=passwd, input_fname=input_fname, \
							phenotype_method_id=phenotype_method_id,
							commit=commit, debug=debug, report=report)
	instance.run()
