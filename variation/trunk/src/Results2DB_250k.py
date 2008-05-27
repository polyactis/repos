#!/usr/bin/env python
"""

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
from variation.src.db import Results, ResultsMethod, Stock_250kDatabase, PhenotypeMethod, CallMethod, SNPs, ResultsMethodType
from pymodule import figureOutDelimiter

import sqlalchemy as sql
"""
2008-04-16 temporarily put here
	-s ...,	short_name*	give a short name of what you did. try to incorporate method, genotype data and phenotype data
	-m ...,	method_description*	a longer description of your method
	-a ...,	data_description*	which data you used
"""

class Results2DB_250k(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('input_fname',1, ): [None, 'i', 1, 'File containing association results'],\
							('phenotype_method_id',1,int): [None, 'e', 1, 'which phenotype you used, check table phenotype_method'],\
							('call_method_id',1,int): [None,],\
							('results_method_type_id',1,int): [None,],\
							('comment',0, ): [None, ],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	"""
	04/28/08 no longer needed
							('results_table',1, ): 'results',\
							('results_method_table',1, ):'results_method',\
							('phenotype_method_table',1, ):'phenotype_method',\	
	"""
	def __init__(self, **keywords):
		"""
		2008-04-28
			use ProcessOptions, newer option handling class
		2008-04-16
		"""
		#this is just for the program to ask user
		more_function_argument_dict = {('short_name',1,): None,\
							('method_description',1,):None,\
							('data_description',1, ):None}
		more_function_argument_dict.update(self.option_default_dict)
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, more_function_argument_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
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
		
	marker_pos2snp_id = None
	is_new_marker_added = False	#2008-05-26 flag whether new markers were generated.
	def reset_marker_pos2snp_id(cls):
		"""
		2008-05-26
			after commit or rollback in plone, session is closed and those new marker objects are gone. need to reset everything.
		"""
		if cls.is_new_marker_added:
			del cls.marker_pos2snp_id
			cls.marker_pos2snp_id = cls.get_marker_pos2snp_id(db)
	
	reset_marker_pos2snp_id = classmethod(reset_marker_pos2snp_id)
	
	def get_marker_pos2snp_id(cls, db):
		"""
		2008-05-24
		"""
		sys.stderr.write("Getting marker_pos2snp_id ...")
		marker_pos2snp_id = {}
		snps_table = db.tables['snps'].alias()
		conn = db.connection
		results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
		for row in results:
			key = (row.chromosome, row.position, row.end_position)
			marker_pos2snp_id[key] = row.id
		sys.stderr.write("Done.\n")
		return marker_pos2snp_id
	get_marker_pos2snp_id = classmethod(get_marker_pos2snp_id)
	
	def submit_results(cls, db, input_fname, rm, user):
		"""
		2008-05-26
			input_fname from plone is not file object although it has file object interface.
		2008-05-26
			csv.Sniffer() can't figure out delimiter if '\n' is in the string, use own dumb function figureOutDelimiter()
		2008-05-25
			save marker(snps) in database if it's not there.
			use marker id in results table
		2008-05-24
			figure out delimiter automatically
			input_fname could be a file object (from plone)
			phenotype method doesn't go with results anymore. it goes with results_method
		2008-04-28
			changed to use Stock_250kDatabase (SQLAlchemy) to do db submission
		"""
		if isinstance(input_fname, str) and os.path.isfile(input_fname):
			sys.stderr.write("Submitting results from %s ..."%(os.path.basename(input_fname)))
			delimiter = figureOutDelimiter(input_fname)
			reader = csv.reader(open(input_fname), delimiter=delimiter)
		elif hasattr(input_fname, 'readline') or hasattr(input_fname, 'read'):	#input_fname is not a file name, but direct file object. it could also be <ZPublisher.HTTPRequest.FileUpload instance at 0xa1774f4c>
			sys.stderr.write("Submitting results from %s on plone ..."%input_fname.filename)
			cs = csv.Sniffer()
			input_fname.seek(0)	#it's already read by plone to put int data['input_fname'], check results2db_250k.py
			if getattr(input_fname, 'readline', None) is not None:
				test_line = input_fname.readline()
				delimiter = cs.sniff(test_line).delimiter
			else:
				test_line = input_fname.read(200)
				delimiter = figureOutDelimiter(test_line)	#counting is a safer solution. if test_line include '\n', cs.sniff() won't figure it out.
			input_fname.seek(0)
			reader = csv.reader(input_fname, delimiter=delimiter)
		else:
			sys.stderr.write("Error: %s is neither a file name nor a file object.\n"%input_fname)
			return
		if cls.marker_pos2snp_id is None:
			cls.marker_pos2snp_id = cls.get_marker_pos2snp_id(db)
		session = db.session
		for row in reader:
			chr = int(row[0])
			start_pos = int(row[1])
			if len(row)==3:
				stop_pos = None
				score = row[2]
				marker_name = '%s_%s'%(chr, start_pos)
			elif len(row)==4:
				stop_pos = int(row[2])
				score = row[3]
				marker_name = '%s_%s_%s'%(chr, start_pos, stop_pos)
			else:
				sys.stderr.write("ERROR: Found %s columns.\n"%(len(row)))
				sys.exit(3)
			
			key = (chr, start_pos, stop_pos)
			if key in cls.marker_pos2snp_id:
				snps_id = cls.marker_pos2snp_id[key]
				if isinstance(snps_id, SNPs):	#it's a new marker object
					r = Results(score=score)
					r.snps = snps_id
				else:	#others are all integer ids
					r = Results(snps_id=snps_id, score=score)
			else:
				#construct a new marker
				marker = SNPs(name=marker_name, chromosome=chr, position=start_pos, end_position=stop_pos, created_by=user)
				#save it in database to get id
				session.save(marker)
				cls.marker_pos2snp_id[key] = marker	#for the next time to encounter same marker
				
				r = Results(score=score)
				r.snps = marker
				del marker
			r.results_method = rm
			session.save(r)
			del r
			#curs.execute("insert into %s(chr, start_pos, stop_pos, score, method_id, phenotype_method_id) values (%s, %s, %s, %s, %s, %s)"%\
			#		(results_table, chr, start_pos, stop_pos, score, results_method_id, phenotype_method_id))
		del reader
		sys.stderr.write("Done.\n")
	
	submit_results = classmethod(submit_results)
	
	def plone_run(cls, db, short_name, phenotype_method_id, call_method_id, data_description, \
				method_description, comment, input_fname, user, results_method_type_id=None, \
				results_method_type_short_name=None, commit=0):
		"""
		2008-05-26
			add results_method_type_id and results_method_type_short_name
		2008-05-24
			to conveniently wrap up all codes so that both this program and plone can call
		"""
		session = db.session
		if getattr(db, 'transaction', None) is None:
			db.transaction = session.create_transaction()
		
		#pm = session.query(PhenotypeMethod).get_by(id=phenotype_method_id)
		#if not pm:	#no corresponding phentype method id
		#	phenotype_method_id = None
		
		#cm = session.query(CallMethod).get_by(id=call_method_id)
		
		rmt = session.query(ResultsMethodType).get_by(id=results_method_type_id)
		if not rmt and results_method_type_short_name is not None:	#create a new results method type
			rmt = ResultsMethodType(short_name=results_method_type_short_name)
			session.save(rmt)
		
		rm = ResultsMethod(short_name=short_name, method_description=method_description, \
						data_description=data_description, comment=comment, phenotype_method_id=phenotype_method_id,\
						call_method_id=call_method_id, created_by=user)
		if rmt:
			rm.results_method_type = rmt
		
		cls.submit_results(db, input_fname, rm, user)
		#session.flush()	#not necessary as no immediate query on the new results after this and commit() would execute this.
		if commit:
			#curs.execute("commit")
			db.transaction.commit()
			session.clear()	#Remove all object instances from this Session
			del session
			db.transaction = None	#delete the transaction
			cls.reset_marker_pos2snp_id()
		#else:	#default is also rollback(). to demonstrate good programming
		#	transaction.rollback()
	plone_run = classmethod(plone_run)
	
	def run(self):
		"""
		2008-04-28
			use Stock_250kDatabase to do database stuff
		2008-04-16
		"""
		#import MySQLdb
		#conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		#curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		
		db = Stock_250kDatabase(username=self.user,
				   password=self.passwd, hostname=self.hostname, database=self.dbname)
		if not self.short_name:
			sys.stderr.write("Error: short_name unspecified.\n"%(self.short_name))
			sys.stderr.exit(2)
		self.plone_run(db, self.short_name, self.phenotype_method_id, self.call_method_id, self.data_description, \
				self.method_description, self.comment, self.input_fname, self.user, results_method_type_id=self.results_method_type_id,\
				commit=self.commit)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Results2DB_250k
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
	"""
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
	"""
