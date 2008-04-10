#!/usr/bin/env python
"""
Usage: Calls2DB_250k.py [OPTIONS] -i input_fname -n strain_name

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-i ...,	input file, snp.RData outputted by write.table() in R
	-l ...,	calls_table, 'calls_250k'(default)
	-a ...,	calls_comment_table, 'calls_250k_duplicate_comment'(default)
	-s ...,	250k snp table, to find out which snpid based on chr+position, 'snps_250k'(default)
	-n ...,	strain name, it'll be used to match nativename and name to get ecotype id.
	-e ...,	comment for this loading
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Calls2DB_250k.py -i /Network/Data/250k/yanli9-11-07/Tamm2B_base-calls.txt -n Tamm  -e "yanli9-11-07/Tamm2B_base-calls.txt" -c
Description:

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, getopt, csv

class Calls2DB_250k_old:
	"""
	2007-12-11
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		input_fname='', calls_table='calls_250k', calls_comment_table='calls_250k_duplicate_comment', snps_table='snps_250k', strain_name=None, \
		comment='', commit=0, debug=0, report=0):
		"""
		2007-12-11
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.calls_table = calls_table
		self.calls_comment_table = 'calls_250k_duplicate_comment'
		self.snps_table = snps_table
		self.strain_name = strain_name
		self.comment = comment
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def create_calls_table(self, curs, calls_table):
		"""
		2007-12-11
			won't be used. 2 tables, calls_250k_duplicate_comment and calls_250k will
			be created separately by mysql.
		"""
		sys.stderr.write("Creating %s ..."%calls_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			ecotypeid	integer,\
			snpid	integer,\
			snpcall	varchar(2),\
			duplicate	integer,\
			date_created	timestamp default current_timestamp,\
			date_modified	timestamp)"%calls_table)
		sys.stderr.write("Done.\n")
	
	def create_calls_comment_table(self, curs, calls_comment_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Creating %s ..."%calls_comment_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			ecotypeid	integer,\
			duplicate	integer,\
			comment	varchar(2000),\
			date_created	timestamp default current_timestamp,\
			date_modified	timestamp)"%calls_comment_table)
		sys.stderr.write("Done.\n")
	
	def get_chr_pos2snpid(self, curs, snps_table):
		"""
		2007-12-11
		2007-12-13
			'chr' in snps_table changed to 'chromosome'
		"""
		sys.stderr.write("Getting chr_pos2snpid from %s ..."%(snps_table))
		chr_pos2snpid = {}
		curs.execute("select id, chromosome, position from %s"%(snps_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, chromosome, position = row
			chr_pos_key = (chromosome, position)
			chr_pos2snpid[chr_pos_key] = snpid
		sys.stderr.write("Done.\n")
		return chr_pos2snpid
	
	def find_out_ecotypeid_given_strain_name(self, curs, strain_name, ecotype_table='ecotype'):
		sys.stderr.write("Trying to find ecotypeid based on name=%s ...\n"%(strain_name))
		curs.execute("select e.id, e.name, e.stockparent, e.nativename, s.name, c.abbr from %s e, site s, address a, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and (e.nativename rlike '%s'or e.name rlike '%s')"%(ecotype_table, strain_name, strain_name))
		rows = curs.fetchall()
		choice_id2ecotypeid = {}
		choice_id = 0
		header_ls = ['\t', 'ecotypeid', 'name', 'stockparent', 'nativename', 'site_name', 'country']
		sys.stderr.write("%s\n"%('\t'.join(header_ls)))
		for row in rows:
			ecotypeid, name, stockparent, nativename, site_name, country = row
			choice_id += 1
			choice_id2ecotypeid[repr(choice_id)] = ecotypeid
			ls = ['%s:'%choice_id]
			ls.append(list(row))
			ls = map(repr, ls)
			sys.stderr.write("%s\n"%'\t'.join(ls))
		if len(choice_id2ecotypeid)==0:
			sys.stderr.write("Found no match for %s.\n"%strain_name)
			sys.exit(3)
		choice_id = raw_input("Enter the row number corresponding to the ecotypeid(q for exit):")
		while 1:
			if choice_id =='q':
				sys.exit(2)
			elif choice_id in choice_id2ecotypeid:
				break
			else:
				sys.stderr.write("%s is not a choice.\n"%(choice_id))
				choice_id = raw_input("Enter the row number corresponding to the ecotypeid(q for exit):")
		return choice_id2ecotypeid[choice_id]
	
	def find_out_duplicate_id_given_ecotypeid(self, curs, ecotypeid, calls_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Finding out the duplicate_id for ecotypeid=%s ..."%(ecotypeid))
		curs.execute("select distinct duplicate from %s where ecotypeid=%s"%(calls_table, ecotypeid))
		rows = curs.fetchall()
		duplicate_ls = []
		for row in rows:
			duplicate_ls.append(row[0])
		if duplicate_ls:
			duplicate_id = max(duplicate_ls) + 1	#1 plus the maximum
		else:
			duplicate_id = 1
		sys.stderr.write("is %s.\n"%duplicate_id)
		return duplicate_id
	
	def get_calls(self, input_fname, chr_pos2snpid):
		"""
		2007-12-11
		"""
		sys.stderr.write("Getting calls ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		calls_ls = []
		reader.next()	#toss out the 1st header row
		for row in reader:
			chr, position, allele1, allele2, antisense1, sense1, antisense2, sense2, snpcall = row
			chr = int(chr)
			position = int(position)
			chr_pos_key = (chr, position)
			if snpcall=='?':	#2007-12-11 '?' was used as NA
				snpcall = 'N'
			snpid = chr_pos2snpid[chr_pos_key]
			calls_ls.append([snpid, snpcall])
		del reader
		sys.stderr.write("Done.\n")
		return calls_ls
	
	def submit_comment(self, curs, ecotypeid, duplicate_id, comment, calls_comment_table):
		"""
		2007-12-11
		"""
		sys.stderr.write("Submitting comment ...")
		curs.execute("insert into %s(ecotypeid, duplicate, comment) values(%s, %s, '%s')"%(calls_comment_table, ecotypeid, duplicate_id, comment))
		sys.stderr.write("Done.\n")
	
	def submit_calls_ls(self, curs, calls_ls, ecotypeid, duplicate_id, calls_table):
		sys.stderr.write("Submitting calls_ls ...")
		for snpid, snpcall in calls_ls:
			curs.execute("insert into %s(ecotypeid, snpid, snpcall, duplicate) values(%s, %s, '%s', %s)"%(calls_table, ecotypeid, snpid, snpcall, duplicate_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		ecotypeid = self.find_out_ecotypeid_given_strain_name(curs, self.strain_name)
		duplicate_id = self.find_out_duplicate_id_given_ecotypeid(curs, ecotypeid, self.calls_table)
		chr_pos2snpid = self.get_chr_pos2snpid(curs, self.snps_table)
		calls_ls = self.get_calls(self.input_fname, chr_pos2snpid)
		if self.commit:
			self.submit_comment(curs, ecotypeid, duplicate_id, self.comment, self.calls_comment_table)
			self.submit_calls_ls(curs, calls_ls, ecotypeid, duplicate_id, self.calls_table)

import getopt, csv
import traceback, gc
from pymodule import process_function_arguments
class Calls2DB_250k(object):
	"""
	2008-04-08
	Calls2DB_250k.py [OPTIONS] -i INPUT_DIR -m METHOD_ID
	
	Argument list:
		-z ..., --hostname=...	the hostname, papaya.usc.edu(default)
		-d ..., --dbname=...	the database name, stock_250k(default)
		-u ..., --user=...	the username, (otherwise it will ask for it).
		-p ..., --password=...	the password, (otherwise it will ask for it).
		-i ...,	input_dir*	The directory containing calling algorithm output files.
		-m ...,	method_id*
		-o ...,	output_dir to store the renamed call files. '/Network/Data/250k/db/calls/'(default)
		-t ...,	call_info_table, 'stock_250k.call_info'(default)
		-a ...,	call_method_table, 'stock_250k.call_method'(default)
		-c,	commit db transaction
		-b,	toggle debug
		-r, toggle report
	Examples:
		Calls2DB_250k.py -i /tmp/simplecalls -m 1 -c
	
	Description:
		Turn calling algorithm's results into db and associated filesystem directory.
		
		Each file in input_dir shall be named like 'array_id'_call.tsv.
		
		The format is 2-column and tab-delimited. example:
			SNP_ID	'array_id'
			1_657_C_T	C
	
	"""
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		argument_default_dict = {('hostname',1, ):'papaya.usc.edu',\
								('dbname',1, ):'stock_250k',\
								('user',1, ):None,\
								('passwd',1, ):None,\
								('input_dir',1, ):None,\
								('output_dir',1, ):'/Network/Data/250k/db/calls/',\
								('method_id',1,):None,\
								('call_info_table',1, ):'stock_250k.call_info',\
								('call_method_table',1, ):'stock_250k.call_method',\
								('array_info_table',1, ):'stock_250k.array_info',\
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
		self.cur_max_call_id = None
	
	def get_cur_max_call_id(self, curs, call_info_table):
		"""
		get current maximum call id in db
		"""
		sys.stderr.write("Getting current maximum call id in db.\n")
		curs.execute("select max(id) from %s"%(call_info_table))
		rows = curs.fetchall()
		cur_max_call_id = rows[0][0]
		if cur_max_call_id!=None:
			return cur_max_call_id
		else:
			return 0
	
	def check_method_id_exists(self, curs, call_method_table, method_id):
		"""
		"""
		curs.execute("select id from %s where id=%s"%(call_method_table, method_id))
		rows = curs.fetchall()
		if len(rows)==0:
			return 0
		else:
			return 1
	
	def get_new_call_id(self, curs, call_info_table, array_id, method_id):
		"""
		"""
		if self.cur_max_call_id==None:
			self.cur_max_call_id = self.get_cur_max_call_id(curs, call_info_table)
		curs.execute("select id from %s where array_id=%s and method_id=%s"%(call_info_table, array_id, method_id))
		rows = curs.fetchall()
		if len(rows)>0:
			sys.stderr.write("\tarray_id=%s and method_id=%s already exists in %s. Ignored.\n"%(array_id, method_id, call_info_table))
			return -1
		else:
			self.cur_max_call_id += 1
			return self.cur_max_call_id
	
	def submit_one_call_entry(self, curs, call_info_table, call_id, filename, array_id, method_id):
		"""
		2008-04-08
			not used right now.
		"""
		pass
	
	def submit_call_dir2db(self, curs, input_dir, call_info_table, output_dir, method_id, user):
		"""
		2008-04-09
			add method_id as sub-directory
			submit user into table as 'created_by'
		2008-04-08
		"""
		output_dir = os.path.join(output_dir, '%s'%method_id)
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		file_ls = os.listdir(input_dir)
		sys.stderr.write("\n\tTotally, %d files to be processed.\n"%len(file_ls))
		file_ls.sort()
		for i in range(len(file_ls)):
			filename = file_ls[i]
			array_id = filename.split('_')[0]
			sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(file_ls),filename))
			new_call_id = self.get_new_call_id(curs, call_info_table, array_id, method_id)
			if new_call_id!=-1:
				output_fname = os.path.join(output_dir, '%s_call.tsv'%new_call_id)
				curs.execute("insert into %s(id, filename, array_id, method_id, created_by) values (%s, '%s', %s, %s, '%s')"%\
						(call_info_table, new_call_id, output_fname, array_id, method_id, user))
				input_fname = os.path.join(input_dir, filename)
				pipe_f = os.popen('cp %s %s'%(input_fname, output_fname))
				pipe_f_out = pipe_f.read()
				if pipe_f_out:
					sys.stderr.write("\tcp output: %s\n"%pipe_f_out)
	
	def run(self):
		"""
		2008-04-08
			-check_method_id_exists()
			-submit_call_dir2db()
				-get_new_call_id()
					-get_cur_max_call_id()
		"""
		
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		if not self.check_method_id_exists(curs, self.call_method_table, self.method_id):
			sys.stderr.write("Error: method_id=%s not in %s. The method has to be put into db beforehand.\n"%\
							(self.method_id, self.call_method_table))
			sys.exit(2)
		if self.commit:
			self.submit_call_dir2db(curs, self.input_dir, self.call_info_table, self.output_dir, self.method_id, self.user)
			curs.execute("commit")
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print Calls2DB_250k.__doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "user=", "passwd=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:u:p:i:m:o:t:a:cbrh", long_options_list)
	except:
		traceback.print_exc()
		print sys.exc_info()
		print Calls2DB_250k.__doc__
		sys.exit(2)
	
	hostname = None
	dbname = None
	user = None
	passwd = None
	input_dir = None
	method_id = None
	output_dir = None
	call_info_table = None
	call_method_table = None
	commit = 0
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
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
		elif opt in ("-i",):
			input_dir = arg
		elif opt in ("-m",):
			method_id = arg
		elif opt in ("-o",):
			output_dir = arg
		elif opt in ("-t",):
			call_info_table = arg
		elif opt in ("-a",):
			call_method_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	instance = Calls2DB_250k(hostname=hostname, dbname=dbname, user=user, passwd=passwd, input_dir=input_dir, output_dir=output_dir,
					method_id=method_id, call_info_table=call_info_table, call_method_table=call_method_table,\
					commit=commit, debug=debug, report=report)
	instance.run()
