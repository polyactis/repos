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

class Calls2DB_250k:
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
		"""
		sys.stderr.write("Getting chr_pos2snpid from %s ..."%(snps_table))
		chr_pos2snpid = {}
		curs.execute("select id, chr, position from %s"%(snps_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, chr, position = row
			chr_pos_key = (chr, position)
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


if __name__ == '__main__':
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:l:a:s:n:e:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	input_fname = ''
	calls_table = 'calls_250k'
	calls_comment_table = 'calls_250k_duplicate_comment'
	snps_table = 'snps_250k'
	strain_name = None
	comment = ''
	new_table = 0
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
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-l",):
			calls_table = arg
		elif opt in ("-a",):
			calls_comment_table = arg
		elif opt in ("-s",):
			snps_table = arg
		elif opt in ("-n",):
			strain_name = arg
		elif opt in ("-e",):
			comment = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if hostname and dbname and schema and input_fname and calls_table and calls_comment_table and snps_table and strain_name:
		instance = Calls2DB_250k(hostname, dbname, schema, input_fname,\
			calls_table, calls_comment_table, snps_table, strain_name, comment, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
