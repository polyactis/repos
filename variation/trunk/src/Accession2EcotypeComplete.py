#!/usr/bin/env python
"""
Usage: Accession2EcotypeComplete.py [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, at(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-l ...,	ecotype2accession_table, 'at.ecotype2accession' (default)
	-e ...,	ecotype table, 'ecotype' (default)
	-a ...,	accession2ecotype_complete_table, 'at.accession2ecotype_complete' (default)
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/Accession2EcotypeComplete.py -h
	
	./src/Accession2EcotypeComplete.py -c

Description:
	find all possible matching ecotype ids for at.accession ids with the help of ecotype2accession. use nativename to identify other ecotype ids although stockparent might be the original one. 
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, getopt

class Accession2EcotypeComplete:
	"""
	2007-10-29
		create a complete table linking accession table to ecotype table in at
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		ecotype2accession_table='at.ecotype2accession', ecotype_table='ecotype',\
		accession2ecotype_complete_table='at.accession2ecotype_complete', commit=0, debug=0, report=0):
		"""
		2007-10-29
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ecotype2accession_table = ecotype2accession_table
		self.ecotype_table = ecotype_table
		self.accession2ecotype_complete_table = accession2ecotype_complete_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def create_table(self, curs, accession2ecotype_complete_table):
		sys.stderr.write("Creating %s ..."%accession2ecotype_complete_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			accession_id	integer	not null,\
			ecotype_id	integer not null,\
			is_stockparent_original	integer)"%accession2ecotype_complete_table)
		sys.stderr.write("Done.\n")
	
	def get_accession_ecotype_ls(self, curs, ecotype2accession_table, ecotype_table):
		"""
		2007-10-29
			get the link via same nativename
		"""
		sys.stderr.write("Getting accession_ecotype_ls ...")
		curs.execute("select ea.accession_id, e.id, e.stockparent, e2.id as eid, e2.stockparent as sp from %s e, %s e2 , %s ea where e.id=ea.ecotype_id and e.nativename=e2.nativename order by id, eid"%(ecotype_table, ecotype_table, ecotype2accession_table))
		rows = curs.fetchall()
		accession_ecotype_ls = []
		for row in rows:
			accession_id, ecotype_id, original_sp, eid, sp = row
			if sp!=original_sp:
				is_stockparent_original = 0
			else:
				is_stockparent_original = 1
			accession_ecotype_ls.append([accession_id, eid, is_stockparent_original])
		sys.stderr.write("Done.\n")
		return accession_ecotype_ls
	
	def submit_accession_ecotype_ls(self, curs, accession2ecotype_complete_table, accession_ecotype_ls):
		sys.stderr.write("Submitting  accession_ecotype_ls...")
		for accession_id, ecotype_id, is_stockparent_original in accession_ecotype_ls:
			curs.execute("insert into %s(accession_id, ecotype_id, is_stockparent_original) values(%s, %s, %s)"%\
			(accession2ecotype_complete_table, accession_id, ecotype_id, is_stockparent_original))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		accession_ecotype_ls = self.get_accession_ecotype_ls(curs, self.ecotype2accession_table, self.ecotype_table)
		if self.commit:
			self.create_table(curs, self.accession2ecotype_complete_table)
			self.submit_accession_ecotype_ls(curs, self.accession2ecotype_complete_table, accession_ecotype_ls)


if __name__ == '__main__':
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:l:e:a:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype2accession_table = 'at.ecotype2accession'
	ecotype_table = 'ecotype'
	accession2ecotype_complete_table = 'at.accession2ecotype_complete'
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
		elif opt in ("-l",):
			ecotype2accession_table = arg
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-a",):
			accession2ecotype_complete_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if hostname and dbname and schema and ecotype2accession_table and ecotype_table and accession2ecotype_complete_table:
		instance = Accession2EcotypeComplete(hostname, dbname, schema, ecotype2accession_table, \
			ecotype_table, accession2ecotype_complete_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
