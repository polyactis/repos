#!/usr/bin/env python
"""
Usage: Accession2EcotypeComplete.py [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, at(default)
	-k ..., --schema=...	which schema in the database, dbsnp (default) IGNORE
	-l ...,	ecotype2accession_table (INPUT), 'at.ecotype2accession' (default)
	-e ...,	ecotype table (INPUT), 'ecotype' (default)
	-f ...,	ecotype_duplicate2tg_ecotypeid_table (INPUT), i.e 'stock20071008.ecotype_duplicate2tg_ecotypeid'
	-a ...,	accession2ecotype_complete_table (OUTPUT), 'at.accession2ecotype_complete' (default)
	-m ...,	accession2ecotype_table (OUTPUT), 'at.accession2ecotype' (default)
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	./src/Accession2EcotypeComplete.py -f stock20071008.ecotype_duplicate2tg_ecotypeid
	
	./src/Accession2EcotypeComplete.py -f stock20071008.ecotype_duplicate2tg_ecotypeid -c

Description:
	find all possible matching ecotype ids for at.accession ids.

	method 1: use nativename to identify other ecotype ids although stockparent might not be the original one. with the help of ecotype2accession.
		accession2ecotype_complete_table is the OUTPUT.
	method 2: use (nativename, stockparent) to identify ecotype ids that could be mapped to accession ids.
		accession2ecotype_table is the OUTPUT.
		would be same to at.accession2ecotype_complete(is_stockparent_original=1), except less coverage. some entries in at.ecotype2accession are not covered by ecotype_duplicate2tg_ecotypeid_table due to no genotyping or no GPS.
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
		ecotype2accession_table='at.ecotype2accession', ecotype_table='ecotype', ecotype_duplicate2tg_ecotypeid_table='', \
		accession2ecotype_complete_table='at.accession2ecotype_complete', accession2ecotype_table='at.accession2ecotype', \
		commit=0, debug=0, report=0):
		"""
		2007-10-29
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ecotype2accession_table = ecotype2accession_table
		self.ecotype_table = ecotype_table
		self.ecotype_duplicate2tg_ecotypeid_table = ecotype_duplicate2tg_ecotypeid_table
		self.accession2ecotype_complete_table = accession2ecotype_complete_table
		self.accession2ecotype_table = accession2ecotype_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def create_table1(self, curs, accession2ecotype_complete_table):
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
	
	def get_tg_ecotypeid2ecotypeid_set(self, curs, ecotype_duplicate2tg_ecotypeid_table='ecotype_duplicate2tg_ecotypeid'):
		"""
		2008-01-07
			take the grouping of duplicated strains in ecotype_duplicate2tg_ecotypeid_table
		"""
		sys.stderr.write("Getting ecotype_gp_ls ...")
		tg_ecotypeid2ecotypeid_set = {}
		from sets import Set
		curs.execute("select ecotypeid, duplicate, tg_ecotypeid from %s"%ecotype_duplicate2tg_ecotypeid_table)
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, duplicate, tg_ecotypeid = row
			if tg_ecotypeid not in tg_ecotypeid2ecotypeid_set:
				tg_ecotypeid2ecotypeid_set[tg_ecotypeid] = Set()
			tg_ecotypeid2ecotypeid_set[tg_ecotypeid].add(ecotypeid)
		sys.stderr.write("Done.\n")
		return tg_ecotypeid2ecotypeid_set
	
	def get_seed_ecotype_id2accession_id(self, curs, ecotype2accession_table='at.ecotype2accession'):
		"""
		2008-01-07
		"""
		sys.stderr.write("Getting seed_ecotype_id2accession_id ...")
		seed_ecotype_id2accession_id = {}
		curs.execute("select ecotype_id, accession_id from %s"%ecotype2accession_table)
		rows = curs.fetchall()
		for row in rows:
			ecotype_id, accession_id = row
			if ecotype_id in seed_ecotype_id2accession_id:
				sys.stderr.write("Error: ecotype_id %s already mapped to accession_id %s.\n"%(ecotype_id, ecotype_id2accession_id[ecotype_id]))
				sys.exit(2)
			seed_ecotype_id2accession_id[ecotype_id] = accession_id
		sys.stderr.write("Done.\n")
		return seed_ecotype_id2accession_id
	
	def map_tg_ecotypeid2ecotypeid_set(self, curs, tg_ecotypeid2ecotypeid_set, seed_ecotype_id2accession_id):
		"""
		2008-01-07
			connect more ecotype id to accession id conditioned on sharing same (nativename, stockparent)
				1. take the grouping of duplicated strains in ecotype_duplicate2tg_ecotypeid_table
				2. if the same-strain group are all(if any) mapped to a single accession id, they are all mapped to that accession.
		"""
		sys.stderr.write("Mapping tg_ecotypeid2ecotypeid_set to accession id ...")
		ecotype_id2accession_id = {}
		from sets import Set
		for tg_ecotypeid, ecotypeid_set in tg_ecotypeid2ecotypeid_set.iteritems():
			accession_id_set_mapped = Set()
			for ecotype_id in ecotypeid_set:
				if ecotype_id in seed_ecotype_id2accession_id:
					accession_id_set_mapped.add(seed_ecotype_id2accession_id[ecotype_id])
			if len(accession_id_set_mapped)==1:
				accession_id = accession_id_set_mapped.pop()
				for ecotype_id in ecotypeid_set:
					ecotype_id2accession_id[ecotype_id] = accession_id
			elif len(accession_id_set_mapped)>1:
				sys.stderr.write("Warning: set %s (tg_ecotypeid=%s) mapped to multiple accession_id: %s. Omitted.\n"%(ecotypeid_set, tg_ecotypeid, accession_id_set_mapped))
		sys.stderr.write("Done.\n")
		return ecotype_id2accession_id
		
	def create_table2(self, curs, accession2ecotype_table):
		"""
		2008-01-07
		"""
		sys.stderr.write("Creating %s ..."%accession2ecotype_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			accession_id	integer	not null,\
			ecotype_id	integer not null)"%accession2ecotype_table)
		sys.stderr.write("Done.\n")
	
	def submit_ecotype_id2accession_id(self, curs, accession2ecotype_table, ecotype_id2accession_id):
		"""
		2008-01-07
		"""
		sys.stderr.write("Submitting ecotype_id2accession_id...")
		for ecotype_id, accession_id in ecotype_id2accession_id.iteritems():
			curs.execute("insert into %s(accession_id, ecotype_id) values(%s, %s)"%\
			(accession2ecotype_table, accession_id, ecotype_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		accession_ecotype_ls = self.get_accession_ecotype_ls(curs, self.ecotype2accession_table, self.ecotype_table)
		if self.commit:
			self.create_table1(curs, self.accession2ecotype_complete_table)
			self.submit_accession_ecotype_ls(curs, self.accession2ecotype_complete_table, accession_ecotype_ls)
		
		#2008-01-07 method 2 of linking via (nativename, stockparent)
		tg_ecotypeid2ecotypeid_set = self.get_tg_ecotypeid2ecotypeid_set(curs, self.ecotype_duplicate2tg_ecotypeid_table)
		seed_ecotype_id2accession_id = self.get_seed_ecotype_id2accession_id(curs, self.ecotype2accession_table)
		ecotype_id2accession_id = self.map_tg_ecotypeid2ecotypeid_set(curs, tg_ecotypeid2ecotypeid_set, seed_ecotype_id2accession_id)
		if self.commit:
			self.create_table2(curs, self.accession2ecotype_table)
			self.submit_ecotype_id2accession_id(curs, self.accession2ecotype_table, ecotype_id2accession_id)

if __name__ == '__main__':
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:l:e:f:a:m:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype2accession_table = 'at.ecotype2accession'
	ecotype_table = 'ecotype'
	ecotype_duplicate2tg_ecotypeid_table = ''
	accession2ecotype_complete_table = 'at.accession2ecotype_complete'
	accession2ecotype_table = 'at.accession2ecotype'
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
		elif opt in ("-f",):
			ecotype_duplicate2tg_ecotypeid_table = arg
		elif opt in ("-a",):
			accession2ecotype_complete_table = arg
		elif opt in ("-m",):
			accession2ecotype_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if hostname and dbname and schema and ecotype2accession_table and ecotype_table and ecotype_duplicate2tg_ecotypeid_table and accession2ecotype_complete_table and accession2ecotype_table:
		instance = Accession2EcotypeComplete(hostname, dbname, schema, ecotype2accession_table, \
			ecotype_table, ecotype_duplicate2tg_ecotypeid_table, accession2ecotype_complete_table, \
			accession2ecotype_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
